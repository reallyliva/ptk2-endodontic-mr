library(data.table)
library(TwoSampleMR)
dir <- "##"

## ============================================================
## 0. 载入 5447 干净工具集
## ============================================================
inst <- as.data.table(readRDS(file.path(dir, "exposure_instruments_FINAL_r0.001_clean.rds")))
stopifnot(nrow(inst) == 5447)
cat("工具集:", nrow(inst), "SNP /", uniqueN(inst$Gene), "基因\n\n")

## ============================================================
## 1. 桥文件 rsID → hg38 pos（FinnGen 无 rsID，用 K11_PULP_PERIAPICAL 回补）
##    坑#8：rsids 可能逗号分隔一格多 ID，先 explode 再匹配
## ============================================================
br <- fread(file.path(dir, "summary_stats_release_finngen_R12_K11_PULP_PERIAPICAL.gz"),
            select = c("rsids", "#chrom", "pos"))
setnames(br, "#chrom", "chrom")
cm     <- grepl(",", br$rsids, fixed = TRUE)
bridge <- rbind(br[!cm],
                br[cm][, .(rsids = unlist(strsplit(rsids, ",", fixed = TRUE))),
                       by = .(chrom, pos)])
bridge <- unique(bridge[rsids %in% inst$SNP], by = "rsids")
bridge[, key38 := paste(chrom, pos, sep = ":")]
cat("桥接到 hg38 的工具 SNP:", nrow(bridge), "/", uniqueN(inst$SNP), "\n\n")
rm(br, cm); gc()

## ============================================================
## 2. 暴露端（TwoSampleMR 格式，沿用 5447 里的原始 beta/se/eaf）
## ============================================================
exp_df <- inst[, .(
  SNP,
  beta.exposure          = beta_exp,
  se.exposure            = se_exp,
  effect_allele.exposure = AssessedAllele,
  other_allele.exposure  = OtherAllele,
  eaf.exposure           = eaf,
  pval.exposure          = Pvalue,
  id.exposure            = Gene,
  exposure               = Gene)]
gene_map <- unique(inst[, .(id.exposure = Gene, GeneSymbol)])

## ============================================================
## 3. 单结局 MR：harmonise(action=2, palindromic 按频率对齐) + IVW/Wald + BH-FDR
## ============================================================
run_mr <- function(gz, oc_name) {
  cat("==========", oc_name, "==========\n")
  o <- fread(file.path(dir, gz)); setnames(o, "#chrom", "chrom", skip_absent = TRUE)
  o[, key38 := paste(chrom, pos, sep = ":")]
  o <- o[key38 %in% bridge$key38]
  o <- merge(o, bridge[, .(key38, rsids)], by = "key38", allow.cartesian = TRUE)
  o <- unique(o, by = "rsids")

  out_df <- o[, .(
    SNP                   = rsids,
    beta.outcome          = beta,
    se.outcome            = sebeta,
    effect_allele.outcome = alt,      # FinnGen 效应量以 alt 计
    other_allele.outcome  = ref,
    eaf.outcome           = af_alt,
    pval.outcome          = pval,
    id.outcome            = oc_name,
    outcome               = oc_name)]
  rm(o); gc()

  harm <- harmonise_data(exp_df, out_df, action = 2)     # ★ 默认 action=2
  res  <- as.data.table(mr(harm, method_list = c("mr_ivw", "mr_wald_ratio")))
  res  <- merge(res, gene_map, by = "id.exposure", all.x = TRUE)
  res[, `:=`(OR = exp(b), fdr = p.adjust(pval, "fdr"))]
  setorder(res, pval)

  cat("检验基因数:", nrow(res), "| FDR<0.05:", sum(res$fdr < 0.05, na.rm = TRUE), "\n\n")
  saveRDS(res, file.path(dir, paste0("MR_results_", oc_name, "_5447.rds")))
  res
}

r_nec <- run_mr("necrosis_or_apical_periodontitis.gz", "necrosis")
r_k04 <- run_mr("pulpal_and_apical_diseases.gz",        "K04")
r_pul <- run_mr("pulpitis.gz",                          "pulpitis")

## ============================================================
## 4. 硬校验：必须复现 2545 / 19 / ROCK1 未过 / PTK2 稳
## ============================================================
cat("########## 保真校验 ##########\n")
cat("necrosis 检验数:", nrow(r_nec), "(应 2545) | FDR<0.05:",
    sum(r_nec$fdr < 0.05, na.rm = TRUE), "(应 19)\n")
print(r_nec[GeneSymbol %in% c("PTK2","ROCK1","DLK1","SIK2"),
            .(GeneSymbol, nsnp, b = round(b,4), pval = signif(pval,3),
              OR = round(OR,3), fdr = signif(fdr,4), method)])

stopifnot(
  nrow(r_nec) == 2545,
  sum(r_nec$fdr < 0.05, na.rm = TRUE) == 19,
  abs(r_nec[GeneSymbol == "PTK2"]$b - 0.1139) < 1e-3,   # PTK2 稳
  r_nec[GeneSymbol == "PTK2"]$fdr < 1e-3,
  r_nec[GeneSymbol == "ROCK1"]$fdr > 0.05,              # ★ ROCK1 临界未过（≈0.0521）
  r_nec[GeneSymbol == "DLK1"]$fdr < 0.05                # DLK1 三表型过（necrosis）
)
cat("\n✅ 复现成功：2545 / 19，ROCK1 临界未过，PTK2 稳\n")

## necrosis 的 19 名单（收口 coloc 用这个）
cat("\n########## necrosis FDR<0.05（19 名单，coloc 依据）##########\n")
print(r_nec[fdr < 0.05, .(GeneSymbol, nsnp, b = round(b,3),
                          OR = round(OR,3), fdr = signif(fdr,3))], row.names = FALSE)

## ============================================================
## 5. 四靶 × 三表型（分层依据，对齐手稿）
## ============================================================
tetra <- c("PTK2","ROCK1","DLK1","SIK2")
show <- function(d, nm) d[GeneSymbol %in% tetra,
  .(outcome = nm, GeneSymbol, nsnp, b = round(b,3),
    OR = round(OR,3), fdr = signif(fdr,3), method)]
cat("\n########## 四靶 × 三表型 ##########\n")
print(rbind(show(r_nec,"necrosis"), show(r_k04,"K04"), show(r_pul,"pulpitis"))[
        order(match(GeneSymbol, tetra))], row.names = FALSE)