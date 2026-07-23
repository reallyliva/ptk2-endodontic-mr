Stage 0 · 原始输入

NIHMS80906-supplement-Table_S1.xlsx(Finan,sheet "Data")
eQTLGen cis-eQTL .gz
2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz
EUR(1000G Phase3 EUR,bed/bim/fam)+ plink.exe(1.9)

Stage 1 · 可成药基因 × eQTLGen 取交集 → eqtl_druggable_cis.rds

r
library(data.table); library(readxl)
dir <- "##"

finan <- as.data.table(read_excel(
  file.path(dir, "NIHMS80906-supplement-Table_S1.xlsx"), sheet = "Data"))

## ← 换成磁盘上 eQTLGen 文件的实际名字
eqtl <- fread(file.path(dir,
  "2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"))

eqtl_drug <- eqtl[Gene %in% finan$ensembl_gene_id]
eqtl_drug <- merge(eqtl_drug,
  finan[, .(ensembl_gene_id, druggability_tier, hgnc_names,
            small_mol_druggable, bio_druggable)],
  by.x = "Gene", by.y = "ensembl_gene_id", all.x = TRUE)

saveRDS(eqtl_drug, file.path(dir, "eqtl_druggable_cis.rds"))
cat("Finan 基因:", nrow(finan),
    "| 交集基因:", uniqueN(eqtl_drug$Gene),
    "| cis-eQTL 行:", nrow(eqtl_drug), "\n")

预期:Finan 4479;交集 2715 基因 / 1,891,727 行。

Stage 2 · 合频率 → Z→beta/se → F>10 → eqtl_exposure_prepped.rds

r
library(data.table)
dir <- "##"
eqtl <- readRDS(file.path(dir, "eqtl_druggable_cis.rds"))

af <- fread(file.path(dir,
  "2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz"),
  select = c("SNP", "AlleleA", "AlleleB", "AlleleB_all"))

## eaf 方向：AssessedAllele==AlleleB 用 AlleleB_all，==AlleleA 用 1-它，都不是→NA
eqtl <- merge(eqtl, af, by = "SNP", all.x = TRUE)
eqtl[, eaf := fifelse(AssessedAllele == AlleleB, AlleleB_all,
              fifelse(AssessedAllele == AlleleA, 1 - AlleleB_all, NA_real_))]
eqtl <- eqtl[!is.na(eaf) & eaf > 0 & eaf < 1]

## Z → beta/se（暴露标度，N=NrSamples）
eqtl[, `:=`(
  beta_exp = Zscore / sqrt(2 * eaf * (1 - eaf) * (NrSamples + Zscore^2)),
  se_exp   = 1      / sqrt(2 * eaf * (1 - eaf) * (NrSamples + Zscore^2)))]

## F = Z^2；eQTLGen 显著 cis 天然 |Z|>3.16
eqtl[, F_stat := (beta_exp / se_exp)^2]
eqtl <- eqtl[F_stat > 10]

saveRDS(eqtl, file.path(dir, "eqtl_exposure_prepped.rds"))
cat("prepped:", nrow(eqtl), "行 /", uniqueN(eqtl$Gene), "基因\n")

预期:AssessedAllele==AlleleB ≈ 78.5%(其余翻转,无"都不是");F>10 丢 0;1,891,727 行 / 2715 基因。

Stage 3a · 第一趟"宽松" clump（gene-body±100kb, r²<0.1） → exposure_instruments_final.rds

r
library(data.table); library(readxl); library(ieugwasr)
dir <- "##"; plink <- "##/plink.exe"; bfile <- "##/EUR"

eqtl  <- readRDS(file.path(dir, "eqtl_exposure_prepped.rds"))
finan <- as.data.table(read_excel(
  file.path(dir, "NIHMS80906-supplement-Table_S1.xlsx"), sheet = "Data"))

## gene-body ±100kb（注意：这一趟用基因体坐标，不是 TSS）
gb   <- finan[, .(Gene = ensembl_gene_id, g_start = start_b37, g_end = end_b37)]
inst <- merge(eqtl, gb, by = "Gene", all.x = TRUE)
inst_body <- inst[SNPPos >= g_start - 1e5 & SNPPos <= g_end + 1e5]

genes <- unique(inst_body$Gene)
keep  <- vector("list", length(genes))
pb <- txtProgressBar(min = 0, max = length(genes), style = 3)
for (i in seq_along(genes)) {
  d  <- inst_body[Gene == genes[i], .(rsid = SNP, pval = Pvalue)]
  cl <- tryCatch(
    ld_clump(d, clump_r2 = 0.1, clump_kb = 100, clump_p = 1,
             plink_bin = plink, bfile = bfile),
    error = function(e) d[which.min(pval)])
  keep[[i]] <- data.table(Gene = genes[i], SNP = cl$rsid)
  setTxtProgressBar(pb, i)
}
close(pb)

kp <- rbindlist(keep)
exposure_instruments_final <- inst_body[kp, on = .(Gene, SNP), nomatch = 0L]
saveRDS(exposure_instruments_final, file.path(dir, "exposure_instruments_final.rds"))
cat("宽松版:", nrow(exposure_instruments_final), "SNP /",
    uniqueN(exposure_instruments_final$Gene), "基因（应 ≈ 40287 / 2650）\n")

预期:40,287 SNP / 2650 基因

Stage 3b+3c · 第二趟"严格" clump（TSS±100kb, r²<0.001）+ 按 (Gene,SNP) 配对聚合（白名单 bug 从根去掉） → exposure_instruments_FINAL_r0.001_clean.rds

r
library(data.table); library(readxl)
dir   <- "##"; plink <- "##/plink.exe"; bfile <- "##/EUR"

instruments <- as.data.table(readRDS(file.path(dir, "exposure_instruments_final.rds")))
finan <- as.data.table(read_excel(
  file.path(dir, "NIHMS80906-supplement-Table_S1.xlsx"), sheet = "Data"))
stopifnot(is.numeric(finan$strand),
          uniqueN(finan$ensembl_gene_id) == nrow(finan))

## TSS：正链 start_b37、负链 end_b37（b37/hg19，与 SNPPos 同系）
gc19 <- finan[, .(Gene = ensembl_gene_id,
                  tss  = fifelse(strand >= 0, start_b37, end_b37))]
inst     <- merge(instruments, gc19, by = "Gene", all.x = TRUE)
inst_tss <- inst[abs(SNPPos - tss) <= 1e5]              # ★ 严格趟：TSS±100kb

genes_multi  <- inst_tss[, .N, by = Gene][N > 1]$Gene
genes_single <- inst_tss[, .N, by = Gene][N == 1]$Gene
cat("窗口内基因:", uniqueN(inst_tss$Gene),
    "| 需 clump:", length(genes_multi),
    "| 单 SNP:", length(genes_single), "（应 191）\n")

## 逐基因严格 clump，结果写进 clump_keep/（文件名=基因）
keepdir <- file.path(dir, "clump_keep"); dir.create(keepdir, showWarnings = FALSE)
assoc   <- file.path(dir, "clump_tmp"); dir.create(assoc, showWarnings = FALSE)
af_file <- file.path(assoc, "assoc.txt"); out <- file.path(assoc, "out")

pb <- txtProgressBar(min = 0, max = length(genes_multi), style = 3)
for (i in seq_along(genes_multi)) {
  g  <- genes_multi[i]
  kf <- file.path(keepdir, paste0(g, ".txt"))
  if (file.exists(kf)) { setTxtProgressBar(pb, i); next }   # 断点续跑
  d  <- inst_tss[Gene == g, .(SNP, P = Pvalue)]
  fwrite(d, af_file, sep = "\t")
  system2(plink, c("--bfile", shQuote(bfile), "--clump", shQuote(af_file),
                   "--clump-snp-field", "SNP", "--clump-field", "P",
                   "--clump-r2", "0.001", "--clump-kb", "100",
                   "--clump-p1", "1", "--clump-p2", "1", "--out", shQuote(out)),
          stdout = FALSE, stderr = FALSE)
  cf <- paste0(out, ".clumped")
  if (file.exists(cf)) { writeLines(fread(cf)$SNP, kf); file.remove(cf) }
  else writeLines(d$SNP[which.min(d$P)], kf)               # 极少无输出，留最强
  setTxtProgressBar(pb, i)
}
close(pb)

## ★ 修复版聚合：按 (Gene, SNP) 配对，绝不拉平成全局白名单
keep_multi <- rbindlist(lapply(list.files(keepdir, full.names = TRUE), function(f) {
  s <- readLines(f); if (!length(s)) return(NULL)
  data.table(Gene = sub("\\.txt$", "", basename(f)), SNP = s)
}))
keep_pairs <- rbind(keep_multi, inst_tss[Gene %in% genes_single, .(Gene, SNP)])
clean <- unique(inst_tss[keep_pairs, on = .(Gene, SNP), nomatch = 0L],
                by = c("Gene", "SNP"))

## ---- 硬校验：必须 5447 / 2602，两主角 SNP 精确一致 ----
cat("\n定稿:", nrow(clean), "工具 /", uniqueN(clean$Gene), "基因（预期 5447 / 2602）\n")
for (s in c("PTK2", "ROCK1"))
  cat(s, ":", paste(clean[GeneSymbol == s]$SNP, collapse = ", "), "\n")
stopifnot(
  nrow(clean) == 5447,
  uniqueN(clean$Gene) == 2602,
  clean[GeneSymbol == "PTK2",  .N] == 2,   # rs10088133, rs13250787
  clean[GeneSymbol == "ROCK1", .N] == 1    # rs112006560
)

saveRDS(clean, file.path(dir, "exposure_instruments_FINAL_r0.001_clean.rds"))
cat("\n✅ 已产出干净工具集 exposure_instruments_FINAL_r0.001_clean.rds（5447）\n")

预期:定稿 5447 / 2602;PTK2 = rs10088133, rs13250787;ROCK1 = rs112006560。stopifnot 静默通过才算复现成功。