suppressMessages({library(data.table); library(coloc); library(readxl)})
dir <- "##"; out_dir <- file.path(dir, "FROZEN")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

## —— 宽原料 + 本次可复现的 MR 结果——
eqtl_sig <- as.data.table(readRDS(file.path(dir, "coloc_eqtl_21genes.rds")))
nec_sub  <- as.data.table(readRDS(file.path(dir, "coloc_nec_21genes.rds")))
mr_nec   <- as.data.table(readRDS(file.path(dir, "MR_results_necrosis_5447.rds")))

## ★ coloc 基因集 = necrosis 的 19 个 FDR 显著 + ROCK1（方向一致、有 FDA 药，显式纳入）
sig19  <- mr_nec[fdr < 0.05, .(id.exposure, GeneSymbol)]
rock1  <- mr_nec[GeneSymbol == "ROCK1", .(id.exposure, GeneSymbol)]
sig_genes <- unique(rbind(sig19, rock1))
setnames(sig_genes, "GeneSymbol", "exposure")
cat("coloc 基因集:", nrow(sig_genes), "（19 FDR + ROCK1 = 20；ROCK1 在否:",
    "ROCK1" %in% sig_genes$exposure, "）\n")

## —— eaf/beta/se（col.md 第 318-326 行，原样）——
af <- fread(file.path(dir, "2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz"),
            select = c("SNP","AlleleA","AlleleB","AlleleB_all"))
eqtl_sig <- merge(eqtl_sig, af, by = "SNP", all.x = TRUE)
eqtl_sig[, eaf := fifelse(AssessedAllele==AlleleB, AlleleB_all,
                   fifelse(AssessedAllele==AlleleA, 1-AlleleB_all, NA_real_))]
eqtl_sig <- eqtl_sig[!is.na(eaf) & eaf>0 & eaf<1]
eqtl_sig[, maf := pmin(eaf, 1-eaf)]
eqtl_sig[, se_exp   := 1/sqrt(2*maf*(1-maf)*(NrSamples+Zscore^2))]
eqtl_sig[, beta_exp := Zscore*se_exp]

## —— Finan hg19 TSS（col.md 第 329-331 行）——
finan <- as.data.table(read_excel(file.path(dir,"NIHMS80906-supplement-Table_S1.xlsx"), sheet="Data"))
gc19 <- finan[ensembl_gene_id %in% sig_genes$id.exposure,
              .(Gene=ensembl_gene_id, tss=fifelse(strand>=0, start_b37, end_b37))]

N_case <- 103832; N_ctrl <- 353106; S_CC <- N_case/(N_case+N_ctrl)

## —— 单基因 coloc ——
coloc_gene <- function(g_ens, p12){
  tss <- gc19[Gene==g_ens]$tss; if (length(tss)==0) return(NULL)
  e <- eqtl_sig[Gene==g_ens & abs(SNPPos - tss) <= 1e5]        # ★ 第341行：TSS±100kb
  o <- nec_sub[rsids %in% e$SNP]                               # rsID 捞，不碰坐标
  common <- intersect(e$SNP, o$rsids)
  if (length(common) < 20) return(list(nsnp=length(common), s=NULL))
  e <- e[SNP %in% common][!duplicated(SNP)]; o <- o[rsids %in% common][!duplicated(rsids)]
  setkey(e,SNP); setkey(o,rsids); o <- o[e$SNP]
  D_exp <- list(beta=e$beta_exp, varbeta=e$se_exp^2, snp=e$SNP,  type="quant", N=max(e$NrSamples), MAF=e$maf)
  D_out <- list(beta=o$beta,     varbeta=o$sebeta^2, snp=o$rsids, type="cc", s=S_CC, N=N_case+N_ctrl, MAF=o$af_alt)
  r <- tryCatch(suppressWarnings(coloc.abf(D_exp, D_out, p12=p12)), error=function(x) NULL)
  list(nsnp=length(common), s=if(is.null(r)) NULL else r$summary)
}

## —— 默认先验 p12=1e-5：出 coloc 表 ——
cat("\n########## coloc（19 + ROCK1，默认 p12=1e-5）##########\n")
def <- rbindlist(lapply(1:nrow(sig_genes), function(i){
  cg <- coloc_gene(sig_genes$id.exposure[i], 1e-5)
  if (is.null(cg$s)) data.table(gene=sig_genes$exposure[i], nsnp=cg$nsnp, PP.H3=NA, PP.H4=NA)
  else data.table(gene=sig_genes$exposure[i], nsnp=cg$nsnp,
                  PP.H3=round(cg$s["PP.H3.abf"],3), PP.H4=round(cg$s["PP.H4.abf"],3))
}), fill=TRUE)[order(-PP.H4)]
print(def, row.names=FALSE)
fwrite(def, file.path(out_dir, paste0("FROZEN_coloc_19plusROCK1_",stamp,".csv")))

## —— 与旧表 S2 核对四靶（值应不变）——
old <- as.data.table(readRDS(file.path(dir,"coloc_results_21genes.rds")))
chk <- merge(old[,.(gene, H4_old=PP.H4, nsnp_old=nsnp)],
             def[,.(gene, PP.H4, nsnp)], by="gene")[gene %in% c("PTK2","ROCK1","DLK1","SIK2")]
chk[, `:=`(dH4=round(PP.H4-H4_old,3), dN=nsnp-nsnp_old)]
cat("\n--- 四靶 vs 旧表 S2（dH4、dN 应为 0）---\n"); print(chk, row.names=FALSE)

## —— 先验敏感性（同一窄输入，p12 = 1e-5 / 5e-6 / 1e-6）——
cat("\n########## 先验敏感性 ##########\n")
grid <- rbindlist(lapply(c(1e-5,5e-6,1e-6), function(pp){
  rbindlist(lapply(1:nrow(sig_genes), function(i){
    cg <- coloc_gene(sig_genes$id.exposure[i], pp)
    data.table(gene=sig_genes$exposure[i], p12=pp,
               PP.H4=if(is.null(cg$s)) NA else round(cg$s["PP.H4.abf"],4))
  }))
}))
w <- dcast(grid, gene~p12, value.var="PP.H4"); setorder(w, -`1e-05`)
print(w, row.names=FALSE)
fwrite(grid, file.path(out_dir, paste0("FROZEN_coloc_p12grid_",stamp,".csv")))

## —— 硬校验：四靶必须精确复现旧表 S2 ——
d <- setNames(def$PP.H4, def$gene)
stopifnot(
  abs(d["PTK2"]  - 0.796) < 0.01,
  abs(d["ROCK1"] - 0.789) < 0.01,
  abs(d["DLK1"]  - 0.661) < 0.02,
  abs(d["SIK2"]  - 0.632) < 0.02
)
cat("\n✅ 四靶 coloc 精确复现（PTK2 0.796 / ROCK1 0.789 / DLK1 0.661 / SIK2 0.632），且与 19 名单自洽\n")
cat("Timestamp:", stamp, "| 基因集 = necrosis 19 FDR + ROCK1 | s =", round(S_CC,4), "\n")