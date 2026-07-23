library(data.table); library(TwoSampleMR)
dir <- "##"
tetra <- c("PTK2","ROCK1","DLK1","SIK2")

## ============================================================
## 0. 重建三表型的 harmonise 对象（和 MR 定稿同口径：5447 + action=2）
##    —— Steiger/LOO 都要 harmonise 后的 exposure+outcome 对齐数据
## ============================================================
inst <- as.data.table(readRDS(file.path(dir,"exposure_instruments_FINAL_r0.001_clean.rds")))
stopifnot(nrow(inst)==5447)

br <- fread(file.path(dir,"summary_stats_release_finngen_R12_K11_PULP_PERIAPICAL.gz"),
            select=c("rsids","#chrom","pos")); setnames(br,"#chrom","chrom")
cm <- grepl(",",br$rsids,fixed=TRUE)
bridge <- rbind(br[!cm], br[cm][,.(rsids=unlist(strsplit(rsids,",",fixed=TRUE))),by=.(chrom,pos)])
bridge <- unique(bridge[rsids %in% inst$SNP], by="rsids"); bridge[,key38:=paste(chrom,pos,sep=":")]

exp_df <- inst[, .(SNP, beta.exposure=beta_exp, se.exposure=se_exp,
  effect_allele.exposure=AssessedAllele, other_allele.exposure=OtherAllele,
  eaf.exposure=eaf, pval.exposure=Pvalue, id.exposure=Gene, exposure=Gene,
  samplesize.exposure=NrSamples)]                       # ← Steiger 需要样本量
gene_map <- unique(inst[, .(id.exposure=Gene, GeneSymbol)])

## necrosis 病例/对照数（三表型共用对照；Steiger 结局侧算 R² 用）
N_meta <- list(
  necrosis = c(nca=103832, nco=353106),
  K04      = c(nca=132124, nco=353106),
  pulpitis = c(nca=48120,  nco=353106))
gz_map <- c(necrosis="necrosis_or_apical_periodontitis.gz",
            K04="pulpal_and_apical_diseases.gz", pulpitis="pulpitis.gz")

build_harm <- function(ph){
  o <- fread(file.path(dir, gz_map[[ph]])); setnames(o,"#chrom","chrom",skip_absent=TRUE)
  o[,key38:=paste(chrom,pos,sep=":")]; o <- o[key38 %in% bridge$key38]
  o <- unique(merge(o, bridge[,.(key38,rsids)], by="key38", allow.cartesian=TRUE), by="rsids")
  out_df <- o[, .(SNP=rsids, beta.outcome=beta, se.outcome=sebeta,
    effect_allele.outcome=alt, other_allele.outcome=ref, eaf.outcome=af_alt,
    pval.outcome=pval, id.outcome=ph, outcome=ph,
    ncase.outcome=N_meta[[ph]]["nca"], ncontrol.outcome=N_meta[[ph]]["nco"],
    samplesize.outcome=N_meta[[ph]]["nca"]+N_meta[[ph]]["nco"], prevalence.outcome=0.05)]
  h <- harmonise_data(exp_df, out_df, action=2)
  as.data.table(merge(h, gene_map, by="id.exposure", all.x=TRUE))
}
harm <- lapply(names(gz_map), build_harm); names(harm) <- names(gz_map)

## ============================================================
## 1. STEIGER —— 四靶 × 三表型方向（expression → disease 应为 TRUE）
## ============================================================
steiger_all <- rbindlist(lapply(names(harm), function(ph){
  h <- harm[[ph]][GeneSymbol %in% tetra & mr_keep==TRUE]
  h[, `:=`(units.exposure="SD", units.outcome="log odds")]
  st <- as.data.table(directionality_test(h))          # 每基因一行
  st <- merge(st, unique(h[,.(exposure=id.exposure, GeneSymbol)]), by="exposure", all.x=TRUE)
  st[, phenotype := ph]
  st[, .(phenotype, GeneSymbol,
         snp_r2.exposure=round(snp_r2.exposure,5),
         snp_r2.outcome=signif(snp_r2.outcome,3),
         correct_causal_direction, steiger_pval=signif(steiger_pval,3))]
}))
cat("########## STEIGER（四靶 × 三表型）##########\n")
print(steiger_all[order(match(GeneSymbol,tetra), phenotype)], row.names=FALSE)
cat("\n方向全部为 expression→disease:",
    all(steiger_all$correct_causal_direction==TRUE), "（应 TRUE）\n\n")

## ============================================================
## 2. LOO / 单SNP —— 四靶（用 mr_singlesnp，2 变异靶避开 mr_leaveoneout 折叠坑）
##    necrosis 为主结局
## ============================================================
hnec <- harm[["necrosis"]]
loo <- rbindlist(lapply(tetra, function(g){
  h <- hnec[GeneSymbol==g & mr_keep==TRUE]
  ss <- as.data.table(mr_singlesnp(h))
  ss <- ss[grepl("^rs", SNP)]                           # 只留单 SNP 行（去掉 "All"/Egger 汇总）
  if (nrow(ss)==0) return(data.table(GeneSymbol=g, SNP=NA, b=NA, se=NA, p=NA, OR=NA))
  ss[, .(GeneSymbol=g, SNP, b=round(b,4), se=round(se,4),
         p=signif(p,3), OR=round(exp(b),3))]
}), fill=TRUE)
cat("########## LOO / 单SNP（necrosis，四靶）##########\n")
print(loo, row.names=FALSE)

## ============================================================
## 3. Cochran's Q —— 仅 ≥2 变异靶（PTK2/DLK1/SIK2；ROCK1 单变异→NA）
## ============================================================
Q_all <- rbindlist(lapply(tetra, function(g){
  h <- hnec[GeneSymbol==g & mr_keep==TRUE]
  if (nrow(h) < 2) return(data.table(GeneSymbol=g, nSNP=nrow(h),
                                     Q=NA, Q_df=NA, Q_pval=NA, note="single variant → NA"))
  het <- as.data.table(mr_heterogeneity(h, method_list="mr_ivw"))
  data.table(GeneSymbol=g, nSNP=nrow(h),
             Q=round(het$Q,3), Q_df=het$Q_df, Q_pval=signif(het$Q_pval,3), note="")
}), fill=TRUE)
cat("\n########## Cochran's Q（necrosis）##########\n")
print(Q_all, row.names=FALSE)

## ============================================================
## 4. 硬校验：对齐 v6.1 §2.4（PTK2 两 SNP 各显著、Q p≈0.79；ROCK1 LOO/Q 不可做）
## ============================================================
ptk2_snp <- loo[GeneSymbol=="PTK2"]
stopifnot(
  all(steiger_all$correct_causal_direction==TRUE),      # Steiger 四靶三表型全对
  nrow(ptk2_snp)==2,                                    # PTK2 两个单 SNP 都在
  all(ptk2_snp$p < 1e-3),                               # 各自 p<1e-3（1.1e-4 / 7.5e-4）
  abs(Q_all[GeneSymbol=="PTK2"]$Q_pval - 0.79) < 0.05,  # PTK2 Q p≈0.79
  is.na(Q_all[GeneSymbol=="ROCK1"]$Q_pval)              # ROCK1 单变异 Q=NA
)
cat("\n✅ 校验通过：Steiger 方向全对；PTK2 两变异各显著、Q p≈0.79；ROCK1 单变异 LOO/Q 如实 NA\n")

## ============================================================
## 5. 固化
## ============================================================
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(dir,"FROZEN"); if(!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)
fwrite(steiger_all, file.path(out_dir, paste0("FROZEN_steiger_5447_",stamp,".csv")))
fwrite(loo,         file.path(out_dir, paste0("FROZEN_loo_singlesnp_5447_",stamp,".csv")))
fwrite(Q_all,       file.path(out_dir, paste0("FROZEN_cochranQ_5447_",stamp,".csv")))
cat("\n已固化 → FROZEN_steiger_5447 / FROZEN_loo_singlesnp_5447 / FROZEN_cochranQ_5447 （", stamp, "）\n")