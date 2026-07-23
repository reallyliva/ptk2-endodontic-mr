## ============================================================
## scRNA 复现 Step 2：QC + 标准化 + 降维聚类（→ 44397 / 21 簇）
## ============================================================
library(Seurat); library(data.table); library(ggplot2)
dir <- "##"

set.seed(42)                                    # 固定随机种子（可复现前提）

## ---- 参数集中定义（改这里，日志自动同步）----
P <- list(
  qc_min_cells=3, qc_min_features=200,
  filt_nFeature_min=200, filt_nFeature_max=6000, filt_nCount_max=50000, filt_percent_mt=15,
  mt_pattern="^MT-", norm_method="LogNormalize", scale_factor=10000,
  n_hvg=2000, vars_to_regress="none",
  pca_npcs=30, cluster_dims=30, cluster_resolution=0.5, seed=42
)

cap <- readRDS(file.path(dir, "cap_merged_raw_repro.rds"))
n0  <- ncol(cap)

## ---- 1) 质控过滤 ----
cap <- subset(cap, subset = nFeature_RNA > P$filt_nFeature_min &
                            nFeature_RNA < P$filt_nFeature_max &
                            nCount_RNA   < P$filt_nCount_max &
                            percent.mt   < P$filt_percent_mt)
n1 <- ncol(cap)
cat("质控：", n0, "→", n1, "细胞（去掉", n0-n1, "）\n"); print(table(cap$sample))

## ---- 2) Seurat v5：JoinLayers → 标准化 → HVG → scale → PCA ----
cap <- JoinLayers(cap)
cap <- NormalizeData(cap, normalization.method=P$norm_method, scale.factor=P$scale_factor, verbose=FALSE)
cap <- FindVariableFeatures(cap, nfeatures=P$n_hvg, verbose=FALSE)
cap <- ScaleData(cap, verbose=FALSE)
cap <- RunPCA(cap, npcs=P$pca_npcs, verbose=FALSE)

## ---- 3) UMAP + 聚类（整合前，看批次）----
cap <- RunUMAP(cap, dims=1:P$cluster_dims, seed.use=P$seed, verbose=FALSE)
cap <- FindNeighbors(cap, dims=1:P$cluster_dims, verbose=FALSE)
cap <- FindClusters(cap, resolution=P$cluster_resolution, random.seed=P$seed, verbose=FALSE)
cat("\n聚类数（res=0.5）:", nlevels(cap$seurat_clusters), "\n"); print(table(cap$seurat_clusters))

## ---- 4) 两张图：批次检查 + 簇（_repro 命名）----
ggsave(file.path(dir,"QC_umap_bysample_repro.png"),
       DimPlot(cap, group.by="sample")+ggtitle("By sample (batch check)"), width=6.5, height=5, dpi=200)
ggsave(file.path(dir,"QC_umap_clusters_repro.png"),
       DimPlot(cap, group.by="seurat_clusters", label=TRUE)+ggtitle("Clusters"), width=6.5, height=5, dpi=200)

## ---- 5) 四靶各簇平均表达（粗看）----
cat("\n四靶各簇平均表达：\n")
print(round(AverageExpression(cap, features=c("PTK2","ROCK1","DLK1","SIK2"),
                              group.by="seurat_clusters")$RNA, 2))

## ---- 参数日志（记录=实际执行）----
log_f <- file.path(dir, "scRNA_params_log_repro.txt")
sink(log_f)
cat("=== CAP scRNA 复现参数日志 ===\n生成时间:", format(Sys.time()), "\n\n")
cat("[样本] 5 CAP: A,B,C(GSE181688) + D,E(GSE197680)\n")
cat("[细胞数] 质控前", n0, "→ 质控后", n1, "\n")
cat("[每样本质控后]\n"); print(table(cap$sample))
cat("[聚类] 数目", nlevels(cap$seurat_clusters), "@ resolution", P$cluster_resolution, "\n\n")
cat("=== 参数 ===\n"); for (k in names(P)) cat(sprintf("  %-20s = %s\n", k, P[[k]]))
cat("\n=== sessionInfo ===\n"); print(sessionInfo())
sink()
cat("\n参数已写入 scRNA_params_log_repro.txt\n")

## ---- 存产出（_repro）----
saveRDS(cap, file.path(dir,"cap_clustered_repro.rds"))
cat("已存 cap_clustered_repro.rds\n")

## ---- 软校验（不 stopifnot，因簇数可能因版本微漂）----
cat("\n【对齐原始记录】质控后应≈44397（原始 44397）；聚类应≈21 簇\n")
cat("  本次：质控后", n1, "细胞 /", nlevels(cap$seurat_clusters), "簇\n")
if (abs(n1-44397)<=50) cat("  ✓ 细胞数对齐\n") else cat("  ⚠️ 细胞数偏离\n")