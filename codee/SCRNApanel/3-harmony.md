## ============================================================
## scRNA 复现 Step 3：Harmony 整合 + 重聚类
## 修正：AverageExpression 返回 S4 稀疏矩阵，取 which.max 前先 as.matrix()
## ============================================================
library(Seurat); library(data.table); library(ggplot2)
if (!requireNamespace("harmony", quietly=TRUE)) install.packages("harmony")
library(harmony)
dir <- "##"
set.seed(42)

cap <- readRDS(file.path(dir,"cap_clustered_repro.rds"))

## ---- Harmony 整合（按 sample 校正批次）----
cap <- RunHarmony(cap, group.by.vars="sample",
                  reduction.use="pca", reduction.save="harmony",
                  dims.use=1:30, verbose=FALSE)

## ---- harmony embedding 重做 UMAP + 聚类 ----
cap <- RunUMAP(cap, reduction="harmony", dims=1:30,
               seed.use=42, reduction.name="umap_h", verbose=FALSE)
cap <- FindNeighbors(cap, reduction="harmony", dims=1:30, verbose=FALSE)
cap <- FindClusters(cap, resolution=0.5, random.seed=42, verbose=FALSE)
cat("整合后聚类数:", nlevels(cap$seurat_clusters), "\n"); print(table(cap$seurat_clusters))

## ---- 出图（_repro）----
ggsave(file.path(dir,"QC_umap_bysample_harmony_repro.png"),
       DimPlot(cap, group.by="sample", reduction="umap_h")+ggtitle("After Harmony (batch)"),
       width=6.5, height=5, dpi=200)
ggsave(file.path(dir,"QC_umap_clusters_harmony_repro.png"),
       DimPlot(cap, group.by="seurat_clusters", reduction="umap_h", label=TRUE)+ggtitle("Clusters (Harmony)"),
       width=6.5, height=5, dpi=200)

## ---- 四靶各簇表达（★ 先 as.matrix 转普通矩阵，避免 S4 问题）----
avg4 <- as.matrix(AverageExpression(cap, features=c("PTK2","ROCK1","DLK1","SIK2"),
                                    group.by="seurat_clusters")$RNA)
cat("\n四靶各簇平均表达（整合后簇号）:\n")
print(round(avg4, 2))

## ---- 追加 Harmony 参数到日志 ----
sink(file.path(dir,"scRNA_params_log_repro.txt"), append=TRUE)
cat("\n\n=== Harmony 整合（追加）===\n生成时间:", format(Sys.time()), "\n")
cat("  integration=Harmony  group.by.vars=sample  dims.use=1:30\n")
cat("  downstream=harmony (umap_h, dims 1:30, res 0.5, seed 42)\n")
cat("  clusters_after =", nlevels(cap$seurat_clusters), "\n")
cat("  harmony_version =", as.character(packageVersion("harmony")), "\n")
sink()

## ---- 存产出（_repro）----
saveRDS(cap, file.path(dir,"cap_harmony_repro.rds"))
cat("\n已存 cap_harmony_repro.rds；参数已追加日志\n")

## ---- 软提示（★ 修正：as.matrix 后再 which.max）----
nc <- nlevels(cap$seurat_clusters)
cat("\n【对齐原始记录】整合后应≈17 簇\n")
cat("  本次整合后簇数:", nc, if(abs(nc-17)<=2) "（✓ 接近 17）" else "（⚠️ 偏离）", "\n")

pk <- as.matrix(AverageExpression(cap, features="PTK2", group.by="seurat_clusters")$RNA)
best <- colnames(pk)[which.max(pk[1,])]           # pk 是 1×簇 矩阵，取第一行
cat("  PTK2 最高簇 =", best, "，值 =", round(max(pk[1,]),2), "（原始：某簇 1.62）\n")

## 顺带把 ROCK1 最高簇也打印，确认和 PTK2 是否同一簇（原始两者都在 g4 最高）
rk <- as.matrix(AverageExpression(cap, features="ROCK1", group.by="seurat_clusters")$RNA)
cat("  ROCK1 最高簇 =", colnames(rk)[which.max(rk[1,])],
    "，值 =", round(max(rk[1,]),2), "（原始：也在同一簇最高，1.15）\n")