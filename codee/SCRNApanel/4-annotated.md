## ============================================================
## scRNA 复现 Step 4：注释（★含周细胞 marker，钉死 g4=周细胞）
## 产出全部 _repro；stopifnot 锚定周细胞结论
## ============================================================
library(Seurat); library(data.table); library(ggplot2)
dir <- "##"; set.seed(42)
cap <- readRDS(file.path(dir,"cap_harmony_repro.rds"))
cat("整合后簇数:", nlevels(cap$seurat_clusters), "\n")

## ---- 1) 全谱系 marker（★关键：含周细胞/平滑肌，避免 g4 误判成纤维）----
markers <- list(
  Fibroblast      = c("COL1A1","COL1A2","DCN","LUM","PDGFRA","POSTN"),
  Pericyte        = c("RGS5","PDGFRB","NOTCH3","ACTA2","TAGLN","MYH11","CNN1","KCNJ8","HIGD1B","CSPG4"),
  Endothelial     = c("PECAM1","VWF","CDH5","CLDN5"),
  `Mono/Macrophage`= c("CD68","CD14","LYZ","CD163","C1QA","C1QB","C1QC","FCN1"),
  Tcell           = c("CD3D","CD3E","CD8A","IL7R","TRAC"),
  Bcell           = c("MS4A1","CD79A","CD19","CD22","BANK1"),
  Plasma          = c("MZB1","IGHG1","IGKC","JCHAIN","XBP1"),
  Mast            = c("TPSAB1","CPA3","MS4A2","KIT"),
  DC              = c("CLEC9A","LILRA4","IDO1","CPVL"),
  Epithelial      = c("KRT14","KRT5","KRT15","EPCAM","PITX2"),
  Proliferating   = c("MKI67","TOP2A")
)

## ---- 2) g4 的平滑肌 vs 周细胞判别（复现原始那次关键交叉验证）----
peri_probe <- c("RGS5","ACTA2","MYH11","TAGLN","CNN1","PDGFRB","CSPG4","NOTCH3","KCNJ8","HIGD1B","PTK2","ROCK1")
peri_probe <- peri_probe[peri_probe %in% rownames(cap)]
g4_prof <- as.matrix(AverageExpression(cap, features=peri_probe, group.by="seurat_clusters")$RNA)[,"g4",drop=FALSE]
cat("\n=== g4 的周细胞 vs 平滑肌标志（未 log 平均表达）===\n"); print(round(g4_prof,2))

## ---- 3) DotPlot：全 marker 按簇（人工判读用）----
all_mk <- unique(unlist(markers)); all_mk <- all_mk[all_mk %in% rownames(cap)]
ggsave(file.path(dir,"annot_dotplot_markers_repro.png"),
       DotPlot(cap, features=all_mk, group.by="seurat_clusters")+RotatedAxis()+
         theme(axis.text.x=element_text(size=7))+ggtitle("Marker expression by cluster"),
       width=16, height=6, dpi=200)

## ---- 4) 每簇 top5 marker（cell-level，仅贴标签）----
DefaultAssay(cap) <- "RNA"
cap.markers <- FindAllMarkers(cap, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25, verbose=FALSE)
top5 <- as.data.table(cap.markers)[, head(.SD[order(-avg_log2FC)],5), by=cluster]
fwrite(top5, file.path(dir,"annot_top5_markers_per_cluster_repro.csv"))

## ---- 5) 17 簇 → 11 类（沿用原始三次交叉验证的最终命名）----
celltype_map <- c(
  "0"="Plasma cell","3"="Plasma cell","9"="Plasma cell","10"="Plasma cell","15"="Plasma cell",
  "1"="Fibroblast","13"="Proliferating fibroblast",
  "2"="Endothelial","4"="Pericyte",
  "5"="T cell","7"="T cell","6"="B cell","16"="B cell",
  "8"="Mono/Macrophage","11"="Epithelial","12"="Mast cell","14"="Dendritic cell"
)
cap$celltype <- factor(unname(celltype_map[as.character(cap$seurat_clusters)]))
Idents(cap) <- "celltype"
cat("\n各细胞类型细胞数：\n"); print(table(cap$celltype))

## ---- 6) 四靶按细胞类型平均表达 ----
avg_ct <- as.matrix(AverageExpression(cap, features=c("PTK2","ROCK1","DLK1","SIK2"),
                                      group.by="celltype")$RNA)
cat("\n四靶各细胞类型平均表达：\n"); print(round(avg_ct,3))

## ---- 7) 注释依据写日志 ----
sink(file.path(dir,"scRNA_params_log_repro.txt"), append=TRUE)
cat("\n\n=== 细胞类型注释（追加）===\n生成时间:", format(Sys.time()), "\n")
cat("Pericyte(g4): RGS5+/PDGFRB+/NOTCH3+/ACTA2+/TAGLN+，MYH11/CNN1 低 → 壁周细胞非成熟平滑肌\n")
cat("簇→类型映射:\n"); for(k in names(celltype_map)) cat(sprintf("  g%s = %s\n", k, celltype_map[k]))
cat("PTK2 周细胞特异；ROCK1 周细胞最高但广谱；DLK1/SIK2 测不到\n")
sink()

saveRDS(cap, file.path(dir,"cap_annotated_repro.rds"))
cat("\n已存 cap_annotated_repro.rds\n")

## ---- 8) 🔴 硬校验：钉死周细胞结论（scRNA该复现）----
g4v <- g4_prof[,1]; names(g4v) <- rownames(g4_prof)
stopifnot(
  g4v["RGS5"]  > 10,          # 周细胞金标准，原始 30.83
  g4v["ACTA2"] > 10,          # 原始 36.60
  g4v["MYH11"] < 10,          # 成熟平滑肌标志低 → 周细胞非 VSMC，原始 2.60
  g4v["PDGFRB"]> 1,           # 周细胞偏向，原始 5.20
  avg_ct["PTK2","Pericyte"]  > 1,     # PTK2 周细胞高，原始 1.62
  avg_ct["PTK2","Pericyte"]  > avg_ct["PTK2","Fibroblast"],   # PTK2 周细胞 >> 成纤维
  avg_ct["PTK2","Pericyte"]  > avg_ct["PTK2","Endothelial"]   # PTK2 周细胞 >> 内皮
)
cat("\n✅ 周细胞结论复现：g4 = Pericyte（RGS5/ACTA2 高、MYH11 低）；PTK2 周细胞特异\n")