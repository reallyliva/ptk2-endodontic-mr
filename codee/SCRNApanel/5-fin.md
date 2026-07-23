## ============================================================
## scRNA 复现 Step 5：正式定位图 + 目标表达表（收尾）
## 产出全部 _repro；
## ============================================================
library(Seurat); library(ggplot2); library(data.table)
dir <- "##"; set.seed(42)
cap <- readRDS(file.path(dir,"cap_annotated_repro.rds"))
cat("细胞类型:", paste(levels(cap$celltype), collapse=", "), "\n")

## ---- 1) 带细胞类型标签的 UMAP ----
p_umap <- DimPlot(cap, reduction="umap_h", group.by="celltype", label=TRUE,
                  repel=TRUE, label.size=3.5) + ggtitle("CAP periapical tissue — cell types")
ggsave(file.path(dir,"FIG_scRNA_umap_celltypes_repro.png"), p_umap, width=8, height=6, dpi=300)
ggsave(file.path(dir,"FIG_scRNA_umap_celltypes_repro.pdf"), p_umap, width=8, height=6)

## ---- 2) PTK2/ROCK1 FeaturePlot ----
p_feat <- FeaturePlot(cap, features=c("PTK2","ROCK1"), reduction="umap_h",
                      order=TRUE, ncol=2, cols=c("lightgrey","#D55E00"))
ggsave(file.path(dir,"FIG_scRNA_featureplot_PTK2_ROCK1_repro.png"), p_feat, width=11, height=5, dpi=300)
ggsave(file.path(dir,"FIG_scRNA_featureplot_PTK2_ROCK1_repro.pdf"), p_feat, width=11, height=5)

## ---- 3) VlnPlot：PTK2/ROCK1 按细胞类型（最能说明特异性）----
p_vln <- VlnPlot(cap, features=c("PTK2","ROCK1"), group.by="celltype",
                 pt.size=0, ncol=1) & theme(axis.text.x=element_text(size=8, angle=45, hjust=1))
ggsave(file.path(dir,"FIG_scRNA_violin_PTK2_ROCK1_repro.png"), p_vln, width=8, height=7, dpi=300)

## ---- 4) DotPlot：两靶 + 周细胞标志（坐实周细胞节点）----
dot_feats <- c("PTK2","ROCK1","RGS5","PDGFRB","NOTCH3","ACTA2","TAGLN")
p_dot <- DotPlot(cap, features=dot_feats, group.by="celltype") + RotatedAxis() +
  ggtitle("Expression of PTK2/ROCK1 and pericyte markers across cell types")
ggsave(file.path(dir,"FIG_scRNA_dotplot_pericyte_repro.png"), p_dot, width=8, height=6, dpi=300)

## ---- 5) 目标表达表（进补充材料）----
avg_ct <- as.matrix(AverageExpression(cap, features=c("PTK2","ROCK1","DLK1","SIK2"),
                                      group.by="celltype")$RNA)
out_tab <- data.table(gene=rownames(avg_ct), as.data.table(round(avg_ct,3)))
fwrite(out_tab, file.path(dir,"TABLE_scRNA_target_expression_bycelltype_repro.csv"))
cat("\n四靶各细胞类型平均表达（已存表）:\n"); print(round(avg_ct,3))

## ---- 6) 收尾日志 ----
sink(file.path(dir,"scRNA_params_log_repro.txt"), append=TRUE)
cat("\n\n=== Step5 出图收尾 ===\n生成时间:", format(Sys.time()), "\n")
cat("产出: FIG_scRNA_{umap_celltypes, featureplot_PTK2_ROCK1, violin_PTK2_ROCK1, dotplot_pericyte}_repro\n")
cat("      TABLE_scRNA_target_expression_bycelltype_repro.csv\n")
cat("结论口径: PTK2 特异富集周细胞；ROCK1 周细胞最高但广谱；两者共汇周细胞（FAK-ROCK 节点）\n")
cat("          DLK1/SIK2 表达低于检测，定位不适用（由 MR 支撑）\n")
sink()
cat("\n✅ scRNA 复现全流程完成：cap_annotated_repro.rds + 5 图 + 1 表 + 日志\n")