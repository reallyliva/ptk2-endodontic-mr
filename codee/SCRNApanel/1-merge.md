## ============================================================
## scRNA 复现 Step 1：解压 + 读入 + 合并（不过滤）
## ============================================================
need <- c("Seurat","data.table")
for (p in need) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
library(Seurat); library(data.table)
cat("Seurat 版本:", as.character(packageVersion("Seurat")), "\n\n")

dir <- "##"

## ---- 解压两个 RAW.tar 到各自子目录 ----
for (g in c("GSE181688","GSE197680")) {
  td <- file.path(dir, g); dir.create(td, showWarnings=FALSE)
  untar(file.path(dir, paste0(g,"_RAW.tar")), exdir=td)
  cat("===", g, "解压完成 ===\n"); print(list.files(td))
}

## ---- 样本清单：样本名 → 数据集 + GSM 前缀 ----
meta <- data.table(
  sample = c("A","B","C","D","E"),
  gse    = c("GSE181688","GSE181688","GSE181688","GSE197680","GSE197680"),
  gsm    = c("GSM5509264","GSM5509265","GSM5509266","GSM5928092","GSM5928093"),
  tag    = c("SampleA","SampleB","SampleC","SampleD","SampleE")
)

## ---- 为每个样本建标准三件套目录（去前缀 → Read10X 才认）----
## 用 std_repro 目录，避免和上次的 std 混
std_root <- file.path(dir, "std_repro"); dir.create(std_root, showWarnings=FALSE)
obj_list <- list()
for (i in seq_len(nrow(meta))) {
  s <- meta$sample[i]; g <- meta$gse[i]
  pfx <- paste0(meta$gsm[i], "_", meta$tag[i], "_")
  src <- file.path(dir, g); dst <- file.path(std_root, s); dir.create(dst, showWarnings=FALSE)
  file.copy(file.path(src, paste0(pfx,"barcodes.tsv.gz")), file.path(dst,"barcodes.tsv.gz"), overwrite=TRUE)
  file.copy(file.path(src, paste0(pfx,"features.tsv.gz")), file.path(dst,"features.tsv.gz"), overwrite=TRUE)
  file.copy(file.path(src, paste0(pfx,"matrix.mtx.gz")),   file.path(dst,"matrix.mtx.gz"),   overwrite=TRUE)

  mat <- Read10X(data.dir = dst)
  obj <- CreateSeuratObject(counts = mat, project = s, min.cells = 3, min.features = 200)
  obj$sample  <- s
  obj$dataset <- g
  obj <- RenameCells(obj, add.cell.id = s)             # 细胞名加样本前缀，防撞名
  obj_list[[s]] <- obj
  cat(sprintf("样本 %s (%s): %d 基因 × %d 细胞\n", s, g, nrow(obj), ncol(obj)))
}

## ---- 合并 5 个样本 ----
cap <- merge(obj_list[[1]], y = obj_list[-1], add.cell.ids = NULL, project = "CAP_scRNA")
cat("\n合并后：", nrow(cap), "基因 ×", ncol(cap), "细胞\n")
print(table(cap$sample))

## ---- 线粒体比例（质控要用，先算好，本条不过滤）----
cap[["percent.mt"]] <- PercentageFeatureSet(cap, pattern = "^MT-")
cat("\nnFeature / nCount / percent.mt 概览：\n")
print(summary(cap$nFeature_RNA))
print(summary(cap$nCount_RNA))
print(summary(cap$percent.mt))

## ---- 四靶是否在特征里（定位分析前提）----
cat("\n四靶存在性检查：\n")
for (gsym in c("PTK2","ROCK1","DLK1","SIK2"))
  cat("  ", gsym, "在基因列表中:", gsym %in% rownames(cap), "\n")

## ---- 硬校验：对齐原始记录的建库结果 ----
stopifnot(
  ncol(cap) == 45483,                       # 合并后细胞数
  nrow(cap) == 24949,                       # 合并后基因数
  all(c("PTK2","ROCK1","DLK1","SIK2") %in% rownames(cap))   # 四靶全在
)
cat("\n✅ 建库校验通过：24949 基因 × 45483 细胞，四靶全在\n")

## ---- 存产出（★ _repro 后缀，不覆盖原 cap_merged_raw.rds）----
saveRDS(cap, file.path(dir, "cap_merged_raw_repro.rds"))
cat("已存 cap_merged_raw_repro.rds（未过滤）\n")