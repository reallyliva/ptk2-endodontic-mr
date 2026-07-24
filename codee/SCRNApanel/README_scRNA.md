# Single-cell localization (scRNA-seq) — PTK2 in periapical disease

Single-cell arm of the study: locating the Mendelian-randomization targets
(PTK2, ROCK1, DLK1, SIK2) in chronic apical periodontitis (CAP) tissue.
**Main finding: PTK2 is specifically enriched in vascular pericytes; ROCK1 is
highest in pericytes but broadly expressed.**

Part of the `ptk2-endodontic-mr` repository (see the top-level README for the
Mendelian-randomization arm).

## Overview

Five CAP periapical samples (no healthy control) are merged, quality-controlled,
integrated across samples with Harmony, clustered, and annotated by canonical
markers. Because there is no matched healthy control, this arm is **descriptive
cell-type localization only** — where the targets are expressed — not a
case-versus-control comparison.

Key result: PTK2 is pericyte-specific (Pericyte 1.62 >> Endothelial 0.37 >>
Fibroblast 0.20; only pericytes express it as a whole population). ROCK1 is
highest in pericytes (1.15) but broadly expressed across several types. DLK1 and
SIK2 are at or below detection, so single-cell localization is not applicable for
them (their causal support comes from the MR arm).

## Repository layout

```
codee/SCRNApanel/
├── 1-merge.md         read 10x triplets, build 5 Seurat objects, merge (-> cap_merged_raw_repro.rds)
├── 2-clustered.md     QC filter + normalize + PCA + UMAP + cluster (-> cap_clustered_repro.rds)
├── 3-harmony.md       Harmony integration by sample + re-cluster (-> cap_harmony_repro.rds)
├── 4-annotated.md     marker annotation incl. pericyte markers (-> cap_annotated_repro.rds)
├── 5-fin.md           final localization figures + target-expression table
├── cap_annotated_repro.rds                         final annotated object
├── cap_harmony_repro.rds                           integrated object (optional intermediate)
└── TABLE_scRNA_target_expression_bycelltype_repro.csv   target expression per cell type
```

Scripts are `.md` files containing R code; copy the code into R and run in order
(1 -> 5). All outputs carry a `_repro` suffix.

## Requirements

R 4.5 with `Seurat` (v5), `harmony`, `data.table`, `ggplot2`. Set the working
directory (`dir <- "X:/##"`) to where the raw data are placed.

## Source data (not included here — obtain from GEO)

Raw files are **not redistributed**. Download and place in the working directory:

- **GSE181688** (samples A/B/C) and **GSE197680** (samples D/E), `_RAW.tar`,
  10x MTX+TSV triplets — NCBI GEO. Together they form five CAP periapical samples
  from the same group.

## Key parameters (reproducibility)

- Random seed 42 (`set.seed(42)`, UMAP `seed.use=42`, `FindClusters random.seed=42`).
- Object construction: `min.cells=3`, `min.features=200`; mitochondrial fraction
  from `^MT-`.
- QC filter: nFeature 200-6000, nCount < 50000, percent.mt < 15 (relaxed for
  inflamed tissue).
- Normalization LogNormalize (scale factor 1e4); 2000 variable features; no
  variables regressed; 30 PCs.
- Clustering on dims 1:30, resolution 0.5.
- Integration: Harmony by `sample` (dims 1:30); UMAP/clustering redone on the
  Harmony embedding (`umap_h`, dims 1:30, resolution 0.5).
- All parameters are also written to `scRNA_params_log_repro.txt` at run time
  (record = actual execution).

## Key reproduction checkpoints

- Merge: 24,949 genes x 45,483 cells (A 8107, B 10982, C 7265, D 9585, E 9544);
  all four targets present.
- QC: 45,483 -> 44,397 cells; ~21 clusters (pre-integration).
- Batch: sample A separates on its own before integration (Harmony required).
- After Harmony: ~17 clusters; PTK2 and ROCK1 both highest in the same pericyte
  cluster.
- Annotation: the PTK2-high cluster shows RGS5 ~31, ACTA2 ~37, TAGLN ~35, PDGFRB
  ~5, NOTCH3 ~5, with MYH11 ~2.6 and CNN1 ~0.5 low -> **wall pericyte, not mature
  vascular smooth muscle**. PTK2 in pericytes 1.62 >> endothelial 0.37 >>
  fibroblast 0.20.

## Notes (important)

- **Pericyte identity requires pericyte markers.** With a fibroblast-only marker
  panel the PTK2-high cluster looks like a fibroblast (it expresses some collagen).
  Its true identity is revealed only by adding RGS5 / PDGFRB / NOTCH3 / ACTA2 /
  TAGLN (high) and MYH11 / CNN1 (low). Script 4 includes these; do not drop them.
- **Cluster numbers may vary by version.** Exact cluster count (e.g. 21 -> 17) and
  cluster numbering depend on the Seurat/uwot version; identify the pericyte
  population by its markers (RGS5/ACTA2 high, PTK2 highest), not by a fixed
  cluster index. The PTK2-pericyte conclusion is robust to this.
- **Aggregation level (pseudobulk).** Localization is shown with FeaturePlot /
  DotPlot / VlnPlot and average expression per cell type; **no between-group
  differential-expression p-values are computed at the cell level.** Any
  between-group differential expression must be done on pseudobulk aggregated by
  donor/sample, never with cell-level tests. Marker detection for annotation
  (cell-level `FindAllMarkers`) is used only to label clusters.
- **Interpretation scope.** Single-cell results are localization/association, not
  mechanism. Wording: "PTK2 specifically enriched in pericytes; ROCK1 highest in
  pericytes but broadly expressed; both converge on the pericyte population." Avoid
  "co-localize" as a blanket term (accurate for PTK2, overstated for ROCK1).

## Outputs

`cap_annotated_repro.rds`, `TABLE_scRNA_target_expression_bycelltype_repro.csv`,
and figures: `FIG_scRNA_umap_celltypes_repro`, `FIG_scRNA_featureplot_PTK2_ROCK1_repro`,
`FIG_scRNA_violin_PTK2_ROCK1_repro`, `FIG_scRNA_dotplot_pericyte_repro`, plus
`scRNA_params_log_repro.txt`.

## License

MIT (see LICENSE). GEO datasets remain under their original terms; cite GSE181688
and GSE197680.
