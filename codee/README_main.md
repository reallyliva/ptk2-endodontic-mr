# PTK2 endodontic drug-target Mendelian randomization + single-cell localization

Reproducibility repository for the study identifying **PTK2 (focal adhesion kinase,
FAK)** as a causal, druggable target for pulpal and apical (endodontic) disease,
localizing it to **vascular pericytes**, with **ROCK1** as a supportive same-pathway
target.

- **Manuscript:** *Druggable-genome Mendelian randomization implicates PTK2 (focal adhesion kinase) in pulpal and apical disease and localizes the target to vascular pericytes*
- **Archived release:** *https://doi.org/10.5281/zenodo.21535754*
- **License:** MIT (see LICENSE)

## What this study does

Cis-eQTL instruments for druggable genes (eQTLGen blood) are tested by Mendelian
randomization against three FinnGen release 12 endodontic phenotypes (necrosis of
pulp or apical periodontitis [primary], the composite category of pulpal and apical
diseases, and pulpitis), then prioritized by colocalization, Steiger directionality,
leave-one-out / Cochran's Q, and phenome-wide pleiotropy screening. The leading
targets are then localized in an integrated single-cell atlas of five chronic apical
periodontitis samples.

**Headline result.** PTK2 is the leading target — significant in all three phenotypes
(necrosis OR 1.12, FDR 7.8e-4), colocalizing with the disease signal (PP.H4 0.796),
correctly oriented from expression to disease, and specifically expressed in
pericytes. ROCK1 is retained as a supportive target (consistent direction, cleanest
colocalization PP.H4 0.789, but a borderline FDR on a single instrument). Both lie on
the FAK–RhoA–ROCK pathway and are druggable (approved FAK and ROCK inhibitors exist).

## Repository structure

The repository has two self-contained modules, each with its own README:

```
.
├── codee/MRpanel/       Mendelian randomization arm  (see codee/MRpanel/README)
│   ├── 1-expose.md … 5-Druggability.md        analysis scripts (R, in .md files)
│   ├── exposure_instruments_*.rds             instrument sets (5447 primary; 6652 sensitivity)
│   ├── MR_results_*_5447.rds                  MR results per phenotype
│   └── FROZEN/                                frozen result tables (source of the reported numbers)
│
└── codee/SCRNApanel/    Single-cell localization arm (see codee/SCRNApanel/README)
    ├── 1-merge.md … 5-fin.md                  analysis scripts (R, in .md files)
    ├── cap_annotated_repro.rds                final annotated Seurat object
    └── TABLE_scRNA_target_expression_bycelltype_repro.csv
```

- **MR arm** (`codee/MRpanel/`): builds the cis-eQTL instrument set, runs the
  three-phenotype MR, and performs colocalization, Steiger, leave-one-out, Cochran's Q
  and druggability annotation. See its README for the full pipeline and checkpoints.
- **Single-cell arm** (`codee/SCRNApanel/`): merges, integrates (Harmony), clusters
  and annotates the five-sample atlas, and locates the targets across cell types. See
  its README for parameters and checkpoints.

Scripts are provided as `.md` files containing R code; copy the code into R and run in
the numbered order within each module.

## Requirements

R ≥ 4.5. MR arm: `TwoSampleMR`, `coloc`, `data.table`, `readxl`, `ieugwasr`, plus
**PLINK v1.9** and the 1000 Genomes Phase 3 EUR panel. Single-cell arm: `Seurat` (v5),
`harmony`, `data.table`, `ggplot2`.

## Source data (not included — obtain from the original sources)

Raw third-party files are **not redistributed**; download them and place them in the
working directory before running the scripts:

- **eQTLGen** cis-eQTL and allele-frequency files — https://www.eqtlgen.org/
- **FinnGen release 12** three-phenotype summary statistics and the
  `K11_PULP_PERIAPICAL` bridge file (used to restore rsIDs) — FinnGen public release
- **Finan et al.** druggable-genome catalogue — *Sci Transl Med* 2017;9:eaag1166
- **1000 Genomes** Phase 3 EUR panel — https://www.internationalgenome.org/
- **Single-cell** datasets — NCBI GEO **GSE181688** and **GSE197680**

## Key reproduction checkpoints

- Instrument set: **5447** instruments / 2602 genes (two-stage clumping; a single-stage
  clump gives 6652, provided as a sensitivity analysis — PTK2 unchanged).
- Necrosis MR: 2545 genes tested, **19** at FDR < 0.05; PTK2 b = 0.114, p = 3.0e-7;
  ROCK1 FDR = 0.052 (borderline, retained as supportive).
- Colocalization: PTK2 650/0.796, ROCK1 218/0.789, DLK1 764/0.661, SIK2 415/0.632.
- Single-cell: 24,949 genes × 45,483 cells → 44,397 after QC → ~17 clusters after
  Harmony; PTK2 pericyte-specific (1.62 >> endothelial 0.37 >> fibroblast 0.20).

Reproducibility notes specific to each arm (two-stage clumping provenance; pericyte
annotation; pseudobulk rule) are documented in the two module READMEs.

## License and citation

Code is released under the MIT License (see LICENSE). Please cite the manuscript and
this archived release (Zenodo DOI). Third-party datasets (eQTLGen, FinnGen, Finan et
al., 1000 Genomes, GEO GSE181688/GSE197680) remain under their original terms; cite
the original sources listed above.
