# Drug-target Mendelian randomization of the druggable genome in pulpal and apical disease

Reproducibility code for the study identifying **PTK2 (FAK)** as a causal, druggable
target for pulpal and apical (endodontic) disease, with single-cell localization to
vascular pericytes and **ROCK1** as a supportive same-pathway target.

## Overview

Cis-eQTL instruments for druggable genes (eQTLGen blood) were tested by Mendelian
randomization against three FinnGen R12 endodontic phenotypes (necrosis of pulp or
apical periodontitis [primary], pulpal and apical diseases [K04], pulpitis), then
prioritized by colocalization, Steiger directionality, leave-one-out / Cochran's Q,
and single-cell localization. **PTK2** is the leading target (necrosis OR 1.12,
FDR 7.8e-4, coloc PP.H4 0.796, pericyte-specific); **ROCK1** is supportive (borderline
FDR, coloc PP.H4 0.789).

## Repository layout

```
codee/MRpanel/
├── 1-expose.md            build the cis-eQTL instrument set (-> ..._clean.rds, 5447 instruments)
├── 2-MR-main.md           harmonize + IVW/Wald + BH-FDR across the three phenotypes
├── 3-coloc.md             colocalization (coloc.abf) for the 19 FDR genes + ROCK1
├── 4-steigr+LOO+CQ.md     Steiger directionality, leave-one-out, Cochran's Q
├── 5-Druggability.md      druggability annotation (Finan / DGIdb / Open Targets)
├── exposure_instruments_final.rds                permissive-clump intermediate (40,287 variants; see Notes)
├── exposure_instruments_FINAL_r0.001_clean.rds   5447 instruments / 2602 genes (primary)
├── exposure_instruments_6652_singlepass.rds      single-stage clumping (sensitivity)
├── MR_results_{necrosis,K04,pulpitis}_5447.rds   MR results per phenotype
└── FROZEN/                frozen result tables (the single source of the reported numbers)
```

Scripts are provided as `.md` files containing R code; copy the code into R to run.
Run them in the numbered order (1 -> 5).

## Requirements

R 4.5 with `TwoSampleMR`, `coloc`, `data.table`, `readxl`, `ieugwasr`; **PLINK v1.9**
for clumping; 1000 Genomes Phase 3 EUR panel as the LD reference.

## Source data (not included here — obtain from the original sources)

Raw third-party files are **not redistributed**. Download them and place them in the
working directory before running the scripts:

- **eQTLGen** cis-eQTL and allele-frequency files — https://www.eqtlgen.org/
- **FinnGen R12** three-phenotype GWAS summary statistics and the `K11_PULP_PERIAPICAL`
  bridge file (used to restore rsIDs) — FinnGen public release
- **Finan et al.** druggable-genome catalogue (Supplementary Table S1),
  Sci Transl Med 2017;9:eaag1166
- **1000 Genomes** Phase 3 EUR panel — https://www.internationalgenome.org/
- Single-cell datasets — NCBI GEO **GSE181688** and **GSE197680**

## Key reproduction checkpoints

- Instrument set: **5447** instruments / 2602 genes
- Necrosis MR: 2545 genes tested, **19** at FDR < 0.05; PTK2 b = 0.114, p = 3.0e-7;
  ROCK1 FDR = 0.052 (borderline, retained as supportive)
- Colocalization: PTK2 650/0.796, ROCK1 218/0.789, DLK1 764/0.661, SIK2 415/0.632

## Notes

- **Two-stage clumping and the intermediate file.** The 5447 instrument set is produced
  in two stages: a permissive pass (r2 < 0.1, gene-body +/- 100 kb) that yields
  `exposure_instruments_final.rds` (40,287 variants), followed by a strict pass
  (r2 < 0.001, TSS +/- 100 kb) with per-(gene,SNP) aggregation that yields
  `..._clean.rds` (5447). The intermediate `exposure_instruments_final.rds` is provided
  directly so that the strict pass and all downstream results can be reproduced exactly
  without re-running the permissive clump. Regenerating that intermediate from scratch
  may differ by a few variants depending on the clumping implementation/version, which
  can shift the borderline gene count slightly; the four reported targets (PTK2, ROCK1,
  DLK1, SIK2) and their instruments are unaffected.
- **Sensitivity.** A single-stage r2 < 0.001 clump gives 6652 instruments
  (`..._6652_singlepass.rds`); the leading target PTK2 is unchanged.
- **Colocalization** uses all variants within TSS +/- 100 kb (not the clumped
  instruments), with rsID-based matching between the hg19 eQTL and hg38 GWAS data.

## License

MIT (see LICENSE). Third-party datasets remain under their original terms; cite the
original sources listed above.
