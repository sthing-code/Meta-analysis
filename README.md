# Intratumoral Microbiome × Cancer Gene Expression — Correlation Analysis

## Overview

This repository contains the full analysis pipeline correlating intratumoral bacterial abundance
(from TCMbio) with the expression of cancer-related genes (from TCGA) across multiple cancer types.
The primary output is a series of Pearson correlation heatmaps, one per cancer type, intended as
manuscript figures.

**Cancer types covered (in order):**
- Colon adenocarcinoma (COAD) — primary analysis
- Breast invasive carcinoma (BRCA) — extension
- Pancreatic adenocarcinoma (PAAD) — extension
- Prostate adenocarcinoma (PRAD) — extension

---

## Priority Genes

| Gene symbol | Protein | Notes |
|---|---|---|
| IGFBP7 | IGF-binding protein 7 | Strong positive correlations flagged |
| LGALS1 | Galectin-1 | Strong positive correlations flagged |
| VEGFA | VEGF-A | Downstream: MAPK pathway |
| MMP9 | MMP-9 | |
| IL6 | IL-6 | Downstream: STAT3 |
| STAT3 | STAT3 | |
| APC | APC | Downstream: WNT/β-catenin |
| BRCA2 | BRCA2 | |
| CTSB | Cathepsin B | Strong negative correlations flagged |
| CTSD | Cathepsin D | Strong negative correlations flagged |
| CTSL | Cathepsin L | Strong negative correlations flagged |
| CTSS | Cathepsin S | Strong negative correlations flagged |
| CTSV | Cathepsin V | Strong negative correlations flagged |

---

## Data Sources

### Gene expression
- **Source:** TCGA via Bioconductor `TCGAbiolinks`
- **Workflow:** STAR — Counts (GDC harmonised, hg38)
- **Normalisation:** TPM (transcripts per million), computed from STAR raw counts + gene lengths
- **Sample filter:** Primary tumour samples only (TCGA barcode suffix `-01`)

### Bacterial abundance
- **Source:** TCMbio (https://microbiomex.sdu.edu.cn/) — cancer-type-specific download
- **Format:** Species × sample abundance matrix (CSV), downloaded manually per cancer type
- **Preprocessing:** log₁₀(x + 1) transformation; prevalence filter ≥ 10% of samples

---

## Project Structure

```
.
├── README.md
├── .gitignore
│
├── R/                          # Shared helper functions
│   ├── utils.R                 # Barcode matching, TPM computation, FDR helpers
│   └── plot_themes.R           # ggplot2 / ComplexHeatmap theme settings
│
├── scripts/
│   ├── 01_download_TCGA.R      # Download TCGA data for all cancer types
│   ├── 02_preprocess_TCGA.R    # Normalise, filter, extract priority genes
│   ├── 03_preprocess_TCMbio.R  # Load and clean TCMbio bacterial abundance data
│   ├── 04_correlate.R          # Pearson correlations + FDR correction
│   ├── 05_heatmap.R            # Main correlation heatmaps (one per cancer type)
│   └── 06_threeway.R           # Three-way correlation analysis (IL6/VEGFA/APC axes)
│
├── data/
│   ├── raw/                    # Downloaded files (not tracked by git — see .gitignore)
│   └── processed/              # RDS files output by preprocessing scripts
│
└── figures/
    ├── colon/
    ├── breast/
    ├── pancreatic/
    └── prostate/
```

---

## Reproduction Steps

Run scripts in order:

```r
source("scripts/01_download_TCGA.R")      # ~15–30 min depending on connection
source("scripts/02_preprocess_TCGA.R")
# Manually download TCMbio data — see instructions in script 03
source("scripts/03_preprocess_TCMbio.R")
source("scripts/04_correlate.R")
source("scripts/05_heatmap.R")
source("scripts/06_threeway.R")
```

---

## Dependencies

```r
# Bioconductor
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "DESeq2",
                       "ComplexHeatmap", "circlize"))

# CRAN
install.packages(c("tidyverse", "data.table", "corrplot",
                   "RColorBrewer", "ggrepel", "patchwork", "here"))
```

R version: ≥ 4.3.0 recommended.

---

## Correlation Thresholds and Flagging

| Threshold | Use |
|---|---|
| \|r\| > 0.15 | Supplementary heatmap (broader signal) |
| \|r\| > 0.20 | Main manuscript heatmap |
| q < 0.05 (BH) | Significance annotation on heatmap |

**Flagging logic:**
- IGFBP7, Galectin-1: stars on cells where r > 0.20 and q < 0.05
- Cathepsins: stars on cells where r < −0.20 and q < 0.05

---

## Literature Benchmarks (Colon)

The following bacteria have published CRC-association evidence and serve as
positive controls. Correlations involving these species are flagged as
**literature-supported**; all others are flagged as **novel observations**.

- *Fusobacterium nucleatum* — most validated CRC intratumoral bacterium
- *Parvimonas micra*
- *Peptostreptococcus stomatis*
- *Solobacterium moorei*
- *Bacteroides fragilis* (enterotoxigenic ETBF)
- *Campylobacter* spp.
- *Leptotrichia* spp.
- *Porphyromonas gingivalis*

---

## Publication Notes

Primary goal: one comparative paper across cancer types.
If a strong cancer-type-specific or gene-specific story emerges (e.g. same gene
behaves directionally opposite across two cancer types), scope for a separate
focused paper exists.
