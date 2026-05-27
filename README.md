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

Correlations involving the species below are flagged as **literature-supported**
in the heatmap annotation; all others are flagged as **novel observations**.
The list is divided into bacteria enriched in CRC (pro-tumorigenic) and bacteria
depleted in CRC (protective commensals). Sources used are listed per species.

### Enriched in CRC

| Species | Mechanism | Citation |
|---|---|---|
| *Fusobacterium nucleatum* | FadA adhesin activates Wnt/β-catenin; NF-κB-driven inflammation | Bautista et al., 2026; Dalal et al., 2021; Kostic et al., 2013; Rubinstein et al., 2013; Tahara et al., 2014; Wang et al., 2024 |
| *Escherichia coli* (pks⁺) | Colibactin genotoxin causes DNA double-strand breaks and chromosomal instability | Bautista et al., 2026; Dalal et al., 2021; Kostic et al., 2013 |
| *Bacteroides fragilis* (ETBF) | BFT (fragilysin) cleaves E-cadherin; activates STAT3 and Wnt signalling | Bautista et al., 2026; Dalal et al., 2021 |
| *Streptococcus gallolyticus* | COX-2 activation; suppresses apoptosis; longstanding clinical CRC association | Dalal et al., 2021 |
| *Enterococcus faecalis* | Extracellular ROS production; DNA damage and chromosomal instability | Dalal et al., 2021 |
| *Peptostreptococcus anaerobius* | TLR2/TLR4 activation; elevated ROS and intracellular cholesterol; promotes dysplasia | Dalal et al., 2021; Liu et al., 2023 |
| *Peptostreptococcus stomatis* | Enriched in CRC tissue; part of validated metagenomic CRC signature | Yu et al., 2017 |
| *Parvimonas micra* | Enriched in CRC tissue; part of validated metagenomic CRC signature | Yu et al., 2017; Zhao et al., 2021 |
| *Solobacterium moorei* | Enriched in CRC tissue; part of validated metagenomic CRC signature | Yu et al., 2017 |
| *Helicobacter pylori* | Pro-inflammatory cytokines (IL-1, IL-6, TNF-α); independent adenoma/CRC risk factor | Dalal et al., 2021 |
| *Clostridium septicum* | α-toxin activates MAPK signalling; associated with occult colonic malignancy | Dalal et al., 2021 |
| *Salmonella enterica* | AvrA effector activates β-catenin and STAT3; flagellin antibodies elevated in CRC | Dalal et al., 2021 |
| *Porphyromonas gingivalis* | Immune suppression; promotes epithelial–mesenchymal transition (EMT) | Bautista et al., 2026; Sheng et al., 2024 |

### Depleted in CRC (protective commensals)

These species are expected to show **negative** correlations with pro-tumorigenic
genes and are flagged as literature-supported negative controls.

| Species | Mechanism | Citation |
|---|---|---|
| *Faecalibacterium prausnitzii* | Major butyrate producer; anti-inflammatory; consistently depleted in CRC | Bautista et al., 2026; Cao et al., 2025 |
| *Roseburia intestinalis* | Butyrate producer; loss characteristic of CRC-associated dysbiosis | Bautista et al., 2026 |
| *Roseburia hominis* | Butyrate producer; loss characteristic of CRC-associated dysbiosis | Bautista et al., 2026 |

---

## Literature Benchmarks (Breast — BRCA)

The BRCA benchmark list differs fundamentally from COAD. The dominant mechanism is **estrogen metabolism via the estrobolome** — bacteria producing β-glucuronidase, sulfatase, and hydroxysteroid dehydrogenase enzymes that reactivate conjugated oestrogens and elevate circulating hormone levels, driving hormone-receptor-positive (HR+) tumour growth. A second distinct mechanism operates specifically in triple-negative breast cancer (TNBC) via **TMAO-producing Clostridiales**, which promote antitumour immunity.

### Enriched in breast cancer (pro-tumorigenic or risk-associated)

| Species/Genus | Mechanism | Notes | Citation |
|---|---|---|---|
| *Escherichia coli* | β-glucuronidase; oestrogen deconjugation → elevated HR+ tumour growth; prognostic risk factor | Strong, consistent | Larnder et al. 2025; Mahno et al. 2024; Sheng et al. 2024 |
| *Bacteroides fragilis* | Broadest estrobolome enzyme activity (β-glucuronidase + sulfatase + 3β-HSD + 17β-HSD) | Strong | Larnder et al. 2025 |
| *Fusobacterium nucleatum* | Enriched in breast tumour tissue; pro-inflammatory; prognostic risk factor | Moderate | Sheng et al. 2024 |
| *Lactobacillus* spp. | Enriched in breast tumour tissue; β-glucuronidase activity | Moderate | Sheng et al. 2024; Larnder et al. 2025 |
| *Staphylococcus* spp. | Enriched in breast tumour tissue !known sequencing contaminant — interpret with caution! | Moderate/caveat | Sheng et al. 2024 |
| *Ruminococcus* spp. | β-glucuronidase activity; **conflicting evidence** across cohorts | Mixed | Larnder et al. 2025; Mahno et al. 2024 |
| *Klebsiella* spp. | β-glucuronidase activity; enriched in cases | Moderate | Mahno et al. 2024 |
| *Citrobacter* spp. | β-glucuronidase activity; enriched in cases | Moderate | Larnder et al. 2025 |
| *Enterobacter* spp. | β-glucuronidase activity; enriched in cases | Moderate | Larnder et al. 2025 |
| *Sphingomonas* spp. | Enriched in breast tumour tissue; **also top IGFBP7 correlator in COAD** !cross-cancer signal of interest! | Moderate | Sheng et al. 2024 |

### Depleted in breast cancer (protective)

| Species/Genus | Mechanism | Citation |
|---|---|---|
| *Faecalibacterium prausnitzii* | Directly suppresses breast cancer cell growth via IL-6/STAT3 inhibition; butyrate producer; depleted in cases | Ruo et al. 2021; Larnder et al. 2025 |
| *Roseburia inulinivorans* | β-glucosidase; activates protective phytoestrogens; depleted in cases | Larnder et al. 2025 |
| *Roseburia hominis / intestinalis* | Butyrate producers; depleted in dysbiotic microenvironment | Larnder et al. 2025 |
| *Bifidobacterium* spp. | Anti-inflammatory; β-glucuronidase depleted; consistently protective | Larnder et al. 2025; Mahno et al. 2024 |
| *Collinsella aerofaciens* | Depleted in cases; protective gut commensal | Larnder et al. 2025 |

### TMAO-producing Clostridiales — antitumour in TNBC

These bacteria are enriched specifically in the immunomodulatory (IM) subtype of TNBC and their
metabolite TMAO directly induces GSDME-mediated pyroptosis in tumour cells, enhancing CD8+ T cell
recruitment. Higher plasma TMAO correlates with better immunotherapy response in TNBC patients.

| Genera | Evidence | Citation |
|---|---|---|
| *Blautia*, *Dorea*, *Ruminococcus*, *Tyzzerella*, *Roseburia* (order Clostridiales) | TMAO → PERK ER stress → GSDME pyroptosis → CD8+ T cell immunity; n=360 TNBC cohort | Wang et al., 2022 |

### Key note on conflicting evidence
Several genera (*Ruminococcus*, *Fusobacterium*, *Prevotella*, *Faecalibacterium*, *Lactobacillus*)
show opposing associations in different breast cancer cohorts (Mahno et al. 2024). In the heatmap
annotation, these are flagged "literature-supported" but the manuscript methods section should
acknowledge this heterogeneity explicitly.

---

## Literature Benchmarks (Pancreatic — PAAD)

PAAD has three defining microbial features distinct from COAD and BRCA:
(1) **oral-to-pancreas bacterial translocation** — oral pathobionts enriched in tumour tissue and early cystic precursors;
(2) **intratumoral gemcitabine inactivation** — *Mycoplasma* and Gammaproteobacteria directly metabolise the standard chemotherapy drug;
(3) **TMAO anti-tumour immunity** — the same Clostridiales–TMAO mechanism from BRCA also operates in PDAC.

### Section A — Pro-tumorigenic: oral-to-pancreas translocation

| Species/Genus | Mechanism | Citation |
|---|---|---|
| *Porphyromonas gingivalis* | Promotes PDAC progression in K-rasG12D transgenic mouse models; induces immunosuppressive TME; enriched in early cystic precursors to invasive PDAC; pro-tumorigenic effects attenuated by *Lactobacillus* | Sheng et al. 2024; Wang et al. 2024 |
| *Fusobacterium nucleatum* | Oral pathobiont enriched in PAAD tumour tissue; pro-inflammatory microenvironment; pan-cancer prognostic risk factor | Sheng et al. 2024 |
| *Prevotella oris* | Oral pathobiont; enriched in early cystic PDAC precursors; shared prognostic risk bacterium across multiple cancer types | Sheng et al. 2024 |
| *Prevotella melaninogenica* | Oral pathobiont; shared pan-cancer prognostic risk bacterium | Sheng et al. 2024 |
| *Veillonella parvula* | Oral commensal; shared pan-cancer prognostic risk bacterium | Sheng et al. 2024 |

### Section B — Pro-tumorigenic: gemcitabine inactivation

| Species/Group | Mechanism | Citation |
|---|---|---|
| *Mycoplasma* spp. | Metabolises and inactivates gemcitabine intracellularly; key determinant of chemotherapy resistance and treatment failure in PAAD | Geller et al., 2017, cited in Sheng et al. 2024 |
| Gammaproteobacteria (class) | Cytidine deaminase activity inactivates gemcitabine within the tumour microenvironment; broad class-level benchmark for drug resistance | Bautista et al. 2026 |

### Section C — Antitumour: TMAO-producing bacteria (cross-cancer bridge with BRCA)

| Species/Genus | Mechanism | Citation |
|---|---|---|
| Clostridiales order (*Blautia*, *Dorea*, *Ruminococcus*, *Roseburia*, *Bacillus*) | Produce TMAO via choline trimethylamine-lyase; TMAO activates type I IFN in macrophages → dendritic cell and cytotoxic T cell activation → sensitisation to immune checkpoint inhibitors in PDAC; high plasma TMAO correlates with better ICI response | Mirji et al., 2022, cited in Liu et al. 2023 and W. Zhang et al. 2024 |

> **Cross-cancer note:** This is the same Clostridiales–TMAO mechanism documented in BRCA
> (Wang et al. *Cell Metabolism* 2022). If these genera show consistent positive gene correlations
> in both BRCA and PAAD heatmaps, this constitutes a strong cross-cancer mechanistic story.

### Section D — Protective benchmark

| Species/Genus | Mechanism | Citation |
|---|---|---|
| *Lactobacillus* spp. | Attenuates pro-tumorigenic effects of *P. gingivalis* in pancreatic cancer models; dietary *Lactobacillus*-derived exopolysaccharides also enhance immune checkpoint blockade efficacy | Wang et al. 2024; Liu et al. 2023 |

### Evidence confidence summary

| Category | Confidence | Basis |
|---|---|---|
| *P. gingivalis* pro-tumorigenic | High | K-rasG12D mouse model + pan-cancer atlas |
| *Mycoplasma* gemcitabine resistance | High | Mechanistic *Science* paper + two independent reviews |
| Gammaproteobacteria gemcitabine resistance | High | Mechanistic evidence, Bautista et al. 2026 |
| Oral microbiota group (translocation) | Moderate | Pan-cancer atlas; no individual species RCT |
| TMAO / Clostridiales antitumour | Moderate–High | PDAC mouse model; cross-validated in BRCA |
| *Lactobacillus* protective | Moderate | Mechanistic mouse model; species unspecified |
