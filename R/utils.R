# R/utils.R
# Shared helper functions for the microbiome-cancer correlation pipeline.
# Source this file at the top of every analysis script.

# ── Gene lists ────────────────────────────────────────────────────────────────

# Priority genes, grouped by biological function.
# These are the column variables in every heatmap.
PRIORITY_GENES <- list(
  igfbp_galectin = c("IGFBP7", "LGALS1"),      # Galectin-1 gene symbol is LGALS1
  angiogenesis   = c("VEGFA", "MMP9"),
  jak_stat       = c("IL6", "STAT3"),
  wnt            = c("APC"),
  dna_repair     = c("BRCA2"),
  cathepsins     = c("CTSB", "CTSD", "CTSL", "CTSS", "CTSV")
)

# Flat vector of all priority genes (for TCGA query)
ALL_PRIORITY_GENES <- unlist(PRIORITY_GENES, use.names = FALSE)

# Downstream pathway genes added for three-way correlation analysis
PATHWAY_GENES <- list(
  il6_stat3  = c("IL6", "STAT3"),
  vegfa_mapk = c("VEGFA", "MAPK1", "MAPK3", "KRAS"),
  apc_wnt    = c("APC", "CTNNB1", "MYC", "CCND1")
)

ALL_PATHWAY_GENES <- unique(unlist(PATHWAY_GENES, use.names = FALSE))

# ── TCGA cancer-type config ───────────────────────────────────────────────────

CANCER_CONFIG <- list(
  colon = list(
    tcga_project  = "TCGA-COAD",
    tcmbio_file   = "data/raw/TCMbio_COAD_bacteria.csv",
    figures_dir   = "figures/colon",
    label         = "Colon adenocarcinoma (COAD)"
  ),
  breast = list(
    tcga_project  = "TCGA-BRCA",
    tcmbio_file   = "data/raw/TCMbio_BRCA_bacteria.csv",
    figures_dir   = "figures/breast",
    label         = "Breast invasive carcinoma (BRCA)"
  ),
  pancreatic = list(
    tcga_project  = "TCGA-PAAD",
    tcmbio_file   = "data/raw/TCMbio_PAAD_bacteria.csv",
    figures_dir   = "figures/pancreatic",
    label         = "Pancreatic adenocarcinoma (PAAD)"
  ),
  prostate = list(
    tcga_project  = "TCGA-PRAD",
    tcmbio_file   = "data/raw/TCMbio_PRAD_bacteria.csv",
    figures_dir   = "figures/prostate",
    label         = "Prostate adenocarcinoma (PRAD)"
  )
)

# ── Cancer-type-specific heatmap thresholds ──────────────────────────────────
# COAD has strong intratumoral microbiome–gene correlations; r > 0.20 is appropriate.
# BRCA signal is weaker overall; r > 0.15 is used as the primary heatmap threshold.
# PAAD and PRAD thresholds set to 0.15 pending data (adjust after running script 04).
HEATMAP_THRESHOLDS <- list(
  colon      = list(primary = 0.20, supplementary = 0.15),
  breast     = list(primary = 0.15, supplementary = 0.10),
  pancreatic = list(primary = 0.15, supplementary = 0.10),
  prostate   = list(primary = 0.15, supplementary = 0.10)
)

# ── Literature-supported CRC bacteria ────────────────────────────────────────
# Used to flag correlations as literature-supported vs novel.
# Extend analogous lists for other cancer types as you add them.

BENCHMARK_BACTERIA <- list(
  colon = c(
    # ── Enriched in CRC (pro-tumorigenic) ──────────────────────────────────
    # Kostic et al. Genome Res 2012; Rubinstein et al. Cell Host Microbe 2013
    "Fusobacterium nucleatum",
    # Yu et al. Gut 2017; Proctor et al. NAR Cancer 2021
    "Parvimonas micra",
    "Peptostreptococcus stomatis",
    "Solobacterium moorei",
    # Dalal et al. 2021; Wu et al. Nat Med 2004
    "Bacteroides fragilis",
    # Dalal et al. 2021
    "Streptococcus gallolyticus",
    "Enterococcus faecalis",
    "Peptostreptococcus anaerobius",
    "Helicobacter pylori",
    "Clostridium septicum",
    "Salmonella enterica",
    # Bautista et al. 2026; Sheng et al. 2024
    "Porphyromonas gingivalis",
    # ── Depleted in CRC (protective commensals) ────────────────────────────
    # Flagged as negative controls — strong negative correlation expected
    # Bautista et al. 2026; Cao et al. 2025
    "Faecalibacterium prausnitzii",
    # Bautista et al. 2026
    "Roseburia intestinalis",
    "Roseburia hominis"
  ),
  breast = c(
    # ── Enriched in breast cancer / estrobolome-active ─────────────────────
    # Larnder et al. 2025; Mahno et al. 2024; Sheng et al. 2024
    "Escherichia coli",
    # Larnder et al. 2025 — broadest estrobolome enzyme activity
    "Bacteroides fragilis",
    # Sheng et al. 2024
    "Fusobacterium nucleatum",
    "Lactobacillus iners",
    "Staphylococcus aureus",       # contamination caveat — flag in figure
    "Staphylococcus epidermidis",  # contamination caveat — flag in figure
    # Mahno et al. 2024; Larnder et al. 2025 (conflicting evidence)
    "Ruminococcus gnavus",
    # Larnder et al. 2025
    "Klebsiella pneumoniae",
    "Citrobacter freundii",
    "Enterobacter hormaechei",
    # Sheng et al. 2024 — also top IGFBP7 correlator in COAD (cross-cancer interest)
    "Sphingomonas panacisoli",
    "Sphingomonas paucimobilis",
    # ── Depleted in breast cancer / protective ──────────────────────────────
    # Ruo et al. 2021; Larnder et al. 2025
    "Faecalibacterium prausnitzii",
    # Larnder et al. 2025
    "Roseburia inulinivorans",
    "Roseburia hominis",
    "Roseburia intestinalis",
    "Bifidobacterium longum",
    "Bifidobacterium adolescentis",
    "Collinsella aerofaciens",
    # ── TMAO-producing Clostridiales — antitumour in TNBC ──────────────────
    # Wang et al. 2022 Cell Metabolism — TMAO → pyroptosis → CD8+ T cells
    "Blautia obeum",
    "Dorea longicatena",
    "Ruminococcus torques"
  ),
  pancreatic = c(
    # ── Pro-tumorigenic: oral-to-pancreas translocation ───────────────────
    # Sheng et al. 2024; Wang et al. 2024 (K-rasG12D mouse model)
    "Porphyromonas gingivalis",
    # Sheng et al. 2024 — oral microbiota enriched in early PDAC precursors
    "Fusobacterium nucleatum",
    "Prevotella oris",
    "Prevotella melaninogenica",
    "Veillonella parvula",
    # ── Pro-tumorigenic: gemcitabine inactivation ─────────────────────────
    # Geller et al. Science 2017, cited in Sheng et al. 2024
    "Mycoplasma hyorhinis",
    # Bautista et al. 2026 — Gammaproteobacteria inactivate gemcitabine
    # via cytidine deaminase; class-level benchmark
    "Pseudomonas aeruginosa",     # Gammaproteobacteria representative
    "Klebsiella pneumoniae",      # Gammaproteobacteria representative
    # ── TMAO-producing bacteria — antitumour in PAAD ──────────────────────
    # Liu et al. 2023; W. Zhang et al. 2024
    # TMAO → type I IFN → DC + cytotoxic T cell activation → ICI sensitisation
    "Blautia obeum",
    "Dorea longicatena",
    "Faecalibacterium prausnitzii",
    # ── Protective benchmark ──────────────────────────────────────────────
    # Wang et al. 2024 — attenuates P. gingivalis pro-tumorigenic effects
    "Lactobacillus acidophilus",
    "Lactobacillus rhamnosus"
  ),
  prostate = c(
    "Fusobacterium nucleatum",
    "Propionibacterium acnes",  # now Cutibacterium acnes
    "Trichomonas vaginalis"     # not a bacterium but often included
  )
)

# ── Barcode utilities ─────────────────────────────────────────────────────────

#' Truncate a TCGA barcode to 12-character patient ID
#' e.g. "TCGA-AA-3977-01A-01R-1635-07" → "TCGA-AA-3977"
shorten_barcode <- function(x) substr(x, 1, 12)

#' Truncate to 15 characters (patient + sample type)
#' e.g. "TCGA-AA-3977-01A-01R-1635-07" → "TCGA-AA-3977-01"
sample_barcode <- function(x) substr(x, 1, 15)

#' Extract sample type code from TCGA barcode (positions 14–15)
#' "01" = Primary tumour, "11" = solid normal tissue
get_sample_type <- function(x) substr(x, 14, 15)

#' Filter a vector of barcodes to primary tumour samples only
filter_primary_tumour <- function(barcodes) {
  barcodes[get_sample_type(barcodes) == "01"]
}

# ── TPM computation ───────────────────────────────────────────────────────────

#' Compute TPM from a raw counts matrix and a gene-length vector.
#' counts : genes × samples integer matrix (STAR unstranded counts)
#' gene_lengths : named numeric vector (gene name → length in bp)
#' Returns a genes × samples TPM matrix.
counts_to_tpm <- function(counts, gene_lengths) {
  stopifnot(all(rownames(counts) %in% names(gene_lengths)))
  lengths  <- gene_lengths[rownames(counts)]
  rpk      <- counts / (lengths / 1e3)           # reads per kilobase
  col_sums <- colSums(rpk, na.rm = TRUE)
  tpm      <- t(t(rpk) / (col_sums / 1e6))       # scale to per million
  return(tpm)
}

# ── Bacterial abundance preprocessing ────────────────────────────────────────

#' Log-transform bacterial abundance: log10(x + 1)
log_transform_abundance <- function(mat) {
  log10(mat + 1)
}

#' Filter bacteria by prevalence across samples.
#' Keeps only species present (abundance > 0) in >= min_prev fraction of samples.
filter_by_prevalence <- function(mat, min_prev = 0.10) {
  prev <- rowMeans(mat > 0)
  mat[prev >= min_prev, , drop = FALSE]
}

# ── Correlation utilities ─────────────────────────────────────────────────────

#' Compute Pearson r and two-sided p-value between two numeric vectors.
#' Returns NA silently if fewer than 10 finite paired observations exist.
pearson_pair <- function(x, y) {
  valid <- is.finite(x) & is.finite(y)
  if (sum(valid) < 10) return(c(r = NA_real_, p = NA_real_))
  ct <- cor.test(x[valid], y[valid], method = "pearson")
  c(r = unname(ct$estimate), p = ct$p.value)
}

#' Compute full Pearson correlation matrix (bacteria × genes) with p-values.
#'
#' bact_mat : bacteria × samples (rows = species, cols = samples)
#' expr_mat : genes × samples   (rows = genes, cols = samples)
#'
#' Returns a list:
#'   $r  — bacteria × genes correlation matrix
#'   $p  — bacteria × genes raw p-value matrix
#'   $q  — bacteria × genes BH-adjusted q-value matrix
compute_pearson_matrix <- function(bact_mat, expr_mat) {

  # Ensure samples are aligned
  common_samples <- intersect(colnames(bact_mat), colnames(expr_mat))
  if (length(common_samples) < 10) {
    stop("Fewer than 10 matched samples — check barcode formatting.")
  }
  message("  Matched samples: ", length(common_samples))

  bact_mat <- bact_mat[, common_samples, drop = FALSE]
  expr_mat <- expr_mat[, common_samples, drop = FALSE]

  n_bact <- nrow(bact_mat)
  n_gene <- nrow(expr_mat)

  r_mat <- matrix(NA_real_, nrow = n_bact, ncol = n_gene,
                  dimnames = list(rownames(bact_mat), rownames(expr_mat)))
  p_mat <- r_mat

  for (i in seq_len(n_bact)) {
    for (j in seq_len(n_gene)) {
      res            <- pearson_pair(bact_mat[i, ], expr_mat[j, ])
      r_mat[i, j]   <- res["r"]
      p_mat[i, j]   <- res["p"]
    }
  }

  # BH FDR correction across all tests in the matrix
  q_mat            <- matrix(p.adjust(as.vector(p_mat), method = "BH"),
                             nrow = n_bact, ncol = n_gene,
                             dimnames = dimnames(p_mat))

  list(r = r_mat, p = p_mat, q = q_mat,
       n_samples = length(common_samples))
}

#' Filter correlation results to retain only |r| >= threshold.
#' Returns the same list structure but with rows (bacteria) that have at least
#' one gene meeting the threshold.
filter_by_r <- function(cor_list, r_threshold = 0.20) {
  keep <- rowSums(abs(cor_list$r) >= r_threshold, na.rm = TRUE) > 0
  list(
    r = cor_list$r[keep, , drop = FALSE],
    p = cor_list$p[keep, , drop = FALSE],
    q = cor_list$q[keep, , drop = FALSE],
    n_samples = cor_list$n_samples
  )
}

# ── Flagging helpers ──────────────────────────────────────────────────────────

#' Given a gene name, return whether it should be flagged for strong POSITIVE r
is_positive_flag_gene <- function(gene) gene %in% c("IGFBP7", "LGALS1")

#' Given a gene name, return whether it should be flagged for strong NEGATIVE r
is_negative_flag_gene <- function(gene) gene %in% PRIORITY_GENES$cathepsins

#' Given a bacterium name and a cancer type, classify as benchmark or novel.
#' Handles TCMbio-format names (s__Genus_species) automatically.
flag_bacterium <- function(bacterium, cancer_type = "colon") {
  benchmarks <- BENCHMARK_BACTERIA[[cancer_type]]
  if (is.null(benchmarks)) return("novel")
  plain <- normalise_species_name(bacterium)
  genus <- strsplit(plain, " ")[[1]][1]
  if (plain %in% benchmarks ||
      any(grepl(paste0("^", genus, " "), benchmarks))) {
    "literature-supported"
  } else {
    "novel"
  }
}

#' Normalise a TCMbio species name to plain text for benchmark matching.
#' Strips "s__" prefix and converts underscores to spaces.
#' e.g. "s__Bacteroides_fragilis" -> "Bacteroides fragilis"
normalise_species_name <- function(x) {
  x <- sub("^s__", "", x)
  x <- gsub("_", " ", x)
  trimws(x)
}
