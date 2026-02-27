#!/usr/bin/env Rscript
# =============================================================================
# beta_stats.R — Beta diversity: CLR-PCoA + PERMANOVA with covariates
# =============================================================================
# Usage: Rscript beta_stats.R <table_tsv> <metadata_tsv> <bray_dm_tsv>
#                             <wunifrac_dm_tsv> <group_col> <covariates_csv>
#                             <out_dir>
#
# Required: vegan, ggplot2, compositions, ape, dplyr
# Install:  mamba install -n qiime2 -c conda-forge r-vegan r-ggplot2 r-compositions r-ape r-dplyr
# =============================================================================

suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(ape)
  library(dplyr)
})

# ── Try to load compositions for CLR (optional) ───────────────────────────────
has_compositions <- requireNamespace("compositions", quietly = TRUE)

# ── Arguments ─────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) stop("Usage: beta_stats.R <table_tsv> <metadata> <bray_dm> <wunifrac_dm> <group_col> <covariates_csv> <out_dir>")

table_tsv    <- args[1]
metadata_file <- args[2]
bray_file    <- args[3]
wunifrac_file <- args[4]
group_col    <- args[5]
covariates   <- if (nchar(args[6]) > 0) strsplit(args[6], ",")[[1]] else character(0)
out_dir      <- args[7]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load data ─────────────────────────────────────────────────────────────────
meta <- read.table(metadata_file, header = TRUE, sep = "\t",
                   comment.char = "", row.names = 1, check.names = FALSE)

# Feature table (skip biom comment lines starting with #)
raw_tab <- read.table(table_tsv, header = TRUE, sep = "\t",
                      skip = 1, row.names = 1, check.names = FALSE)
# Transpose: samples as rows
tab <- t(raw_tab)

# Keep only samples present in both metadata and feature table
common_samples <- intersect(rownames(meta), rownames(tab))
if (length(common_samples) < 2) stop("Fewer than 2 common samples between table and metadata.")
meta <- meta[common_samples, , drop = FALSE]
tab  <- tab[common_samples, , drop = FALSE]

# ── CLR transformation for Aitchison PCoA ────────────────────────────────────
clr_tab <- NULL
if (has_compositions) {
  library(compositions)
  # Replace zeros with small pseudocount before CLR
  tab_pseudo <- tab + 0.5
  clr_tab <- as.matrix(compositions::clr(tab_pseudo))
  message("CLR transformation applied (compositions package)")
} else {
  # Fallback: manual CLR
  tab_pseudo <- tab + 0.5
  gm <- apply(tab_pseudo, 1, function(x) exp(mean(log(x))))
  clr_tab <- log(tab_pseudo / gm)
  message("CLR transformation applied (manual fallback)")
}

# ── Aitchison distance (Euclidean on CLR) PCoA ───────────────────────────────
ait_dist <- dist(clr_tab)
pcoa_ait <- pcoa(ait_dist)
pcoa_df  <- as.data.frame(pcoa_ait$vectors[, 1:2])
colnames(pcoa_df) <- c("PC1", "PC2")
pcoa_df  <- merge(pcoa_df, meta, by = "row.names", all.x = TRUE)
colnames(pcoa_df)[1] <- "SampleID"

pct_var <- round(pcoa_ait$values$Relative_eig[1:2] * 100, 1)

# ── PERMANOVA (adonis2) using distance matrices ───────────────────────────────
load_dm <- function(f) {
  dm <- read.table(f, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  dm <- dm[common_samples, common_samples]
  as.dist(dm)
}

permanova_results <- list()

run_permanova <- function(dist_obj, label) {
  all_terms <- c(group_col, covariates[covariates %in% colnames(meta)])
  formula_str <- paste("dist_obj ~", paste(all_terms, collapse = " + "))
  tryCatch({
    perm <- adonis2(as.formula(formula_str), data = meta, permutations = 999)
    perm_df <- as.data.frame(perm)
    perm_df$term   <- rownames(perm_df)
    perm_df$metric <- label
    message(label, " PERMANOVA:\n", paste(capture.output(perm), collapse = "\n"))
    perm_df
  }, error = function(e) {
    message("PERMANOVA failed for ", label, ": ", e$message)
    NULL
  })
}

if (file.exists(bray_file)) {
  bray_dist <- load_dm(bray_file)
  permanova_results[["Bray-Curtis"]] <- run_permanova(bray_dist, "Bray-Curtis")
}
if (file.exists(wunifrac_file)) {
  wu_dist <- load_dm(wunifrac_file)
  permanova_results[["Weighted-UniFrac"]] <- run_permanova(wu_dist, "Weighted-UniFrac")
}

perm_df_all <- do.call(rbind, Filter(Negate(is.null), permanova_results))
if (!is.null(perm_df_all) && nrow(perm_df_all) > 0) {
  write.table(perm_df_all, file.path(out_dir, "permanova_results.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message("Saved: permanova_results.tsv")
} else {
  # Create empty file so Snakemake target is satisfied
  write.table(data.frame(), file.path(out_dir, "permanova_results.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# ── PCoA plots ────────────────────────────────────────────────────────────────
pdf(file.path(out_dir, "pcoa_plots.pdf"), width = 9, height = 7)

# Aitchison (CLR) PCoA
if (nrow(pcoa_df) >= 2) {
  p_ait <- ggplot(pcoa_df, aes_string(x = "PC1", y = "PC2",
                                       color = group_col, label = "SampleID")) +
    geom_point(size = 4, alpha = 0.85) +
    ggrepel_labels(pcoa_df) +
    xlab(paste0("PC1 (", pct_var[1], "%)")) +
    ylab(paste0("PC2 (", pct_var[2], "%)")) +
    stat_ellipse(aes_string(group = group_col), type = "t", level = 0.8, linetype = 2) +
    theme_bw(base_size = 13) +
    theme(legend.position = "right") +
    labs(title = paste("Aitchison PCoA (CLR) by", group_col),
         subtitle = "Ellipses: 80% confidence interval")
  print(p_ait)
}

# Bray-Curtis PCoA (from distance matrix)
if (file.exists(bray_file)) {
  bray_pcoa <- pcoa(load_dm(bray_file))
  bc_df <- as.data.frame(bray_pcoa$vectors[, 1:2])
  colnames(bc_df) <- c("PC1", "PC2")
  bc_df <- merge(bc_df, meta, by = "row.names", all.x = TRUE)
  bc_pct <- round(bray_pcoa$values$Relative_eig[1:2] * 100, 1)

  p_bc <- ggplot(bc_df, aes_string(x = "PC1", y = "PC2", color = group_col)) +
    geom_point(size = 4, alpha = 0.85) +
    xlab(paste0("PC1 (", bc_pct[1], "%)")) +
    ylab(paste0("PC2 (", bc_pct[2], "%)")) +
    stat_ellipse(aes_string(group = group_col), type = "t", level = 0.8, linetype = 2) +
    theme_bw(base_size = 13) +
    labs(title = paste("Bray-Curtis PCoA by", group_col))
  print(p_bc)
}

dev.off()
message("Saved: pcoa_plots.pdf")

# ── Helper: ggrepel labels if available ──────────────────────────────────────
ggrepel_labels <- function(df) {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_text_repel(size = 3, max.overlaps = Inf)
  } else {
    geom_text(size = 3, vjust = -0.8)
  }
}
