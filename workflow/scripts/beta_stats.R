#!/usr/bin/env Rscript
# =============================================================================
# beta_stats.R — Robust Beta diversity: CLR-PCoA + PERMANOVA
# =============================================================================

suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(ape)
  library(dplyr)
})

# ── Arguments ─────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) stop("Usage: beta_stats.R <table_tsv> <metadata> <bray_dm> <wunifrac_dm> <group_col> <covariates_csv> <out_dir>")

table_tsv     <- args[1]
metadata_file <- args[2]
bray_file     <- args[3]
wunifrac_file <- args[4]
group_col     <- args[5]
covariates    <- if (nchar(args[6]) > 0) strsplit(args[6], ",")[[1]] else character(0)
out_dir       <- args[7]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load data (Defensive Loading) ─────────────────────────────────────────────
meta <- read.table(metadata_file, header = TRUE, sep = "\t",
                   comment.char = "", row.names = 1, check.names = FALSE, quote = "")

raw_tab <- read.table(table_tsv, header = TRUE, sep = "\t",
                      skip = 1, row.names = 1, check.names = FALSE, 
                      comment.char = "", quote = "")
tab <- t(raw_tab) 

common_samples <- intersect(rownames(meta), rownames(tab))
if (length(common_samples) < 2) {
    stop("Fewer than 2 common samples. Verify Sample IDs match exactly.")
}

meta <- meta[common_samples, , drop = FALSE]
tab  <- tab[common_samples, , drop = FALSE]

# ── Helper for Distance Matrices ──────────────────────────────────────────────
load_dm <- function(f, samples) {
  if (!file.exists(f)) return(NULL)
  dm <- read.table(f, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  common_dm <- intersect(samples, rownames(dm))
  if (length(common_dm) < 2) return(NULL)
  as.dist(dm[common_dm, common_dm])
}

# ── CLR Transformation & Aitchison PCoA ──────────────────────────────────────
tab_pseudo <- tab + 0.5
gm <- apply(tab_pseudo, 1, function(x) exp(mean(log(x))))
clr_tab <- log(tab_pseudo / gm)

ait_dist <- dist(clr_tab)
pcoa_ait <- pcoa(ait_dist)
pcoa_df  <- as.data.frame(pcoa_ait$vectors[, 1:2])
colnames(pcoa_df) <- c("PC1", "PC2")
pcoa_df$SampleID <- rownames(pcoa_df)
pcoa_df  <- merge(pcoa_df, tibble::rownames_to_column(meta, "SampleID"), by = "SampleID")

pct_var <- round(pcoa_ait$values$Relative_eig[1:2] * 100, 1)

# ── PERMANOVA (adonis2) ───────────────────────────────────────────────────────
permanova_results <- list()

run_permanova <- function(dist_obj, label, metadata) {
  valid_covs <- covariates[covariates %in% colnames(metadata)]
  all_terms <- c(group_col, valid_covs)
  formula_str <- paste("dist_obj ~", paste(all_terms, collapse = " + "))
  
  tryCatch({
    perm <- adonis2(as.formula(formula_str), data = metadata, permutations = 999)
    perm_df <- as.data.frame(perm)
    perm_df$term   <- rownames(perm_df)
    perm_df$metric <- label
    perm_df
  }, error = function(e) {
    message("PERMANOVA failed for ", label, ": ", e$message)
    NULL
  })
}

permanova_results[["Aitchison"]] <- run_permanova(ait_dist, "Aitchison", meta)

bray_dist <- load_dm(bray_file, common_samples)
if (!is.null(bray_dist)) {
  permanova_results[["Bray-Curtis"]] <- run_permanova(bray_dist, "Bray-Curtis", meta)
}

wu_dist <- load_dm(wunifrac_file, common_samples)
if (!is.null(wu_dist)) {
  permanova_results[["Weighted-UniFrac"]] <- run_permanova(wu_dist, "Weighted-UniFrac", meta)
}

# ── FIXED: Combine and Save Stats ─────────────────────────────────────────────
# We actually need to create the perm_df_all variable from the list
perm_df_all <- do.call(rbind, Filter(Negate(is.null), permanova_results))

if (!is.null(perm_df_all) && nrow(perm_df_all) > 0) {
  write.table(perm_df_all, file.path(out_dir, "permanova_results.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message("Saved: permanova_results.tsv")
} else {
  write.table(data.frame(Status="No_Results"), file.path(out_dir, "permanova_results.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# ── Plotting ──────────────────────────────────────────────────────────────────
pdf(file.path(out_dir, "pcoa_plots.pdf"), width = 10, height = 8)

p_ait <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = !!sym(group_col))) +
  geom_point(size = 4, alpha = 0.7) +
  stat_ellipse(level = 0.8, linetype = 2) +
  theme_bw() +
  labs(title = "Aitchison PCoA (CLR-Euclidean)",
       x = paste0("PC1 (", pct_var[1], "%)"),
       y = paste0("PC2 (", pct_var[2], "%)"))

if (requireNamespace("ggrepel", quietly = TRUE)) {
  p_ait <- p_ait + ggrepel::geom_text_repel(aes(label = SampleID), size = 2, max.overlaps = 10)
}
print(p_ait)

if (!is.null(bray_dist)) {
  bc_pcoa <- pcoa(bray_dist)
  bc_df <- as.data.frame(bc_pcoa$vectors[, 1:2])
  colnames(bc_df) <- c("PC1", "PC2")
  bc_df$SampleID <- rownames(bc_df)
  bc_df <- merge(bc_df, tibble::rownames_to_column(meta, "SampleID"), by = "SampleID")
  bc_pct <- round(bc_pcoa$values$Relative_eig[1:2] * 100, 1)
  
  p_bc <- ggplot(bc_df, aes(x = PC1, y = PC2, color = !!sym(group_col))) +
    geom_point(size = 4, alpha = 0.7) +
    stat_ellipse(level = 0.8, linetype = 2) +
    theme_bw() +
    labs(title = "Bray-Curtis PCoA",
         x = paste0("PC1 (", bc_pct[1], "%)"),
         y = paste0("PC2 (", bc_pct[2], "%)"))
  print(p_bc)
}

dev.off()