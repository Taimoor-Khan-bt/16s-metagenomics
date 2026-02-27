#!/usr/bin/env Rscript
# =============================================================================
# lefse_analysis.R — LEfSe differential analysis via microbiomeMarker + R
# =============================================================================
# Usage: Rscript lefse_analysis.R <genus_tsv> <taxonomy_tsv> <metadata_tsv>
#                                 <group_col> <out_dir>
#
# Required: microbiomeMarker, phyloseq, ggplot2, dplyr
# Install (conda):  mamba install -n qiime2 -c conda-forge -c bioconda bioconductor-microbiomemarker
# Install (R):      BiocManager::install("microbiomeMarker")
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# ── Arguments ─────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("Usage: lefse_analysis.R <genus_tsv> <taxonomy_tsv> <metadata> <group_col> <out_dir>")

genus_tsv  <- args[1]
tax_tsv    <- args[2]
meta_file  <- args[3]
group_col  <- args[4]
out_dir    <- args[5]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load feature table ────────────────────────────────────────────────────────
tab <- read.table(genus_tsv, header = TRUE, sep = "\t",
                  skip = 1, row.names = 1, check.names = FALSE)

# ── Load metadata ─────────────────────────────────────────────────────────────
meta <- read.table(meta_file, header = TRUE, sep = "\t",
                   comment.char = "", row.names = 1, check.names = FALSE)

# ── Check for microbiomeMarker ────────────────────────────────────────────────
if (!requireNamespace("microbiomeMarker", quietly = TRUE)) {
  message("microbiomeMarker not installed. Falling back to manual Wilcoxon LEfSe approximation.")

  # Manual fallback: Wilcoxon + effect size (log2FC)
  common_samples <- intersect(colnames(tab), rownames(meta))
  if (length(common_samples) < 2) stop("No common samples between table and metadata.")
  tab_sub  <- tab[, common_samples, drop = FALSE]
  meta_sub <- meta[common_samples, , drop = FALSE]
  grp      <- as.factor(meta_sub[[group_col]])
  lvls     <- levels(grp)
  if (length(lvls) < 2) stop("Group column must have at least 2 levels.")

  g1 <- which(grp == lvls[1])
  g2 <- which(grp == lvls[2])

  results <- lapply(rownames(tab_sub), function(taxon) {
    x1 <- as.numeric(tab_sub[taxon, g1])
    x2 <- as.numeric(tab_sub[taxon, g2])
    if (sum(x1) == 0 && sum(x2) == 0) return(NULL)
    p  <- tryCatch(wilcox.test(x1, x2)$p.value, error = function(e) NA)
    m1 <- mean(x1 + 1e-10); m2 <- mean(x2 + 1e-10)
    data.frame(taxon    = taxon,
               group1   = lvls[1],
               group2   = lvls[2],
               lda_score = log2(m1 / m2),
               p_value  = p,
               stringsAsFactors = FALSE)
  })

  res_df <- do.call(rbind, Filter(Negate(is.null), results))
  res_df$p_adj <- p.adjust(res_df$p_value, method = "BH")
  res_df       <- res_df[order(res_df$p_adj), ]
  write.table(res_df, file.path(out_dir, "lefse_results.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # LDA score barplot
  sig <- res_df[!is.na(res_df$p_adj) & res_df$p_adj < 0.2, ]
  if (nrow(sig) == 0) {
    message("No significant taxa found (p_adj < 0.2). Lowering threshold to top 20.")
    sig <- head(res_df, 20)
  }
  sig$direction <- ifelse(sig$lda_score > 0, lvls[1], lvls[2])

  pdf(file.path(out_dir, "lefse_plots.pdf"), width = 10, height = max(6, nrow(sig) * 0.35))
  p <- ggplot(sig, aes(x = reorder(taxon, lda_score), y = lda_score, fill = direction)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("steelblue", "tomato")) +
    theme_bw(base_size = 12) +
    labs(title = "LEfSe-style LDA Scores (Log2FC)",
         x = "Taxon", y = "Log2(Fold-Change)", fill = lvls[1])
  print(p)
  dev.off()
  message("Fallback LEfSe complete. Saved to: ", out_dir)
  quit(save = "no", status = 0)
}

# ── microbiomeMarker LEfSe ─────────────────────────────────────────────────────
library(microbiomeMarker)
library(phyloseq)

# Build phyloseq object
common_samples <- intersect(colnames(tab), rownames(meta))
tab_sub  <- tab[, common_samples]
meta_sub <- meta[common_samples, , drop = FALSE]
otu      <- otu_table(as.matrix(tab_sub), taxa_are_rows = TRUE)
sam      <- sample_data(meta_sub)
ps       <- phyloseq(otu, sam)

# Run LEfSe
message("Running LEfSe via microbiomeMarker...")
mm <- tryCatch(
  run_lefse(ps, group = group_col,
             multigrp_strat = TRUE,
             lda_cutoff = 2.0,
             wilcoxon_cutoff = 0.05),
  error = function(e) { message("LEfSe failed: ", e$message); NULL }
)

if (!is.null(mm) && ntaxa(mm) > 0) {
  res_df <- marker_table(mm) %>% as.data.frame()
  write.table(res_df, file.path(out_dir, "lefse_results.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  pdf(file.path(out_dir, "lefse_plots.pdf"), width = 10, height = 8)
  print(plot_ef_bar(mm) +
        theme_bw(base_size = 12) +
        labs(title = paste("LEfSe — Differentially abundant taxa by", group_col)))
  print(plot_cladogram(mm, color = setNames(
    c("steelblue","tomato"), sort(unique(meta_sub[[group_col]]))),
    clade_label_level = 5))
  dev.off()
  message("LEfSe complete. Saved to: ", out_dir)
} else {
  message("No significant taxa found by LEfSe. Writing empty results.")
  write.table(data.frame(), file.path(out_dir, "lefse_results.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  pdf(file.path(out_dir, "lefse_plots.pdf")); dev.off()
}
