#!/usr/bin/env Rscript
# =============================================================================
# core_microbiome.R — Core microbiome chi-square / Fisher exact test
# =============================================================================
# Usage: Rscript core_microbiome.R <table_tsv> <metadata_tsv>
#                                  <group_col> <prevalence_threshold> <out_dir>
#
# Required: ggplot2, dplyr, tidyr
# Optional: UpSetR (for UpSet plots)
# Install:  mamba install -n qiime2 -c conda-forge r-ggplot2 r-dplyr r-tidyr r-upsetr
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# ── Arguments ─────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("Usage: core_microbiome.R <table_tsv> <metadata> <group_col> <prevalence> <out_dir>")

table_tsv <- args[1]
meta_file <- args[2]
group_col <- args[3]
prevalence <- as.numeric(args[4])
out_dir   <- args[5]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load data ─────────────────────────────────────────────────────────────────
tab <- read.table(table_tsv, header = TRUE, sep = "\t",
                  skip = 1, row.names = 1, check.names = FALSE)

meta <- read.table(meta_file, header = TRUE, sep = "\t",
                   comment.char = "", row.names = 1, check.names = FALSE)

common_samples <- intersect(colnames(tab), rownames(meta))
if (length(common_samples) < 2) stop("Fewer than 2 common samples.")
tab  <- tab[, common_samples, drop = FALSE]
meta <- meta[common_samples, , drop = FALSE]

grp <- as.factor(meta[[group_col]])
grp_levels <- levels(grp)

# ── Define core taxa per group: present in ≥prevalence fraction ──────────────
core_taxa <- list()
for (g in grp_levels) {
  samp_g    <- common_samples[grp == g]
  tab_g     <- tab[, samp_g, drop = FALSE]
  prev_g    <- rowSums(tab_g > 0) / ncol(tab_g)
  core_taxa[[g]] <- names(prev_g[prev_g >= prevalence])
  message("Group ", g, ": ", length(core_taxa[[g]]), " core taxa at prevalence ≥", prevalence)
}

# ── Chi-square / Fisher exact test: each taxon core vs non-core ────────────
all_taxa <- unique(unlist(core_taxa))
if (length(grp_levels) == 2) {
  g1_core <- core_taxa[[grp_levels[1]]]
  g2_core <- core_taxa[[grp_levels[2]]]

  stats_rows <- lapply(all_taxa, function(tx) {
    in_g1 <- tx %in% g1_core
    in_g2 <- tx %in% g2_core
    if (!in_g1 && !in_g2) return(NULL)
    contingency <- matrix(c(
      sum(grp == grp_levels[1] & as.logical(tab[tx, grp == grp_levels[1]] > 0)),
      sum(grp == grp_levels[1] & as.logical(tab[tx, grp == grp_levels[1]] == 0)),
      sum(grp == grp_levels[2] & as.logical(tab[tx, grp == grp_levels[2]] > 0)),
      sum(grp == grp_levels[2] & as.logical(tab[tx, grp == grp_levels[2]] == 0))
    ), nrow = 2, dimnames = list(
      c("Present", "Absent"),
      c(grp_levels[1], grp_levels[2])
    ))
    p <- tryCatch(fisher.test(contingency)$p.value, error = function(e) NA)
    data.frame(taxon = tx,
               core_group1 = in_g1,
               core_group2 = in_g2,
               p_value = p,
               stringsAsFactors = FALSE)
  })

  stats_df <- do.call(rbind, Filter(Negate(is.null), stats_rows))
  if (!is.null(stats_df) && nrow(stats_df) > 0) {
    stats_df$p_adj <- p.adjust(stats_df$p_value, method = "BH")
    stats_df <- stats_df[order(stats_df$p_adj), ]
    colnames(stats_df)[2:3] <- paste0("core_", grp_levels)
    write.table(stats_df, file.path(out_dir, "core_stats.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    message("Saved: core_stats.tsv (", nrow(stats_df), " taxa)")
  }
}

# ── Plots ─────────────────────────────────────────────────────────────────────
pdf(file.path(out_dir, "core_plots.pdf"), width = 10, height = 8)

# UpSet plot if available
if (requireNamespace("UpSetR", quietly = TRUE)) {
  library(UpSetR)
  upset_input <- setNames(
    lapply(grp_levels, function(g) which(all_taxa %in% core_taxa[[g]])),
    grp_levels
  )
  UpSetR::upset(fromList(upset_input),
                nsets = length(grp_levels),
                order.by = "freq",
                main.bar.color = "#2196F3",
                sets.bar.color = "#4CAF50",
                text.scale = 1.3)
  title(paste("Core Microbiome Overlap (prevalence ≥", prevalence, ")"))
} else {
  # Fallback: Venn-style summary barplot
  union_taxa <- unique(c(unlist(core_taxa)))
  counts_df <- data.frame(
    group = grp_levels,
    n_core = sapply(grp_levels, function(g) length(core_taxa[[g]]))
  )
  p <- ggplot(counts_df, aes(x = group, y = n_core, fill = group)) +
    geom_bar(stat = "identity", width = 0.5) +
    geom_text(aes(label = n_core), vjust = -0.4, size = 5) +
    theme_bw(base_size = 13) +
    labs(title = paste("Core taxa count per group (prevalence ≥", prevalence, ")"),
         x = group_col, y = "Number of core taxa") +
    theme(legend.position = "none")
  print(p)
}

dev.off()
message("Saved: core_plots.pdf")
