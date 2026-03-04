#!/usr/bin/env Rscript
# =============================================================================
# core_microbiome.R — Core microbiome analysis (Multi-Group Support)
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: core_microbiome.R <table_tsv> <taxonomy_tsv> <metadata> <group_col> <prevalence> <out_dir>")
}

table_tsv  <- args[1]
tax_tsv    <- args[2]
meta_file  <- args[3]
group_col  <- args[4]
prevalence <- as.numeric(args[5])
out_dir    <- args[6]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── 1. Load Data ─────────────────────────────────────────────────────────────
tab <- read.table(table_tsv, header = TRUE, sep = "\t", skip = 1, 
                  row.names = 1, check.names = FALSE, comment.char = "", quote = "")

tax_df <- read.table(tax_tsv, header = TRUE, sep = "\t", 
                     row.names = 1, check.names = FALSE, quote = "")

meta <- read.table(meta_file, header = TRUE, sep = "\t", 
                   row.names = 1, check.names = FALSE, comment.char = "", quote = "")

# ── 2. Sync Samples ──────────────────────────────────────────────────────────
common_samples <- intersect(colnames(tab), rownames(meta))
if (length(common_samples) < 2) stop("Fewer than 2 common samples.")

tab  <- tab[, common_samples, drop = FALSE]
meta <- meta[common_samples, , drop = FALSE]
grp  <- as.factor(meta[[group_col]])
grp_levels <- levels(grp)

# ── 3. Identify Core Taxa per Group ──────────────────────────────────────────
core_ids <- list()
for (g in grp_levels) {
  samp_g    <- common_samples[grp == g]
  tab_g     <- tab[, samp_g, drop = FALSE]
  prev_g    <- rowSums(tab_g > 0) / ncol(tab_g)
  core_ids[[g]] <- names(prev_g[prev_g >= prevalence])
  message("Group ", g, ": ", length(core_ids[[g]]), " core features.")
}

all_core_ids <- unique(unlist(core_ids))

# ── 4. Helper: Map Taxonomy ──────────────────────────────────────────────────
get_tax_name <- function(id) {
  if (id %in% rownames(tax_df)) {
    return(gsub("[a-z]__", "", tax_df[id, "Taxon"]))
  }
  return(id)
}

# ── 5. Multi-Group Statistical Test (Omnibus) ────────────────────────────────
if (length(all_core_ids) > 0) {
  stats_rows <- lapply(all_core_ids, function(tx_id) {
    # Create 2xN Contingency Table (Present/Absent vs Groups)
    count_matrix <- matrix(NA, nrow = 2, ncol = length(grp_levels))
    for (i in seq_along(grp_levels)) {
      g <- grp_levels[i]
      count_matrix[1, i] <- sum(tab[tx_id, grp == g] > 0)  # Present
      count_matrix[2, i] <- sum(tab[tx_id, grp == g] == 0) # Absent
    }
    
    # Use Chi-square for the omnibus difference across all groups
    p_val <- tryCatch(chisq.test(count_matrix)$p.value, error = function(e) NA)
    
    res <- data.frame(FeatureID = tx_id, Taxon = get_tax_name(tx_id), p_value = p_val)
    # Add Boolean flags for each group
    for (g in grp_levels) { res[[paste0("is_core_", g)]] <- tx_id %in% core_ids[[g]] }
    return(res)
  })
  
  stats_df <- do.call(rbind, stats_rows)
  stats_df$p_adj <- p.adjust(stats_df$p_value, method = "BH")
  write.table(stats_df, file.path(out_dir, "core_stats.tsv"), 
              sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  # Write an empty file to satisfy Snakemake if no core taxa exist
  write.table(data.frame(Note="No core taxa found"), 
              file.path(out_dir, "core_stats.tsv"), sep="\t")
}

# ── 6. Plotting ─────────────────────────────────────────────────────────────
pdf(file.path(out_dir, "core_plots.pdf"), width = 12, height = 8)
if (requireNamespace("UpSetR", quietly = TRUE) && length(grp_levels) > 1) {
  upset_list <- lapply(core_ids, function(ids) if(length(ids)>0) ids else character(0))
  # Ensure we have data for the plot
  if(sum(sapply(upset_list, length)) > 0) {
     print(UpSetR::upset(UpSetR::fromList(upset_list), order.by = "freq", 
                         main.bar.color = "steelblue", sets.bar.color = "darkgreen"))
  }
}
dev.off()