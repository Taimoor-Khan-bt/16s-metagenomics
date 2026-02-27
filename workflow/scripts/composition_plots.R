#!/usr/bin/env Rscript
# =============================================================================
# composition_plots.R — Stacked barplots at phylum / class / genus level
# =============================================================================
# Usage: Rscript composition_plots.R <phylum_tsv> <class_tsv> <genus_tsv>
#                                    <metadata_tsv> <group_col> <out_pdf>
#
# Required: ggplot2, dplyr, tidyr, RColorBrewer, scales
# Install:  mamba install -n qiime2 -c conda-forge r-ggplot2 r-dplyr r-tidyr r-rcolorbrewer r-scales
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(scales)
})

# ── Arguments ─────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) stop("Usage: composition_plots.R <phylum_tsv> <class_tsv> <genus_tsv> <metadata> <group_col> <out_pdf>")

phylum_tsv <- args[1]
class_tsv  <- args[2]
genus_tsv  <- args[3]
meta_file  <- args[4]
group_col  <- args[5]
out_pdf    <- args[6]

# ── Load metadata ─────────────────────────────────────────────────────────────
meta <- read.table(meta_file, header = TRUE, sep = "\t",
                   comment.char = "", row.names = 1, check.names = FALSE)

# ── Function: load relative frequency table ───────────────────────────────────
load_relfreq <- function(tsv_file, level_name) {
  df <- read.table(tsv_file, header = TRUE, sep = "\t",
                   skip = 1, row.names = 1, check.names = FALSE)
  # Transpose: samples as rows, taxa as columns
  df_t <- as.data.frame(t(df))
  df_t$SampleID <- rownames(df_t)
  df_long <- pivot_longer(df_t, -SampleID, names_to = "taxon", values_to = "rel_abund")
  df_long$level <- level_name
  df_long
}

# ── Function: collapse low-abundance taxa into "Other" ───────────────────────
collapse_other <- function(df, top_n = 15) {
  top_taxa <- df %>%
    group_by(taxon) %>%
    summarise(mean_abund = mean(rel_abund)) %>%
    arrange(desc(mean_abund)) %>%
    slice_head(n = top_n) %>%
    pull(taxon)
  df %>%
    mutate(taxon = if_else(taxon %in% top_taxa, taxon, "Other")) %>%
    group_by(SampleID, taxon, level) %>%
    summarise(rel_abund = sum(rel_abund), .groups = "drop")
}

# ── Function: stacked barplot ─────────────────────────────────────────────────
make_barplot <- function(df, meta, group_col, level_name, palette_n = 16) {
  # Merge group info
  meta_sub <- meta[, group_col, drop = FALSE]
  meta_sub$SampleID <- rownames(meta_sub)
  df <- left_join(df, meta_sub, by = "SampleID")

  # Color palette
  n_taxa <- length(unique(df$taxon))
  cols <- if (n_taxa <= 12) brewer.pal(max(3, n_taxa), "Set3") else
          colorRampPalette(brewer.pal(12, "Set3"))(n_taxa)
  # Put "Other" last (grey)
  cols <- c(cols[seq_len(n_taxa - 1)], "grey70")
  names(cols) <- c(setdiff(unique(df$taxon), "Other"), "Other")

  ggplot(df, aes(x = SampleID, y = rel_abund, fill = taxon)) +
    geom_bar(stat = "identity", width = 0.85) +
    facet_grid(~ !!sym(group_col), scales = "free_x", space = "free_x") +
    scale_fill_manual(values = cols) +
    scale_y_continuous(labels = percent_format()) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
          legend.key.size = unit(0.4, "cm"),
          legend.text = element_text(size = 8),
          strip.background = element_rect(fill = "#E3F2FD")) +
    labs(title = paste("Relative Abundance —", level_name, "by", group_col),
         x = "Sample", y = "Relative abundance", fill = level_name)
}

# ── Load and process all levels ───────────────────────────────────────────────
levels_data <- list(
  Phylum = load_relfreq(phylum_tsv, "Phylum"),
  Class  = load_relfreq(class_tsv,  "Class"),
  Genus  = load_relfreq(genus_tsv,  "Genus")
)

# ── Write PDF ─────────────────────────────────────────────────────────────────
pdf(out_pdf, width = 14, height = 7)

for (lvl_name in names(levels_data)) {
  df <- collapse_other(levels_data[[lvl_name]], top_n = 15)
  p  <- make_barplot(df, meta, group_col, lvl_name)
  print(p)
  message("Plotted: ", lvl_name)
}

# ── Mean relative abundance heatmap summary per group ────────────────────────
for (lvl_name in names(levels_data)) {
  df <- collapse_other(levels_data[[lvl_name]], top_n = 15)
  meta_sub <- meta[, group_col, drop = FALSE]
  meta_sub$SampleID <- rownames(meta_sub)
  df <- left_join(df, meta_sub, by = "SampleID")

  hmap_df <- df %>%
    group_by(!!sym(group_col), taxon) %>%
    summarise(mean_abund = mean(rel_abund), .groups = "drop")

  p_heat <- ggplot(hmap_df, aes_string(x = group_col, y = "taxon", fill = "mean_abund")) +
    geom_tile(color = "white") +
    geom_text(aes(label = scales::percent(mean_abund, accuracy = 0.1)), size = 3) +
    scale_fill_gradient(low = "#EFF3FF", high = "#2171B5", labels = percent_format()) +
    theme_bw(base_size = 11) +
    labs(title = paste("Mean Relative Abundance —", lvl_name),
         x = group_col, y = lvl_name, fill = "Mean rel. abund.")
  print(p_heat)
}

dev.off()
message("Saved: ", out_pdf)
