#!/usr/bin/env Rscript
# =============================================================================
# composition_plots.R — Robust Stacked barplots
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) stop("Usage: composition_plots.R <phylum_tsv> <class_tsv> <genus_tsv> <metadata> <group_col> <out_pdf>")

phylum_tsv <- args[1]
class_tsv  <- args[2]
genus_tsv  <- args[3]
meta_file  <- args[4]
group_col  <- args[5]
out_pdf    <- args[6]

# ── 1. Load metadata (Defensive) ──────────────────────────────────────────────
meta <- read.table(meta_file, header = TRUE, sep = "\t", 
                   row.names = 1, check.names = FALSE, 
                   comment.char = "", quote = "", fill = TRUE)

# ── 2. Function: load relative frequency table ────────────────────────────────
load_relfreq <- function(tsv_file, level_name) {
  # QIIME 2 exports have a comment line first, then the #OTU ID line
  df <- read.table(tsv_file, header = TRUE, sep = "\t",
                   skip = 1, row.names = 1, check.names = FALSE, 
                   comment.char = "", quote = "")
  
  # Check if we have data
  if(ncol(df) == 0) stop(paste("File is empty or incorrectly formatted:", tsv_file))

  # Transpose: samples as rows, taxa as columns
  df_t <- as.data.frame(t(df))
  df_t$SampleID <- rownames(df_t)
  
  # Pivot to long format for ggplot
  df_long <- pivot_longer(df_t, -SampleID, names_to = "taxon", values_to = "rel_abund")
  
  # Clean up taxon names (optional: removes the d__ p__ prefixes for cleaner plots)
  df_long$taxon <- gsub(".*__[a-z]__", "", df_long$taxon)
  df_long$taxon <- gsub(";__", "", df_long$taxon) # clean empty levels
  
  df_long$level <- level_name
  return(df_long)
}

# ── 3. Function: collapse low-abundance taxa into "Other" ─────────────────────
collapse_other <- function(df, top_n = 12) {
  # Calculate mean abundance across all samples
  top_taxa <- df %>%
    group_by(taxon) %>%
    summarise(mean_abund = mean(rel_abund, na.rm = TRUE)) %>%
    arrange(desc(mean_abund)) %>%
    slice_head(n = top_n) %>%
    pull(taxon)
    
  df %>%
    mutate(taxon = if_else(taxon %in% top_taxa, taxon, "Other")) %>%
    group_by(SampleID, taxon, level) %>%
    summarise(rel_abund = sum(rel_abund), .groups = "drop")
}

# ── 4. Function: stacked barplot ──────────────────────────────────────────────
make_barplot <- function(df, meta, group_col, level_name) {
  # Merge metadata
  df <- inner_join(df, tibble::rownames_to_column(meta, "SampleID"), by = "SampleID")
  
  if(nrow(df) == 0) {
      stop(paste("No matching Sample IDs found between table and metadata in", level_name))
  }

  # Arrange taxons so "Other" is always at the bottom/end
  taxa_levels <- c(setdiff(unique(df$taxon), "Other"), "Other")
  df$taxon <- factor(df$taxon, levels = taxa_levels)

  # Color palette (Set3 + manual extension for many taxa)
  n_taxa <- length(taxa_levels)
  getPalette = colorRampPalette(brewer.pal(min(n_taxa, 12), "Set3"))
  cols <- getPalette(n_taxa)
  names(cols) <- taxa_levels
  cols["Other"] <- "#D3D3D3" # Force 'Other' to be light grey

  ggplot(df, aes(x = SampleID, y = rel_abund, fill = taxon)) +
    geom_bar(stat = "identity", width = 0.9) +
    facet_grid(cols = vars(!!sym(group_col)), scales = "free_x", space = "free_x") +
    scale_fill_manual(values = cols) +
    scale_y_continuous(labels = percent_format(), expand = c(0,0)) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
          panel.spacing = unit(0.2, "lines"),
          strip.background = element_rect(fill = "grey90", color = NA)) +
    labs(title = paste("Composition:", level_name), 
         x = "Samples", y = "Relative Abundance", fill = level_name)
}

# ── 5. Run Execution ──────────────────────────────────────────────────────────
levels_data <- list(
  Phylum = load_relfreq(phylum_tsv, "Phylum"),
  Class  = load_relfreq(class_tsv,  "Class"),
  Genus  = load_relfreq(genus_tsv,  "Genus")
)

pdf(out_pdf, width = 16, height = 8)

for (lvl in names(levels_data)) {
  # Process and Plot
  df_plot <- collapse_other(levels_data[[lvl]], top_n = 15)
  print(make_barplot(df_plot, meta, group_col, lvl))
  
  # Summary Heatmap
  hmap_df <- inner_join(df_plot, tibble::rownames_to_column(meta, "SampleID"), by = "SampleID") %>%
    group_by(!!sym(group_col), taxon) %>%
    summarise(mean_abund = mean(rel_abund), .groups = "drop")

  print(ggplot(hmap_df, aes(x = !!sym(group_col), y = taxon, fill = mean_abund)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "white", high = "darkblue", labels = percent_format()) +
    geom_text(aes(label = percent(mean_abund, accuracy = 0.1)), size = 3) +
    theme_minimal() + labs(title = paste("Heatmap Mean Abundance:", lvl)))
}

dev.off()