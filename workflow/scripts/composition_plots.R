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
if (length(args) < 6) stop("Usage: composition_plots.R <phylum_tsv> <class_tsv> <genus_tsv> <metadata> <group_col> <out_dir> [strategy] [dual_plots]")

phylum_tsv <- args[1]
class_tsv  <- args[2]
genus_tsv  <- args[3]
meta_file  <- args[4]
group_col  <- args[5]
out_dir    <- args[6]
strategy   <- if (length(args) >= 7) args[7] else "rename"
dual_plots <- isTRUE(tolower(if (length(args) >= 8) args[8] else "false") == "true")

# =============================================================================
# format_taxon_label()
# Converts raw SILVA taxonomy strings to publication-ready labels.
#   In:  "k__Bacteria;p__Firmicutes;...;f__Lachnospiraceae;g__"
#   Out: "Unclassified Lachnospiraceae"
# =============================================================================
format_taxon_label <- function(x) {
  if (is.na(x) || nchar(trimws(x)) == 0 || x == "Other") return(x)
  parts <- strsplit(gsub("\\|", ";", x), ";")[[1]]
  parts <- trimws(parts)
  clean <- sub("^[dpcofgskDPCOFGSK]__", "", parts)
  non_empty <- clean[nchar(clean) > 0 & tolower(clean) != "uncultured"]
  if (length(non_empty) == 0) return("Unclassified")
  last_raw_clean <- sub("^[dpcofgskDPCOFGSK]__", "", parts[length(parts)])
  if (nchar(trimws(last_raw_clean)) == 0) {
    return(paste("Unclassified", non_empty[length(non_empty)]))
  }
  return(non_empty[length(non_empty)])
}

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
make_barplot <- function(df, meta, group_col, level_name, cleaned = TRUE) {
  # Apply taxonomic label cleaning when requested (preserves "Other" as-is)
  if (cleaned && strategy != "none") {
    df <- df %>% mutate(taxon = sapply(taxon, format_taxon_label))
  }
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
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

levels_data <- list(
  Phylum = load_relfreq(phylum_tsv, "Phylum"),
  Class  = load_relfreq(class_tsv,  "Class"),
  Genus  = load_relfreq(genus_tsv,  "Genus")
)

make_heatmap <- function(df_plot, meta, group_col, level_name) {
  hmap_df <- inner_join(df_plot, tibble::rownames_to_column(meta, "SampleID"), by = "SampleID") %>%
    group_by(!!sym(group_col), taxon) %>%
    summarise(mean_abund = mean(rel_abund), .groups = "drop")
  ggplot(hmap_df, aes(x = !!sym(group_col), y = taxon, fill = mean_abund)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "white", high = "darkblue", labels = percent_format()) +
    geom_text(aes(label = percent(mean_abund, accuracy = 0.1)), size = 3) +
    theme_minimal() +
    labs(title = paste("Heatmap Mean Abundance:", level_name))
}

write_composition_pdf <- function(out_pdf, cleaned) {
  pdf(out_pdf, width = 16, height = 8)
  for (lvl in names(levels_data)) {
    df_plot <- collapse_other(levels_data[[lvl]], top_n = 15)
    print(make_barplot(df_plot, meta, group_col, lvl, cleaned = cleaned))
    print(make_heatmap(df_plot, meta, group_col, lvl))
  }
  dev.off()
}

# Raw PDF (SILVA strings intact — no prefix stripping)
write_composition_pdf(file.path(out_dir, "composition_plots_raw.pdf"), cleaned = FALSE)

# Cleaned PDF (format_taxon_label applied)
write_composition_pdf(file.path(out_dir, "composition_plots.pdf"), cleaned = TRUE)

# Save cleaned barplots as 600 DPI PNG (one per taxonomy level)
for (lvl in names(levels_data)) {
  df_plot <- collapse_other(levels_data[[lvl]], top_n = 15)
  p_png   <- make_barplot(df_plot, meta, group_col, lvl, cleaned = TRUE)
  ggsave(file.path(out_dir, paste0("composition_", tolower(lvl), ".png")),
         plot = p_png, width = 16, height = 8, units = "in", dpi = 600)
  message("Saved: composition_", tolower(lvl), ".png (600 DPI)")
}
message("Saved: composition_plots.pdf, composition_plots_raw.pdf, and per-level PNGs (600 DPI)")