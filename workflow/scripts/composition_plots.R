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
  if (!requireNamespace("patchwork", quietly = TRUE))
    install.packages("patchwork", repos = "https://cloud.r-project.org", quiet = TRUE)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) stop("Usage: composition_plots.R <phylum_tsv> <class_tsv> <genus_tsv> <species_tsv> <metadata> <group_col> <out_dir> [strategy] [dual_plots]")

phylum_tsv  <- args[1]
class_tsv   <- args[2]
genus_tsv   <- args[3]
species_tsv <- args[4]
meta_file   <- args[5]
group_col   <- args[6]
out_dir     <- args[7]
strategy    <- if (length(args) >= 8) args[8] else "rename"
dual_plots  <- isTRUE(tolower(if (length(args) >= 9) args[9] else "false") == "true")

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

# ── shared helpers ────────────────────────────────────────────────────────────
taxa_palette <- function(taxa_levels) {
  n <- length(taxa_levels)
  getPalette <- colorRampPalette(brewer.pal(min(n, 12), "Set3"))
  cols <- getPalette(n)
  names(cols) <- taxa_levels
  cols["Other"] <- "#D3D3D3"
  cols
}

merge_meta <- function(df, meta) {
  inner_join(df, tibble::rownames_to_column(meta, "SampleID"), by = "SampleID")
}

# ── 4. Plot 1: Group-averaged stacked bar ─────────────────────────────────────
# Samples within each group are averaged into a single bar per group.
make_barplot <- function(df, meta, group_col, level_name, cleaned = TRUE) {
  if (cleaned && strategy != "none")
    df <- df %>% mutate(taxon = sapply(taxon, format_taxon_label))
  df <- merge_meta(df, meta)
  if (nrow(df) == 0) stop(paste("No matching SampleIDs in", level_name))

  df <- df %>%
    group_by(!!sym(group_col), taxon) %>%
    summarise(rel_abund = mean(rel_abund, na.rm = TRUE), .groups = "drop")

  taxa_levels <- c(setdiff(unique(df$taxon), "Other"), "Other")
  df$taxon <- factor(df$taxon, levels = taxa_levels)
  cols <- taxa_palette(taxa_levels)

  ggplot(df, aes(x = !!sym(group_col), y = rel_abund, fill = taxon)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = cols) +
    scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x        = element_text(angle = 30, hjust = 1, face = "bold", size = 11),
      axis.ticks.x       = element_line(),
      panel.grid.major.x = element_blank()
    ) +
    labs(
      title    = paste("Mean Composition by Group:", level_name),
      subtitle = "Each bar = mean relative abundance across all samples within the group",
      x = group_col, y = "Mean Relative Abundance", fill = level_name
    )
}

# ── 5. Plot 2: Heatmap (log-scaled, taxa ordered by abundance) ─────────────────
make_heatmap <- function(df, meta, group_col, level_name, cleaned = TRUE) {
  if (cleaned && strategy != "none")
    df <- df %>% mutate(taxon = sapply(taxon, format_taxon_label))

  hmap_df <- merge_meta(df, meta) %>%
    group_by(!!sym(group_col), taxon) %>%
    summarise(mean_abund = mean(rel_abund, na.rm = TRUE), .groups = "drop") %>%
    filter(taxon != "Other")

  taxon_order <- hmap_df %>%
    group_by(taxon) %>%
    summarise(tot = mean(mean_abund), .groups = "drop") %>%
    arrange(desc(tot)) %>%
    pull(taxon)
  hmap_df$taxon <- factor(hmap_df$taxon, levels = rev(taxon_order))

  ggplot(hmap_df, aes(x = !!sym(group_col), y = taxon, fill = log1p(mean_abund * 100))) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(mean_abund >= 0.005, percent(mean_abund, accuracy = 0.1), "")),
              size = 3, colour = "black") +
    scale_fill_distiller(palette = "YlOrRd", direction = 1,
                         name = "log1p(Mean %)", labels = function(x) round(x, 1)) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x     = element_text(angle = 30, hjust = 1, face = "bold"),
      panel.grid      = element_blank(),
      legend.position = "right"
    ) +
    labs(
      title    = paste("Mean Relative Abundance Heatmap:", level_name),
      subtitle = "Fill: log1p(mean %) — values \u2265 0.5% labelled",
      x = group_col, y = NULL
    )
}

# ── 6. Plot 3: Bubble plot (size = mean abundance, colour = prevalence) ─────────
make_bubble_plot <- function(df, meta, group_col, level_name, cleaned = TRUE, top_n = 25) {
  if (cleaned && strategy != "none")
    df <- df %>% mutate(taxon = sapply(taxon, format_taxon_label))

  bubble_df <- merge_meta(df, meta) %>%
    group_by(!!sym(group_col), taxon) %>%
    summarise(
      mean_abund = mean(rel_abund, na.rm = TRUE),
      prevalence = mean(rel_abund > 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(taxon != "Other")

  top_taxa <- bubble_df %>%
    group_by(taxon) %>%
    summarise(tot = mean(mean_abund), .groups = "drop") %>%
    slice_max(tot, n = top_n) %>%
    pull(taxon)

  bubble_df <- bubble_df %>%
    filter(taxon %in% top_taxa) %>%
    mutate(taxon = factor(taxon, levels = rev(top_taxa)))

  ggplot(bubble_df, aes(x = !!sym(group_col), y = taxon,
                        size = mean_abund * 100, colour = prevalence)) +
    geom_point(alpha = 0.85) +
    scale_size_area(max_size = 14, name = "Mean Abund. (%)") +
    scale_colour_distiller(palette = "RdYlBu", direction = -1,
                           name = "Prevalence", labels = percent_format(), limits = c(0, 1)) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x      = element_text(angle = 30, hjust = 1, face = "bold"),
      panel.grid.major = element_line(colour = "grey92")
    ) +
    labs(
      title    = paste("Abundance & Prevalence by Group:", level_name),
      subtitle = paste0("Top ", top_n, " taxa — dot size = mean %, colour = fraction of samples detected"),
      x = group_col, y = NULL
    )
}

# ── 7. Plot 4: Ranked abundance curve (Whittaker plot) ────────────────────────
make_ranked_abundance_plot <- function(df, meta, group_col, level_name, cleaned = TRUE) {
  if (cleaned && strategy != "none")
    df <- df %>% mutate(taxon = sapply(taxon, format_taxon_label))

  ranked_df <- merge_meta(df, meta) %>%
    group_by(!!sym(group_col), taxon) %>%
    summarise(mean_abund = mean(rel_abund, na.rm = TRUE), .groups = "drop") %>%
    filter(mean_abund > 0, taxon != "Other") %>%
    group_by(!!sym(group_col)) %>%
    arrange(desc(mean_abund), .by_group = TRUE) %>%
    mutate(rank = row_number()) %>%
    ungroup()

  n_groups   <- length(unique(ranked_df[[group_col]]))
  group_cols <- brewer.pal(max(3, min(n_groups, 12)), "Set1")[seq_len(n_groups)]

  ggplot(ranked_df, aes(x = rank, y = log10(mean_abund),
                        colour = !!sym(group_col), group = !!sym(group_col))) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 1.8, alpha = 0.7) +
    scale_colour_manual(values = group_cols, name = group_col) +
    scale_y_continuous(name = "log\u2081\u2080(Mean Relative Abundance)") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank()) +
    labs(
      title    = paste("Ranked Abundance Curve:", level_name),
      subtitle = "Steeper decline = community dominated by few taxa; flatter = more even",
      x = "Rank (most \u2192 least abundant)"
    )
}

# ── 8. Plot 5: Prevalence\u2013abundance scatter (genus level) ─────────────────────
make_prevalence_scatter <- function(df, meta, group_col, level_name = "Genus", cleaned = TRUE) {
  if (cleaned && strategy != "none")
    df <- df %>% mutate(taxon = sapply(taxon, format_taxon_label))

  prev_df <- merge_meta(df, meta) %>%
    filter(taxon != "Other") %>%
    group_by(taxon) %>%
    summarise(
      prevalence = mean(rel_abund > 0, na.rm = TRUE),
      mean_abund = mean(rel_abund, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(mean_abund > 0)

  top20   <- prev_df %>% slice_max(mean_abund, n = 20) %>% pull(taxon)
  prev_df <- prev_df %>% mutate(label = ifelse(taxon %in% top20, taxon, NA_character_))

  p <- ggplot(prev_df, aes(x = prevalence, y = log10(mean_abund), label = label)) +
    geom_point(aes(colour = prevalence), alpha = 0.75, size = 2.5) +
    geom_vline(xintercept = 0.5,    linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = log10(0.01), linetype = "dashed", colour = "grey50") +
    scale_colour_distiller(palette = "RdYlBu", direction = -1,
                           labels = percent_format(), name = "Prevalence") +
    scale_x_continuous(labels = percent_format(), limits = c(0, 1)) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank()) +
    labs(
      title    = paste("Prevalence vs. Abundance:", level_name),
      subtitle = "Dashed: 50% prevalence / 1% mean abundance \u2014 top 20 taxa labelled",
      x = "Prevalence (fraction of samples)",
      y = "log\u2081\u2080(Mean Relative Abundance)"
    )

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(na.rm = TRUE, size = 3, max.overlaps = 15)
  } else {
    p <- p + geom_text(na.rm = TRUE, size = 2.8, hjust = -0.1)
  }
  p
}

# =============================================================================
# Execution
# =============================================================================
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

levels_data <- list(
  Phylum  = load_relfreq(phylum_tsv,  "Phylum"),
  Class   = load_relfreq(class_tsv,   "Class"),
  Genus   = load_relfreq(genus_tsv,   "Genus"),
  Species = load_relfreq(species_tsv, "Species")
)

save_png <- function(p, filename, w = 12, h = 8) {
  ggsave(file.path(out_dir, filename), plot = p,
         width = w, height = h, units = "in", dpi = 600)
  message("Saved: ", filename, " (600 DPI)")
}

# ── Combined PDFs (all plot types, all levels) ────────────────────────────────
write_composition_pdf <- function(out_pdf, cleaned) {
  pdf(out_pdf, width = 14, height = 9)
  for (lvl in names(levels_data)) {
    df_bar  <- collapse_other(levels_data[[lvl]], top_n = 15)
    df_full <- collapse_other(levels_data[[lvl]], top_n = 50)
    print(make_barplot(df_bar,  meta, group_col, lvl, cleaned = cleaned))
    print(make_heatmap(df_bar,  meta, group_col, lvl, cleaned = cleaned))
    print(make_bubble_plot(df_full, meta, group_col, lvl, cleaned = cleaned))
    print(make_ranked_abundance_plot(df_full, meta, group_col, lvl, cleaned = cleaned))
  }
  df_genus_full   <- collapse_other(levels_data[["Genus"]],   top_n = 50)
  df_species_full <- collapse_other(levels_data[["Species"]], top_n = 50)
  print(make_prevalence_scatter(df_genus_full,   meta, group_col, "Genus",   cleaned = cleaned))
  print(make_prevalence_scatter(df_species_full, meta, group_col, "Species", cleaned = cleaned))
  dev.off()
}

write_composition_pdf(file.path(out_dir, "composition_plots_raw.pdf"), cleaned = FALSE)
write_composition_pdf(file.path(out_dir, "composition_plots.pdf"),     cleaned = TRUE)

# ── Multipanel 600 DPI PNGs — one per plot type, all levels in a 2×2 grid ────

# Helper: keep only the level name as the panel title
panel_ready <- function(p, lvl) {
  p + labs(title = lvl, subtitle = NULL) +
    theme(plot.title    = element_text(face = "bold", size = 12, hjust = 0.5),
          plot.subtitle = element_blank())
}

save_panel_png <- function(plots, title, filename, w = 22, h = 16) {
  panel <- patchwork::wrap_plots(plots, ncol = 2) +
    patchwork::plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5))
    )
  ggsave(file.path(out_dir, filename), plot = panel,
         width = w, height = h, units = "in", dpi = 600)
  message("Saved: ", filename, " (600 DPI)")
}

# Barplot panel (2×2)
bar_plots <- lapply(names(levels_data), function(lvl) {
  df <- collapse_other(levels_data[[lvl]], top_n = 15)
  panel_ready(make_barplot(df, meta, group_col, lvl), lvl)
})
save_panel_png(bar_plots, "Mean Taxonomic Composition by Group",
               "composition_barplot_panel.png", w = 22, h = 16)

# Heatmap panel (2×2)
heat_plots <- lapply(names(levels_data), function(lvl) {
  df <- collapse_other(levels_data[[lvl]], top_n = 15)
  panel_ready(make_heatmap(df, meta, group_col, lvl), lvl)
})
save_panel_png(heat_plots, "Mean Relative Abundance Heatmaps",
               "composition_heatmap_panel.png", w = 22, h = 20)

# Bubble panel (2×2)
bubble_plots <- lapply(names(levels_data), function(lvl) {
  df <- collapse_other(levels_data[[lvl]], top_n = 50)
  panel_ready(make_bubble_plot(df, meta, group_col, lvl, top_n = 25), lvl)
})
save_panel_png(bubble_plots, "Abundance & Prevalence Bubble Plots",
               "composition_bubble_panel.png", w = 22, h = 22)

# Ranked abundance panel (2×2)
ranked_plots <- lapply(names(levels_data), function(lvl) {
  df <- collapse_other(levels_data[[lvl]], top_n = 50)
  panel_ready(make_ranked_abundance_plot(df, meta, group_col, lvl), lvl)
})
save_panel_png(ranked_plots, "Ranked Abundance Curves",
               "composition_ranked_panel.png", w = 22, h = 14)

# Prevalence scatter panel (1×2 — Genus + Species)
prev_plots <- lapply(c("Genus", "Species"), function(lvl) {
  df <- collapse_other(levels_data[[lvl]], top_n = 50)
  panel_ready(make_prevalence_scatter(df, meta, group_col, lvl), lvl)
})
save_panel_png(prev_plots, "Prevalence vs. Abundance — Genus & Species",
               "composition_prevalence_panel.png", w = 22, h = 10)

message("All composition plots saved: 2 PDFs + 5 multipanel PNGs (600 DPI) in ", out_dir)