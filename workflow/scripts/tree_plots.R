#!/usr/bin/env Rscript
# =============================================================================
# tree_plots.R -- Publication-quality phylogenetic tree visualizations
# =============================================================================
# Usage:
#   Rscript tree_plots.R <tree_nwk> <taxonomy_tsv> <table_tsv> \
#                        <metadata_tsv> <group_col> <lefse_tsv> <out_dir>
#
# Readability strategy:
#   The raw tree can have 2000+ ASV tips -- completely unreadable.
#   This script picks ONE representative ASV per genus (highest mean abundance),
#   then caps at MAX_TIPS total. Only the most abundant genera are shown.
# =============================================================================

# -- 0. Packages ---------------------------------------------------------------
suppressPackageStartupMessages({
  for (pkg in c("ape", "treeio", "ggtree", "ggtreeExtra",
                "ggplot2", "dplyr", "tidyr", "tibble", "scales", "RColorBrewer")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (pkg %in% c("treeio", "ggtree", "ggtreeExtra")) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager", repos = "https://cloud.r-project.org")
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
      } else {
        install.packages(pkg, repos = "https://cloud.r-project.org")
      }
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
})

# -- 1. Args -------------------------------------------------------------------
args       <- commandArgs(trailingOnly = TRUE)
tree_file  <- args[1]
tax_file   <- args[2]
feat_file  <- args[3]
meta_file  <- args[4]
group_col  <- args[5]
lefse_file <- args[6]
out_dir    <- args[7]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat("[Tree] Output directory:", out_dir, "\n")

# -- 2. Tunable parameters -----------------------------------------------------
MAX_TIPS <- 80   # max representative tips shown in every plot

# -- 3. Helpers ----------------------------------------------------------------
save_plot <- function(p, path, w = 14, h = 20) {
  tryCatch(
    ggsave(path, plot = p, width = w, height = h,
           units = "in", dpi = 300, limitsize = FALSE),
    error = function(e)
      message("[Tree] WARN: could not save ", basename(path), ": ", e$message)
  )
  png_path <- sub("\\.pdf$", ".png", path)
  tryCatch(
    ggsave(png_path, plot = p, width = w, height = h,
           units = "in", dpi = 600, limitsize = FALSE),
    error = function(e)
      message("[Tree] WARN: could not save ", basename(png_path), ": ", e$message)
  )
  cat("[Tree] Saved:", path, "\n")
}

strip_prefix <- function(x) gsub("^[dpcofgs]__", "", trimws(as.character(x)))

# Return the deepest resolved rank name from a SILVA lineage string
format_taxon_label <- function(taxon_str) {
  if (is.na(taxon_str) || nchar(trimws(taxon_str)) == 0) return("Unclassified")
  parts  <- trimws(strsplit(taxon_str, ";")[[1]])
  parts  <- parts[nchar(parts) > 0]
  filled <- parts[!grepl("^[dpcofgs]__$", parts)]
  if (length(filled) == 0) return("Unclassified")
  clean <- strip_prefix(filled[length(filled)])
  if (nchar(clean) == 0) return("Unclassified")
  clean
}

parse_rank <- function(taxon_str, prefix) {
  parts <- trimws(strsplit(as.character(taxon_str), ";")[[1]])
  hit   <- parts[grepl(paste0("^", prefix, ".+"), parts)]
  if (length(hit) > 0) strip_prefix(hit[1]) else NA_character_
}

# -- 4. Load data --------------------------------------------------------------
cat("[Tree] Loading tree ...\n")
tree <- ape::read.tree(tree_file)
cat("[Tree] ", length(tree$tip.label), "tips\n")

cat("[Tree] Loading taxonomy ...\n")
tax_raw <- read.delim(tax_file, sep = "\t", header = TRUE,
                      comment.char = "", stringsAsFactors = FALSE)
colnames(tax_raw)[1] <- sub("^#", "", colnames(tax_raw)[1])
tax_raw <- tax_raw[!grepl("^#q2:types", tax_raw[[1]]), , drop = FALSE]
tax_df  <- tax_raw %>%
  rename(feature_id = 1, Taxon = 2) %>%
  mutate(
    phylum    = sapply(Taxon, parse_rank, prefix = "p__"),
    genus     = sapply(Taxon, parse_rank, prefix = "g__"),
    tax_label = sapply(Taxon, format_taxon_label),
    phylum    = ifelse(is.na(phylum) | phylum == "", "Other", phylum),
    genus     = ifelse(is.na(genus)  | genus  == "", phylum, genus)
  )

cat("[Tree] Loading feature table ...\n")
feat_raw  <- read.delim(feat_file, sep = "\t", header = TRUE,
                        skip = 1,          # skip "# Constructed from biom file"
                        comment.char = "", stringsAsFactors = FALSE,
                        check.names = FALSE)
colnames(feat_raw)[1] <- sub("^#", "", colnames(feat_raw)[1])
feat_raw  <- feat_raw[!grepl("^#q2:types", feat_raw[[1]]), , drop = FALSE]
feat_mat  <- as.data.frame(
  lapply(feat_raw[, -1, drop = FALSE], as.numeric),
  row.names = feat_raw[[1]], check.names = FALSE
)
col_tot   <- colSums(feat_mat, na.rm = TRUE)
col_tot[col_tot == 0] <- 1
rel_abund   <- sweep(feat_mat, 2, col_tot, "/")
mean_abund  <- setNames(rowMeans(rel_abund, na.rm = TRUE), rownames(rel_abund))

cat("[Tree] Loading metadata ...\n")
meta <- read.delim(meta_file, sep = "\t", header = TRUE,
                   comment.char = "", stringsAsFactors = FALSE)
colnames(meta)[1] <- sub("^#", "", colnames(meta)[1])
meta <- meta[!grepl("^#q2:types", meta[[1]]), , drop = FALSE]
if (!group_col %in% colnames(meta)) {
  cat("[Tree] WARN: group_col '", group_col, "' not found in metadata.\n", sep = "")
  meta$group__ <- "All"; group_col <- "group__"
}

cat("[Tree] Loading LEfSe results ...\n")
lefse_ok <- file.exists(lefse_file) && file.size(lefse_file) > 0
lefse_df <- NULL
if (lefse_ok) {
  lefse_df <- tryCatch(
    read.delim(lefse_file, sep = "\t", header = TRUE,
               stringsAsFactors = FALSE, comment.char = "") %>%
      filter(!is.na(lda_score), lda_score != "") %>%
      mutate(lda_score = as.numeric(lda_score),
             direction = as.character(direction)),
    error = function(e) NULL
  )
  lefse_ok <- !is.null(lefse_df) && nrow(lefse_df) > 0
}

# -- 5. Select representative tips ---------------------------------------------
# Strategy: 1 tip per genus (highest mean abundance) -> sort -> take top MAX_TIPS
cat("[Tree] Selecting representative tips (max", MAX_TIPS, ") ...\n")

tip_base <- data.frame(tip = tree$tip.label, stringsAsFactors = FALSE) %>%
  left_join(tax_df %>% select(feature_id, Taxon, genus, phylum, tax_label),
            by = c("tip" = "feature_id")) %>%
  mutate(
    mean_abund = mean_abund[tip],
    mean_abund = ifelse(is.na(mean_abund), 0, mean_abund),
    phylum     = ifelse(is.na(phylum), "Other", phylum),
    genus      = ifelse(is.na(genus) | genus == "", phylum, genus),
    tax_label  = ifelse(is.na(tax_label) | tax_label == "", "Unclassified", tax_label)
  ) %>%
  group_by(genus) %>%
  slice_max(mean_abund, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(mean_abund)) %>%
  head(MAX_TIPS)

keep_tips <- tip_base$tip
cat("[Tree]  Retained", length(keep_tips), "representative tips\n")

tree_p   <- ape::keep.tip(tree, keep_tips)
tip_data <- tip_base

# LEfSe annotation: match by genus name
if (lefse_ok) {
  lda_best <- lefse_df %>%
    group_by(taxon) %>%
    slice_max(abs(lda_score), n = 1, with_ties = FALSE) %>%
    ungroup()
  tip_data <- tip_data %>%
    left_join(lda_best %>% select(taxon, lda_score, direction),
              by = c("genus" = "taxon")) %>%
    mutate(lda_score = as.numeric(lda_score),
           direction = as.character(direction))
} else {
  tip_data <- tip_data %>%
    mutate(lda_score = NA_real_, direction = NA_character_)
}

n_kept <- length(keep_tips)

# -- 6. Colour palettes --------------------------------------------------------
all_phyla <- sort(unique(tip_data$phylum[!is.na(tip_data$phylum)]))
n_phyla   <- length(all_phyla)
base_pal  <- if (n_phyla <= 8) {
  RColorBrewer::brewer.pal(max(3, n_phyla), "Set2")[seq_len(n_phyla)]
} else if (n_phyla <= 12) {
  RColorBrewer::brewer.pal(12, "Set3")[seq_len(n_phyla)]
} else {
  colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_phyla)
}
phy_pal          <- setNames(base_pal, all_phyla)
phy_pal["Other"] <- "#CCCCCC"

base_theme <- theme_tree2() +
  theme(
    legend.position   = "right",
    legend.text       = element_text(size = 7),
    legend.title      = element_text(size = 8, face = "bold"),
    legend.key.size   = unit(0.4, "cm"),
    plot.title        = element_text(size = 11, face = "bold", hjust = 0.5),
    plot.margin       = margin(5, 130, 5, 5)
  )

# Shared data frames needed by multiple plots (must come before plot blocks)
strip_df <- tip_data %>%
  select(tip, fr_phylum = phylum) %>%
  mutate(fr_phylum = ifelse(is.na(fr_phylum), "Other", fr_phylum))

diff_df <- tip_data %>%
  mutate(
    fr_lda = as.numeric(lda_score),
    fr_dir = as.character(case_when(
      is.na(lda_score) ~ "Not tested",
      lda_score > 0    ~ "Enriched",
      lda_score < 0    ~ "Depleted",
      TRUE             ~ "Not tested"
    ))
  ) %>%
  select(tip, fr_lda, fr_dir)

lda_present <- lefse_ok && any(!is.na(diff_df$fr_lda))

# -- 7. Plot 01: basic rectangular cladogram ----------------------------------
cat("[Tree] Building plot 01: basic cladogram ...\n")

p1 <- ggtree(tree_p, layout = "rectangular", size = 0.35, color = "grey40") %<+% tip_data +
  geom_tippoint(aes(size = mean_abund, color = phylum), alpha = 0.85) +
  geom_tiplab(aes(label = tax_label), size = 2.2, hjust = -0.08) +
  scale_color_manual(values = phy_pal, na.value = "#AAAAAA", name = "Phylum") +
  scale_size_continuous(range = c(1.5, 5), name = "Mean rel. abund.",
                        labels = scales::percent_format(accuracy = 0.1)) +
  base_theme +
  labs(title = "Phylogenetic Tree - Top genera (size = mean abundance)")

save_plot(p1, file.path(out_dir, "01_tree_basic.pdf"), w = 14, h = max(10, n_kept * 0.22))

# -- 8. Plot 02: abundance heatmap ring ----------------------------------------
cat("[Tree] Building plot 02: abundance heatmap ...\n")

heat_cols <- colnames(rel_abund)
if (length(heat_cols) > 24) heat_cols <- heat_cols[round(seq(1, length(heat_cols), length.out = 24))]

heat_df <- rel_abund[keep_tips, heat_cols, drop = FALSE] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("tip") %>%
  tidyr::pivot_longer(-tip, names_to = "fr_sample", values_to = "fr_abund") %>%
  mutate(fr_abund = ifelse(is.na(fr_abund), 0, fr_abund))

p2 <- ggtree(tree_p, layout = "rectangular", size = 0.3, color = "grey40") %<+% tip_data +
  geom_tiplab(aes(label = tax_label), size = 2.0, hjust = -0.05) +
  geom_fruit(
    data    = heat_df,
    geom    = geom_tile,
    mapping = aes(y = tip, x = fr_sample, fill = fr_abund),
    offset  = 0.05, pwidth = 0.5
  ) +
  scale_fill_gradientn(
    colours  = c("white", "#FFF5EB", "#FD8D3C", "#D94701"),
    name     = "Rel. abund.",
    labels   = scales::percent_format(accuracy = 0.01),
    na.value = "white"
  ) +
  theme_tree2() +
  theme(
    legend.position = "right",
    legend.text     = element_text(size = 7),
    legend.title    = element_text(size = 8, face = "bold"),
    plot.title      = element_text(size = 11, face = "bold", hjust = 0.5),
    plot.margin     = margin(5, 20, 40, 5),
    axis.text.x     = element_blank()
  ) +
  labs(title = "Phylogenetic Tree - Per-sample relative abundance heatmap")

save_plot(p2, file.path(out_dir, "02_tree_abundance_heatmap.pdf"),
          w = 16 + length(heat_cols) * 0.25, h = max(10, n_kept * 0.22))

# -- 9. Plot 03: phylum colour strip -------------------------------------------
cat("[Tree] Building plot 03: phylum colour strip ...\n")

p3 <- ggtree(tree_p, layout = "rectangular", size = 0.3, color = "grey40") %<+% tip_data +
  geom_tippoint(aes(color = phylum), size = 1.8, alpha = 0.8) +
  geom_tiplab(aes(label = tax_label), size = 2.2, hjust = -0.08) +
  geom_fruit(
    data    = strip_df,
    geom    = geom_tile,
    mapping = aes(y = tip, x = 1, fill = fr_phylum),
    offset  = 0.04, pwidth = 0.04
  ) +
  scale_color_manual(values = phy_pal, na.value = "#AAAAAA", name = "Phylum") +
  scale_fill_manual(values  = phy_pal, na.value = "#CCCCCC",  guide  = "none") +
  base_theme +
  labs(title = "Phylogenetic Tree - Phylum colour strip")

save_plot(p3, file.path(out_dir, "03_tree_phylum_colorstrip.pdf"), w = 14, h = max(10, n_kept * 0.22))

# -- 10. Plot 04: LEfSe differential LDA bars ----------------------------------
cat("[Tree] Building plot 04: differential tree ...\n")

if (lda_present) {
  p4 <- ggtree(tree_p, layout = "rectangular", size = 0.3, color = "grey40") %<+% tip_data +
    geom_tippoint(aes(color = phylum), size = 1.5, alpha = 0.75) +
    geom_tiplab(aes(label = tax_label), size = 2.0, hjust = -0.08) +
    geom_fruit(
      data    = diff_df,
      geom    = geom_col,
      mapping = aes(y = tip, x = fr_lda, fill = fr_dir),
      offset = 0.05, pwidth = 0.35, orientation = "y"
    ) +
    scale_color_manual(values = phy_pal, na.value = "#AAAAAA", name = "Phylum") +
    scale_fill_manual(
      values   = c("Enriched" = "#E41A1C", "Depleted" = "#377EB8", "Not tested" = "#DDDDDD"),
      na.value = "#DDDDDD", name = "Differential"
    ) +
    base_theme +
    labs(title = "Phylogenetic Tree - LEfSe differential abundance (LDA score)")
} else {
  p4 <- ggtree(tree_p, layout = "rectangular", size = 0.3, color = "grey40") %<+% tip_data +
    geom_tippoint(aes(color = phylum), size = 1.8) +
    geom_tiplab(aes(label = tax_label), size = 2.2, hjust = -0.08) +
    scale_color_manual(values = phy_pal, na.value = "#AAAAAA", name = "Phylum") +
    base_theme +
    labs(title = "Phylogenetic Tree - LEfSe (no data available)")
}

save_plot(p4, file.path(out_dir, "04_tree_differential.pdf"), w = 14, h = max(10, n_kept * 0.22))

# -- 11. Plot 05: circular layout ----------------------------------------------
cat("[Tree] Building plot 05: circular tree ...\n")

p5 <- ggtree(tree_p, layout = "fan", open.angle = 25,
             size = 0.3, color = "grey50") %<+% tip_data +
  geom_tippoint(aes(size = mean_abund, color = phylum), alpha = 0.85) +
  geom_tiplab(aes(label = tax_label), size = 1.8, offset = 0.002, align = FALSE) +
  geom_fruit(
    data    = strip_df,
    geom    = geom_tile,
    mapping = aes(y = tip, x = 1, fill = fr_phylum),
    offset  = 0.05, pwidth = 0.05
  ) +
  scale_color_manual(values = phy_pal, na.value = "#AAAAAA", name = "Phylum") +
  scale_fill_manual(values  = phy_pal, na.value = "#CCCCCC",  guide  = "none") +
  scale_size_continuous(range = c(1, 5), name = "Mean rel. abund.",
                        labels = scales::percent_format(accuracy = 0.1)) +
  theme_tree2() +
  theme(
    legend.position  = "right",
    legend.text      = element_text(size = 7),
    legend.title     = element_text(size = 8, face = "bold"),
    plot.title       = element_text(size = 11, face = "bold", hjust = 0.5),
    plot.margin      = margin(10, 10, 10, 10)
  ) +
  labs(title = "Phylogenetic Tree - Circular layout (top genera by abundance)")

save_plot(p5, file.path(out_dir, "05_tree_circular.pdf"), w = 18, h = 18)

cat("[Tree] All 5 plots complete.\n")
