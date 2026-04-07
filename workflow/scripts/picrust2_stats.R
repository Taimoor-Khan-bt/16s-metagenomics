#!/usr/bin/env Rscript
# =============================================================================
# picrust2_stats.R — Robust Multi-group Differential Pathway Analysis
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# ── Arguments ─────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: picrust2_stats.R <pathways_tsv> <metadata> <group_col> <out_dir> [strategy] [dual_plots] [nsti_file]")

pathways_file <- args[1]
meta_file     <- args[2]
group_col     <- args[3]
out_dir       <- args[4]
strategy      <- if (length(args) >= 5) args[5] else "rename"
dual_plots    <- isTRUE(tolower(if (length(args) >= 6) args[6] else "false") == "true")
nsti_file     <- if (length(args) >= 7 && nchar(trimws(args[7])) > 0) args[7] else NULL

# =============================================================================
# clean_pathway_label()
# Converts MetaCyc pathway IDs/names to readable title-cased strings.
#   In:  "PWY-7111: pyruvate_fermentation_to_isobutanol"
#   Out: "Pyruvate Fermentation To Isobutanol"
# =============================================================================
clean_pathway_label <- function(x) {
  if (is.na(x) || nchar(trimws(x)) == 0) return(x)
  # Strip pathway ID prefix (e.g. "PWY-1234: ")
  lbl <- sub("^[A-Z0-9.-]+[:-] ", "", x)
  # Convert underscores and hyphens to spaces
  lbl <- gsub("[_-]", " ", lbl)
  # Title-case each word
  lbl <- paste(sapply(strsplit(lbl, " ")[[1]], function(w) {
    if (nchar(w) == 0) return(w)
    paste0(toupper(substr(w, 1, 1)), tolower(substr(w, 2, nchar(w))))
  }), collapse = " ")
  # Truncate overly long labels
  if (nchar(lbl) > 60) lbl <- paste0(substr(lbl, 1, 57), "...")
  return(lbl)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── 1. Load Data ─────────────────────────────────────────────────────────────
tab <- read.table(pathways_file, header = TRUE, sep = "\t",
                  row.names = 1, check.names = FALSE, comment.char = "#")

meta <- read.table(meta_file, header = TRUE, sep = "\t",
                   comment.char = "", row.names = 1, check.names = FALSE)

# Sync samples
common_samples <- intersect(colnames(tab), rownames(meta))
if (length(common_samples) < 3) stop("Fewer than 3 common samples. Check Sample IDs.")
tab  <- tab[, common_samples, drop = FALSE]
meta <- meta[common_samples, , drop = FALSE]

# Ensure group is a factor and remove levels with zero samples
meta[[group_col]] <- as.factor(meta[[group_col]])
grp  <- meta[[group_col]]
lvls <- levels(grp)

# ── NSTI filtering ─────────────────────────────────────────────────────────
nsti_max <- 2.0  # samples with mean NSTI >= this are excluded
if (!is.null(nsti_file) && file.exists(nsti_file)) {
  nsti_tab <- tryCatch(
    read.table(nsti_file, header = TRUE, sep = "\t",
               row.names = 1, check.names = FALSE, comment.char = "#"),
    error = function(e) { message("NSTI file unreadable: ", e$message); NULL }
  )
  if (!is.null(nsti_tab)) {
    # PICRUSt2 NSTI column is named "metadata_NSTI" in marker_predicted_and_nsti.tsv
    nsti_col <- if ("metadata_NSTI" %in% colnames(nsti_tab)) "metadata_NSTI"
               else colnames(nsti_tab)[1]
    nsti_vals <- setNames(as.numeric(nsti_tab[[nsti_col]]), rownames(nsti_tab))
    nsti_summary <- data.frame(
      sample    = names(nsti_vals),
      mean_NSTI = nsti_vals,
      status    = ifelse(nsti_vals < nsti_max, "kept", "excluded"),
      stringsAsFactors = FALSE
    )
    write.table(nsti_summary, file.path(out_dir, "nsti_summary.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    excluded <- nsti_summary$sample[nsti_summary$status == "excluded"]
    if (length(excluded) > 0) {
      message("NSTI filter: excluding ", length(excluded), " sample(s) with mean NSTI >= ",
              nsti_max, ": ", paste(excluded, collapse = ", "))
      common_samples <- setdiff(common_samples, excluded)
      tab  <- tab[, common_samples, drop = FALSE]
      meta <- meta[common_samples, , drop = FALSE]
      grp  <- meta[[group_col]]
    } else {
      message("NSTI filter: all ", length(common_samples), " samples retained (NSTI < ", nsti_max, ")")
    }
  }
} else {
  write.table(data.frame(Note = if (is.null(nsti_file)) "NSTI file not provided"
                                else paste0("NSTI file not found: ", nsti_file)),
              file.path(out_dir, "nsti_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message("NSTI filtering skipped — nsti_summary.tsv placeholder written.")
}

message("Analysis Groups: ", paste(lvls, collapse=", "))

# Dark palette (consistent with alpha/beta scripts)
grp_levels_p  <- levels(meta[[group_col]])
dark_pal_p    <- c("#1B4F72","#922B21","#1D8348","#6C3483","#784212",
                   "#0E6655","#4A235A","#1A5276","#7D6608","#212F3D")
dark_colors_p <- setNames(dark_pal_p[seq_along(grp_levels_p)], grp_levels_p)

bold_theme_p <- theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.title    = element_text(face = "bold", size = 12),
    legend.text     = element_text(face = "bold", size = 11),
    axis.title      = element_text(face = "bold", size = 12),
    axis.text       = element_text(face = "bold", size = 11, color = "black"),
    plot.title      = element_text(face = "bold", size = 13),
    panel.grid.minor = element_blank()
  )

# ── 2. Kruskal-Wallis Test (Omnibus) ─────────────────────────────────────────
# Correct for 3+ groups instead of just Wilcoxon
results <- lapply(rownames(tab), function(pw) {
  vals <- as.numeric(tab[pw, ])
  if (sum(vals) == 0) return(NULL)
  
  # Omnibus test: Is there any difference across the 3 groups?
  p_val <- tryCatch(kruskal.test(vals ~ grp)$p.value, error = function(e) NA)
  
  # Calculate means per group
  group_means <- tapply(vals, grp, mean)
  
  res <- data.frame(
    pathway = pw,
    p_value = p_val,
    stringsAsFactors = FALSE
  )
  # Append means dynamically for any number of groups
  for(l in lvls) { res[[paste0("mean_", l)]] <- group_means[l] }
  return(res)
})

res_df <- do.call(rbind, Filter(Negate(is.null), results))
res_df$p_adj <- p.adjust(res_df$p_value, method = "BH")
res_df <- res_df[order(res_df$p_adj), ]

write.table(res_df, file.path(out_dir, "pathway_differential.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ── 3. Visualization ──────────────────────────────────────────────────────────
# Always take the top 30 pathways by ascending q-value (no significance filter).
# Significance is communicated via transparency rather than exclusion.
sig <- head(res_df[!is.na(res_df$p_adj), ], 30)

# Prepare data for plotting (Pivot to long format)
plot_df <- sig %>%
  select(pathway, p_adj, starts_with("mean_")) %>%
  pivot_longer(cols = starts_with("mean_"), names_to = "Group", values_to = "Abundance") %>%
  mutate(
    Group    = gsub("mean_", "", Group),
    sig_tier = case_when(
      p_adj < 0.05  ~ "q < 0.05",
      p_adj < 0.20  ~ "q < 0.20",
      TRUE          ~ "NS"
    ),
    sig_alpha = case_when(
      sig_tier == "q < 0.05" ~ 1.0,
      sig_tier == "q < 0.20" ~ 0.65,
      TRUE                   ~ 0.35
    ),
    sig_tier = factor(sig_tier, levels = c("NS", "q < 0.20", "q < 0.05"))
  )

# ── Helper to build bar + heatmap panels ──────────────────────────────────
make_bar_panel <- function(df, title_str) {
  ggplot(df, aes(x = reorder(pathway, Abundance), y = Abundance,
                 fill = Group, alpha = sig_tier)) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip() +
    scale_fill_manual(values = dark_colors_p, name = grp_display_p) +
    scale_alpha_manual(values = c("NS" = 0.35, "q < 0.20" = 0.65, "q < 0.05" = 1.0),
                       name = "FDR tier") +
    bold_theme_p +
    labs(title    = title_str,
         subtitle = "Omnibus Kruskal-Wallis Test  |  Shading = FDR significance tier",
         x = "MetaCyc Pathway", y = "Relative Abundance")
}

make_heat_panel <- function(df, title_str) {
  ggplot(df, aes(x = Group, y = pathway, fill = log10(Abundance + 1))) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_viridis_c(option = "plasma", name = "log10(Abund+1)") +
    bold_theme_p +
    labs(title = title_str, x = "Group", y = "Pathway")
}

# ── Category summary panel ──────────────────────────────────────────────
cat_keywords <- c("Biosynthesis", "Degradation", "Fermentation",
                   "Carbohydrate", "Amino.Acid|Amino Acid", "Nucleotide",
                   "Lipid", "Energy", "Cofactor", "Fatty.Acid|Fatty Acid",
                   "Sugar", "Vitamin", "Glycolysis", "TCA", "Respiration")

cat_labels <- c("Biosynthesis", "Degradation", "Fermentation",
                "Carbohydrate", "Amino Acid", "Nucleotide",
                "Lipid", "Energy", "Cofactor", "Fatty Acid",
                "Sugar", "Vitamin", "Glycolysis", "TCA", "Respiration")

assign_category <- function(pw) {
  for (i in seq_along(cat_keywords)) {
    if (grepl(cat_keywords[i], pw, ignore.case = TRUE)) return(cat_labels[i])
  }
  "Other"
}

build_category_panel <- function(df_long, all_pathways) {
  # df_long: long-format data with columns pathway, Group, Abundance
  full_long <- all_pathways %>%
    select(pathway, starts_with("mean_")) %>%
    pivot_longer(cols = starts_with("mean_"), names_to = "Group", values_to = "Abundance") %>%
    mutate(Group    = gsub("mean_", "", Group),
           Category = vapply(pathway, assign_category, character(1)))
  cat_df <- full_long %>%
    dplyr::group_by(Category, Group) %>%
    dplyr::summarise(MeanAbundance = mean(Abundance, na.rm = TRUE), .groups = "drop")
  ggplot(cat_df, aes(x = MeanAbundance, y = reorder(Category, MeanAbundance),
                     fill = Group)) +
    geom_col(position = "dodge", alpha = 0.88) +
    scale_fill_manual(values = dark_colors_p, name = grp_display_p) +
    bold_theme_p +
    labs(title = "Pathway Category Summary",
         x = "Mean Relative Abundance", y = "Category")
}

grp_display_p <- tools::toTitleCase(gsub("_", " ", group_col))

pdf(file.path(out_dir, "pathway_plots_raw.pdf"), width = 18, height = 9)

plot_df_raw <- plot_df  # preserve raw pathway names

# Page 1: multipanel bar + heatmap side-by-side
p_bar_raw  <- make_bar_panel(plot_df_raw  %>% mutate(pathway = pathway),
                              "Top Predicted Pathways (raw labels)")
p_heat_raw <- make_heat_panel(plot_df_raw %>% mutate(pathway = pathway),
                               "Pathway Heatmap — raw labels (Log10)")

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  print(p_bar_raw + p_heat_raw + patchwork::plot_layout(widths = c(1.4, 1)) +
          patchwork::plot_annotation(tag_levels = "A"))
} else {
  print(p_bar_raw); print(p_heat_raw)
}

# Page 2: category summary
tryCatch({
  p_cat_raw <- build_category_panel(plot_df_raw, res_df)
  print(p_cat_raw)
}, error = function(e) {
  message("Category summary failed (raw): ", e$message)
  plot.new(); title("Pathway Category Summary unavailable")
})
dev.off()

ggsave(file.path(out_dir, "pathway_plots_raw.png"),
       plot = p_bar_raw, width = 14, height = 9, units = "in", dpi = 600)
ggsave(file.path(out_dir, "pathway_plots_raw_heatmap.png"),
       plot = p_heat_raw, width = 10, height = 9, units = "in", dpi = 600)

# Cleaned labels
if (strategy != "none") {
  plot_df_clean <- plot_df %>% mutate(pathway = sapply(pathway, clean_pathway_label))
} else {
  plot_df_clean <- plot_df
}

pdf(file.path(out_dir, "pathway_plots.pdf"), width = 18, height = 9)

# Page 1: multipanel bar + heatmap
p_bar  <- make_bar_panel(plot_df_clean, "Top Predicted Pathways (Abundance by Group)")
p_heat <- make_heat_panel(plot_df_clean, "Pathway Abundance Heatmap (Log10)")

if (requireNamespace("patchwork", quietly = TRUE)) {
  print(p_bar + p_heat + patchwork::plot_layout(widths = c(1.4, 1)) +
          patchwork::plot_annotation(tag_levels = "A"))
} else {
  print(p_bar); print(p_heat)
}

# Page 2: category summary
tryCatch({
  p_cat <- build_category_panel(plot_df_clean, res_df)
  print(p_cat)
}, error = function(e) {
  message("Category summary failed: ", e$message)
  plot.new(); title("Pathway Category Summary unavailable")
})
dev.off()

ggsave(file.path(out_dir, "pathway_plots.png"),
       plot = p_bar, width = 14, height = 9, units = "in", dpi = 600)
ggsave(file.path(out_dir, "pathway_plots_heatmap.png"),
       plot = p_heat, width = 10, height = 9, units = "in", dpi = 600)
message("Saved: pathway_differential.tsv, pathway_plots.pdf/.png, pathway_plots_heatmap.png, pathway_plots_raw.pdf/.png, pathway_plots_raw_heatmap.png")