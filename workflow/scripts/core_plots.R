#!/usr/bin/env Rscript
# =============================================================================
# core_plots.R — Publication-quality PNGs from core pipeline exported TSVs
#
# All outputs at 600 DPI in <viz_root>/{qc,diversity,taxonomy}/
#
# PNGs generated:
#   qc/dada2_stats.png                    — read-retention waterfall
#   qc/read_depth.png                     — per-sample total reads barplot
#   diversity/alpha_diversity_overview.png — faceted alpha-metric boxplots
#   diversity/pcoa_<metric>.png           — PCoA scatter per beta metric
#   diversity/pcoa_overview.png           — 3-panel: Aitchison / BC PCoA / BC boxplot
#   taxonomy/taxonomy_barplot.png         — phylum stacked bars per sample
#   taxonomy/top_taxa_abundance.png       — top-20 phyla lollipop chart
#
# Usage:
#   Rscript core_plots.R \
#       <dada2_stats_tsv>     \   # exported/dada2_stats/stats.tsv
#       <alpha_dir>           \   # exported/alpha_diversity/
#       <metadata_tsv>        \   # config/metadata.tsv
#       <group_col>           \   # e.g. Stunting
#       <viz_root>            \   # visualizations/
#       <feature_table_tsv>   \   # exported/feature_table/feature-table.tsv
#       <taxonomy_tsv>        \   # exported/taxonomy/taxonomy.tsv
#       <beta_dir>            \   # exported/beta_diversity/
#       [palette]                 # color palette name (default: okabe_ito)
# ============================================================================= 

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(ape)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 9) {
  stop(paste(
    "Usage: core_plots.R",
    "<dada2_stats_tsv> <alpha_dir> <metadata_tsv> <group_col> <viz_root>",
    "<feature_table_tsv> <taxonomy_tsv> <beta_dir> <bray_dm_tsv>"
  ))
}

dada2_tsv      <- args[1]
alpha_dir      <- args[2]
meta_file      <- args[3]
group_col      <- args[4]
viz_root       <- args[5]
feat_table_tsv <- args[6]
taxonomy_tsv   <- args[7]
beta_dir       <- args[8]
bray_dm_tsv    <- args[9]
palette_name   <- if (length(args) >= 10) args[10] else "okabe_ito"
grp_display    <- tools::toTitleCase(gsub("_", " ", group_col))

# ── Output subdirectories ─────────────────────────────────────────────────────
qc_dir  <- file.path(viz_root, "qc")
div_dir <- file.path(viz_root, "diversity")
tax_dir <- file.path(viz_root, "taxonomy")
for (d in c(qc_dir, div_dir, tax_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ── Shared helpers ─────────────────────────────────────────────────────────────
save_png <- function(p, path, w = 12, h = 7) {
  tryCatch({
    ggsave(path, plot = p, width = w, height = h, units = "in", dpi = 600,
           bg = "white")
    message("Saved: ", basename(path))
  }, error = function(e) message("ERROR saving ", basename(path), ": ", e$message))
}

read_tsv_q2 <- function(path, ...) {
  df <- read.table(path, header = TRUE, sep = "\t",
                   comment.char = "", quote = "", check.names = FALSE, ...)
  names(df)[1] <- "sample_id"
  df[!grepl("^#", df$sample_id), ]
}

read_metadata <- function(path) {
  tryCatch(read_tsv_q2(path), error = function(e) {
    message("Cannot read metadata: ", e$message); NULL
  })
}

# =============================================================================
# Palette dispatcher — resolves a named palette to a hex color vector.
# All palettes are defined as hex vectors (zero extra R package dependencies).
# To switch palette: change plots.color_palette in config/config.yaml
# =============================================================================
get_palette <- function(name, n) {
  palettes <- list(
    # Wong 2011, Nature Methods — explicit colorblind-safe recommendation
    okabe_ito    = c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                     "#0072B2", "#D55E00", "#CC79A7", "#000000"),
    # Nature Publishing Group (ggsci::pal_npg compatible hex values)
    npg          = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488",
                     "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
                     "#7E6148", "#B09C85"),
    # The Lancet journal (ggsci::pal_lancet compatible hex values)
    lancet       = c("#00468B", "#ED0000", "#42B540", "#0099B4",
                     "#925E9F", "#FDAF91", "#AD002A", "#ADB6B6", "#1B1919"),
    # Journal of Clinical Oncology (ggsci::pal_jco compatible hex values)
    jco          = c("#0073C2", "#EFC000", "#868686", "#CD534C",
                     "#7AA6DC", "#003C67", "#8F7700", "#3B3B3B",
                     "#A73030", "#4A6990"),
    # Science / AAAS (ggsci::pal_aaas compatible hex values)
    aaas         = c("#3B4992", "#EE0000", "#008B45", "#631879",
                     "#008280", "#BB0021", "#5F559B", "#A20056",
                     "#808180", "#1B1919"),
    # Tableau 10 — high-contrast categorical palette
    tableau10    = c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2",
                     "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7",
                     "#9C755F", "#BAB0AC"),
    # RColorBrewer Dark2 — print-safe, up to 8 groups
    brewer_dark2 = RColorBrewer::brewer.pal(8,  "Dark2"),
    # RColorBrewer Set1 — high-contrast, up to 9 groups
    brewer_set1  = RColorBrewer::brewer.pal(9,  "Set1")
  )
  base <- palettes[[name]]
  if (is.null(base)) {
    message("[core_plots] Unknown palette '", name, "', falling back to okabe_ito")
    base <- palettes[["okabe_ito"]]
  }
  if (n <= length(base)) return(base[seq_len(n)])
  # Interpolate when more groups than palette entries
  colorRampPalette(base)(n)
}

get_group_colors <- function(levels_vec) {
  setNames(get_palette(palette_name, length(levels_vec)), levels_vec)
}

# =============================================================================
# Plot 1: DADA2 read-retention waterfall
# =============================================================================
plot_dada2 <- function() {
  dada2 <- tryCatch(read_tsv_q2(dada2_tsv),
                    error = function(e) { message("Skipping DADA2 plot: ", e$message); NULL })
  if (is.null(dada2)) return(invisible(NULL))

  count_cols <- intersect(c("input", "filtered", "denoised", "merged", "non-chimeric"),
                          names(dada2))
  dada2_long <- dada2 %>%
    select(sample_id, all_of(count_cols)) %>%
    mutate(across(all_of(count_cols), as.numeric)) %>%
    pivot_longer(all_of(count_cols), names_to = "step", values_to = "reads") %>%
    mutate(step = factor(step, levels = count_cols))

  step_means <- dada2_long %>%
    group_by(step) %>%
    summarise(mean_reads = mean(reads, na.rm = TRUE), .groups = "drop")

  p <- ggplot(dada2_long, aes(x = step, y = reads, group = sample_id)) +
    geom_line(alpha = 0.25, colour = "steelblue", linewidth = 0.5) +
    geom_point(alpha = 0.45, colour = "steelblue", size = 2) +
    geom_line(data = step_means, aes(x = step, y = mean_reads, group = 1),
              colour = "firebrick", linewidth = 1.5, inherit.aes = FALSE) +
    geom_point(data = step_means, aes(x = step, y = mean_reads),
               colour = "firebrick", size = 4, inherit.aes = FALSE) +
    scale_y_continuous(labels = comma_format()) +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(face = "bold"),
          panel.grid.minor = element_blank()) +
    labs(title    = "DADA2 Read Retention",
         subtitle = "Blue = per-sample  |  Red = mean across samples",
         x = "Denoising Step", y = "Read Count")

  save_png(p, file.path(qc_dir, "dada2_stats.png"), w = 10, h = 6)
}

# =============================================================================
# Plot 2: Per-sample read depth barplot
# =============================================================================
plot_read_depth <- function(meta) {
  ft <- tryCatch({
    df <- read.table(feat_table_tsv, header = TRUE, sep = "\t",
                     comment.char = "", quote = "", check.names = FALSE,
                     skip = 1)  # skip the '# Constructed from biom file' line
    # First column is '#OTU ID' (feature IDs); samples start at column 2
    sample_cols <- names(df)[-1]
    colSums(df[, sample_cols, drop = FALSE])
  }, error = function(e) { message("Skipping read_depth: ", e$message); NULL })
  if (is.null(ft)) return(invisible(NULL))

  depth_df <- data.frame(sample_id = names(ft), reads = as.numeric(ft),
                          stringsAsFactors = FALSE)
  if (!is.null(meta) && group_col %in% names(meta)) {
    depth_df <- left_join(depth_df, meta[, c("sample_id", group_col)], by = "sample_id")
    fill_col <- group_col
  } else {
    depth_df[[".group"]] <- "All"
    fill_col <- ".group"
  }

  p <- ggplot(depth_df,
              aes(x = reorder(sample_id, -reads), y = reads,
                  fill = !!sym(fill_col))) +
    geom_col(colour = "white", linewidth = 0.3) +
    scale_y_continuous(labels = comma_format()) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x          = element_text(angle = 55, hjust = 1, size = 9),
          panel.grid.major.x   = element_blank(),
          legend.title         = element_text(face = "bold")) +
    labs(title = "Per-Sample Read Depth (raw feature table)",
         x = "Sample", y = "Total Reads", fill = grp_display)

  save_png(p, file.path(qc_dir, "read_depth.png"), w = 12, h = 6)
}

# =============================================================================
# Plot 3: Alpha diversity — bar chart (mean ± SE + jitter), patchwork A/B/C
# =============================================================================
plot_alpha <- function(meta) {
  if (is.null(meta) || !(group_col %in% names(meta))) {
    message("Skipping alpha plot: group '", group_col, "' not in metadata")
    return(invisible(NULL))
  }

  metric_dirs <- list.dirs(alpha_dir, recursive = FALSE, full.names = TRUE)
  alpha_list  <- lapply(metric_dirs, function(mdir) {
    tsv <- file.path(mdir, "alpha-diversity.tsv")
    if (!file.exists(tsv)) return(NULL)
    df <- tryCatch(read_tsv_q2(tsv), error = function(e) NULL)
    if (is.null(df) || ncol(df) < 2) return(NULL)
    df$metric <- basename(mdir)
    df$value  <- suppressWarnings(as.numeric(df[[2]]))
    df[, c("sample_id", "metric", "value")]
  })
  alpha_df <- bind_rows(Filter(Negate(is.null), alpha_list))
  if (nrow(alpha_df) == 0) { message("No alpha data found"); return(invisible(NULL)) }

  alpha_df <- inner_join(alpha_df,
                          meta[, c("sample_id", group_col)], by = "sample_id") %>%
    filter(!is.na(value))

  metric_labels <- c(
    shannon           = "Shannon Index",
    simpson           = "Simpson Index",
    chao1             = "Chao1 Richness",
    faith_pd          = "Faith's PD",
    evenness          = "Pielou's Evenness",
    observed_features = "Observed ASVs"
  )

  all_metrics <- unique(alpha_df$metric)
  pub_metrics <- intersect(c("shannon", "simpson", "chao1"), all_metrics)
  if (length(pub_metrics) == 0) pub_metrics <- all_metrics

  grp_levels  <- sort(unique(as.character(meta[[group_col]])))
  dark_colors <- get_group_colors(grp_levels)

  panel_plots <- list()
  for (m in pub_metrics) {
    df_m <- alpha_df[alpha_df$metric == m, ]
    df_m[[group_col]] <- factor(df_m[[group_col]], levels = grp_levels)

    summ <- df_m %>%
      group_by(!!sym(group_col)) %>%
      summarise(mean_val = mean(value, na.rm = TRUE),
                se_val   = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
                .groups  = "drop")

    p_val <- tryCatch({
      y   <- df_m$value
      grp <- as.factor(df_m[[group_col]])
      if (length(levels(grp)) == 2) wilcox.test(y ~ grp)$p.value
      else                          kruskal.test(y ~ grp)$p.value
    }, error = function(e) NA_real_)

    p_lab <- if (!is.na(p_val)) {
      if      (p_val < 0.001) "p < 0.001 ***"
      else if (p_val < 0.01)  paste0("p = ", signif(p_val, 2), " **")
      else if (p_val < 0.05)  paste0("p = ", signif(p_val, 2), " *")
      else                    paste0("p = ", signif(p_val, 2), " (NS)")
    } else ""

    m_label <- if (!is.na(metric_labels[m])) metric_labels[m] else m

    p <- ggplot(summ, aes(x = !!sym(group_col), y = mean_val,
                          fill = !!sym(group_col))) +
      geom_col(width = 0.55, alpha = 0.90, color = "grey20", linewidth = 0.45) +
      geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val),
                    width = 0.18, linewidth = 0.9, color = "grey20") +
      geom_jitter(data = df_m,
                  aes(x = !!sym(group_col), y = value, fill = !!sym(group_col)),
                  width = 0.12, size = 2.4, shape = 21,
                  color = "grey20", alpha = 0.75, inherit.aes = FALSE) +
      scale_fill_manual(values = dark_colors, name = grp_display) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
      theme_classic(base_size = 13) +
      theme(
        legend.position = "bottom",
        legend.title    = element_text(face = "bold", size = 12),
        legend.text     = element_text(face = "bold", size = 11),
        axis.title      = element_text(face = "bold", size = 12),
        axis.text       = element_text(face = "bold", size = 11, color = "black"),
        plot.title      = element_text(face = "bold", size = 13),
        axis.line       = element_line(linewidth = 0.6)
      ) +
      labs(x = NULL, y = m_label, caption = p_lab)
    panel_plots[[m]] <- p
  }

  n_panels <- length(panel_plots)
  if (n_panels == 0) { message("No alpha panels created"); return(invisible(NULL)) }

  n_cols <- min(n_panels, 3L)           # max 3 columns; rows grow automatically
  n_rows <- ceiling(n_panels / n_cols)
  fig_w  <- 4.5 * n_cols               # 4.5 in per column
  fig_h  <- 5.5 * n_rows               # 5.5 in per row
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined_plot <- patchwork::wrap_plots(panel_plots, ncol = n_cols, nrow = n_rows) +
      patchwork::plot_annotation(tag_levels = "A") +
      patchwork::plot_layout(guides = "collect") &
      theme(legend.position = "bottom")
  } else {
    combined_plot <- panel_plots[[1]]
  }

  save_png(combined_plot, file.path(div_dir, "alpha_diversity_overview.png"),
           w = fig_w, h = fig_h)
}

# =============================================================================
# Shared PCoA parser (scikit-bio / QIIME 2 ordination.txt format)
# =============================================================================
parse_ordination <- function(ord_file) {
  lines <- readLines(ord_file, warn = FALSE)

  # Proportion explained — line starts with 'Proportion explained'; values on next line
  prop_row  <- which(startsWith(trimws(lines), "Proportion explained"))
  prop_vals <- if (length(prop_row) > 0 && length(lines) >= prop_row + 1)
    suppressWarnings(as.numeric(strsplit(trimws(lines[prop_row + 1]), "\t")[[1]]))
  else rep(NA_real_, 10)

  # Site block — line starts with 'Site'; coordinate rows follow immediately
  site_row <- which(startsWith(trimws(lines), "Site") &
                    !startsWith(trimws(lines), "Site constraints"))
  if (length(site_row) == 0) {
    message("No 'Site' block in ", ord_file); return(NULL)
  }
  keywords    <- c("Biplot", "Site constraints", "Species", "Eigvals", "Proportion explained")
  after_site  <- which(Reduce(`|`, lapply(keywords, function(k) startsWith(trimws(lines), k))) &
                       seq_along(lines) > site_row[1])
  end_row     <- if (length(after_site) > 0) after_site[1] - 2 else length(lines)

  data_lines <- lines[(site_row[1] + 1):end_row]
  data_lines <- data_lines[nchar(trimws(data_lines)) > 0]
  if (length(data_lines) == 0) { message("Empty Site block in ", ord_file); return(NULL) }

  parsed <- lapply(data_lines, function(l) {
    pts <- strsplit(trimws(l), "\t")[[1]]
    if (length(pts) < 3) return(NULL)
    data.frame(sample_id = pts[1], PC1 = as.numeric(pts[2]), PC2 = as.numeric(pts[3]),
               stringsAsFactors = FALSE)
  })
  df <- bind_rows(Filter(Negate(is.null), parsed))
  list(df = df, pct_var = round(prop_vals * 100, 1))
}

make_pcoa_gg <- function(ord, label, meta) {
  if (is.null(ord)) return(NULL)
  df <- ord$df
  pv <- ord$pct_var

  if (!is.null(meta) && group_col %in% names(meta))
    df <- left_join(df, meta[, c("sample_id", group_col)], by = "sample_id")
  else
    df[[group_col]] <- "All"

  df[[group_col]] <- factor(df[[group_col]])
  grp_levels  <- levels(df[[group_col]])
  dark_colors <- get_group_colors(grp_levels)

  xlab <- if (!is.na(pv[1])) paste0("PC1 (", pv[1], "% var)") else "PC1"
  ylab <- if (!is.na(pv[2])) paste0("PC2 (", pv[2], "% var)") else "PC2"

  p <- ggplot(df, aes(x = PC1, y = PC2, color = !!sym(group_col))) +
    geom_point(size = 3.8, alpha = 0.82) +
    stat_ellipse(level = 0.80, linetype = 2, linewidth = 0.7) +
    scale_color_manual(values = dark_colors, name = grp_display) +
    theme_bw(base_size = 13) +
    theme(
      legend.position  = "bottom",
      legend.title     = element_text(face = "bold", size = 12),
      legend.text      = element_text(face = "bold", size = 11),
      axis.title       = element_text(face = "bold", size = 12),
      axis.text        = element_text(face = "bold", size = 11, color = "black"),
      plot.title       = element_text(face = "bold", size = 13),
      panel.grid.minor = element_blank()
    ) +
    labs(title = paste("PCoA —", label), x = xlab, y = ylab)
  p
}

# =============================================================================
# Plots 4 & 5: Beta diversity — individual PCoA per metric + 3-panel overview
#   Panel A: Aitchison PCoA (CLR-transformed feature table, ape::pcoa)
#   Panel B: Bray-Curtis PCoA (from QIIME2 ordination.txt)
#   Panel C: Bray-Curtis within- vs between-group dissimilarity boxplot
# =============================================================================
plot_beta <- function(meta) {
  metric_dirs <- list.dirs(beta_dir, recursive = FALSE, full.names = TRUE)
  if (length(metric_dirs) == 0) {
    message("No beta PCoA data in ", beta_dir); return(invisible(NULL))
  }

  label_map <- c(
    bray_curtis        = "Bray-Curtis",
    jaccard            = "Jaccard",
    weighted_unifrac   = "Weighted UniFrac",
    unweighted_unifrac = "Unweighted UniFrac"
  )

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    message("Install 'patchwork' for pcoa_overview.png (skipped)")
    return(invisible(NULL))
  }

  grp_levels  <- sort(unique(as.character(meta[[group_col]])))
  dark_colors <- get_group_colors(grp_levels)

  bold_theme <- theme_bw(base_size = 13) +
    theme(
      legend.position  = "bottom",
      legend.title     = element_text(face = "bold", size = 12),
      legend.text      = element_text(face = "bold", size = 11),
      axis.title       = element_text(face = "bold", size = 12),
      axis.text        = element_text(face = "bold", size = 11, color = "black"),
      plot.title       = element_text(face = "bold", size = 13),
      panel.grid.minor = element_blank()
    )

  # ── Panel A: Aitchison PCoA (CLR from feature table) ─────────────────────
  p_ait <- tryCatch({
    raw_tab <- read.table(feat_table_tsv, header = TRUE, sep = "\t",
                          comment.char = "", quote = "", check.names = FALSE,
                          skip = 1)
    # rows = features, cols = samples → transpose so rows = samples
    mat <- as.matrix(raw_tab[, -1, drop = FALSE])
    storage.mode(mat) <- "numeric"
    tab <- t(mat)

    # Align with metadata
    if (!is.null(meta)) {
      common <- intersect(rownames(tab), meta$sample_id)
      tab    <- tab[common, , drop = FALSE]
    }
    if (nrow(tab) < 3) stop("Too few samples for Aitchison PCoA")

    # CLR transform (pseudo-count 0.5)
    tab_ps  <- tab + 0.5
    gm      <- apply(tab_ps, 1, function(x) exp(mean(log(x))))
    clr_tab <- log(tab_ps / gm)

    ait_dist <- dist(clr_tab)
    pcoa_ait <- ape::pcoa(ait_dist)
    pcoa_df  <- as.data.frame(pcoa_ait$vectors[, 1:2])
    colnames(pcoa_df) <- c("PC1", "PC2")
    pcoa_df$sample_id <- rownames(pcoa_df)
    pct_var  <- round(pcoa_ait$values$Relative_eig[1:2] * 100, 1)

    if (!is.null(meta) && group_col %in% names(meta))
      pcoa_df <- left_join(pcoa_df, meta[, c("sample_id", group_col)], by = "sample_id")
    else
      pcoa_df[[group_col]] <- "All"

    pcoa_df[[group_col]] <- factor(pcoa_df[[group_col]], levels = grp_levels)

    ggplot(pcoa_df, aes(x = PC1, y = PC2, color = !!sym(group_col))) +
      geom_point(size = 3.8, alpha = 0.82) +
      stat_ellipse(level = 0.80, linetype = 2, linewidth = 0.7) +
      scale_color_manual(values = dark_colors, name = grp_display) +
      bold_theme +
      labs(title = "Aitchison PCoA (CLR)",
           x = paste0("PC1 (", pct_var[1], "%)"),
           y = paste0("PC2 (", pct_var[2], "%)")
      )
  }, error = function(e) {
    message("[core_plots] Aitchison PCoA failed: ", e$message); NULL
  })

  # ── Panel B: Bray-Curtis PCoA (from ordination.txt) ──────────────────────
  bc_ord_file <- file.path(beta_dir, "bray_curtis", "ordination.txt")
  p_bc <- tryCatch({
    if (!file.exists(bc_ord_file)) stop("Missing ", bc_ord_file)
    ord <- parse_ordination(bc_ord_file)
    if (is.null(ord)) stop("parse_ordination returned NULL")
    df  <- ord$df
    pv  <- ord$pct_var

    if (!is.null(meta) && group_col %in% names(meta))
      df <- left_join(df, meta[, c("sample_id", group_col)], by = "sample_id")
    else
      df[[group_col]] <- "All"

    df[[group_col]] <- factor(df[[group_col]], levels = grp_levels)

    xlab <- if (!is.na(pv[1])) paste0("PC1 (", pv[1], "% var)") else "PC1"
    ylab <- if (!is.na(pv[2])) paste0("PC2 (", pv[2], "% var)") else "PC2"

    ggplot(df, aes(x = PC1, y = PC2, color = !!sym(group_col))) +
      geom_point(size = 3.8, alpha = 0.82) +
      stat_ellipse(level = 0.80, linetype = 2, linewidth = 0.7) +
      scale_color_manual(values = dark_colors, name = grp_display) +
      bold_theme +
      labs(title = "Bray-Curtis PCoA", x = xlab, y = ylab)
  }, error = function(e) {
    message("[core_plots] Bray-Curtis PCoA failed: ", e$message); NULL
  })

  # ── Panel C: within- vs between-group Bray-Curtis dissimilarity boxplot ───
  p_bc_box <- tryCatch({
    if (!file.exists(bray_dm_tsv)) stop("Missing ", bray_dm_tsv)
    dm_raw  <- read.table(bray_dm_tsv, header = TRUE, sep = "\t",
                          row.names = 1, check.names = FALSE)
    bc_mat  <- as.matrix(dm_raw)

    common_samp <- if (!is.null(meta)) intersect(rownames(bc_mat), meta$sample_id)
                   else rownames(bc_mat)
    bc_mat <- bc_mat[common_samp, common_samp, drop = FALSE]

    grp_vec <- if (!is.null(meta) && group_col %in% names(meta)) {
      setNames(as.character(meta[[group_col]]), meta$sample_id)[common_samp]
    } else setNames(rep("All", length(common_samp)), common_samp)

    n_s <- length(common_samp)
    box_rows <- vector("list", n_s * (n_s - 1L) / 2L)
    idx <- 1L
    for (i in seq_len(n_s - 1L)) {
      for (j in seq(i + 1L, n_s)) {
        si <- common_samp[i]; sj <- common_samp[j]
        box_rows[[idx]] <- data.frame(
          Dissimilarity = bc_mat[si, sj],
          Comparison    = if (grp_vec[si] == grp_vec[sj]) "Within-group" else "Between-group",
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
    bc_box_df <- do.call(rbind, box_rows)
    bc_box_df$Comparison <- factor(bc_box_df$Comparison,
                                    levels = c("Within-group", "Between-group"))

    box_p <- tryCatch(
      wilcox.test(Dissimilarity ~ Comparison, data = bc_box_df)$p.value,
      error = function(e) NA_real_
    )
    box_lab <- if (!is.na(box_p)) {
      if      (box_p < 0.001) "p < 0.001 ***"
      else if (box_p < 0.01)  paste0("p = ", signif(box_p, 2), " **")
      else if (box_p < 0.05)  paste0("p = ", signif(box_p, 2), " *")
      else                    paste0("p = ", signif(box_p, 2), " (NS)")
    } else ""

    box_pal <- setNames(get_palette(palette_name, 2), c("Within-group", "Between-group"))

    ggplot(bc_box_df, aes(x = Comparison, y = Dissimilarity, fill = Comparison)) +
      geom_boxplot(alpha = 0.85, outlier.shape = 21, outlier.size = 1.5,
                   color = "grey20", linewidth = 0.6) +
      geom_jitter(width = 0.12, size = 1.2, shape = 21,
                  color = "grey20", alpha = 0.50) +
      scale_fill_manual(values = box_pal, name = "Comparison") +
      bold_theme +
      labs(title   = "Bray-Curtis Dissimilarity",
           x       = NULL,
           y       = "Bray-Curtis Dissimilarity",
           caption = box_lab)
  }, error = function(e) {
    message("[core_plots] BC boxplot failed: ", e$message); NULL
  })

  # ── Assemble 3-panel overview (A / B / C) ────────────────────────────────
  panels_ov <- Filter(Negate(is.null), list(p_ait, p_bc, p_bc_box))
  if (length(panels_ov) == 0) {
    message("[core_plots] No beta overview panels produced")
    return(invisible(NULL))
  }

  overview <- patchwork::wrap_plots(panels_ov, nrow = 1) +
    patchwork::plot_annotation(tag_levels = "A") +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  save_png(overview, file.path(div_dir, "pcoa_overview.png"), w = 16, h = 6)
}

# =============================================================================
# Plot 6: Taxonomy phylum stacked barplot per sample (grouped by group_col)
# =============================================================================
plot_taxonomy_bar <- function(meta) {
  ft <- tryCatch({
    df <- read.table(feat_table_tsv, header = TRUE, sep = "\t",
                     comment.char = "", quote = "", check.names = FALSE,
                     skip = 1)  # skip '# Constructed from biom file'
    row.names(df) <- df[[1]]; df[, -1, drop = FALSE]
  }, error = function(e) { message("Skipping taxonomy bar: ", e$message); NULL })

  tax <- tryCatch({
    df <- read.table(taxonomy_tsv, header = TRUE, sep = "\t",
                     comment.char = "", quote = "", check.names = FALSE)
    names(df)[1:2] <- c("feature_id", "Taxon"); df
  }, error = function(e) { message("Skipping taxonomy bar: ", e$message); NULL })

  if (is.null(ft) || is.null(tax)) return(invisible(NULL))

  tax$Phylum <- sub(".*p__([^;]+).*", "\\1", tax$Taxon)
  tax$Phylum[!grepl("p__", tax$Taxon)] <- "Unclassified"
  tax$Phylum <- trimws(tax$Phylum)

  ft_long <- as.data.frame(as.table(as.matrix(ft))) %>%
    rename(feature_id = Var1, sample_id = Var2, count = Freq) %>%
    mutate(count = as.numeric(as.character(count)))

  phylum_df <- ft_long %>%
    left_join(tax[, c("feature_id", "Phylum")], by = "feature_id") %>%
    group_by(sample_id, Phylum) %>%
    summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
    group_by(sample_id) %>%
    mutate(rel_abund = count / sum(count, na.rm = TRUE)) %>%
    ungroup()

  # Keep top 12 phyla; lump rest as "Other"
  top12 <- phylum_df %>%
    group_by(Phylum) %>% summarise(m = mean(rel_abund), .groups = "drop") %>%
    arrange(desc(m)) %>% slice_head(n = 12) %>% pull(Phylum)

  phylum_df <- phylum_df %>%
    mutate(Phylum = ifelse(Phylum %in% top12, Phylum, "Other")) %>%
    group_by(sample_id, Phylum) %>%
    summarise(rel_abund = sum(rel_abund), .groups = "drop")

  # Aggregate to group level — mean relative abundance per group per phylum
  if (!is.null(meta) && group_col %in% names(meta)) {
    phylum_df <- left_join(phylum_df, meta[, c("sample_id", group_col)], by = "sample_id")
    plot_df   <- phylum_df %>%
      group_by(!!sym(group_col), Phylum) %>%
      summarise(rel_abund = mean(rel_abund, na.rm = TRUE), .groups = "drop")
    x_col   <- group_col
    x_label <- grp_display
  } else {
    plot_df <- phylum_df %>%
      group_by(sample_id, Phylum) %>%
      summarise(rel_abund = sum(rel_abund, na.rm = TRUE), .groups = "drop")
    x_col   <- "sample_id"
    x_label <- "Sample"
  }

  n_phy <- length(unique(plot_df$Phylum))
  pal   <- get_palette(palette_name, n_phy)

  p <- ggplot(plot_df, aes(x = !!sym(x_col), y = rel_abund, fill = Phylum)) +
    geom_col(position = "stack", width = 0.65) +
    scale_y_continuous(labels = percent_format()) +
    scale_fill_manual(values = setNames(pal, sort(unique(plot_df$Phylum)))) +
    theme_classic(base_size = 13) +
    theme(axis.text.x    = element_text(face = "bold", size = 12, color = "black"),
          axis.title      = element_text(face = "bold", size = 12),
          legend.position = "right",
          panel.grid      = element_blank()) +
    labs(x = x_label, y = "Mean Relative Abundance", fill = "Phylum")

  save_png(p, file.path(tax_dir, "taxonomy_barplot.png"), w = 10, h = 7)
}

# =============================================================================
# Plot 7: Top-20 phyla lollipop (mean % + prevalence)
# =============================================================================
plot_top_taxa <- function() {
  ft <- tryCatch({
    df <- read.table(feat_table_tsv, header = TRUE, sep = "\t",
                     comment.char = "", quote = "", check.names = FALSE,
                     skip = 1)  # skip '# Constructed from biom file'
    row.names(df) <- df[[1]]; df[, -1, drop = FALSE]
  }, error = function(e) NULL)

  tax <- tryCatch({
    df <- read.table(taxonomy_tsv, header = TRUE, sep = "\t",
                     comment.char = "", quote = "", check.names = FALSE)
    names(df)[1:2] <- c("feature_id", "Taxon"); df
  }, error = function(e) NULL)

  if (is.null(ft) || is.null(tax)) {
    message("Skipping top_taxa plot: missing ft or taxonomy")
    return(invisible(NULL))
  }

  tax$Phylum <- trimws(sub(".*p__([^;]+).*", "\\1", tax$Taxon))
  tax$Phylum[!grepl("p__", tax$Taxon)] <- "Unclassified"

  mat     <- as.matrix(ft)
  ft_rel  <- sweep(mat, 2, colSums(mat), "/")
  ft_long <- as.data.frame(as.table(ft_rel)) %>%
    rename(feature_id = Var1, sample_id = Var2, rel_abund = Freq) %>%
    mutate(rel_abund = as.numeric(as.character(rel_abund)))

  top_df <- ft_long %>%
    left_join(tax[, c("feature_id", "Phylum")], by = "feature_id") %>%
    group_by(Phylum) %>%
    summarise(mean_abund  = mean(rel_abund, na.rm = TRUE),
              prevalence  = mean(rel_abund > 0, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(desc(mean_abund)) %>% slice_head(n = 20)

  p <- ggplot(top_df, aes(y = reorder(Phylum, mean_abund))) +
    geom_segment(aes(x = 0, xend = mean_abund,
                     yend = reorder(Phylum, mean_abund)),
                 colour = "grey70", linewidth = 0.9) +
    geom_point(aes(x = mean_abund, size = prevalence, colour = mean_abund)) +
    scale_colour_gradient(low = "#2166ac", high = "#d73027",
                          labels = percent_format()) +
    scale_size_continuous(range = c(3, 10), labels = percent_format()) +
    scale_x_continuous(labels = percent_format()) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank()) +
    labs(title  = "Top-20 Phyla — Mean Relative Abundance",
         x      = "Mean Relative Abundance", y = NULL,
         colour = "Mean Abund.", size = "Prevalence")

  save_png(p, file.path(tax_dir, "top_taxa_abundance.png"), w = 11, h = 8)
}

# =============================================================================
# Execute all plots
# =============================================================================
meta <- read_metadata(meta_file)

message("── [1/4] DADA2 read-retention waterfall ──────────────────────────────")
plot_dada2()

message("── [2/4] Per-sample read depth ───────────────────────────────────────")
plot_read_depth(meta)

message("── [3/4] Taxonomy phylum barplot ─────────────────────────────────────")
plot_taxonomy_bar(meta)

message("── [4/4] Top-20 phyla lollipop ───────────────────────────────────────")
plot_top_taxa()

message("core_plots.R complete.")

