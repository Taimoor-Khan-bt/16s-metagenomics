#!/usr/bin/env Rscript
# =============================================================================
# alpha_stats.R — Alpha diversity statistical analysis
# =============================================================================
# Usage: Rscript alpha_stats.R <metadata_tsv> <alpha_dir> <group_col>
#                              <covariates_csv> <out_dir>
#
# Required R packages: ggplot2, ggpubr, dplyr, tidyr, broom
# Install: mamba install -n qiime2 -c conda-forge r-ggplot2 r-ggpubr r-dplyr r-tidyr r-broom
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# ── Arguments ──────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("Usage: alpha_stats.R <metadata> <alpha_dir> <group_col> <covariates_csv> <out_dir>")

metadata_file <- args[1]
alpha_dir     <- args[2]
group_col     <- args[3]
covariates    <- if (nchar(args[4]) > 0) strsplit(args[4], ",")[[1]] else character(0)
out_dir       <- args[5]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load data ─────────────────────────────────────────────────────────────────
meta <- read.table(metadata_file, header = TRUE, sep = "\t",
                   comment.char = "", row.names = 1, check.names = FALSE)

metrics <- c("faith_pd", "shannon", "evenness", "observed_features")
alpha_list <- list()

for (m in metrics) {
  f <- file.path(alpha_dir, m, "alpha-diversity.tsv")
  if (file.exists(f)) {
    df <- read.table(f, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
    alpha_list[[m]] <- df[, 1, drop = FALSE]
    colnames(alpha_list[[m]]) <- m
    message("Loaded: ", m, " (", nrow(df), " samples)")
  } else {
    message("Skipping (not found): ", f)
  }
}

if (length(alpha_list) == 0) stop("No alpha diversity files found in: ", alpha_dir)

# Merge all metrics with metadata
alpha_df <- do.call(cbind, alpha_list)
combined <- merge(alpha_df, meta, by = "row.names", all.x = TRUE)
rownames(combined) <- combined$Row.names
combined$Row.names <- NULL

# ── Statistical tests ─────────────────────────────────────────────────────────
stats_rows <- list()

for (m in names(alpha_list)) {
  if (!m %in% colnames(combined)) next
  y   <- combined[[m]]
  grp <- as.factor(combined[[group_col]])
  lvls <- levels(grp)
  n_per_grp <- table(grp)

  # Wilcoxon rank-sum (2 groups) or Kruskal-Wallis (>2 groups)
  test_name <- if (length(lvls) == 2) "Wilcoxon" else "Kruskal-Wallis"
  p_val <- tryCatch({
    if (length(lvls) == 2) {
      wilcox.test(y ~ grp)$p.value
    } else {
      kruskal.test(y ~ grp)$p.value
    }
  }, error = function(e) { message("Test failed for ", m, ": ", e$message); NA })

  stats_rows[[m]] <- data.frame(
    metric    = m,
    test      = test_name,
    p_value   = p_val,
    n_samples = nrow(combined),
    stringsAsFactors = FALSE
  )

  # GLM with covariates
  cov_cols <- covariates[covariates %in% colnames(combined)]
  all_terms <- c(group_col, cov_cols)
  formula_str <- paste(m, "~", paste(all_terms, collapse = " + "))

  tryCatch({
    # Convert group col to factor in combined for GLM
    combined[[group_col]] <- as.factor(combined[[group_col]])
    glm_fit <- lm(as.formula(formula_str), data = combined)
    coef_df <- as.data.frame(summary(glm_fit)$coefficients)
    coef_df$term   <- rownames(coef_df)
    coef_df$metric <- m
    colnames(coef_df) <- c("estimate", "std_error", "t_value", "p_value", "term", "metric")
    write.table(coef_df[, c("metric","term","estimate","std_error","t_value","p_value")],
                file.path(out_dir, paste0(m, "_glm.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    message("GLM saved: ", m)
  }, error = function(e) message("GLM failed for ", m, ": ", e$message))
}

# FDR correction across metrics
stats_df <- do.call(rbind, stats_rows)
stats_df$p_adj_BH <- p.adjust(stats_df$p_value, method = "BH")
write.table(stats_df, file.path(out_dir, "alpha_statistics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Saved: alpha_statistics.tsv")

# ── Plots ─────────────────────────────────────────────────────────────────────
# Publication metrics to display as bars: Shannon, Simpson, Chao1 (if loaded).
# Falls back to all available metrics if none of the preferred three are present.
pub_metrics <- intersect(c("shannon", "simpson", "chao1"), names(alpha_list))
if (length(pub_metrics) == 0) pub_metrics <- names(alpha_list)

metric_labels <- c(
  shannon           = "Shannon Index",
  simpson           = "Simpson Index",
  chao1             = "Chao1 Richness",
  faith_pd          = "Faith's PD",
  evenness          = "Pielou's Evenness",
  observed_features = "Observed ASVs"
)

# Dark, publication-quality color palette (colorblind-accessible)
grp_levels  <- sort(unique(as.character(combined[[group_col]])))
dark_pal    <- c("#1B4F72", "#922B21", "#1D8348", "#6C3483", "#784212",
                 "#0E6655", "#4A235A", "#1A5276", "#7D6608", "#212F3D")
dark_colors <- setNames(dark_pal[seq_along(grp_levels)], grp_levels)

# Build one bar chart per metric (mean ± SE + individual sample dots)
panel_plots <- list()
for (m in pub_metrics) {
  if (!m %in% colnames(combined)) {
    message("Metric not found in data, skipping: ", m)
    next
  }
  df_m         <- combined[, c(m, group_col), drop = FALSE]
  colnames(df_m)[1] <- "value"
  df_m[[group_col]] <- factor(df_m[[group_col]], levels = grp_levels)

  # Mean ± SE summary
  summ <- df_m %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::summarise(
      mean_val = mean(value, na.rm = TRUE),
      se_val   = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
      .groups  = "drop"
    )

  # Adjusted p-value label from stats table
  p_row <- stats_df[stats_df$metric == m, , drop = FALSE]
  p_lab <- if (nrow(p_row) > 0 && !is.na(p_row$p_adj_BH[1])) {
    pv <- p_row$p_adj_BH[1]
    if      (pv < 0.001) "q < 0.001 ***"
    else if (pv < 0.01)  paste0("q = ", signif(pv, 2), " **")
    else if (pv < 0.05)  paste0("q = ", signif(pv, 2), " *")
    else                 paste0("q = ", signif(pv, 2), " (NS)")
  } else ""

  p <- ggplot(summ, aes(x = .data[[group_col]], y = mean_val,
                        fill = .data[[group_col]])) +
    geom_col(width = 0.55, alpha = 0.90, color = "grey20", linewidth = 0.45) +
    geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val),
                  width = 0.18, linewidth = 0.9, color = "grey20") +
    geom_jitter(data = df_m,
                aes(x = .data[[group_col]], y = value, fill = .data[[group_col]]),
                width = 0.12, size = 2.4, shape = 21,
                color = "grey20", alpha = 0.75, inherit.aes = FALSE) +
    scale_fill_manual(values = dark_colors, name = group_col) +
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
    labs(
      x       = NULL,
      y       = if (!is.na(metric_labels[m])) metric_labels[m] else m,
      title   = if (!is.na(metric_labels[m])) metric_labels[m] else m,
      caption = p_lab
    )
  panel_plots[[m]] <- p
}

# Assemble multi-panel (A / B / C labels) with patchwork; facet fallback if absent
n_panels <- length(panel_plots)
fig_w    <- max(4.5 * n_panels, 8)

if (n_panels > 0) {
  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    combined_plot <- patchwork::wrap_plots(panel_plots, nrow = 1) +
      patchwork::plot_annotation(
        tag_levels = "A",
        title      = paste("Alpha Diversity by", group_col),
        theme      = theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5))
      ) +
      patchwork::plot_layout(guides = "collect") &
      theme(legend.position = "bottom")
  } else {
    message("patchwork not installed — using facet fallback for alpha plots")
    plot_long_pub <- pivot_longer(
      combined[, c(pub_metrics[pub_metrics %in% colnames(combined)], group_col)],
      cols = pub_metrics[pub_metrics %in% colnames(combined)],
      names_to = "metric", values_to = "value"
    )
    plot_long_pub$metric_label <- vapply(
      plot_long_pub$metric,
      function(x) if (!is.na(metric_labels[x])) metric_labels[x] else x,
      character(1)
    )
    combined_plot <- ggplot(plot_long_pub,
                            aes(x = .data[[group_col]], y = value,
                                fill = .data[[group_col]])) +
      stat_summary(fun = mean, geom = "bar", width = 0.55, alpha = 0.90,
                   color = "grey20") +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.18,
                   linewidth = 0.9) +
      geom_jitter(width = 0.12, size = 2.0, shape = 21,
                  color = "grey20", alpha = 0.75) +
      facet_wrap(~ metric_label, scales = "free_y", nrow = 1) +
      scale_fill_manual(values = dark_colors, name = group_col) +
      theme_classic(base_size = 13) +
      theme(
        legend.position = "bottom",
        legend.title    = element_text(face = "bold", size = 12),
        legend.text     = element_text(face = "bold", size = 11),
        axis.title      = element_text(face = "bold", size = 12),
        axis.text       = element_text(face = "bold", size = 11, color = "black"),
        strip.text      = element_text(face = "bold", size = 12)
      ) +
      labs(x = NULL, y = "Alpha Diversity",
           title = paste("Alpha Diversity by", group_col))
  }

  # ── PDF: Page 1 = bar chart, Page 2 = GLM forest plot ─────────────────────
  pdf(file.path(out_dir, "alpha_plots.pdf"), width = fig_w, height = 6)
  print(combined_plot)

  # ── Page 2: GLM Forest Plot (retained from original) ──────────────────────
  glm_rows <- list()
  for (m in names(alpha_list)) {
    glm_file <- file.path(out_dir, paste0(m, "_glm.tsv"))
    if (file.exists(glm_file)) {
      gdf <- read.table(glm_file, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE, check.names = FALSE)
      gdf <- gdf[grepl(paste0("^", group_col), gdf$term), , drop = FALSE]
      if (nrow(gdf) > 0) glm_rows[[m]] <- gdf
    }
  }

  if (length(glm_rows) > 0) {
    forest_df            <- do.call(rbind, glm_rows)
    forest_df$ci_lo      <- forest_df$estimate - 1.96 * forest_df$std_error
    forest_df$ci_hi      <- forest_df$estimate + 1.96 * forest_df$std_error
    forest_df$sig        <- ifelse(!is.na(forest_df$p_value) & forest_df$p_value < 0.05,
                                   "p < 0.05", "p \u2265 0.05")
    forest_df$label      <- paste0("\u03b2=", round(forest_df$estimate, 3),
                                    "  p=", signif(forest_df$p_value, 3))
    forest_df$term_short <- sub(paste0("^", group_col), "", forest_df$term)
    forest_df$row_id     <- paste0(forest_df$metric, "  (", forest_df$term_short, ")")

    pf <- ggplot(forest_df, aes(x = estimate, y = row_id, color = sig)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.25, linewidth = 0.8) +
      geom_point(size = 3.5) +
      geom_text(aes(label = label), nudge_y = 0.35, size = 3, show.legend = FALSE) +
      scale_color_manual(values = c("p < 0.05" = "#D32F2F", "p \u2265 0.05" = "#1565C0"),
                         name = "Significance") +
      theme_bw(base_size = 12) +
      theme(
        legend.position  = "bottom",
        legend.title     = element_text(face = "bold", size = 12),
        legend.text      = element_text(face = "bold", size = 11),
        axis.title       = element_text(face = "bold", size = 12),
        axis.text.y      = element_text(size = 11),
        panel.grid.minor = element_blank()
      ) +
      labs(x        = paste0("GLM Estimate (ref = reference level of ", group_col, ")"),
           y        = "Alpha metric  (group contrast)",
           title    = paste("Multivariable GLM \u2014 Effect of", group_col,
                            "on Alpha Diversity"),
           subtitle = paste("Adjusted for:",
                            if (length(covariates) > 0) paste(covariates, collapse = ", ")
                            else "no covariates"))
    print(pf)
    message("GLM forest plot added to alpha_plots.pdf")
  } else {
    message("No GLM files found — skipping forest plot page")
  }

  dev.off()
  message("Saved: alpha_plots.pdf")

  # ── PNG (600 DPI) — bar chart panel only ──────────────────────────────────
  ggsave(file.path(out_dir, "alpha_plots.png"),
         plot = combined_plot,
         width = fig_w, height = 6, units = "in", dpi = 600)
  message("Saved: alpha_plots.png (600 DPI)")

} else {
  message("No metrics available for plotting — writing placeholder files")
  pdf(file.path(out_dir, "alpha_plots.pdf"), width = 6, height = 4)
  plot.new(); title("No alpha diversity metrics available")
  dev.off()
  grDevices::png(file.path(out_dir, "alpha_plots.png"), width = 600, height = 400)
  plot.new(); title("No alpha diversity metrics available")
  dev.off()
}
