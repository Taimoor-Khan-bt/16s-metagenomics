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
pdf(file.path(out_dir, "alpha_plots.pdf"), width = 5 * length(names(alpha_list)), height = 6)

# Reshape for faceting
plot_df <- combined[, c(names(alpha_list), group_col)]
plot_long <- pivot_longer(plot_df, cols = names(alpha_list),
                          names_to = "metric", values_to = "value")

p <- ggplot(plot_long, aes_string(x = group_col, y = "value", fill = group_col)) +
  geom_violin(alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 1.5, alpha = 0.8) +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5, shape = 21) +
  facet_wrap(~ metric, scales = "free_y", ncol = 4) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", strip.background = element_rect(fill = "#E8EAF6")) +
  labs(x = group_col, y = "Alpha Diversity Value",
       title = paste("Alpha Diversity by", group_col))

# Add p-values from stats_df as subtitle annotations
annot <- stats_df[, c("metric", "p_value", "p_adj_BH")]
annot$label <- paste0("p=", signif(annot$p_value, 3),
                      " (adj=", signif(annot$p_adj_BH, 3), ")")
p <- p + geom_text(data = annot,
                   aes(x = -Inf, y = Inf, label = label, group = NULL, fill = NULL),
                   hjust = -0.05, vjust = 1.5, size = 3.5, color = "black",
                   inherit.aes = FALSE)

print(p)

# ── Page 2: GLM Forest Plot ───────────────────────────────────────────────────
# Collect GLM coefficients from per-metric TSV files and plot estimate ± 95% CI
# for the primary group term across all alpha diversity metrics.
glm_rows <- list()
for (m in names(alpha_list)) {
  glm_file <- file.path(out_dir, paste0(m, "_glm.tsv"))
  if (file.exists(glm_file)) {
    gdf <- read.table(glm_file, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, check.names = FALSE)
    # Keep only rows for the primary group variable (exclude Intercept, covariates)
    gdf <- gdf[grepl(paste0("^", group_col), gdf$term), , drop = FALSE]
    if (nrow(gdf) > 0) glm_rows[[m]] <- gdf
  }
}

if (length(glm_rows) > 0) {
  forest_df <- do.call(rbind, glm_rows)
  forest_df$ci_lo   <- forest_df$estimate - 1.96 * forest_df$std_error
  forest_df$ci_hi   <- forest_df$estimate + 1.96 * forest_df$std_error
  forest_df$sig     <- ifelse(!is.na(forest_df$p_value) & forest_df$p_value < 0.05,
                              "p < 0.05", "p ≥ 0.05")
  forest_df$label   <- paste0("β=", round(forest_df$estimate, 3),
                               "  p=", signif(forest_df$p_value, 3))
  # Shorten term label: strip group_col prefix for display
  forest_df$term_short <- sub(paste0("^", group_col), "", forest_df$term)
  forest_df$row_id  <- paste0(forest_df$metric, "  (", forest_df$term_short, ")")

  pf <- ggplot(forest_df, aes(x = estimate, y = row_id, color = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.25, linewidth = 0.8) +
    geom_point(size = 3.5) +
    geom_text(aes(label = label), nudge_y = 0.35, size = 3, show.legend = FALSE) +
    scale_color_manual(values = c("p < 0.05" = "#D32F2F", "p ≥ 0.05" = "#1565C0"),
                       name = "Significance") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          axis.text.y   = element_text(size = 11),
          panel.grid.minor = element_blank()) +
    labs(x     = paste0("GLM Estimate (ref = reference level of ", group_col, ")"),
         y     = "Alpha metric  (group contrast)",
         title = paste("Multivariable GLM — Effect of", group_col, "on Alpha Diversity"),
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
