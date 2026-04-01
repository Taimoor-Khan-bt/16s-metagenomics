#!/usr/bin/env Rscript
# =============================================================================
# beta_stats.R — Robust Beta diversity: CLR-PCoA + PERMANOVA
# =============================================================================

suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(ape)
  library(dplyr)
})

# ── Arguments ─────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) stop("Usage: beta_stats.R <table_tsv> <metadata> <bray_dm> <wunifrac_dm> <group_col> <covariates_csv> <out_dir>")

table_tsv     <- args[1]
metadata_file <- args[2]
bray_file     <- args[3]
wunifrac_file <- args[4]
group_col     <- args[5]
covariates    <- if (nchar(args[6]) > 0) strsplit(args[6], ",")[[1]] else character(0)
out_dir       <- args[7]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load data (Defensive Loading) ─────────────────────────────────────────────
meta <- read.table(metadata_file, header = TRUE, sep = "\t",
                   comment.char = "", row.names = 1, check.names = FALSE, quote = "")

raw_tab <- read.table(table_tsv, header = TRUE, sep = "\t",
                      skip = 1, row.names = 1, check.names = FALSE, 
                      comment.char = "", quote = "")
tab <- t(raw_tab) 

common_samples <- intersect(rownames(meta), rownames(tab))
if (length(common_samples) < 2) {
    stop("Fewer than 2 common samples. Verify Sample IDs match exactly.")
}

meta <- meta[common_samples, , drop = FALSE]
tab  <- tab[common_samples, , drop = FALSE]

# ── Helper for Distance Matrices ──────────────────────────────────────────────
load_dm <- function(f, samples) {
  if (!file.exists(f)) return(NULL)
  dm <- read.table(f, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  common_dm <- intersect(samples, rownames(dm))
  if (length(common_dm) < 2) return(NULL)
  as.dist(dm[common_dm, common_dm])
}

# ── CLR Transformation & Aitchison PCoA ──────────────────────────────────────
tab_pseudo <- tab + 0.5
gm <- apply(tab_pseudo, 1, function(x) exp(mean(log(x))))
clr_tab <- log(tab_pseudo / gm)

ait_dist <- dist(clr_tab)
pcoa_ait <- pcoa(ait_dist)
pcoa_df  <- as.data.frame(pcoa_ait$vectors[, 1:2])
colnames(pcoa_df) <- c("PC1", "PC2")
pcoa_df$SampleID <- rownames(pcoa_df)
pcoa_df  <- merge(pcoa_df, tibble::rownames_to_column(meta, "SampleID"), by = "SampleID")

pct_var <- round(pcoa_ait$values$Relative_eig[1:2] * 100, 1)

# ── PERMANOVA (adonis2) ───────────────────────────────────────────────────────
permanova_results <- list()

run_permanova <- function(dist_obj, label, metadata) {
  valid_covs <- covariates[covariates %in% colnames(metadata)]
  all_terms <- c(group_col, valid_covs)
  formula_str <- paste("dist_obj ~", paste(all_terms, collapse = " + "))
  
  tryCatch({
    perm <- adonis2(as.formula(formula_str), data = metadata, permutations = 999)
    perm_df <- as.data.frame(perm)
    perm_df$term   <- rownames(perm_df)
    perm_df$metric <- label
    perm_df
  }, error = function(e) {
    message("PERMANOVA failed for ", label, ": ", e$message)
    NULL
  })
}

permanova_results[["Aitchison"]] <- run_permanova(ait_dist, "Aitchison", meta)

bray_dist <- load_dm(bray_file, common_samples)
if (!is.null(bray_dist)) {
  permanova_results[["Bray-Curtis"]] <- run_permanova(bray_dist, "Bray-Curtis", meta)
}

wu_dist <- load_dm(wunifrac_file, common_samples)
if (!is.null(wu_dist)) {
  permanova_results[["Weighted-UniFrac"]] <- run_permanova(wu_dist, "Weighted-UniFrac", meta)
}

# ── FIXED: Combine and Save Stats ─────────────────────────────────────────────
# We actually need to create the perm_df_all variable from the list
perm_df_all <- do.call(rbind, Filter(Negate(is.null), permanova_results))

if (!is.null(perm_df_all) && nrow(perm_df_all) > 0) {
  write.table(perm_df_all, file.path(out_dir, "permanova_results.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message("Saved: permanova_results.tsv")
} else {
  write.table(data.frame(Status="No_Results"), file.path(out_dir, "permanova_results.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# ── Plotting ──────────────────────────────────────────────────────────────────
# Dark color palette (consistent with alpha diversity plots)
grp_levels_b <- sort(unique(as.character(meta[[group_col]])))
dark_pal_b   <- c("#1B4F72", "#922B21", "#1D8348", "#6C3483", "#784212",
                  "#0E6655", "#4A235A", "#1A5276", "#7D6608", "#212F3D")
dark_colors_b <- setNames(dark_pal_b[seq_along(grp_levels_b)], grp_levels_b)

bold_theme_b <- theme_bw(base_size = 13) +
  theme(
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold", size = 12),
    legend.text      = element_text(face = "bold", size = 11),
    axis.title       = element_text(face = "bold", size = 12),
    axis.text        = element_text(face = "bold", size = 11, color = "black"),
    plot.title       = element_text(face = "bold", size = 13),
    panel.grid.minor = element_blank()
  )

# ── Panel A: Aitchison PCoA (no sample labels) ────────────────────────────
p_ait <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = !!sym(group_col))) +
  geom_point(size = 3.8, alpha = 0.82) +
  stat_ellipse(level = 0.80, linetype = 2, linewidth = 0.7) +
  scale_color_manual(values = dark_colors_b, name = group_col) +
  bold_theme_b +
  labs(title = "Aitchison PCoA (CLR)",
       x = paste0("PC1 (", pct_var[1], "%)"),
       y = paste0("PC2 (", pct_var[2], "%)"))

# ── Panels B & C: Bray-Curtis PCoA + within/between dissimilarity boxplot ──
p_bc     <- NULL
p_bc_box <- NULL

if (!is.null(bray_dist)) {
  bc_pcoa <- pcoa(bray_dist)
  bc_df   <- as.data.frame(bc_pcoa$vectors[, 1:2])
  colnames(bc_df) <- c("PC1", "PC2")
  bc_df$SampleID <- rownames(bc_df)
  bc_df <- merge(bc_df, tibble::rownames_to_column(meta, "SampleID"), by = "SampleID")
  bc_pct <- round(bc_pcoa$values$Relative_eig[1:2] * 100, 1)

  # Panel B: Bray-Curtis PCoA (no sample labels)
  p_bc <- ggplot(bc_df, aes(x = PC1, y = PC2, color = !!sym(group_col))) +
    geom_point(size = 3.8, alpha = 0.82) +
    stat_ellipse(level = 0.80, linetype = 2, linewidth = 0.7) +
    scale_color_manual(values = dark_colors_b, name = group_col) +
    bold_theme_b +
    labs(title = "Bray-Curtis PCoA",
         x = paste0("PC1 (", bc_pct[1], "%)"),
         y = paste0("PC2 (", bc_pct[2], "%)"))

  # Panel C: within-group vs. between-group BC dissimilarity boxplot
  bc_mat  <- as.matrix(bray_dist)
  grp_vec <- meta[[group_col]]
  names(grp_vec) <- rownames(meta)
  grp_vec <- grp_vec[rownames(bc_mat)]   # align sample order
  n_samp  <- nrow(bc_mat)

  box_rows <- vector("list", n_samp * (n_samp - 1) / 2)
  idx <- 1L
  for (i in seq_len(n_samp - 1)) {
    for (j in seq(i + 1, n_samp)) {
      si <- rownames(bc_mat)[i]; sj <- rownames(bc_mat)[j]
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

  box_pal <- c("Within-group" = "#1B4F72", "Between-group" = "#922B21")

  p_bc_box <- ggplot(bc_box_df, aes(x = Comparison, y = Dissimilarity,
                                     fill = Comparison)) +
    geom_boxplot(alpha = 0.85, outlier.shape = 21, outlier.size = 1.5,
                 color = "grey20", linewidth = 0.6) +
    geom_jitter(width = 0.12, size = 1.2, shape = 21,
                color = "grey20", alpha = 0.50) +
    scale_fill_manual(values = box_pal, name = "Comparison") +
    bold_theme_b +
    labs(title   = "Bray-Curtis Dissimilarity",
         x       = NULL,
         y       = "Bray-Curtis Dissimilarity",
         caption = box_lab)
}

# ── Assemble multi-panel with patchwork (A / B / C) ───────────────────────
panels <- Filter(Negate(is.null), list(p_ait, p_bc, p_bc_box))
n_p    <- length(panels)

if (n_p > 1 && requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  combined_beta <- patchwork::wrap_plots(panels, nrow = 1) +
    patchwork::plot_annotation(
      tag_levels = "A",
      title      = paste("Beta Diversity by", group_col),
      theme      = theme(plot.title = element_text(face = "bold", size = 15,
                                                   hjust = 0.5))
    ) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
} else {
  combined_beta <- NULL
}

# ── Save outputs ──────────────────────────────────────────────────────────
pdf(file.path(out_dir, "pcoa_plots.pdf"), width = 16, height = 6)
if (!is.null(combined_beta)) {
  print(combined_beta)
} else {
  for (pl in panels) print(pl)
}
dev.off()
message("Saved: pcoa_plots.pdf")

png_plot <- if (!is.null(combined_beta)) combined_beta else p_ait
ggsave(file.path(out_dir, "pcoa_plots.png"),
       plot = png_plot,
       width = if (!is.null(combined_beta)) 16 else 8,
       height = 6, units = "in", dpi = 600)
message("Saved: pcoa_plots.png (600 DPI)")