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
grp_display   <- tools::toTitleCase(gsub("_", " ", group_col))

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

run_permanova_full <- function(dist_obj, label, metadata) {
  valid_covs  <- covariates[covariates %in% colnames(metadata)]
  all_terms   <- c(group_col, valid_covs)
  formula_str <- paste("dist_obj ~", paste(all_terms, collapse = " + "))

  result <- tryCatch({
    perm <- adonis2(as.formula(formula_str), data = metadata, permutations = 999)
    perm_df         <- as.data.frame(perm)
    perm_df$term    <- rownames(perm_df)
    perm_df$metric  <- label
    perm_df
  }, error = function(e) {
    message("PERMANOVA failed for ", label, ": ", e$message)
    NULL
  })

  # PERMDISP: test homogeneity of group dispersions
  betadisper_F <- NA_real_; betadisper_p <- NA_real_; interpretation <- NA_character_
  bd <- tryCatch(
    betadisper(dist_obj, metadata[[group_col]]),
    error = function(e) { message("betadisper failed for ", label); NULL }
  )
  if (!is.null(bd)) {
    pt <- tryCatch(permutest(bd, permutations = 999), error = function(e) NULL)
    if (!is.null(pt)) {
      betadisper_F <- as.numeric(pt$statistic[1])
      betadisper_p <- as.numeric(pt$tab["Groups", "Pr(>F)"])
      if (!is.null(result) && !is.null(betadisper_p) && !is.na(betadisper_p)) {
        perm_row <- result[!is.na(result$term) & result$term == group_col, , drop = FALSE]
        perm_p   <- if (nrow(perm_row) > 0) perm_row[["Pr(>F)"]][1] else NA
        interpretation <- if (!is.na(perm_p) && perm_p < 0.05 && betadisper_p < 0.05) "Both"
          else if (!is.na(perm_p) && perm_p < 0.05)  "Location"
          else if (betadisper_p < 0.05)              "Dispersion"
          else                                        "Insufficient data"
      }
    }
  }

  # Attach PERMDISP results to returned data.frame and cache bd object
  if (!is.null(result)) {
    result$betadisper_F  <- betadisper_F
    result$betadisper_p  <- betadisper_p
    result$interpretation <- interpretation
  }
  attr(result, "bd") <- bd
  result
}
# Alias so get_perm_label still works without change
run_permanova <- run_permanova_full

permanova_results[["Aitchison"]] <- run_permanova_full(ait_dist, "Aitchison", meta)

bray_dist <- load_dm(bray_file, common_samples)
if (!is.null(bray_dist)) {
  permanova_results[["Bray-Curtis"]] <- run_permanova_full(bray_dist, "Bray-Curtis", meta)
}

wu_dist <- load_dm(wunifrac_file, common_samples)
if (!is.null(wu_dist)) {
  permanova_results[["Weighted-UniFrac"]] <- run_permanova_full(wu_dist, "Weighted-UniFrac", meta)
}

# ── FIXED: Combine and Save Stats ─────────────────────────────────────────────dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)# We actually need to create the perm_df_all variable from the list
perm_df_all <- do.call(rbind, Filter(Negate(is.null), permanova_results))

if (!is.null(perm_df_all) && nrow(perm_df_all) > 0) {
  write.table(perm_df_all, file.path(out_dir, "permanova_results.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message("Saved: permanova_results.tsv")
} else {
  write.table(data.frame(Status="No_Results"), file.path(out_dir, "permanova_results.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
# ── Pairwise PERMANOVA (3+ groups only) ──────────────────────────────────
# (saved as a separate TSV; does not change any permanova_results list entries)
grp_levels_pw <- levels(as.factor(meta[[group_col]]))
if (length(grp_levels_pw) > 2) {
  pair_combns <- combn(grp_levels_pw, 2, simplify = FALSE)
  # Build list of dist matrices to test pairwise
  pw_dists <- Filter(Negate(is.null), list(
    Aitchison       = ait_dist,
    `Bray-Curtis`   = bray_dist,
    `Weighted-UniFrac` = wu_dist
  ))
  pw_rows <- lapply(names(pw_dists), function(dm_name) {
    dm <- pw_dists[[dm_name]]
    lapply(pair_combns, function(pair) {
      idx <- rownames(meta)[meta[[group_col]] %in% pair]
      if (length(idx) < 2) return(NULL)
      meta_sub <- meta[idx, , drop = FALSE]
      dm_sub   <- as.dist(as.matrix(dm)[idx, idx])
      tryCatch({
        perm_pw <- adonis2(dm_sub ~ meta_sub[[group_col]], permutations = 999)
        row     <- as.data.frame(perm_pw)[rownames(as.data.frame(perm_pw)) == "meta_sub[[group_col]]", , drop = FALSE]
        data.frame(
          metric     = dm_name,
          comparison = paste(pair, collapse = " vs "),
          R2         = row$R2[1],
          p_value    = row[["Pr(>F)"]][1],
          stringsAsFactors = FALSE
        )
      }, error = function(e) NULL)
    })
  })
  pw_flat <- do.call(rbind, Filter(Negate(is.null), unlist(pw_rows, recursive = FALSE)))
  if (!is.null(pw_flat) && nrow(pw_flat) > 0) {
    # BH correction within each distance metric
    pw_flat$p_adj_BH <- ave(pw_flat$p_value, pw_flat$metric,
                            FUN = function(x) p.adjust(x, method = "BH"))
    write.table(pw_flat, file.path(out_dir, "permanova_pairwise.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    message("Saved: permanova_pairwise.tsv")
  } else {
    write.table(data.frame(Note = "Pairwise PERMANOVA returned no results"),
                file.path(out_dir, "permanova_pairwise.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
} else {
  write.table(data.frame(Note = "2-group comparison — pairwise PERMANOVA not applicable"),
              file.path(out_dir, "permanova_pairwise.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message("Saved: permanova_pairwise.tsv (placeholder — 2 groups)")
}
# ── PERMANOVA label helper ────────────────────────────────────────────────
get_perm_label <- function(perm_list, label) {
  df <- perm_list[[label]]
  if (is.null(df)) return(NULL)
  row <- df[!is.na(df$term) & df$term == group_col, , drop = FALSE]
  if (nrow(row) == 0) return(NULL)
  r2  <- if ("R2" %in% names(row)) signif(row$R2[1], 2) else "?"
  pv  <- if ("Pr(>F)" %in% names(row)) row[["Pr(>F)"]][1] else NA
  sig <- if (!is.na(pv)) {
    if (pv < 0.001) "p < 0.001" else paste0("p = ", signif(pv, 2))
  } else ""
  paste0("PERMANOVA: R² = ", r2, ", ", sig)
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
  scale_color_manual(values = dark_colors_b, name = grp_display) +
  bold_theme_b +
  labs(title   = "Aitchison PCoA (CLR)",
       x       = paste0("PC1 (", pct_var[1], "%)"),
       y       = paste0("PC2 (", pct_var[2], "%)"),
       caption = get_perm_label(permanova_results, "Aitchison"))

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
    scale_color_manual(values = dark_colors_b, name = grp_display) +
    bold_theme_b +
    labs(title   = "Bray-Curtis PCoA",
         x       = paste0("PC1 (", bc_pct[1], "%)"),
         y       = paste0("PC2 (", bc_pct[2], "%)"),
         caption = get_perm_label(permanova_results, "Bray-Curtis"))

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

  # ── Panel D: PERMDISP centroid-distance boxplot ────────────────────────
  p_permdisp <- NULL
  bd_bc <- attr(permanova_results[["Bray-Curtis"]], "bd")
  if (!is.null(bd_bc)) {
    disp_df <- data.frame(
      SampleID     = names(bd_bc$distances),
      CentroidDist = as.numeric(bd_bc$distances),
      stringsAsFactors = FALSE
    )
    disp_df <- merge(disp_df,
                     data.frame(SampleID = rownames(meta),
                                Group    = meta[[group_col]],
                                stringsAsFactors = FALSE),
                     by = "SampleID")
    disp_df$Group <- factor(disp_df$Group, levels = grp_levels_b)

    bd_bc_perm <- tryCatch(permutest(bd_bc, permutations = 999), error = function(e) NULL)
    disp_cap <- if (!is.null(bd_bc_perm)) {
      pv <- as.numeric(bd_bc_perm$tab["Groups", "Pr(>F)"])
      sig <- if (pv < 0.001) "p < 0.001 ***"
             else if (pv < 0.01)  paste0("p = ", signif(pv, 2), " **")
             else if (pv < 0.05)  paste0("p = ", signif(pv, 2), " *")
             else                 paste0("p = ", signif(pv, 2), " (NS)")
      paste0("PERMDISP: ", sig)
    } else ""

    p_permdisp <- ggplot(disp_df, aes(x = Group, y = CentroidDist, fill = Group)) +
      geom_boxplot(alpha = 0.85, outlier.shape = 21, outlier.size = 1.5,
                   color = "grey20", linewidth = 0.6) +
      geom_jitter(width = 0.12, size = 1.2, shape = 21,
                  color = "grey20", alpha = 0.50) +
      scale_fill_manual(values = dark_colors_b, name = grp_display) +
      bold_theme_b +
      labs(title   = "Distance to Group Centroid (Bray-Curtis)",
           x       = NULL,
           y       = "Distance to Centroid",
           caption = disp_cap)
  }
}

# ── Panel C: Weighted UniFrac PCoA ────────────────────────────────────
p_wu <- NULL
if (!is.null(wu_dist)) {
  wu_pcoa <- pcoa(wu_dist)
  wu_df   <- as.data.frame(wu_pcoa$vectors[, 1:2])
  colnames(wu_df) <- c("PC1", "PC2")
  wu_df$SampleID <- rownames(wu_df)
  wu_df <- merge(wu_df, tibble::rownames_to_column(meta, "SampleID"), by = "SampleID")
  wu_pct <- round(wu_pcoa$values$Relative_eig[1:2] * 100, 1)

  p_wu <- ggplot(wu_df, aes(x = PC1, y = PC2, color = !!sym(group_col))) +
    geom_point(size = 3.8, alpha = 0.82) +
    stat_ellipse(level = 0.80, linetype = 2, linewidth = 0.7) +
    scale_color_manual(values = dark_colors_b, name = grp_display) +
    bold_theme_b +
    labs(title   = "Weighted UniFrac PCoA",
         x       = paste0("PC1 (", wu_pct[1], "%)"),
         y       = paste0("PC2 (", wu_pct[2], "%)"),
         caption = get_perm_label(permanova_results, "Weighted-UniFrac"))
}

# ── Assemble multi-panel with patchwork (A / B / C / D) ─────────────────
panels <- Filter(Negate(is.null), list(p_ait, p_bc, p_wu, p_bc_box, p_permdisp))
n_p    <- length(panels)

if (n_p > 1 && requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  n_rows         <- if (n_p <= 3) 1L else 2L
  combined_beta  <- patchwork::wrap_plots(panels, nrow = n_rows) +
    patchwork::plot_annotation(tag_levels = "A") +
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
png_w    <- if (n_p <= 3) 16 else 12
ggsave(file.path(out_dir, "pcoa_plots.png"),
       plot   = png_plot,
       width  = png_w,
       height = if (n_p <= 3) 6 else 12, units = "in", dpi = 600)
message("Saved: pcoa_plots.png (600 DPI)")

# Save PERMDISP centroid-distance boxplot as its own PNG
if (!is.null(p_permdisp)) {
  ggsave(file.path(out_dir, "permdisp_plots.png"),
         plot = p_permdisp, width = 8, height = 6, units = "in", dpi = 600)
  message("Saved: permdisp_plots.png (600 DPI)")
} else {
  # Write placeholder PNG so Snakemake output is always satisfied
  grDevices::png(file.path(out_dir, "permdisp_plots.png"),
                 width = 400, height = 200, res = 96)
  graphics::plot.new()
  graphics::title("PERMDISP not available (betadisper failed or distance matrix missing)")
  grDevices::dev.off()
  message("Saved: permdisp_plots.png (placeholder)")
}