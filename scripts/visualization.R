#!/usr/bin/env Rscript
# Unified visualization script for metagenomics (currently 16S)
# Plots are controlled via cfg$plots$enable toggles.
# This script generates publication-ready visualizations using best practices:
# - Consistent theming with customizable base size, fonts, and palettes.
# - High-resolution outputs (TIFF/PDF) with configurable DPI.
# - Improved error handling, metadata validation, and fallback mechanisms.
# - Enhanced plot aesthetics: subtitles, captions, better legends, axis formatting.
# - Refactored into modular functions for clarity and maintainability.
# - Expanded tree plotting with annotations and multiple layouts.
# - Upgraded per-sample multipanels using cowplot for better composition and TIFF output.
# - Added support for additional ordination methods and statistical annotations (e.g., ellipses with checks).
# - Configurable top taxa, ranks, and other parameters via cfg.

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(cowplot)  # Replaced gridExtra for better multipanel handling
  library(reshape2)
  library(viridisLite)
  library(scales)
  library(ggtree)   # For enhanced tree visualizations
  # ggtext is optional; use if available for richer text annotations
  if (requireNamespace("ggtext", quietly = TRUE)) {
    library(ggtext)
  } else {
    message("[viz] Optional package 'ggtext' not installed; continuing without enhanced text rendering.")
  }
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ------------------- Helper Functions -------------------
# Theme function with enhanced customization for publication-ready plots
plot_theme <- function(cfg) {
  theme_sel <- tolower(cfg$plots$theme %||% "classic")
  base_size <- cfg$plots$base_size %||% 12  # Increased default for better readability
  font_family <- cfg$plots$font_family %||% "sans"
  base <- switch(
    theme_sel,
    minimal = theme_minimal(base_size = base_size, base_family = font_family),
    bw = theme_bw(base_size = base_size, base_family = font_family),
    classic = theme_classic(base_size = base_size, base_family = font_family),
    theme_classic(base_size = base_size, base_family = font_family)
  )
  base + theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold", size = base_size * 1.1),
    axis.text = element_text(size = base_size * 0.9),
    legend.title = element_text(face = "bold", size = base_size),
    legend.text = element_text(size = base_size * 0.8),
    plot.title = element_text(face = "bold", size = base_size * 1.2, hjust = 0.5),
    plot.subtitle = element_text(size = base_size * 0.9, hjust = 0.5),
    plot.caption = element_text(size = base_size * 0.7, hjust = 1)
  )
}

# Color palette generator with more options and recycling
palette_vals <- function(n, cfg) {
  pal <- tolower(cfg$plots$color_palette %||% "viridis")
  if (pal == "okabe-ito") {
    cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
    return(rep(cols, length.out = n))
  } else if (pal == "tableau10") {
    cols <- c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F","#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC")
    return(rep(cols, length.out = n))
  } else if (pal == "viridis") {
    return(viridisLite::viridis(n))
  } else if (pal == "plasma") {
    return(viridisLite::plasma(n))
  } else {
    warning("[viz] Unknown palette; falling back to viridis.")
    viridisLite::viridis(n)
  }
}

# Save plot with enhanced options (e.g., PDF support, compression)
save_plot <- function(p, path, cfg, width = NULL, height = NULL, format = "tiff") {
  w <- width %||% cfg$plots$width %||% 8
  h <- height %||% cfg$plots$height %||% 6
  dpi <- cfg$plots$dpi %||% 600  # Higher default for publication quality
  if (format == "tiff") {
    ggsave(filename = path, plot = p, width = w, height = h, dpi = dpi, compression = "lzw")
  } else if (format == "pdf") {
    ggsave(filename = path, plot = p, width = w, height = h, device = "pdf")
  } else {
    ggsave(filename = path, plot = p, width = w, height = h, dpi = dpi)
  }
}

# Validate and prepare metadata
prepare_metadata <- function(ps, cfg) {
  meta <- tryCatch(as(sample_data(ps), "data.frame"), error = function(e) NULL)
  if (is.null(meta) || ncol(meta) == 0) {
    # Use SampleID as the canonical sample identifier column name
    meta <- data.frame(SampleID = sample_names(ps), Group = "Unknown", row.names = sample_names(ps), check.names = FALSE)
  }
  id_col <- cfg$metadata$id_column %||% "SampleID"
  group_col <- cfg$metadata$group_column %||% "Group"
  missing_cols <- setdiff(c(id_col, group_col), colnames(meta))
  for (mc in missing_cols) {
    meta[[mc]] <- if (mc == group_col) "Unknown" else rownames(meta)
  }
  sample_data(ps) <- sample_data(meta)
  return(ps)
}

# ------------------------- Main Visualization Entrypoint -------------------------
run_plots <- function(cfg) {
  mode <- tolower(cfg$project$sequencing_type %||% "16s")
  if (mode != "16s") {
    warning("[viz] Only 16S visualizations are implemented.")
    return(invisible(NULL))
  }

  base_out <- cfg$project$output_dir %||% "output"
  cohort <- cfg$io$cohort %||% basename(cfg$io$input_dir)
  out_base <- file.path(base_out, cohort)
  analysis_dir <- file.path(out_base, "analysis")
  outdir <- file.path(out_base, "visualizations")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  ps_rds <- file.path(analysis_dir, "phyloseq_rarefied.rds")
  ps_raw_rds <- file.path(analysis_dir, "phyloseq_object_raw.rds")
  if (!file.exists(ps_rds) && !file.exists(ps_raw_rds)) {
    warning("[viz] No phyloseq RDS found. Skipping.")
    return(invisible(NULL))
  }
  ps <- if (file.exists(ps_rds)) readRDS(ps_rds) else readRDS(ps_raw_rds)
  ps <- prepare_metadata(ps, cfg)  # Ensure metadata is robust

  en <- cfg$plots$enable %||% list(
    alpha = TRUE, composition = TRUE, heatmap = TRUE, ordination = TRUE,
    tree_rectangular = TRUE, tree_circular = TRUE, per_sample_panels = TRUE  # Enabled by default now
  )

  # Phylogeny and tree plots (enhanced with checks and annotations)
  if (isTRUE(cfg$amplicon$phylogeny$build_tree) && (isTRUE(en$tree_rectangular) || isTRUE(en$tree_circular))) {
    try({
      ps_tree <- build_and_attach_tree(ps, cfg)
      if (inherits(ps_tree, "phyloseq") && !is.null(phy_tree(ps_tree, errorIfNULL = FALSE))) {
        ps <- ps_tree
        message("[viz] Phylogenetic tree attached.")
        plot_phylo_overview(ps, cfg, outdir, en)
      } else {
        warning("[viz] Tree build failed; check dependencies (phangorn/DECIPHER).")
      }
    })
  }

  # Alpha diversity plots
  if (isTRUE(en$alpha)) plot_alpha_diversity(ps, cfg, outdir)

  # Composition stacked bars per rank
  if (isTRUE(en$composition)) plot_composition(ps, cfg, outdir)

  # Heatmap of top taxa
  if (isTRUE(en$heatmap)) plot_heatmap(ps, cfg, outdir)

  # Ordination (configurable method)
  if (isTRUE(en$ordination)) plot_ordination(ps, cfg, outdir)

  # Per-sample multipanels (enhanced with cowplot and TIFF output)
  if (isTRUE(en$per_sample_panels)) plot_per_sample_panels(ps, cfg, outdir)

  message("[viz] Completed plots.")
}

# ------------------------- Modular Plot Functions -------------------------
# Alpha diversity
plot_alpha_diversity <- function(ps, cfg, outdir) {
  try({
    # Only use measures that phyloseq supports directly
    requested_measures <- cfg$plots$alpha_measures %||% c("Shannon", "Observed", "Simpson")
    available_measures <- c("Shannon", "Observed", "Simpson", "InvSimpson", "Chao1", "ACE")
    measures <- intersect(requested_measures, available_measures)
    if (length(measures) == 0) measures <- c("Shannon", "Observed")
    
    df_alpha <- estimate_richness(ps, measures = measures)
    id_col <- cfg$metadata$id_column %||% "SampleID"
    df_alpha[[id_col]] <- rownames(df_alpha)
    
    # Get metadata safely
    meta <- as(sample_data(ps), "data.frame")
    if (!id_col %in% colnames(meta)) {
      meta[[id_col]] <- rownames(meta)
    }
    
    group_col <- cfg$metadata$group_column %||% "Group"
    # Merge by id_col
    df_alpha <- merge(df_alpha, meta, by = id_col, all.x = TRUE)

    for (meas in measures) {
      # Create transformed column with valid name
      if (meas == "Observed") {
        df_alpha$Observed_log1p <- log1p(df_alpha$Observed)
        y_col <- "Observed_log1p"
        y_label <- "log1p(Observed)"
      } else {
        y_col <- meas
        y_label <- meas
      }
      
      p <- ggplot(df_alpha, aes(x = .data[[group_col]], y = .data[[y_col]], fill = .data[[group_col]])) +
        geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +
        geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
        scale_fill_manual(values = palette_vals(length(unique(df_alpha[[group_col]])), cfg)) +
        plot_theme(cfg) +
        labs(title = paste(meas, "Diversity"),
             subtitle = paste("N =", nrow(df_alpha), "samples"),
             x = group_col, y = y_label,
             caption = "Boxplots show median and IQR; points are individual samples") +
        theme(legend.position = "none")  # Remove redundant legend when x-axis shows groups
      save_plot(p, file.path(outdir, paste0("alpha_", tolower(meas), "_boxplot.tiff")), cfg)
    }
  }, silent = FALSE)  # Show errors for debugging
}

# Composition
plot_composition <- function(ps, cfg, outdir) {
  ranks <- cfg$plots$ranks %||% c("Phylum", "Class", "Order", "Family", "Genus")
  topN <- cfg$plots$top_taxa %||% 10  # Reduced from 15 for clarity
  group_col <- cfg$metadata$group_column %||% "Group"
  for (rk in ranks) {
    try({
      ps_r <- tax_glom(ps, taxrank = rk, NArm = TRUE)
      prop <- transform_sample_counts(ps_r, function(x) x / sum(x))
      taxa_sums_rank <- sort(taxa_sums(prop), TRUE)
      keep <- names(taxa_sums_rank)[seq_len(min(topN, length(taxa_sums_rank)))]
      prop_f <- prune_taxa(keep, prop)
      df <- psmelt(prop_f)
      # Determine sample column name produced by psmelt
      id_col <- cfg$metadata$id_column %||% "SampleID"
      sample_col <- if (id_col %in% colnames(df)) id_col else if ("Sample" %in% colnames(df)) "Sample" else grep("^sample", colnames(df), ignore.case = TRUE, value = TRUE)[1]
      df$Abundance <- pmax(df$Abundance, 0)
      df[[rk]] <- forcats::fct_lump_n(as.character(df[[rk]]), n = topN, other_level = "Other")
      
      # Count samples per group for annotation
      meta <- as(sample_data(ps), "data.frame")
      group_counts <- table(meta[[group_col]])
      subtitle_text <- paste0("Top ", topN, " taxa; N per group: ", paste(names(group_counts), "=", group_counts, collapse = ", "))
      
      p <- ggplot(df, aes(x = .data[[sample_col]], y = Abundance, fill = .data[[rk]])) +
        geom_col(width = 0.9) +
        facet_grid(cols = vars(.data[[group_col]]), scales = "free_x", space = "free_x") +
        scale_y_continuous(labels = scales::label_percent(accuracy = 1), expand = c(0, 0)) +
        scale_fill_manual(values = palette_vals(length(levels(df[[rk]])), cfg)) +
        plot_theme(cfg) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              legend.position = "bottom", legend.direction = "horizontal",
              legend.box = "vertical", legend.key.size = unit(0.7, "cm")) +  # Increased from 0.4
        labs(title = paste("Taxonomic Composition by", rk),
             subtitle = subtitle_text,
             y = "Relative Abundance (%)", x = "Samples (grouped)",
             caption = "Stacked bars show proportions; less abundant taxa grouped as 'Other'")
      save_plot(p, file.path(outdir, paste0("composition_", tolower(rk), ".tiff")), cfg, width = 12, height = 6)
    })
  }
}

# Heatmap
plot_heatmap <- function(ps, cfg, outdir) {
  try({
    tax_rank <- cfg$plots$heatmap_rank %||% "Genus"
    topN <- cfg$plots$heatmap_top %||% 25  # Reduced from 50 for readability
    ps_g <- tax_glom(ps, taxrank = tax_rank, NArm = TRUE)
    prop <- transform_sample_counts(ps_g, function(x) x / sum(x))
    taxa_ord <- names(sort(taxa_sums(prop), decreasing = TRUE))[1:min(topN, ntaxa(prop))]
    prop_t <- prune_taxa(taxa_ord, prop)
    mat <- as(otu_table(prop_t), "matrix")
    mat <- t(mat)  # Taxa as rows, samples as columns
    df_m <- melt(mat, varnames = c("Taxon", "Sample"), value.name = "RA")
    df_m$Abundance <- log10(df_m$RA + 1e-4)
    p <- ggplot(df_m, aes(x = Sample, y = Taxon, fill = Abundance)) +
      geom_tile(color = "white", linewidth = 0.1) +
      scale_fill_viridis_c(option = "plasma", name = "log10(RA)", direction = -1, oob = squish) +
      plot_theme(cfg) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = rel(0.8)),
            axis.text.y = element_text(size = rel(0.75))) +
      labs(title = paste("Heatmap of Top", topN, tax_rank),
           subtitle = "Log10-transformed relative abundances",
           x = "Samples", y = tax_rank,
           caption = "Color intensity represents abundance; white = low/absent")
    save_plot(p, file.path(outdir, paste0("heatmap_top_", tolower(tax_rank), ".tiff")), cfg, width = 12, height = 8)
  })
}

# Ordination
plot_ordination <- function(ps, cfg, outdir) {
  try({
    set.seed(cfg$project$random_seed %||% 1234)
    dist_method <- cfg$analysis$distance %||% "bray"
    ord_method <- cfg$analysis$ordination %||% "PCoA"
    dist <- phyloseq::distance(ps, method = dist_method)
    ord <- ordinate(ps, method = ord_method, distance = dist)
    group_col <- cfg$metadata$group_column %||% "Group"
    meta <- as(sample_data(ps), "data.frame")
    # Use phyloseq's plot_ordination explicitly to avoid naming conflicts with this script
    p <- phyloseq::plot_ordination(ps, ord, color = group_col) +
      geom_point(size = 3, alpha = 0.9) +
      scale_color_manual(values = palette_vals(length(unique(meta[[group_col]])), cfg)) +
      plot_theme(cfg) +
      labs(title = paste(ord_method, "Ordination (", dist_method, "distance)"),
           subtitle = "Beta diversity visualization",
           caption = "Points colored by group; ellipses at 68% confidence if applicable")
    sizes <- table(meta[[group_col]])
    if (all(sizes >= 3)) {
      p <- p + stat_ellipse(type = "t", level = 0.68, linetype = 2, linewidth = 0.5)
    }
    save_plot(p, file.path(outdir, paste0("beta_", tolower(ord_method), "_", dist_method, ".tiff")), cfg)
  })
}

# Per-sample multipanels
plot_per_sample_panels <- function(ps, cfg, outdir) {
  try({
    top_rank <- cfg$plots$per_sample_rank %||% "Genus"
    topN <- cfg$plots$per_sample_top %||% 10
    depth <- data.frame(Sample = sample_names(ps), Depth = sample_sums(ps))
    df_alpha <- estimate_richness(ps, measures = c("Shannon"))
    df_alpha$Sample <- rownames(df_alpha)
    prop <- transform_sample_counts(tax_glom(ps, taxrank = top_rank), function(x) x / sum(x))
    taxa_ord <- names(sort(taxa_sums(prop), decreasing = TRUE))[1:min(topN, ntaxa(prop))]
    prop_t <- prune_taxa(taxa_ord, prop)
    df_top <- psmelt(prop_t) %>% group_by(Sample, .data[[top_rank]]) %>% summarise(Abundance = sum(Abundance), .groups = "drop")

    for (s in unique(depth$Sample)) {
      d1 <- ggplot(depth[depth$Sample == s, ], aes(x = Sample, y = Depth)) +
        geom_col(fill = "grey40") + coord_flip() + plot_theme(cfg) +
        labs(title = paste("Sequencing Depth:", s), y = "Reads", x = NULL)
      d2 <- ggplot(df_alpha[df_alpha$Sample == s, ], aes(x = Sample, y = Shannon)) +
        geom_col(fill = "steelblue") + coord_flip() + plot_theme(cfg) +
        labs(title = "Shannon Diversity", y = "Index", x = NULL)
      d3 <- ggplot(df_top[df_top$Sample == s, ], aes(x = reorder(.data[[top_rank]], Abundance), y = Abundance, fill = .data[[top_rank]])) +
        geom_col() + coord_flip() +
        scale_y_continuous(labels = scales::label_percent(accuracy = 1), expand = c(0, 0)) +
        scale_fill_manual(values = palette_vals(length(unique(df_top[[top_rank]])), cfg)) +
        plot_theme(cfg) + theme(legend.position = "none") +
        labs(title = paste("Top", topN, top_rank), y = "Relative Abundance (%)", x = NULL)
      multi_p <- cowplot::plot_grid(d1, d2, d3, ncol = 3, rel_widths = c(1, 1, 2),
                                    labels = "AUTO", label_size = cfg$plots$base_size %||% 12)
      save_plot(multi_p, file.path(outdir, paste0("per_sample_", s, ".tiff")), cfg, width = 12, height = 4)
    }
  })
}

# ------------------------------ Phylogeny Helpers ------------------------------
build_and_attach_tree <- function(ps, cfg) {
  if (!requireNamespace("phangorn", quietly = TRUE) || !requireNamespace("DECIPHER", quietly = TRUE)) {
    warning("[viz] Missing phangorn or DECIPHER; cannot build tree.")
    return(ps)
  }
  max_tips <- cfg$amplicon$phylogeny$max_tips %||% 500
  keep_taxa <- names(sort(taxa_sums(ps), decreasing = TRUE))[seq_len(min(max_tips, ntaxa(ps)))]
  ps_sub <- prune_taxa(keep_taxa, ps)

  taxa_ids <- taxa_names(ps_sub)
  # Look for ASV fasta in cohort-scoped analysis directory
  base_out <- cfg$project$output_dir %||% "output"
  cohort <- cfg$io$cohort %||% basename(cfg$io$input_dir)
  fasta_file <- file.path(base_out, cohort, "analysis", "asv_sequences.fasta")
  seqs <- NULL
  if (file.exists(fasta_file)) {
    fas <- Biostrings::readDNAStringSet(fasta_file)
    names(fas) <- gsub("\t.*$", "", names(fas))
    have <- intersect(names(fas), taxa_ids)
    if (length(have) == 0) {
      warning("[viz] No matching ASV IDs in fasta; cannot build tree.")
      return(ps)
    }
    if (length(have) < length(taxa_ids)) {
      warning("[viz] Partial ASV fasta found; tree will be built for subset of taxa present in fasta.")
      # prune ps_sub to only those taxa that have sequences
      ps_sub <- prune_taxa(have, ps_sub)
      taxa_ids <- taxa_names(ps_sub)
    }
    seqs <- as.character(fas[have])
    names(seqs) <- have
  } else {
    warning("[viz] asv_sequences.fasta not found in cohort analysis dir; cannot build tree.")
    return(ps)
  }
  names(seqs) <- taxa_ids
  seqs_ds <- Biostrings::DNAStringSet(seqs)
  alignment <- DECIPHER::AlignSeqs(seqs_ds, processors = cfg$project$threads %||% 2)
  phang_align <- phangorn::phyDat(as.matrix(alignment), type = "DNA")
  dm <- phangorn::dist.ml(phang_align)
  treeNJ <- phangorn::NJ(dm)
  fit <- phangorn::pml(treeNJ, data = phang_align)
  fitGTR <- phangorn::optim.pml(fit, model = "GTR", optInv = TRUE, optGamma = TRUE,
                                rearrangement = "stochastic", control = phangorn::pml.control(trace = 0))
  phy_tree(ps_sub) <- fitGTR$tree
  return(ps_sub)
}

# Enhanced tree plotting with better scaling and annotations
plot_phylo_overview <- function(ps, cfg, outdir, en) {
  tree <- phy_tree(ps)
  if (length(tree$tip.label) < 2) return(warning("[viz] Tree too small; skipping."))
  taxa_abund <- taxa_sums(ps)
  tt <- as(tax_table(ps), "matrix")
  df <- data.frame(label = names(taxa_abund),
                   Abundance = taxa_abund,
                   Phylum = tt[names(taxa_abund), "Phylum"] %||% "Unassigned",
                   stringsAsFactors = FALSE)
  df$Phylum[is.na(df$Phylum)] <- "Unassigned"
  df$AbundanceScaled <- scales::rescale(log1p(df$Abundance), to = c(1, 5))  # Log-scaled for better visualization
  df <- df[df$label %in% tree$tip.label, ]
  
  n_phyla <- length(unique(df$Phylum))
  subtitle_text <- paste(nrow(df), "ASVs from", n_phyla, "phyla; sized by log abundance")

  if (isTRUE(en$tree_rectangular)) {
    p1 <- ggtree(tree, layout = "rectangular", size = 0.4) %<+% df +
      geom_tippoint(aes(size = AbundanceScaled, color = Phylum), alpha = 0.8) +
      scale_color_manual(values = palette_vals(n_phyla, cfg)) +
      scale_size_continuous(range = c(1, 4)) +
      guides(size = "none", color = guide_legend(override.aes = list(size = 3))) +
      plot_theme(cfg) +
      theme(legend.position = "right") +
      labs(title = "Phylogenetic Tree (Rectangular)",
           subtitle = subtitle_text,
           caption = "ML tree (GTR+G+I model); tips show ASVs colored by Phylum")
    save_plot(p1, file.path(outdir, "phylo_tree_rectangular.tiff"), cfg, width = 12, height = 10)
    save_plot(p1, file.path(outdir, "phylo_tree_rectangular.pdf"), cfg, width = 12, height = 10, format = "pdf")
  }

  if (isTRUE(en$tree_circular)) {
    p2 <- ggtree(tree, layout = "circular", size = 0.3) %<+% df +
      geom_tippoint(aes(size = AbundanceScaled, color = Phylum), alpha = 0.75) +
      scale_color_manual(values = palette_vals(n_phyla, cfg)) +
      scale_size_continuous(range = c(1, 4)) +
      guides(size = "none", color = guide_legend(override.aes = list(size = 3))) +
      plot_theme(cfg) +
      theme(legend.position = "right") +
      labs(title = "Phylogenetic Tree (Circular)",
           subtitle = subtitle_text,
           caption = "ML tree (GTR+G+I model); tips show ASVs colored by Phylum")
    save_plot(p2, file.path(outdir, "phylo_tree_circular.tiff"), cfg, width = 10, height = 10)
    save_plot(p2, file.path(outdir, "phylo_tree_circular.pdf"), cfg, width = 10, height = 10, format = "pdf")
  }
}