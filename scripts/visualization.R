#!/usr/bin/env Rscript
# Unified visualization script for metagenomics (currently 16S)
# Plots are controlled via cfg$plots$enable toggles.

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(gridExtra)
  library(reshape2)
  library(viridisLite)
  library(scales)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ------------------- helpers (formerly visualization_common.R) -------------------
plot_theme <- function(cfg) {
  theme_sel <- cfg$plots$theme %||% "classic"
  base <- switch(
    tolower(theme_sel),
    minimal = theme_minimal(base_size = 11),
    bw = theme_bw(base_size = 11),
    classic = theme_classic(base_size = 11),
    theme_classic(base_size = 11)
  )
  base + theme(panel.grid = element_blank())
}

palette_vals <- function(n, cfg) {
  pal <- tolower(cfg$plots$color_palette %||% "viridis")
  if (pal == "okabe-ito") {
    cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
    return(rep(cols, length.out = n))
  }
  if (pal == "tableau10") {
    cols <- c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F","#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC")
    return(rep(cols, length.out = n))
  }
  viridisLite::viridis(n)
}

save_plot <- function(p, path, cfg, width = NULL, height = NULL) {
  w <- width %||% cfg$plots$width %||% 8
  h <- height %||% cfg$plots$height %||% 6
  dpi <- cfg$plots$dpi %||% 300
  ggsave(filename = path, plot = p, width = w, height = h, dpi = dpi)
}

# ------------------------- main visualization entrypoint -------------------------
run_plots <- function(cfg) {
  mode <- tolower(cfg$project$sequencing_type %||% "16s")
  if (mode != "16s") {
    warning("[viz] Only 16S visualizations are implemented in this unified script right now.")
  }

  base_out <- cfg$project$output_dir %||% "output"
  cohort <- cfg$io$cohort
  if (is.null(cohort) || is.na(cohort) || cohort == "") cohort <- basename(cfg$io$input_dir)
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

  # Metadata guard
  meta <- tryCatch({ as(sample_data(ps), "data.frame") }, error = function(e) NULL)
  if (is.null(meta) || ncol(meta) == 0) {
    meta <- data.frame(Sample = sample_names(ps), Group = "Unknown", row.names = sample_names(ps), check.names = FALSE)
    sample_data(ps) <- sample_data(meta)
  }
  id_col <- cfg$metadata$id_column %||% "Sample"
  group_col <- cfg$metadata$group_column %||% "Group"
  if (!all(c(id_col, group_col) %in% colnames(meta))) {
    missing_cols <- setdiff(c(id_col, group_col), colnames(meta))
    for (mc in missing_cols) meta[[mc]] <- if (mc == group_col) "Unknown" else NA
    sample_data(ps) <- sample_data(meta)
  }

  en <- cfg$plots$enable %||% list(
    alpha = TRUE, composition = TRUE, heatmap = TRUE, ordination = TRUE,
    tree_rectangular = TRUE, tree_circular = TRUE, per_sample_panels = FALSE
  )

  # Phylogeny and tree plots
  if (isTRUE(cfg$amplicon$phylogeny$build_tree) && (isTRUE(en$tree_rectangular) || isTRUE(en$tree_circular))) {
    try({
      ps_tree <- try(build_and_attach_tree(ps, cfg), silent = TRUE)
      if (inherits(ps_tree, "phyloseq")) {
        ps <- ps_tree
        if (!is.null(phy_tree(ps, errorIfNULL = FALSE))) {
          message("[viz] Phylogenetic tree attached to phyloseq object.")
          plot_phylo_overview(ps, cfg, outdir, en)
        } else {
          warning("[viz] Tree not present after build; dependencies may be missing (phangorn/DECIPHER).")
        }
      }
    }, silent = TRUE)
  }

  # Alpha diversity
  if (isTRUE(en$alpha)) {
    try({
      df_alpha <- estimate_richness(ps, measures = c("Shannon", "Observed"))
      df_alpha$Sample <- rownames(df_alpha)
      md <- meta %>% tibble::rownames_to_column(var = "Sample")
      df_alpha <- left_join(df_alpha, md, by = "Sample")

      p_alpha <- ggplot(df_alpha, aes(x = .data[[group_col]], y = Shannon, fill = .data[[group_col]])) +
        geom_boxplot(outlier.shape = NA, alpha = 0.8) +
        geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
        scale_fill_manual(values = palette_vals(length(unique(df_alpha[[group_col]])), cfg)) +
        plot_theme(cfg) +
        labs(title = "Alpha diversity (Shannon)", x = group_col, y = "Shannon")
      save_plot(p_alpha, file.path(outdir, "alpha_shannon_boxplot.tiff"), cfg)

      p_obs <- ggplot(df_alpha, aes(x = .data[[group_col]], y = log1p(Observed), fill = .data[[group_col]])) +
        geom_boxplot(outlier.shape = NA, alpha = 0.8) +
        geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
        scale_fill_manual(values = palette_vals(length(unique(df_alpha[[group_col]])), cfg)) +
        plot_theme(cfg) +
        labs(title = "Observed richness (log1p)", x = group_col, y = "log1p(Observed)")
      save_plot(p_obs, file.path(outdir, "alpha_observed_log1p_boxplot.tiff"), cfg)
    }, silent = TRUE)
  }

  # Composition stacked bars per rank
  if (isTRUE(en$composition)) {
    ranks <- cfg$plots$ranks %||% c("Phylum", "Class", "Order", "Family", "Genus")
    for (rk in ranks) {
      try({
        ps_r <- suppressWarnings(tax_glom(ps, taxrank = rk, NArm = TRUE))
        prop <- transform_sample_counts(ps_r, function(x) x / sum(x))
        topN <- cfg$plots$top_taxa %||% 15
        taxa_sums_rank <- sort(taxa_sums(prop), TRUE)
        keep <- names(taxa_sums_rank)[seq_len(min(topN, length(taxa_sums_rank)))]
        prop_f <- prune_taxa(keep, prop)
        df <- suppressWarnings(psmelt(prop_f))
        df$Abundance <- pmax(df$Abundance, 0)
        df[[rk]] <- as.character(df[[rk]])
        df[[rk]] <- forcats::fct_lump_n(df[[rk]], n = topN, other_level = "Other")
        p <- ggplot(df, aes(x = Sample, y = Abundance, fill = .data[[rk]])) +
          geom_col(width = 0.9) +
          facet_grid(cols = vars(.data[[group_col]]), scales = "free_x", space = "free_x") +
          scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
          scale_fill_manual(values = palette_vals(length(levels(df[[rk]])), cfg)) +
          plot_theme(cfg) +
          theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom") +
          labs(title = paste("Composition by", rk), y = "Relative abundance", x = "Samples")
        save_plot(p, file.path(outdir, paste0("composition_", tolower(rk), ".tiff")), cfg, width = 12, height = 6)
      }, silent = TRUE)
    }
  }

  # Heatmap of top genera
  if (isTRUE(en$heatmap)) {
    try({
      ps_g <- suppressWarnings(tax_glom(ps, taxrank = "Genus", NArm = TRUE))
      prop <- transform_sample_counts(ps_g, function(x) x / sum(x))
      taxa_ord <- names(sort(taxa_sums(prop), decreasing = TRUE))[1:min(50, ntaxa(prop))]
      prop_t <- prune_taxa(taxa_ord, prop)
      mat <- as(otu_table(prop_t), "matrix")
      mat <- t(mat)
      df_m <- reshape2::melt(mat, varnames = c("Sample", "Taxon"), value.name = "RA")
      df_m$Abundance <- log10(df_m$RA + 1e-4)
      p <- ggplot(df_m, aes(x = Taxon, y = Sample, fill = Abundance)) +
        geom_tile() +
        scale_fill_viridis_c(option = "C", name = "log10(RA+1e-4)", oob = scales::squish) +
        plot_theme(cfg) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        labs(title = "Top genera heatmap (log-scaled)", x = "Genus", y = "Sample")
      save_plot(p, file.path(outdir, "heatmap_top_genera.tiff"), cfg, width = 12, height = 8)
    }, silent = TRUE)
  }

  # Ordination (Bray)
  if (isTRUE(en$ordination)) {
    try({
      set.seed(cfg$project$random_seed %||% 1234)
      dist <- phyloseq::distance(ps, method = "bray")
      ord <- ordinate(ps, method = cfg$analysis$ordination %||% "PCoA", distance = dist)
      p <- plot_ordination(ps, ord, color = group_col) +
        geom_point(size = 2, alpha = 0.9) +
        scale_color_manual(values = palette_vals(length(unique(meta[[group_col]])), cfg)) +
        plot_theme(cfg) +
        labs(title = "PCoA (Bray)")
      sizes <- table(meta[[group_col]])
      if (length(sizes) > 0 && min(sizes, na.rm = TRUE) >= 3) {
        p <- p + stat_ellipse(type = "t", level = 0.68, linetype = 2, alpha = 0.6, linewidth = 0.5)
      }
      save_plot(p, file.path(outdir, "beta_pcoa_bray.tiff"), cfg)
    }, silent = TRUE)
  }

  # Per-sample panels
  if (isTRUE(en$per_sample_panels)) {
    try({
      depth <- data.frame(Sample = sample_names(ps), Depth = sample_sums(ps), check.names = FALSE)
      df_alpha <- estimate_richness(ps, measures = c("Shannon"))
      df_alpha$Sample <- rownames(df_alpha)
      top_rank <- suppressWarnings(tax_glom(ps, taxrank = "Genus"))
      prop <- transform_sample_counts(top_rank, function(x) x / sum(x))
      taxa_ord <- names(sort(taxa_sums(prop), decreasing = TRUE))[1:min(10, ntaxa(prop))]
      prop_t <- prune_taxa(taxa_ord, prop)
      df_top <- psmelt(prop_t) %>% group_by(Sample, Genus) %>% summarise(Abundance = sum(Abundance), .groups = "drop")

      pdf(file.path(outdir, "per_sample_panels.pdf"), width = 11, height = 8.5)
      for (s in unique(depth$Sample)) {
        d1 <- ggplot(depth[depth$Sample == s,], aes(x = Sample, y = Depth)) +
          geom_col(fill = "grey40") + coord_flip() + plot_theme(cfg) + labs(title = paste("Depth:", s))
        d2 <- ggplot(df_alpha[df_alpha$Sample == s,], aes(x = Sample, y = Shannon)) +
          geom_col(fill = "steelblue") + coord_flip() + plot_theme(cfg) + labs(title = "Shannon")
        d3 <- ggplot(df_top[df_top$Sample == s,], aes(x = Genus, y = Abundance, fill = Genus)) +
          geom_col() + coord_flip() +
          scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
          scale_fill_manual(values = palette_vals(length(unique(df_top$Genus)), cfg)) +
          plot_theme(cfg) + labs(title = "Top genera")
        gridExtra::grid.arrange(d1, d2, d3, ncol = 3)
      }
      dev.off()
    }, silent = TRUE)
  }

  message("[viz] Completed plots.")
}

# ------------------------------ phylogeny helpers ------------------------------
build_and_attach_tree <- function(ps, cfg) {
  has_phangorn <- requireNamespace("phangorn", quietly = TRUE)
  has_decipher <- requireNamespace("DECIPHER", quietly = TRUE)
  if (!has_phangorn || !has_decipher) {
    warning("[viz] phangorn and/or DECIPHER not installed. Cannot build tree.")
    return(ps)
  }
  max_tips <- cfg$amplicon$phylogeny$max_tips %||% 500
  keep_taxa <- names(sort(taxa_sums(ps), decreasing = TRUE))[seq_len(min(max_tips, ntaxa(ps)))]
  ps_sub <- prune_taxa(keep_taxa, ps)

  taxa_ids <- taxa_names(ps_sub)
  seqs <- taxa_ids
  if (!all(grepl("^[ACGTNacgtn]+$", seqs))) {
    fasta_file <- file.path(cfg$project$output_dir %||% "output", "asv_sequences.fasta")
    if (file.exists(fasta_file)) {
      fas <- Biostrings::readDNAStringSet(fasta_file)
      names(fas) <- gsub("\t.*$", "", names(fas))
      have <- intersect(names(fas), taxa_ids)
      if (length(have) > 0) {
        fas <- fas[have]
        seqs <- as.character(fas)
        names(seqs) <- names(fas)
      } else {
        warning("[viz] No matching ASV IDs in asv_sequences.fasta; using taxa_names as proxy.")
      }
    } else {
      warning("[viz] asv_sequences.fasta not found; using taxa_names as proxy.")
    }
  }
  names(seqs) <- taxa_ids
  seqs_ds <- Biostrings::DNAStringSet(seqs)

  alignment <- DECIPHER::AlignSeqs(seqs_ds, processors = cfg$project$threads %||% 2)
  aln_mat <- as.matrix(alignment)
  phang_align <- phangorn::phyDat(aln_mat, type = "DNA")
  dm <- phangorn::dist.ml(phang_align)
  treeNJ <- phangorn::NJ(dm)
  fit <- phangorn::pml(treeNJ, data = phang_align)
  fitGTR <- phangorn::optim.pml(
    fit,
    model = "GTR", optInv = TRUE, optGamma = TRUE, k = 4, inv = 0.2, rearrangement = "stochastic",
    control = phangorn::pml.control(trace = 0)
  )
  tree <- fitGTR$tree
  phyloseq::phy_tree(ps_sub) <- tree
  return(ps_sub)
}

plot_phylo_overview <- function(ps, cfg, outdir, en) {
  if (is.null(phy_tree(ps, errorIfNULL = FALSE))) return(invisible(NULL))
  tree <- phy_tree(ps)
  if (length(tree$tip.label) < 2) {
    warning("[viz] Tree has fewer than 2 tips; skipping tree plots.")
    return(invisible(NULL))
  }
  taxa_abund <- taxa_sums(ps)
  tt <- as(tax_table(ps), "matrix")
  df <- data.frame(label = names(taxa_abund),
                   Abundance = as.numeric(taxa_abund),
                   Phylum = tt[names(taxa_abund), "Phylum"],
                   stringsAsFactors = FALSE)
  df$Phylum[is.na(df$Phylum)] <- "Unassigned"
  df$AbundanceScaled <- scales::rescale(df$Abundance, to = c(0.8, 3.5))
  df <- df[df$label %in% tree$tip.label, , drop = FALSE]

  if (isTRUE(en$tree_rectangular)) {
    p1 <- ggtree::ggtree(tree, layout = "rectangular", size = 0.3) %<+% df +
      ggtree::geom_tippoint(aes(size = AbundanceScaled, color = Phylum), alpha = 0.8, stroke = 0) +
      scale_color_manual(values = palette_vals(length(unique(df$Phylum)), cfg)) +
      guides(size = "none") +
      plot_theme(cfg)
    save_plot(p1, file.path(outdir, "phylo_tree_rectangular.tiff"), cfg, width = 10, height = 10)
  }

  if (isTRUE(en$tree_circular)) {
    p2 <- ggtree::ggtree(tree, layout = "circular", size = 0.25) %<+% df +
      ggtree::geom_tippoint(aes(size = AbundanceScaled, color = Phylum), alpha = 0.75, stroke = 0) +
      scale_color_manual(values = palette_vals(length(unique(df$Phylum)), cfg)) +
      guides(size = "none") +
      plot_theme(cfg)
    save_plot(p2, file.path(outdir, "phylo_tree_circular.tiff"), cfg, width = 10, height = 10)
    save_plot(p2, file.path(outdir, "phylo_tree_circular.png"), cfg, width = 10, height = 10)
    ggsave(filename = file.path(outdir, "phylo_tree_circular.pdf"), plot = p2, width = 10, height = 10)
  }
}
