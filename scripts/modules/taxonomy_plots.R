# Taxonomy Plotting Module
# Functions for taxonomic composition visualizations

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(forcats)
})

# Source utilities - handle both direct execution and sourcing from other scripts
if (!exists("plot_theme")) {
  script_dir <- tryCatch(
    dirname(sys.frame(1)$ofile),
    error = function(e) "scripts/modules"
  )
  if (file.exists(file.path(script_dir, "plot_utils.R"))) {
    source(file.path(script_dir, "plot_utils.R"))
  } else if (file.exists("scripts/modules/plot_utils.R")) {
    source("scripts/modules/plot_utils.R")
  }
}

#' Plot taxonomic composition barplot
#'
#' @param ps Phyloseq object
#' @param cfg Configuration list
#' @param outdir Output directory
#' @param rank Taxonomic rank to plot
#' @param top_n Number of top taxa to show
#' @export
plot_composition <- function(ps, cfg, outdir, rank = "Phylum", top_n = 10) {
  tryCatch({
    # Agglomerate to rank
    ps_glom <- tax_glom(ps, taxrank = rank, NArm = FALSE)
    
    # Get relative abundance
    ps_rel <- transform_sample_counts(ps_glom, function(x) x / sum(x))
    
    # Get top taxa
    top_taxa <- get_top_taxa(ps_rel, n = top_n, rank = rank)
    
    # Prune to top taxa
    ps_top <- prune_taxa(top_taxa, ps_rel)
    
    # Convert to data frame
    df <- psmelt(ps_top)
    
    # Get grouping variable
    group_col <- cfg$metadata$group_column %||% "Group"
    if (!group_col %in% colnames(df)) {
      message("[taxonomy] Skipping: '", group_col, "' not found")
      return(invisible(NULL))
    }
    
    # Create plot
    p <- ggplot(df, aes(x = Sample, y = Abundance, fill = .data[[rank]])) +
      geom_bar(stat = "identity", position = "stack") +
      facet_wrap(as.formula(paste("~", group_col)), scales = "free_x", nrow = 1) +
      scale_fill_manual(values = palette_vals(length(unique(df[[rank]])), cfg)) +
      plot_theme(cfg) +
      labs(
        title = paste(rank, "Composition (Top", top_n, ")"),
        subtitle = paste("Relative abundance by", group_col),
        y = "Relative Abundance",
        x = "Samples"
      ) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right"
      )
    
    # Save
    filename <- paste0("composition_", tolower(rank), "_top", top_n, ".tiff")
    save_plot(p, file.path(outdir, filename), cfg, width = 12)
    
    message("[taxonomy] Generated composition plot for ", rank)
    
  }, error = function(e) {
    warning("[taxonomy] Failed to generate composition plot: ", e$message)
  })
}

#' Plot heatmap of top taxa
#'
#' @param ps Phyloseq object
#' @param cfg Configuration list
#' @param outdir Output directory
#' @param rank Taxonomic rank
#' @param top_n Number of top taxa
#' @export
plot_taxa_heatmap <- function(ps, cfg, outdir, rank = "Genus", top_n = 30) {
  message("[taxonomy] Generating heatmap for ", rank, " (top ", top_n, ")")
  
  tryCatch({
    # Agglomerate and get top taxa
    ps_glom <- tax_glom(ps, taxrank = rank, NArm = FALSE)
    top_taxa <- get_top_taxa(ps_glom, n = top_n, rank = rank)
    ps_top <- prune_taxa(top_taxa, ps_glom)
    
    message("[taxonomy] Selected ", ntaxa(ps_top), " taxa")
    
    # Transform to relative abundance
    ps_rel <- transform_sample_counts(ps_top, function(x) x / sum(x))
    
    # Create heatmap using phyloseq function
    p <- phyloseq::plot_heatmap(
      ps_rel,
      taxa.label = rank,
      sample.order = sample_names(ps_rel),
      low = "white",
      high = "darkred",
      na.value = "white"
    ) +
      plot_theme(cfg) +
      labs(
        title = paste(rank, "Heatmap (Top", top_n, ")"),
        subtitle = "Relative abundance across samples"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
    
    # Save
    filename <- paste0("heatmap_", tolower(rank), "_top", top_n, ".tiff")
    save_plot(p, file.path(outdir, filename), cfg, width = 14, height = 10)
    
    message("[taxonomy] Generated heatmap for ", rank)
    
  }, error = function(e) {
    warning("[taxonomy] Failed to generate heatmap: ", e$message)
  })
}

#' Plot top abundant taxa as barplot
#'
#' @param ps Phyloseq object
#' @param cfg Configuration list
#' @param outdir Output directory
#' @param rank Taxonomic rank
#' @param top_n Number of taxa to show
#' @export
plot_top_taxa_barplot <- function(ps, cfg, outdir, rank = "Genus", top_n = 20) {
  tryCatch({
    # Agglomerate
    ps_glom <- tax_glom(ps, taxrank = rank, NArm = FALSE)
    
    # Get abundances
    abund <- taxa_sums(ps_glom)
    top_taxa <- names(sort(abund, decreasing = TRUE)[1:min(top_n, length(abund))])
    
    # Get taxa names
    tax_table_df <- as.data.frame(tax_table(ps_glom))
    taxa_names <- tax_table_df[top_taxa, rank]
    
    # Create data frame
    df <- data.frame(
      Taxa = taxa_names,
      Abundance = abund[top_taxa]
    ) %>%
      arrange(desc(Abundance)) %>%
      mutate(Taxa = fct_reorder(Taxa, Abundance))
    
    # Plot
    p <- ggplot(df, aes(x = Taxa, y = Abundance, fill = Taxa)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_fill_manual(values = palette_vals(nrow(df), cfg)) +
      plot_theme(cfg) +
      labs(
        title = paste("Top", top_n, rank, "by Total Abundance"),
        subtitle = paste("Across all samples (N =", nsamples(ps), ")"),
        x = rank,
        y = "Total Abundance"
      ) +
      theme(legend.position = "none")
    
    # Save
    filename <- paste0("top_", tolower(rank), "_abundance.tiff")
    save_plot(p, file.path(outdir, filename), cfg)
    
    message("[taxonomy] Generated top taxa barplot for ", rank)
    
  }, error = function(e) {
    warning("[taxonomy] Failed to generate top taxa plot: ", e$message)
  })
}
