# Refactored Visualization Script (using modular approach)
# This is the new main visualization script that orchestrates the modules

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
})

# Source all modules
module_dir <- file.path("scripts", "modules")
source(file.path(module_dir, "plot_utils.R"))
source(file.path(module_dir, "alpha_diversity_plots.R"))
source(file.path(module_dir, "beta_diversity_plots.R"))
source(file.path(module_dir, "taxonomy_plots.R"))
source(file.path(module_dir, "tree_plots.R"))
source(file.path(module_dir, "group_comparison_plots.R"))
source(file.path(module_dir, "differential_abundance.R"))
source(file.path(module_dir, "correlation_plots.R"))

#' Main visualization entry point
#'
#' @param cfg Configuration list
#' @export
run_plots <- function(cfg) {
  mode <- tolower(cfg$project$sequencing_type %||% "16s")
  if (mode != "16s") {
    message("[viz] Skipping: only 16S mode supported currently")
    return(invisible(NULL))
  }
  
  message("========================================")
  message("Starting Visualization Pipeline")
  message("========================================")
  
  # Setup output directory
  output_base <- cfg$project$output_dir %||% "output"
  cohort <- cfg$io$cohort %||% basename(cfg$io$input_dir)
  outdir <- file.path(output_base, cohort, "visualizations")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  message("[viz] Output directory: ", outdir)
  
  # Load phyloseq objects
  analysis_dir <- file.path(output_base, cohort, "analysis")
  ps_raw_path <- file.path(analysis_dir, "phyloseq_object_raw.rds")
  ps_rare_path <- file.path(analysis_dir, "phyloseq_rarefied.rds")
  
  if (!file.exists(ps_rare_path) && !file.exists(ps_raw_path)) {
    stop("[viz] No phyloseq objects found. Run analysis first.")
  }
  
  # Use rarefied if available, otherwise raw
  ps_path <- if (file.exists(ps_rare_path)) ps_rare_path else ps_raw_path
  ps <- readRDS(ps_path)
  message("[viz] Loaded: ", basename(ps_path))
  message("[viz] Samples: ", nsamples(ps), " | Taxa: ", ntaxa(ps))
  
  # Prepare metadata
  ps <- prepare_metadata(ps, cfg)
  
  # Determine what plots to generate
  plots_cfg <- cfg$plots %||% list()
  enable_cfg <- plots_cfg$enable %||% list()
  
  enable_alpha <- enable_cfg$alpha %||% TRUE
  enable_beta <- enable_cfg$beta %||% TRUE
  enable_taxonomy <- enable_cfg$composition %||% TRUE
  enable_heatmap <- enable_cfg$heatmap %||% FALSE
  enable_tree_rect <- enable_cfg$tree_rectangular %||% FALSE
  enable_tree_circ <- enable_cfg$tree_circular %||% FALSE
  enable_group_comp <- enable_cfg$group_comparisons %||% FALSE
  enable_da_plots <- enable_cfg$differential_abundance %||% FALSE
  
  # Alpha diversity plots
  if (enable_alpha) {
    message("\n--- Alpha Diversity Plots ---")
    plot_alpha_diversity(ps, cfg, outdir)
    
    # Check for continuous variables
    meta <- as(sample_data(ps), "data.frame")
    numeric_cols <- names(meta)[sapply(meta, is.numeric)]
    if (length(numeric_cols) > 0) {
      plot_alpha_continuous(ps, cfg, outdir, numeric_cols)
    }
  }
  
  # Beta diversity plots
  if (enable_beta) {
    message("\n--- Beta Diversity Plots ---")
    plot_ordination(ps, cfg, outdir)
    
    # Run PERMANOVA if configured
    group_col <- cfg$metadata$group_column %||% "Group"
    formula_str <- paste0("dist_mat ~ ", group_col)
    run_permanova(ps, cfg, formula_str, distance = "bray")
  }
  
  # Taxonomy plots
  if (enable_taxonomy) {
    message("\n--- Taxonomy Plots ---")
    
    # Composition plots for multiple ranks
    ranks <- c("Phylum", "Class", "Order", "Family", "Genus")
    for (rank in ranks) {
      if (rank %in% rank_names(ps)) {
        plot_composition(ps, cfg, outdir, rank = rank, top_n = 10)
        
        # Top taxa barplot
        plot_top_taxa_barplot(ps, cfg, outdir, rank = rank, top_n = 20)
      }
    }
  }
  
  # Heatmaps
  if (enable_heatmap) {
    message("\n--- Taxonomy Heatmaps ---")
    heatmap_ranks <- plots_cfg$heatmap$taxonomic_ranks %||% c("Genus", "Family")
    heatmap_top_n <- plots_cfg$heatmap$top_n %||% 15
    
    message(sprintf("[viz] Generating heatmaps for %d ranks (top %d taxa each)", 
                    length(heatmap_ranks), heatmap_top_n))
    
    for (rank in heatmap_ranks) {
      if (rank %in% rank_names(ps)) {
        plot_taxa_heatmap(ps, cfg, outdir, rank = rank, top_n = heatmap_top_n)
      } else {
        message("[viz] Rank '", rank, "' not available in phyloseq object")
      }
    }
  }
  
  # Phylogenetic trees
  if (enable_tree_rect || enable_tree_circ) {
    message("\n--- Phylogenetic Trees ---")
    
    # Check if tree exists, if not try to build it
    tree <- phy_tree(ps, errorIfNULL = FALSE)
    if (is.null(tree) && isTRUE(cfg$amplicon$phylogeny$build_tree)) {
      message("[viz] No tree found, attempting to build...")
      ps <- build_phylogenetic_tree(ps, cfg)
    }
    
    # Plot trees if available
    if (!is.null(phy_tree(ps, errorIfNULL = FALSE))) {
      if (enable_tree_rect) {
        plot_tree_rectangular(ps, cfg, outdir)
      }
      if (enable_tree_circ) {
        plot_tree_circular(ps, cfg, outdir)
      }
    } else {
      message("[viz] Phylogenetic tree not available, skipping tree plots")
    }
  }
  
  # Group comparison plots
  if (enable_group_comp) {
    message("\n--- Group Comparison Plots ---")
    run_all_group_comparisons(ps, cfg, outdir)
  }
  
  # Differential abundance visualizations
  if (enable_da_plots) {
    message("\n--- Differential Abundance Plots ---")
    
    # Check if DA results exist
    da_files <- list.files(analysis_dir, pattern = "^deseq2_.*\\.csv$", full.names = TRUE)
    
    if (length(da_files) > 0) {
      message("[viz] Found ", length(da_files), " DESeq2 result files")
      
      # Plot volcano and top taxa for each comparison
      for (da_file in da_files) {
        if (grepl("_all_comparisons\\.csv$", da_file) || grepl("_summary\\.csv$", da_file)) {
          next  # Skip summary files
        }
        
        results <- read.csv(da_file)
        comp_name <- gsub("^deseq2_", "", gsub("\\.csv$", "", basename(da_file)))
        rank <- cfg$analysis$differential_abundance$taxonomic_rank %||% "Genus"
        
        # Volcano plot
        plot_volcano(results, cfg, outdir, comp_name, rank)
        
        # Top DA taxa
        plot_top_da_taxa(results, cfg, outdir, comp_name, rank, top_n = 20)
      }
    } else {
      message("[viz] No DESeq2 results found, skipping DA plots")
    }
  }
  
  # Correlation analysis
  message("\n--- Correlation Analysis ---")
  correlation_results <- run_all_correlations(ps, cfg, outdir)
  
  message("\n========================================")
  message("Visualization Complete!")
  message("Output: ", outdir)
  message("========================================")
}

# Allow running as standalone script (only when directly executed, not sourced)
if (!interactive() && length(sys.frames()) <= 1 && exists("cfg")) {
  run_plots(cfg)
}
