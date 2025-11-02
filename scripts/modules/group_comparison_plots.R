# Group Comparison Plots Module
# Functions for comparing taxonomic abundance across metadata grouping variables

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(forcats)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Compare taxonomic abundance across groups for a specific variable
#'
#' @param ps Phyloseq object
#' @param cfg Configuration list
#' @param outdir Output directory
#' @param group_var Grouping variable name
#' @param rank Taxonomic rank
#' @param top_n Number of top taxa to show
#' @export
plot_group_comparison <- function(ps, cfg, outdir, group_var, rank = "Genus", top_n = 20) {
  
  message("[group-comp] Comparing ", rank, " across ", group_var)
  
  tryCatch({
    # Check if variable exists
    meta <- as(sample_data(ps), "data.frame")
    if (!group_var %in% colnames(meta)) {
      message("[group-comp] Variable '", group_var, "' not found in metadata")
      return(invisible(NULL))
    }
    
    # Remove NA values
    meta_clean <- meta[!is.na(meta[[group_var]]), ]
    if (nrow(meta_clean) == 0) {
      message("[group-comp] All samples have NA for '", group_var, "'")
      return(invisible(NULL))
    }
    
    if (nrow(meta_clean) < nrow(meta)) {
      message("[group-comp] Removing ", nrow(meta) - nrow(meta_clean), " samples with NA")
      ps <- prune_samples(rownames(meta_clean), ps)
    }
    
    # Agglomerate to rank
    ps_glom <- tax_glom(ps, taxrank = rank, NArm = FALSE)
    
    # Convert to relative abundance
    ps_rel <- transform_sample_counts(ps_glom, function(x) x / sum(x) * 100)
    
    # Get top taxa by mean abundance
    mean_abund <- apply(otu_table(ps_rel), 1, mean)
    top_taxa <- names(sort(mean_abund, decreasing = TRUE)[1:min(top_n, length(mean_abund))])
    
    # Subset to top taxa
    ps_top <- prune_taxa(top_taxa, ps_rel)
    
    # Melt to long format
    df <- psmelt(ps_top)
    
    # Get taxa names
    df$taxa_label <- df[[rank]]
    df$taxa_label[is.na(df$taxa_label)] <- "Unclassified"
    
    # Reorder taxa by mean abundance
    taxa_order <- df %>%
      group_by(taxa_label) %>%
      summarize(mean_abund = mean(Abundance), .groups = 'drop') %>%
      arrange(desc(mean_abund)) %>%
      pull(taxa_label)
    
    df$taxa_label <- factor(df$taxa_label, levels = rev(taxa_order))
    df[[group_var]] <- as.factor(df[[group_var]])
    
    # Determine plot type
    plot_type <- cfg$plots$group_comparison$plot_type %||% "boxplot"
    
    # Create base plot
    p <- ggplot(df, aes(x = taxa_label, y = Abundance, fill = .data[[group_var]]))
    
    if (plot_type == "boxplot") {
      p <- p + geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.8))
    } else if (plot_type == "violin") {
      p <- p + geom_violin(trim = FALSE, position = position_dodge(width = 0.8)) +
        geom_boxplot(width = 0.1, outlier.size = 0.5, 
                     position = position_dodge(width = 0.8),
                     fill = "white", alpha = 0.5)
    } else {  # barplot with error bars
      df_summary <- df %>%
        group_by(taxa_label, .data[[group_var]]) %>%
        summarize(
          mean_abund = mean(Abundance),
          se = sd(Abundance) / sqrt(n()),
          .groups = 'drop'
        )
      
      p <- ggplot(df_summary, aes(x = taxa_label, y = mean_abund, fill = .data[[group_var]])) +
        geom_col(position = position_dodge(width = 0.8)) +
        geom_errorbar(
          aes(ymin = mean_abund - se, ymax = mean_abund + se),
          width = 0.2,
          position = position_dodge(width = 0.8)
        )
    }
    
    # Add statistical comparisons if requested and ggpubr available
    show_stats <- cfg$plots$group_comparison$show_stats %||% TRUE
    if (show_stats && requireNamespace("ggpubr", quietly = TRUE)) {
      groups <- levels(df[[group_var]])
      
      tryCatch({
        if (length(groups) == 2) {
          # Pairwise comparison for 2 groups
          p <- p + ggpubr::stat_compare_means(
            aes(group = .data[[group_var]]),
            method = "wilcox.test",
            label = "p.format",
            label.x.npc = 0.5,
            size = 3.5,
            hide.ns = FALSE  # Show all p-values
          )
        } else if (length(groups) > 2) {
          # Kruskal-Wallis for 3+ groups
          p <- p + ggpubr::stat_compare_means(
            aes(group = .data[[group_var]]),
            method = "kruskal.test",
            label = "p.format",
            label.x.npc = 0.5,
            label.y.npc = 0.95,
            size = 3.5
          )
        }
      }, error = function(e) {
        message("[group-comp] Could not add statistical annotations: ", e$message)
      })
    }
    
    # Finalize plot
    p <- p +
      coord_flip() +
      scale_fill_manual(
        values = get_palette(length(unique(df[[group_var]])), cfg),
        name = group_var
      ) +
      labs(
        title = paste(rank, "Abundance by", group_var),
        subtitle = paste("Top", min(top_n, ntaxa(ps_top)), "most abundant taxa"),
        x = rank,
        y = "Relative Abundance (%)"
      ) +
      theme_classic(base_size = 14) +
      theme(
        legend.position = "right",
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.y = element_text(face = "italic")
      )
    
    # Save
    filename <- file.path(outdir, paste0("group_comparison_", tolower(rank), "_by_", 
                                         tolower(group_var), ".tiff"))
    ggsave(filename, p, width = 12, height = max(8, top_n * 0.3), 
           dpi = 600, compression = "lzw")
    message("[group-comp] Saved: ", basename(filename))
    
    # Also save PDF
    filename_pdf <- gsub("\\.tiff$", ".pdf", filename)
    ggsave(filename_pdf, p, width = 12, height = max(8, top_n * 0.3))
    
    return(p)
    
  }, error = function(e) {
    message("[group-comp] Error comparing ", rank, " by ", group_var, ": ", e$message)
  })
}

#' Create mean abundance comparison plot
#'
#' @param ps Phyloseq object
#' @param cfg Configuration list
#' @param outdir Output directory
#' @param group_var Grouping variable
#' @param rank Taxonomic rank
#' @param top_n Number of taxa to show
#' @export
plot_mean_abundance_comparison <- function(ps, cfg, outdir, group_var, rank = "Genus", top_n = 15) {
  
  message("[group-comp] Creating mean abundance comparison for ", group_var)
  
  tryCatch({
    # Check variable
    meta <- as(sample_data(ps), "data.frame")
    if (!group_var %in% colnames(meta) || all(is.na(meta[[group_var]]))) {
      message("[group-comp] Skipping: ", group_var, " not available")
      return(invisible(NULL))
    }
    
    # Clean and agglomerate
    ps_clean <- prune_samples(!is.na(meta[[group_var]]), ps)
    ps_glom <- tax_glom(ps_clean, taxrank = rank, NArm = FALSE)
    ps_rel <- transform_sample_counts(ps_glom, function(x) x / sum(x) * 100)
    
    # Calculate mean abundance per group
    df <- psmelt(ps_rel)
    df$taxa_label <- df[[rank]]
    df$taxa_label[is.na(df$taxa_label)] <- "Unclassified"
    
    # Get top taxa overall
    top_taxa <- df %>%
      group_by(taxa_label) %>%
      summarize(total_abund = sum(Abundance), .groups = 'drop') %>%
      arrange(desc(total_abund)) %>%
      head(top_n) %>%
      pull(taxa_label)
    
    # Calculate mean and SE per group
    df_summary <- df %>%
      filter(taxa_label %in% top_taxa) %>%
      group_by(taxa_label, .data[[group_var]]) %>%
      summarize(
        mean_abund = mean(Abundance),
        se = sd(Abundance) / sqrt(n()),
        n_samples = n(),
        .groups = 'drop'
      )
    
    # Order taxa by overall abundance
    taxa_order <- df %>%
      filter(taxa_label %in% top_taxa) %>%
      group_by(taxa_label) %>%
      summarize(mean_abund = mean(Abundance), .groups = 'drop') %>%
      arrange(mean_abund) %>%
      pull(taxa_label)
    
    df_summary$taxa_label <- factor(df_summary$taxa_label, levels = taxa_order)
    
    # Create plot
    p <- ggplot(df_summary, aes(x = mean_abund, y = taxa_label, fill = .data[[group_var]])) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      geom_errorbar(
        aes(xmin = mean_abund - se, xmax = mean_abund + se),
        width = 0.3,
        position = position_dodge(width = 0.8)
      ) +
      scale_fill_manual(
        values = get_palette(length(unique(df_summary[[group_var]])), cfg),
        name = group_var
      ) +
      labs(
        title = paste("Mean", rank, "Abundance by", group_var),
        subtitle = paste("Top", length(top_taxa), "taxa (mean ± SE)"),
        x = "Mean Relative Abundance (%)",
        y = rank
      ) +
      theme_classic(base_size = 14) +
      theme(
        legend.position = "right",
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.y = element_text(face = "italic")
      )
    
    # Save
    filename <- file.path(outdir, paste0("mean_abundance_", tolower(rank), "_by_", 
                                         tolower(group_var), ".tiff"))
    ggsave(filename, p, width = 12, height = max(8, top_n * 0.4), 
           dpi = 600, compression = "lzw")
    message("[group-comp] Saved mean abundance plot: ", basename(filename))
    
    # Also save PDF
    filename_pdf <- gsub("\\.tiff$", ".pdf", filename)
    ggsave(filename_pdf, p, width = 12, height = max(8, top_n * 0.4))
    
    return(p)
    
  }, error = function(e) {
    message("[group-comp] Error in mean abundance plot: ", e$message)
  })
}

#' Run all group comparisons configured in the analysis
#'
#' @param ps Phyloseq object
#' @param cfg Configuration list
#' @param outdir Output directory
#' @export
run_all_group_comparisons <- function(ps, cfg, outdir) {
  
  if (!isTRUE(cfg$analysis$group_comparisons$enabled)) {
    message("[group-comp] Group comparisons not enabled")
    return(invisible(NULL))
  }
  
  message("\n=== TAXONOMIC GROUP COMPARISONS ===")
  
  # Get variables to compare
  variables <- cfg$analysis$group_comparisons$variables %||% 
               c(cfg$metadata$primary_comparison$group_column)
  
  # Get ranks to analyze
  ranks <- cfg$analysis$group_comparisons$taxonomic_ranks %||% 
           c("Phylum", "Family", "Genus")
  
  top_n <- cfg$analysis$group_comparisons$top_n %||% 20
  
  # Filter to available ranks
  available_ranks <- rank_names(ps)
  ranks <- ranks[ranks %in% available_ranks]
  
  if (length(ranks) == 0) {
    message("[group-comp] No specified ranks available in data")
    return(invisible(NULL))
  }
  
  message("[group-comp] Variables: ", paste(variables, collapse = ", "))
  message("[group-comp] Ranks: ", paste(ranks, collapse = ", "))
  
  # Generate plots for each combination
  plot_count <- 0
  for (var in variables) {
    for (rank in ranks) {
      # Boxplot/violin comparison
      plot_group_comparison(ps, cfg, outdir, var, rank, top_n)
      
      # Mean abundance comparison
      plot_mean_abundance_comparison(ps, cfg, outdir, var, rank, 
                                     min(15, top_n))
      
      plot_count <- plot_count + 2
    }
  }
  
  message("[group-comp] Generated ", plot_count, " group comparison plots")
  message("===================================\n")
}

#' Get color palette helper
#' @keywords internal
get_palette <- function(n, cfg) {
  pal <- tolower(cfg$visualization$color_palette %||% "okabe-ito")
  
  if (pal == "okabe-ito") {
    cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
              "#0072B2", "#D55E00", "#CC79A7", "#999999")
    return(rep(cols, length.out = n))
  } else if (pal == "viridis") {
    return(viridisLite::viridis(n))
  } else {
    return(scales::hue_pal()(n))
  }
}
