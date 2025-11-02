# Correlation Analysis Module
# Functions for correlating taxa abundance with continuous variables

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggrepel)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Correlate taxa abundance with continuous variable
#'
#' @param ps Phyloseq object
#' @param cfg Configuration list
#' @param outdir Output directory
#' @param cont_var Continuous variable name
#' @param rank Taxonomic rank
#' @param top_n Number of top taxa to analyze
#' @param method Correlation method (spearman, pearson)
#' @export
correlate_taxa_continuous <- function(ps, cfg, outdir, 
                                      cont_var = "dmft",
                                      rank = "Genus", 
                                      top_n = 15,
                                      method = "spearman") {
  
  message(sprintf("[corr] Correlating %s with %s (%s)", rank, cont_var, method))
  
  tryCatch({
    # Check if variable exists
    meta <- as(sample_data(ps), "data.frame")
    if (!cont_var %in% colnames(meta)) {
      message("[corr] Variable '", cont_var, "' not found in metadata")
      return(invisible(NULL))
    }
    
    # Remove NA values
    meta_clean <- meta[!is.na(meta[[cont_var]]), ]
    if (nrow(meta_clean) == 0) {
      message("[corr] All samples have NA for '", cont_var, "'")
      return(invisible(NULL))
    }
    
    if (nrow(meta_clean) < nrow(meta)) {
      message("[corr] Removing ", nrow(meta) - nrow(meta_clean), " samples with NA")
      ps <- prune_samples(rownames(meta_clean), ps)
    }
    
    # Agglomerate to rank
    ps_glom <- tax_glom(ps, taxrank = rank, NArm = FALSE)
    
    # Convert to relative abundance
    ps_rel <- transform_sample_counts(ps_glom, function(x) x / sum(x) * 100)
    
    # Get OTU table in correct orientation
    otu_mat_check <- as(otu_table(ps_rel), "matrix")
    if (!taxa_are_rows(ps_rel)) {
      otu_mat_check <- t(otu_mat_check)
    }
    
    # Get top taxa by mean abundance
    mean_abund <- rowMeans(otu_mat_check)
    top_taxa <- names(sort(mean_abund, decreasing = TRUE)[1:min(top_n, length(mean_abund))])
    
    # Subset to top taxa
    ps_top <- prune_taxa(top_taxa, ps_rel)
    
    message(sprintf("[corr] Analyzing %d taxa", ntaxa(ps_top)))
    
    # Extract data
    meta <- as(sample_data(ps_top), "data.frame")
    otu_mat <- as(otu_table(ps_top), "matrix")
    if (!taxa_are_rows(ps_top)) {
      otu_mat <- t(otu_mat)
    }
    
    tax_df <- as(tax_table(ps_top), "matrix")
    
    # Calculate correlations for each taxon
    results_list <- list()
    
    for (i in 1:nrow(otu_mat)) {
      taxon_abund <- otu_mat[i, ]
      taxon_name <- tax_df[i, rank]
      if (is.na(taxon_name)) taxon_name <- "Unclassified"
      
      # Only correlate if there's variation
      if (sd(taxon_abund) == 0) {
        next
      }
      
      # Correlation test
      cor_test <- cor.test(taxon_abund, meta[[cont_var]], 
                          method = method, exact = FALSE)
      
      results_list[[taxon_name]] <- data.frame(
        taxon = taxon_name,
        rho = cor_test$estimate,
        pvalue = cor_test$p.value,
        n = length(taxon_abund),
        mean_abundance = mean(taxon_abund),
        stringsAsFactors = FALSE
      )
    }
    
    if (length(results_list) == 0) {
      message("[corr] No taxa with sufficient variation")
      return(invisible(NULL))
    }
    
    # Combine results
    results_df <- do.call(rbind, results_list)
    results_df$padj <- p.adjust(results_df$pvalue, method = "fdr")
    results_df <- results_df[order(abs(results_df$rho), decreasing = TRUE), ]
    
    # Add significance categories
    results_df$significance <- "ns"
    results_df$significance[results_df$padj < 0.1] <- "p<0.1"
    results_df$significance[results_df$padj < 0.05] <- "p<0.05"
    results_df$significance[results_df$padj < 0.01] <- "p<0.01"
    
    # Mark significant for exploratory threshold
    results_df$significant_01 <- results_df$padj < 0.1
    results_df$significant_05 <- results_df$padj < 0.05
    
    # Save results
    outfile <- file.path(outdir, paste0("correlation_", tolower(rank), "_vs_", 
                                        tolower(cont_var), ".csv"))
    write.csv(results_df, outfile, row.names = FALSE)
    message(sprintf("[corr] Saved: %s", basename(outfile)))
    
    # Report findings
    sig_01 <- sum(results_df$significant_01, na.rm = TRUE)
    sig_05 <- sum(results_df$significant_05, na.rm = TRUE)
    
    if (sig_01 > 0) {
      message(sprintf("[corr] Found %d correlations (p<0.1), %d strong (p<0.05)", 
                     sig_01, sig_05))
      
      # Show top correlations
      top_corr <- head(results_df[results_df$significant_01, ], 5)
      for (j in 1:min(5, nrow(top_corr))) {
        message(sprintf("[corr]   %s: rho=%.3f, padj=%.3f", 
                       top_corr$taxon[j], 
                       top_corr$rho[j], 
                       top_corr$padj[j]))
      }
    } else {
      message("[corr] No significant correlations found")
    }
    
    # Create scatter plot matrix for significant taxa
    if (sig_01 > 0) {
      plot_correlation_scatter(ps_top, meta, cont_var, results_df, 
                              cfg, outdir, rank)
    }
    
    # Create correlation heatmap
    plot_correlation_heatmap(results_df, cont_var, cfg, outdir, rank)
    
    return(results_df)
    
  }, error = function(e) {
    message("[corr] Error in correlation analysis: ", e$message)
    return(invisible(NULL))
  })
}

#' Plot scatter plots for significant correlations
#'
#' @export
plot_correlation_scatter <- function(ps, meta, cont_var, results, 
                                    cfg, outdir, rank) {
  
  tryCatch({
    # Get significant taxa (p<0.1)
    sig_taxa <- results$taxon[results$significant_01]
    
    if (length(sig_taxa) == 0) return(invisible(NULL))
    
    # Limit to top 9 for 3x3 grid
    sig_taxa <- sig_taxa[1:min(9, length(sig_taxa))]
    
    # Extract data
    otu_mat <- as(otu_table(ps), "matrix")
    if (!taxa_are_rows(ps)) {
      otu_mat <- t(otu_mat)
    }
    
    tax_df <- as(tax_table(ps), "matrix")
    
    # Prepare data for plotting
    plot_data <- data.frame()
    
    for (taxon in sig_taxa) {
      idx <- which(tax_df[, rank] == taxon)
      if (length(idx) == 0) next
      
      abund <- otu_mat[idx[1], ]
      rho <- results$rho[results$taxon == taxon]
      padj <- results$padj[results$taxon == taxon]
      
      temp_df <- data.frame(
        taxon = taxon,
        abundance = abund,
        variable = meta[[cont_var]],
        rho = rho,
        padj = padj,
        label = sprintf("ρ=%.2f, p=%.3f", rho, padj),
        stringsAsFactors = FALSE
      )
      
      plot_data <- rbind(plot_data, temp_df)
    }
    
    # Create faceted plot
    p <- ggplot(plot_data, aes(x = variable, y = abundance)) +
      geom_point(alpha = 0.6, size = 2) +
      geom_smooth(method = "lm", se = TRUE, color = "#E69F00", linewidth = 1) +
      geom_text(
        data = plot_data %>% group_by(taxon) %>% slice(1),
        aes(label = label),
        x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
        size = 3, fontface = "bold"
      ) +
      facet_wrap(~ taxon, scales = "free_y", ncol = 3) +
      labs(
        title = paste("Correlation:", rank, "vs", cont_var),
        subtitle = paste("Significant correlations (p<0.1) shown"),
        x = cont_var,
        y = "Relative Abundance (%)"
      ) +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        strip.text = element_text(face = "italic", size = 10),
        strip.background = element_rect(fill = "grey90", color = "grey60")
      )
    
    # Save
    filename <- file.path(outdir, paste0("correlation_scatter_", tolower(rank), 
                                         "_vs_", tolower(cont_var), ".tiff"))
    ggsave(filename, p, width = 12, height = 10, dpi = 600, compression = "lzw")
    message("[corr] Saved scatter plot: ", basename(filename))
    
    # Also save PDF
    filename_pdf <- gsub("\\.tiff$", ".pdf", filename)
    ggsave(filename_pdf, p, width = 12, height = 10)
    
  }, error = function(e) {
    message("[corr] Error in scatter plot: ", e$message)
  })
}

#' Plot correlation heatmap
#'
#' @export
plot_correlation_heatmap <- function(results, cont_var, cfg, outdir, rank) {
  
  tryCatch({
    # Sort by correlation strength
    results <- results[order(results$rho), ]
    results$taxon <- factor(results$taxon, levels = results$taxon)
    
    # Create plot
    p <- ggplot(results, aes(x = cont_var, y = taxon, fill = rho)) +
      geom_tile(color = "white", linewidth = 1) +
      geom_text(aes(label = ifelse(padj < 0.05, "**", 
                                   ifelse(padj < 0.1, "*", ""))),
               size = 6, vjust = 0.75) +
      scale_fill_gradient2(
        low = "#3B9AB2", mid = "white", high = "#F21A00",
        midpoint = 0,
        limits = c(-1, 1),
        name = "Spearman's ρ"
      ) +
      labs(
        title = paste("Correlation Heatmap:", rank, "vs", cont_var),
        subtitle = "* p<0.1, ** p<0.05 (FDR-adjusted)",
        x = "",
        y = rank
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        axis.text.y = element_text(face = "italic"),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "right",
        panel.grid = element_blank()
      )
    
    # Save
    filename <- file.path(outdir, paste0("correlation_heatmap_", tolower(rank), 
                                         "_vs_", tolower(cont_var), ".tiff"))
    ggsave(filename, p, width = 8, height = max(6, nrow(results) * 0.3), 
           dpi = 600, compression = "lzw")
    message("[corr] Saved heatmap: ", basename(filename))
    
    # Also save PDF
    filename_pdf <- gsub("\\.tiff$", ".pdf", filename)
    ggsave(filename_pdf, p, width = 8, height = max(6, nrow(results) * 0.3))
    
  }, error = function(e) {
    message("[corr] Error in heatmap: ", e$message)
  })
}

#' Run all correlation analyses
#'
#' @param ps Phyloseq object
#' @param cfg Configuration list
#' @param outdir Output directory
#' @export
run_all_correlations <- function(ps, cfg, outdir) {
  
  message("\n=== CORRELATION ANALYSIS ===")
  
  # Get continuous variables
  cont_vars <- cfg$metadata$continuous_variables
  if (is.null(cont_vars)) {
    message("[corr] No continuous variables configured")
    return(invisible(NULL))
  }
  
  # Extract variable names
  var_names <- sapply(cont_vars, function(x) {
    if (is.list(x) && isTRUE(x$enabled)) x$name else NULL
  })
  var_names <- unlist(var_names[!sapply(var_names, is.null)])
  
  if (length(var_names) == 0) {
    message("[corr] No enabled continuous variables")
    return(invisible(NULL))
  }
  
  # Get ranks to analyze
  ranks <- c("Phylum", "Class", "Order", "Family", "Genus")
  available_ranks <- rank_names(ps)
  ranks <- ranks[ranks %in% available_ranks]
  
  message("[corr] Variables: ", paste(var_names, collapse = ", "))
  message("[corr] Ranks: ", paste(ranks, collapse = ", "))
  
  # Run correlations for each combination
  all_results <- list()
  count <- 0
  
  for (var in var_names) {
    for (rank in ranks) {
      result <- correlate_taxa_continuous(ps, cfg, outdir, var, rank, top_n = 15)
      if (!is.null(result)) {
        all_results[[paste(rank, var, sep = "_")]] <- result
        count <- count + 1
      }
    }
  }
  
  message(sprintf("[corr] Completed %d correlation analyses", count))
  
  return(all_results)
}
