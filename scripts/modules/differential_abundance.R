# Differential Abundance Analysis Module
# Functions for identifying differentially abundant taxa between groups

suppressPackageStartupMessages({
  library(phyloseq)
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Run DESeq2 differential abundance analysis
#'
#' @param ps Phyloseq object
#' @param cfg Configuration list
#' @param outdir Output directory for results
#' @param group_col Grouping variable
#' @param rank Taxonomic rank to test
#' @export
run_deseq2_analysis <- function(ps, cfg, outdir, 
                                group_col = "Group", 
                                rank = "Genus") {
  
  message("[DA] ======================================")
  message("[DA] Starting DESeq2 analysis...")
  message("[DA] Grouping variable: ", group_col)
  message("[DA] Taxonomic rank: ", rank)
  message("[DA] Output directory: ", outdir)
  message("[DA] ======================================")
  
  tryCatch({
    # Check inputs
    if (is.null(ps) || !inherits(ps, "phyloseq")) {
      stop("[DA] ERROR: Invalid phyloseq object")
    }
    message(sprintf("[DA] Input: %d samples, %d taxa", nsamples(ps), ntaxa(ps)))
    
    # Agglomerate to rank
    message("[DA] Agglomerating to rank: ", rank)
    ps_glom <- tax_glom(ps, taxrank = rank, NArm = FALSE)
    
    # Filter low abundance/prevalence taxa
    min_prev <- cfg$analysis$differential_abundance$min_prevalence %||% 0.1
    min_abund <- cfg$analysis$differential_abundance$min_abundance %||% 0.0001
    
    message(sprintf("[DA] Filtering criteria: prevalence >= %.1f%%, abundance >= %.4f%%", 
                    min_prev * 100, min_abund * 100))
    
    # Convert to relative abundance for filtering
    ps_rel <- transform_sample_counts(ps_glom, function(x) x / sum(x))
    
    # Get OTU table in correct orientation (taxa as rows)
    otu_mat <- as(otu_table(ps_rel), "matrix")
    if (!taxa_are_rows(ps_rel)) {
      otu_mat <- t(otu_mat)
    }
    
    # Keep taxa present in >min_prev samples with >min_abund abundance
    prev_threshold <- min_prev * nsamples(ps_rel)
    keep_taxa <- apply(otu_mat, 1, function(x) {
      sum(x > min_abund) >= prev_threshold
    })
    
    message(sprintf("[DA] Taxa passing filter: %d / %d", sum(keep_taxa), length(keep_taxa)))
    
    ps_filt <- prune_taxa(keep_taxa, ps_glom)
    message("[DA] Filtered to ", ntaxa(ps_filt), " taxa (from ", ntaxa(ps_glom), ")")
    
    if (ntaxa(ps_filt) == 0) {
      message("[DA] No taxa passed filtering criteria")
      return(invisible(NULL))
    }
    
    # Check if group_col exists and has enough samples
    meta <- as(sample_data(ps_filt), "data.frame")
    if (!group_col %in% colnames(meta)) {
      message("[DA] Grouping variable '", group_col, "' not found in metadata")
      return(invisible(NULL))
    }
    
    # Remove NA values
    meta[[group_col]] <- as.factor(meta[[group_col]])
    meta_clean <- meta[!is.na(meta[[group_col]]), ]
    
    if (nrow(meta_clean) < nrow(meta)) {
      message("[DA] Removed ", nrow(meta) - nrow(meta_clean), " samples with NA in ", group_col)
      ps_filt <- prune_samples(rownames(meta_clean), ps_filt)
    }
    
    # Check we have at least 2 groups with multiple samples
    group_counts <- table(sample_data(ps_filt)[[group_col]])
    if (length(group_counts) < 2) {
      message("[DA] Need at least 2 groups for comparison")
      return(invisible(NULL))
    }
    if (any(group_counts < 3)) {
      message("[DA] Warning: some groups have <3 samples, results may be unreliable")
    }
    
    # Convert to DESeq2 format
    message("[DA] Converting to DESeq2 format...")
    dds <- phyloseq_to_deseq2(ps_filt, as.formula(paste("~", group_col)))
    
    # Calculate geometric means for size factors
    gm_mean <- function(x, na.rm = TRUE) {
      exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
    }
    geoMeans <- apply(counts(dds), 1, gm_mean)
    dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
    
    # Run DESeq2
    message("[DA] Running DESeq2...")
    dds <- DESeq(dds, fitType = "local", quiet = TRUE)
    
    # Get all pairwise comparisons
    groups <- levels(sample_data(ps_filt)[[group_col]])
    message("[DA] Found ", length(groups), " groups: ", paste(groups, collapse = ", "))
    
    all_results <- list()
    comparison_names <- c()
    
    for (i in 1:(length(groups) - 1)) {
      for (j in (i + 1):length(groups)) {
        group1 <- groups[i]
        group2 <- groups[j]
        comp_name <- paste0(group1, "_vs_", group2)
        comparison_names <- c(comparison_names, comp_name)
        
        message("[DA] Testing: ", group1, " vs ", group2)
        
        res <- results(dds, contrast = c(group_col, group1, group2))
        res <- res[order(res$padj, na.last = TRUE), ]
        
        # Convert to data frame with taxonomy
        res_df <- as.data.frame(res)
        res_df$ASV <- rownames(res_df)
        
        # Add taxonomy
        tax_table_df <- as.data.frame(tax_table(ps_filt))
        res_df <- merge(res_df, tax_table_df, by.x = "ASV", by.y = "row.names", all.x = TRUE)
        
        # Add comparison info
        res_df$comparison <- comp_name
        res_df$group1 <- group1
        res_df$group2 <- group2
        
        # Add significance category
        res_df$significance <- "ns"
        log2fc_threshold <- cfg$analysis$differential_abundance$log2fc_threshold %||% 1.0
        res_df$significance[!is.na(res_df$padj) & res_df$padj < 0.05] <- "p<0.05"
        res_df$significance[!is.na(res_df$padj) & res_df$padj < 0.01] <- "p<0.01"
        res_df$significance[!is.na(res_df$padj) & res_df$padj < 0.001] <- "p<0.001"
        
        # Mark biologically significant
        res_df$significant <- !is.na(res_df$padj) & 
                              res_df$padj < 0.05 & 
                              abs(res_df$log2FoldChange) > log2fc_threshold
        
        all_results[[comp_name]] <- res_df
        
        # Save individual comparison
        outfile <- file.path(outdir, paste0("deseq2_", comp_name, ".csv"))
        write.csv(res_df, outfile, row.names = FALSE)
        message("[DA] Saved: ", basename(outfile))
        
        # Report significant taxa
        sig_taxa <- sum(res_df$significant, na.rm = TRUE)
        if (sig_taxa > 0) {
          message("[DA]   Found ", sig_taxa, " significantly different taxa")
          top5 <- head(res_df[res_df$significant & !is.na(res_df$significant), ], 5)
          for (k in 1:min(5, nrow(top5))) {
            taxa_name <- top5[[rank]][k]
            if (is.na(taxa_name)) taxa_name <- "Unclassified"
            message(sprintf("[DA]     %s: log2FC=%.2f, padj=%.2e", 
                           taxa_name, 
                           top5$log2FoldChange[k], 
                           top5$padj[k]))
          }
        } else {
          message("[DA]   No significantly different taxa (padj<0.05, |log2FC|>", log2fc_threshold, ")")
        }
      }
    }
    
    # Combine all results
    combined_results <- do.call(rbind, all_results)
    outfile_combined <- file.path(outdir, "deseq2_all_comparisons.csv")
    write.csv(combined_results, outfile_combined, row.names = FALSE)
    message("[DA] Saved combined results: ", basename(outfile_combined))
    
    # Generate summary
    summary_stats <- combined_results %>%
      group_by(comparison) %>%
      summarize(
        total_taxa = n(),
        significant_taxa = sum(significant, na.rm = TRUE),
        upregulated = sum(significant & log2FoldChange > 0, na.rm = TRUE),
        downregulated = sum(significant & log2FoldChange < 0, na.rm = TRUE),
        .groups = 'drop'
      )
    
    summary_file <- file.path(outdir, "deseq2_summary.csv")
    write.csv(summary_stats, summary_file, row.names = FALSE)
    message("[DA] Saved summary: ", basename(summary_file))
    
    message("[DA] DESeq2 analysis complete!")
    
    return(list(
      results = all_results,
      combined = combined_results,
      summary = summary_stats,
      dds = dds
    ))
    
  }, error = function(e) {
    message("[DA] Error in DESeq2 analysis: ", e$message)
    message("[DA] This may occur with low sample sizes or convergence issues")
    return(invisible(NULL))
  })
}

#' Plot volcano plot for differential abundance
#'
#' @param results DESeq2 results data frame
#' @param cfg Configuration list
#' @param outdir Output directory
#' @param comparison_name Name of comparison
#' @export
plot_volcano <- function(results, cfg, outdir, comparison_name, rank = "Genus") {
  tryCatch({
    # Remove NA values
    res_clean <- results[!is.na(results$padj) & !is.na(results$log2FoldChange), ]
    
    if (nrow(res_clean) == 0) {
      message("[DA] No valid data for volcano plot")
      return(invisible(NULL))
    }
    
    # Add labels for significant taxa
    res_clean$label <- ""
    sig_taxa <- res_clean[res_clean$significant & !is.na(res_clean$significant), ]
    if (nrow(sig_taxa) > 0) {
      top_sig <- head(sig_taxa[order(sig_taxa$padj), ], 10)
      res_clean$label[match(top_sig$ASV, res_clean$ASV)] <- top_sig[[rank]]
    }
    
    # Create plot
    log2fc_threshold <- cfg$analysis$differential_abundance$log2fc_threshold %||% 1.0
    
    p <- ggplot(res_clean, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
      geom_point(alpha = 0.6, size = 2) +
      geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), 
                 linetype = "dashed", color = "gray40") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
      scale_color_manual(
        values = c("ns" = "gray60", "p<0.05" = "skyblue", 
                   "p<0.01" = "orange", "p<0.001" = "red"),
        name = "Significance"
      ) +
      labs(
        title = paste("Differential Abundance:", gsub("_", " ", comparison_name)),
        subtitle = paste0("DESeq2 analysis at ", rank, " level"),
        x = "Log2 Fold Change",
        y = "-Log10(adjusted p-value)"
      ) +
      theme_classic(base_size = 14) +
      theme(
        legend.position = "right",
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)
      )
    
    # Add labels if ggrepel available
    if (requireNamespace("ggrepel", quietly = TRUE) && any(res_clean$label != "")) {
      p <- p + ggrepel::geom_text_repel(
        data = res_clean[res_clean$label != "", ],
        aes(label = label),
        size = 3,
        max.overlaps = 15,
        box.padding = 0.5
      )
    }
    
    # Save
    filename <- file.path(outdir, paste0("volcano_", comparison_name, ".tiff"))
    ggsave(filename, p, width = 10, height = 8, dpi = 600, compression = "lzw")
    message("[DA] Saved volcano plot: ", basename(filename))
    
    # Also save PDF
    filename_pdf <- file.path(outdir, paste0("volcano_", comparison_name, ".pdf"))
    ggsave(filename_pdf, p, width = 10, height = 8)
    
    return(p)
    
  }, error = function(e) {
    message("[DA] Failed to create volcano plot: ", e$message)
  })
}

#' Plot top differentially abundant taxa
#'
#' @param results DESeq2 results data frame
#' @param cfg Configuration list
#' @param outdir Output directory
#' @param comparison_name Name of comparison
#' @param top_n Number of top taxa to show
#' @export
plot_top_da_taxa <- function(results, cfg, outdir, comparison_name, rank = "Genus", top_n = 20) {
  tryCatch({
    # Get significant taxa
    sig_taxa <- results[results$significant & !is.na(results$significant), ]
    
    if (nrow(sig_taxa) == 0) {
      message("[DA] No significant taxa to plot for ", comparison_name)
      return(invisible(NULL))
    }
    
    # Order by absolute log2FC and take top N
    sig_taxa <- sig_taxa[order(-abs(sig_taxa$log2FoldChange)), ]
    top_taxa <- head(sig_taxa, top_n)
    
    # Prepare labels
    top_taxa$taxa_label <- top_taxa[[rank]]
    top_taxa$taxa_label[is.na(top_taxa$taxa_label)] <- "Unclassified"
    top_taxa$taxa_label <- factor(top_taxa$taxa_label, 
                                   levels = top_taxa$taxa_label[order(top_taxa$log2FoldChange)])
    
    # Create plot
    p <- ggplot(top_taxa, aes(x = taxa_label, y = log2FoldChange, fill = log2FoldChange > 0)) +
      geom_col() +
      coord_flip() +
      scale_fill_manual(
        values = c("TRUE" = "#D55E00", "FALSE" = "#0072B2"),
        labels = c("TRUE" = paste("Enriched in", top_taxa$group1[1]), 
                   "FALSE" = paste("Enriched in", top_taxa$group2[1])),
        name = ""
      ) +
      labs(
        title = paste("Top Differentially Abundant Taxa:", gsub("_", " ", comparison_name)),
        subtitle = paste0("DESeq2 analysis at ", rank, " level (padj < 0.05)"),
        x = rank,
        y = "Log2 Fold Change"
      ) +
      theme_classic(base_size = 14) +
      theme(
        legend.position = "top",
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)
      )
    
    # Save
    filename <- file.path(outdir, paste0("top_da_taxa_", comparison_name, ".tiff"))
    ggsave(filename, p, width = 10, height = max(8, nrow(top_taxa) * 0.3), 
           dpi = 600, compression = "lzw")
    message("[DA] Saved top DA taxa plot: ", basename(filename))
    
    # Also save PDF
    filename_pdf <- file.path(outdir, paste0("top_da_taxa_", comparison_name, ".pdf"))
    ggsave(filename_pdf, p, width = 10, height = max(8, nrow(top_taxa) * 0.3))
    
    return(p)
    
  }, error = function(e) {
    message("[DA] Failed to create top DA taxa plot: ", e$message)
  })
}
