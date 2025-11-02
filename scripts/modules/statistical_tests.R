#!/usr/bin/env Rscript
# Statistical Testing Module with Multiple Testing Corrections
# Provides corrected p-values for multiple comparisons

suppressPackageStartupMessages({
  library(phyloseq)
})

#' Perform alpha diversity statistical tests with multiple testing correction
#'
#' @param ps Phyloseq object
#' @param group_column Column name for grouping variable
#' @param measures Vector of diversity measures to test
#' @param correction Method for p-value correction ("fdr", "bonferroni", "holm", "none")
#' @return Data frame with test results and corrected p-values
#' @export
test_alpha_diversity <- function(ps, group_column, 
                                 measures = c("Shannon", "Observed", "Simpson"),
                                 correction = "fdr") {
  tryCatch({
    # Calculate diversity
    df_alpha <- estimate_richness(ps, measures = measures)
    
    # Get metadata and add group column
    meta <- as(sample_data(ps), "data.frame")
    alpha_ids <- gsub("\\.", "-", rownames(df_alpha))
    meta_ids <- rownames(meta)
    df_alpha[[group_column]] <- meta[[group_column]][match(alpha_ids, meta_ids)]
    
    # Remove NA groups
    df_alpha <- df_alpha[!is.na(df_alpha[[group_column]]), ]
    
    # Initialize results
    results <- data.frame(
      measure = character(),
      test = character(),
      statistic = numeric(),
      p_value = numeric(),
      p_adjusted = numeric(),
      n_groups = integer(),
      correction_method = character(),
      stringsAsFactors = FALSE
    )
    
    # Determine number of groups
    n_groups <- length(unique(df_alpha[[group_column]]))
    
    if (n_groups < 2) {
      warning("Need at least 2 groups for statistical testing")
      return(NULL)
    }
    
    # Test each measure
    for (meas in measures) {
      if (!meas %in% colnames(df_alpha)) next
      
      # Choose appropriate test
      if (n_groups == 2) {
        # Wilcoxon rank-sum test for 2 groups
        test_result <- wilcox.test(
          as.formula(paste(meas, "~", group_column)),
          data = df_alpha
        )
        test_name <- "Wilcoxon"
        statistic <- test_result$statistic
      } else {
        # Kruskal-Wallis test for > 2 groups
        test_result <- kruskal.test(
          as.formula(paste(meas, "~", group_column)),
          data = df_alpha
        )
        test_name <- "Kruskal-Wallis"
        statistic <- test_result$statistic
      }
      
      # Store result
      results <- rbind(results, data.frame(
        measure = meas,
        test = test_name,
        statistic = as.numeric(statistic),
        p_value = test_result$p.value,
        p_adjusted = NA,  # Will be filled after correction
        n_groups = n_groups,
        correction_method = correction,
        stringsAsFactors = FALSE
      ))
    }
    
    # Apply multiple testing correction
    if (nrow(results) > 0 && correction != "none") {
      results$p_adjusted <- p.adjust(results$p_value, method = correction)
    } else {
      results$p_adjusted <- results$p_value
    }
    
    # Add significance flags
    results$significant_raw <- results$p_value < 0.05
    results$significant_adjusted <- results$p_adjusted < 0.05
    
    return(results)
    
  }, error = function(e) {
    warning("Alpha diversity testing failed: ", e$message)
    return(NULL)
  })
}

#' Perform pairwise PERMANOVA tests with multiple testing correction
#'
#' @param ps Phyloseq object
#' @param group_column Column name for grouping variable
#' @param distance Distance method
#' @param correction Method for p-value correction ("fdr", "bonferroni", "holm", "none")
#' @return Data frame with pairwise test results and corrected p-values
#' @export
pairwise_permanova <- function(ps, group_column, distance = "bray", correction = "fdr") {
  tryCatch({
    # Get metadata
    meta <- as(sample_data(ps), "data.frame")
    
    # Get groups
    groups <- unique(meta[[group_column]])
    groups <- groups[!is.na(groups)]
    n_groups <- length(groups)
    
    if (n_groups < 2) {
      warning("Need at least 2 groups for pairwise testing")
      return(NULL)
    }
    
    # Calculate distance matrix once
    dist_mat <- phyloseq::distance(ps, method = distance)
    
    # Initialize results
    results <- data.frame(
      group1 = character(),
      group2 = character(),
      F_statistic = numeric(),
      R2 = numeric(),
      p_value = numeric(),
      p_adjusted = numeric(),
      n1 = integer(),
      n2 = integer(),
      correction_method = character(),
      stringsAsFactors = FALSE
    )
    
    # Perform pairwise comparisons
    for (i in 1:(n_groups - 1)) {
      for (j in (i + 1):n_groups) {
        g1 <- groups[i]
        g2 <- groups[j]
        
        # Subset to these two groups
        keep_samples <- meta[[group_column]] %in% c(g1, g2)
        ps_sub <- prune_samples(keep_samples, ps)
        meta_sub <- as(sample_data(ps_sub), "data.frame")
        
        # Subset distance matrix
        sample_names_sub <- sample_names(ps_sub)
        dist_sub <- as.dist(as.matrix(dist_mat)[sample_names_sub, sample_names_sub])
        
        # Run PERMANOVA
        perm <- vegan::adonis2(
          as.formula(paste("dist_sub ~", group_column)),
          data = meta_sub,
          permutations = 999
        )
        
        # Extract results
        results <- rbind(results, data.frame(
          group1 = as.character(g1),
          group2 = as.character(g2),
          F_statistic = perm$F[1],
          R2 = perm$R2[1],
          p_value = perm$`Pr(>F)`[1],
          p_adjusted = NA,  # Will be filled after correction
          n1 = sum(meta_sub[[group_column]] == g1),
          n2 = sum(meta_sub[[group_column]] == g2),
          correction_method = correction,
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # Apply multiple testing correction
    if (nrow(results) > 0 && correction != "none") {
      results$p_adjusted <- p.adjust(results$p_value, method = correction)
    } else {
      results$p_adjusted <- results$p_value
    }
    
    # Add significance flags
    results$significant_raw <- results$p_value < 0.05
    results$significant_adjusted <- results$p_adjusted < 0.05
    
    return(results)
    
  }, error = function(e) {
    warning("Pairwise PERMANOVA failed: ", e$message)
    return(NULL)
  })
}

#' Print statistical test results with formatting
#'
#' @param results Data frame from test_alpha_diversity() or pairwise_permanova()
#' @param test_type Type of test ("alpha" or "pairwise_permanova")
#' @export
print_statistical_results <- function(results, test_type = "alpha") {
  if (is.null(results) || nrow(results) == 0) {
    cat("No results to display\n")
    return(invisible(NULL))
  }
  
  cat("\n")
  cat("══════════════════════════════════════════════════════════\n")
  
  if (test_type == "alpha") {
    cat("Alpha Diversity Statistical Tests\n")
    cat("══════════════════════════════════════════════════════════\n")
    cat("Correction method:", unique(results$correction_method), "\n")
    cat("Number of tests:", nrow(results), "\n\n")
    
    for (i in 1:nrow(results)) {
      row <- results[i, ]
      cat(sprintf("%-15s %s\n", row$measure, row$test))
      cat(sprintf("  Statistic: %.3f\n", row$statistic))
      cat(sprintf("  p-value (raw): %.4f %s\n", 
                 row$p_value,
                 if (!is.na(row$significant_raw) && row$significant_raw) "*" else ""))
      cat(sprintf("  p-value (adj): %.4f %s\n",
                 row$p_adjusted,
                 if (!is.na(row$significant_adjusted) && row$significant_adjusted) "*" else ""))
      cat("\n")
    }
    
  } else if (test_type == "pairwise_permanova") {
    cat("Pairwise PERMANOVA Tests\n")
    cat("══════════════════════════════════════════════════════════\n")
    cat("Correction method:", unique(results$correction_method), "\n")
    cat("Number of comparisons:", nrow(results), "\n\n")
    
    for (i in 1:nrow(results)) {
      row <- results[i, ]
      cat(sprintf("%s vs %s (n=%d, n=%d)\n", 
                 row$group1, row$group2, row$n1, row$n2))
      cat(sprintf("  F = %.3f, R² = %.4f\n", row$F_statistic, row$R2))
      cat(sprintf("  p-value (raw): %.4f %s\n",
                 row$p_value,
                 if (!is.na(row$significant_raw) && row$significant_raw) "*" else ""))
      cat(sprintf("  p-value (adj): %.4f %s\n",
                 row$p_adjusted,
                 if (!is.na(row$significant_adjusted) && row$significant_adjusted) "*" else ""))
      cat("\n")
    }
  }
  
  cat("══════════════════════════════════════════════════════════\n")
  cat("* p < 0.05\n\n")
  
  # Summary of significant results
  n_sig_raw <- sum(results$significant_raw, na.rm = TRUE)
  n_sig_adj <- sum(results$significant_adjusted, na.rm = TRUE)
  
  cat("Summary:\n")
  cat(sprintf("  Significant (uncorrected): %d/%d (%.1f%%)\n",
             n_sig_raw, nrow(results), 100*n_sig_raw/nrow(results)))
  cat(sprintf("  Significant (corrected):   %d/%d (%.1f%%)\n",
             n_sig_adj, nrow(results), 100*n_sig_adj/nrow(results)))
  
  if (n_sig_raw > n_sig_adj) {
    cat(sprintf("\n⚠ Note: %d test(s) became non-significant after correction\n",
               n_sig_raw - n_sig_adj))
    cat("  → This is expected when correcting for multiple comparisons\n")
    cat("  → Use adjusted p-values for formal conclusions\n")
  }
  
  cat("\n")
}

#' Save statistical test results to CSV
#'
#' @param results Data frame from test_alpha_diversity() or pairwise_permanova()
#' @param output_file Path to output CSV file
#' @export
save_statistical_results <- function(results, output_file) {
  if (is.null(results) || nrow(results) == 0) {
    warning("No results to save")
    return(invisible(NULL))
  }
  
  tryCatch({
    write.csv(results, output_file, row.names = FALSE)
    message("Saved statistical results to: ", output_file)
  }, error = function(e) {
    warning("Failed to save results: ", e$message)
  })
}
