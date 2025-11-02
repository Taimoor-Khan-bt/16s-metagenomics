# Alpha Diversity Analysis Module
# Flexible group comparisons, stratified analyses, and continuous correlations
# Author: Taimoor Khan
# Date: 2025-11-01

library(phyloseq)
library(vegan)

#' Run Primary Alpha Diversity Analysis
#' 
#' Compare alpha diversity between groups defined in primary comparison
#' 
#' @param ps phyloseq object
#' @param cfg Configuration list
#' @param plan Analysis plan from get_analysis_plan()
#' @param metrics Vector of diversity metrics (default: Shannon, Observed, Simpson)
#' @return Data frame with diversity values and statistics
#' 
#' @export
run_primary_alpha_analysis <- function(ps, cfg, plan, 
                                       metrics = c("Shannon", "Observed", "Simpson")) {
  
  group_col <- plan$primary$group_column
  message(sprintf("[alpha] Running primary analysis: %s", group_col))
  
  # Calculate diversity
  div_df <- calculate_alpha_diversity(ps, metrics)
  
  # Add grouping variable
  meta <- as(sample_data(ps), "data.frame")
  div_df[[group_col]] <- meta[[group_col]][match(rownames(div_df), rownames(meta))]
  
  # Statistical tests
  results <- list(
    diversity = div_df,
    statistics = test_alpha_by_group(div_df, group_col, metrics)
  )
  
  # Save results
  outdir <- cfg$output$directory
  write.csv(div_df, 
           file.path(outdir, "alpha_diversity_primary.csv"),
           row.names = TRUE)
  
  message(sprintf("[alpha] Results saved: alpha_diversity_primary.csv"))
  
  return(results)
}

#' Run Secondary Alpha Diversity Analysis
#' 
#' Compare alpha diversity by secondary variables (e.g., Sex, Age category)
#' 
#' @param ps phyloseq object
#' @param cfg Configuration list
#' @param plan Analysis plan
#' @param metrics Diversity metrics
#' @return List of results for each secondary variable
#' 
#' @export
run_secondary_alpha_analysis <- function(ps, cfg, plan, 
                                         metrics = c("Shannon", "Observed", "Simpson")) {
  
  if (length(plan$secondary_vars) == 0) {
    message("[alpha] No secondary variables defined")
    return(NULL)
  }
  
  message(sprintf("[alpha] Running secondary analyses: %s", 
                 paste(plan$secondary_vars, collapse = ", ")))
  
  # Calculate diversity once
  div_df <- calculate_alpha_diversity(ps, metrics)
  meta <- as(sample_data(ps), "data.frame")
  
  results <- list()
  
  for (sec_var in plan$secondary_vars) {
    message(sprintf("[alpha] Testing: %s", sec_var))
    
    # Add secondary variable
    div_sec <- div_df
    div_sec[[sec_var]] <- meta[[sec_var]][match(rownames(div_df), rownames(meta))]
    
    # Statistical tests
    stats <- test_alpha_by_group(div_sec, sec_var, metrics)
    
    results[[sec_var]] <- list(
      diversity = div_sec,
      statistics = stats
    )
    
    # Save
    outdir <- cfg$output$directory
    write.csv(div_sec, 
             file.path(outdir, sprintf("alpha_diversity_%s.csv", sec_var)),
             row.names = TRUE)
  }
  
  message(sprintf("[alpha] Secondary analyses complete"))
  
  return(results)
}

#' Run Stratified Alpha Diversity Analysis
#' 
#' Compare primary groups WITHIN each level of secondary variable
#' Example: Compare Caries groups separately for Male and Female
#' 
#' @param ps phyloseq object
#' @param cfg Configuration list
#' @param plan Analysis plan
#' @param metrics Diversity metrics
#' @return Nested list of results
#' 
#' @export
run_stratified_alpha_analysis <- function(ps, cfg, plan,
                                          metrics = c("Shannon", "Observed", "Simpson")) {
  
  if (length(plan$secondary_vars) == 0 || !plan$run_alpha_stratified) {
    message("[alpha] Stratified analysis not enabled")
    return(NULL)
  }
  
  group_col <- plan$primary$group_column
  message(sprintf("[alpha] Running stratified analyses: %s by %s", 
                 group_col, paste(plan$secondary_vars, collapse = ", ")))
  
  div_df <- calculate_alpha_diversity(ps, metrics)
  meta <- as(sample_data(ps), "data.frame")
  
  results <- list()
  
  for (strata_var in plan$secondary_vars) {
    message(sprintf("[alpha] Stratifying by: %s", strata_var))
    
    div_strat <- div_df
    div_strat[[group_col]] <- meta[[group_col]][match(rownames(div_df), rownames(meta))]
    div_strat[[strata_var]] <- meta[[strata_var]][match(rownames(div_df), rownames(meta))]
    
    # Get unique strata levels
    strata_levels <- unique(div_strat[[strata_var]])
    strata_levels <- strata_levels[!is.na(strata_levels)]
    
    strata_results <- list()
    
    for (level in strata_levels) {
      message(sprintf("[alpha]   %s = %s", strata_var, level))
      
      # Subset to this stratum
      div_subset <- div_strat[div_strat[[strata_var]] == level, ]
      
      # Check if stratum has enough samples
      n_in_stratum <- nrow(div_subset)
      n_groups_in_stratum <- length(unique(div_subset[[group_col]]))
      
      if (n_in_stratum < 3 || n_groups_in_stratum < 2) {
        message(sprintf("[alpha]   Skipping: insufficient samples (n=%d, groups=%d)", 
                       n_in_stratum, n_groups_in_stratum))
        next
      }
      
      # Test group differences within stratum
      stats <- test_alpha_by_group(div_subset, group_col, metrics)
      
      strata_results[[as.character(level)]] <- list(
        n = nrow(div_subset),
        diversity = div_subset,
        statistics = stats
      )
    }
    
    results[[strata_var]] <- strata_results
    
    # Save combined results
    outdir <- cfg$output$directory
    write.csv(div_strat, 
             file.path(outdir, sprintf("alpha_diversity_stratified_%s.csv", strata_var)),
             row.names = TRUE)
  }
  
  return(results)
}

#' Correlate Alpha Diversity with Continuous Variables
#' 
#' @param ps phyloseq object
#' @param cfg Configuration list
#' @param plan Analysis plan
#' @param metrics Diversity metrics
#' @return Data frame with correlation results
#' 
#' @export
run_alpha_continuous_correlations <- function(ps, cfg, plan,
                                              metrics = c("Shannon", "Observed", "Simpson")) {
  
  if (length(plan$continuous_vars) == 0 || !plan$run_alpha_continuous) {
    message("[alpha] Continuous correlations not enabled")
    return(NULL)
  }
  
  message(sprintf("[alpha] Correlating with continuous variables: %s", 
                 paste(plan$continuous_vars, collapse = ", ")))
  
  div_df <- calculate_alpha_diversity(ps, metrics)
  meta <- as(sample_data(ps), "data.frame")
  
  # Add continuous variables
  for (cont_var in plan$continuous_vars) {
    div_df[[cont_var]] <- as.numeric(meta[[cont_var]][match(rownames(div_df), rownames(meta))])
  }
  
  # Calculate correlations
  cor_results <- data.frame()
  
  for (metric in metrics) {
    for (cont_var in plan$continuous_vars) {
      # Remove NAs
      valid_idx <- !is.na(div_df[[metric]]) & !is.na(div_df[[cont_var]])
      
      if (sum(valid_idx) < 3) {
        warning(sprintf("Insufficient data for %s vs %s", metric, cont_var))
        next
      }
      
      # Spearman correlation (non-parametric)
      cor_test <- cor.test(div_df[[metric]][valid_idx], 
                          div_df[[cont_var]][valid_idx],
                          method = "spearman")
      
      cor_results <- rbind(cor_results, data.frame(
        metric = metric,
        variable = cont_var,
        n = sum(valid_idx),
        rho = cor_test$estimate,
        p_value = cor_test$p.value,
        significance = ifelse(cor_test$p.value < 0.001, "***",
                             ifelse(cor_test$p.value < 0.01, "**",
                                   ifelse(cor_test$p.value < 0.05, "*", "ns")))
      ))
      
      message(sprintf("[alpha]   %s ~ %s: rho=%.3f, p=%.4f", 
                     metric, cont_var, cor_test$estimate, cor_test$p.value))
    }
  }
  
  # Save
  outdir <- cfg$output$directory
  write.csv(cor_results, 
           file.path(outdir, "alpha_diversity_continuous_correlations.csv"),
           row.names = FALSE)
  
  return(list(
    diversity_with_continuous = div_df,
    correlations = cor_results
  ))
}

#' Calculate Alpha Diversity Metrics
#' 
#' @param ps phyloseq object
#' @param metrics Vector of metric names
#' @return Data frame with diversity values
calculate_alpha_diversity <- function(ps, metrics = c("Shannon", "Observed", "Simpson")) {
  div_df <- data.frame(row.names = sample_names(ps))
  
  # Get OTU table with correct orientation (samples as rows)
  otu_mat <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) {
    # Already samples as rows, taxa as columns - this is DADA2 default
    otu_mat <- otu_mat
  } else {
    # Transpose if taxa are rows
    otu_mat <- t(otu_mat)
  }
  
  for (metric in metrics) {
    if (metric == "Shannon") {
      div_df$Shannon <- vegan::diversity(otu_mat, index = "shannon", MARGIN = 1)
    } else if (metric == "Observed") {
      div_df$Observed <- rowSums(otu_mat > 0)
    } else if (metric == "Simpson") {
      div_df$Simpson <- vegan::diversity(otu_mat, index = "simpson", MARGIN = 1)
    } else if (metric == "Chao1") {
      # Requires incidence data - calculate per sample (row)
      chao <- apply(otu_mat, 1, function(x) {
        S_obs <- sum(x > 0)
        f1 <- sum(x == 1)
        f2 <- sum(x == 2)
        if (f2 == 0) return(S_obs)
        S_obs + (f1 * (f1 - 1)) / (2 * (f2 + 1))
      })
      div_df$Chao1 <- chao
    } else {
      warning(sprintf("Unknown metric: %s", metric))
    }
  }
  
  return(div_df)
}

#' Test Alpha Diversity by Group
#' 
#' Automatically selects appropriate test based on number of groups
#' 
#' @param div_df Data frame with diversity and grouping variable
#' @param group_col Name of grouping column
#' @param metrics Diversity metrics to test
#' @return Data frame with test statistics
test_alpha_by_group <- function(div_df, group_col, metrics) {
  
  # Remove NA groups
  div_clean <- div_df[!is.na(div_df[[group_col]]), ]
  
  n_groups <- length(unique(div_clean[[group_col]]))
  
  if (n_groups < 2) {
    warning(sprintf("Only %d group found for %s", n_groups, group_col))
    return(NULL)
  }
  
  results <- data.frame()
  
  for (metric in metrics) {
    # Remove NAs in metric
    valid_idx <- !is.na(div_clean[[metric]])
    test_data <- div_clean[valid_idx, ]
    
    if (nrow(test_data) < 2) {
      warning(sprintf("Insufficient data for %s", metric))
      next
    }
    
    # Select test
    if (n_groups == 2) {
      # Wilcoxon rank-sum
      groups <- split(test_data[[metric]], test_data[[group_col]])
      
      # Check if both groups have data
      if (length(groups[[1]]) < 1 || length(groups[[2]]) < 1) {
        warning(sprintf("Insufficient samples in groups for %s", metric))
        next
      }
      
      test_res <- tryCatch({
        wilcox.test(groups[[1]], groups[[2]])
      }, error = function(e) {
        warning(sprintf("Wilcoxon test failed for %s: %s", metric, e$message))
        return(NULL)
      })
      
      if (is.null(test_res)) next
      
      method <- "Wilcoxon"
    } else {
      # Kruskal-Wallis
      formula_str <- sprintf("%s ~ %s", metric, group_col)
      test_res <- tryCatch({
        kruskal.test(as.formula(formula_str), data = test_data)
      }, error = function(e) {
        warning(sprintf("Kruskal-Wallis test failed for %s: %s", metric, e$message))
        return(NULL)
      })
      
      if (is.null(test_res)) next
      
      method <- "Kruskal-Wallis"
    }
    
    results <- rbind(results, data.frame(
      metric = metric,
      grouping = group_col,
      n_groups = n_groups,
      n_samples = nrow(test_data),
      method = method,
      statistic = test_res$statistic,
      p_value = test_res$p.value,
      significance = ifelse(test_res$p.value < 0.001, "***",
                           ifelse(test_res$p.value < 0.01, "**",
                                 ifelse(test_res$p.value < 0.05, "*", "ns")))
    ))
  }
  
  return(results)
}
