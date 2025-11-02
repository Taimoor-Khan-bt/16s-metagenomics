# Beta Diversity Analysis Module
# PERMANOVA tests, ordinations, and multifactor models
# Author: Taimoor Khan
# Date: 2025-11-01

library(phyloseq)
library(vegan)
library(ape)

#' Run Primary Beta Diversity Analysis
#' 
#' PERMANOVA test for primary group comparison
#' 
#' @param ps phyloseq object
#' @param cfg Configuration list
#' @param plan Analysis plan
#' @param distance Distance metric (default: bray)
#' @return List with distance matrix, PERMANOVA results, and ordination
#' 
#' @export
run_primary_beta_analysis <- function(ps, cfg, plan, distance = "bray") {
  
  group_col <- plan$primary$group_column
  message(sprintf("[beta] Running primary analysis: %s (%s distance)", 
                 group_col, distance))
  
  # Check if UniFrac requires tree
  if (grepl("unifrac", distance, ignore.case = TRUE)) {
    if (is.null(phy_tree(ps, errorIfNULL = FALSE))) {
      message(sprintf("[beta] Skipping %s: phylogenetic tree not available", distance))
      return(NULL)
    }
  }
  
  # Calculate distance
  dist_matrix <- tryCatch({
    phyloseq::distance(ps, method = distance)
  }, error = function(e) {
    warning(sprintf("Distance calculation failed for %s: %s", distance, e$message))
    return(NULL)
  })
  
  if (is.null(dist_matrix)) return(NULL)
  
  # PERMANOVA
  meta <- as(sample_data(ps), "data.frame")
  formula_str <- sprintf("dist_matrix ~ %s", group_col)
  
  set.seed(cfg$analysis$seed %||% 42)
  permanova <- adonis2(as.formula(formula_str), 
                      data = meta, 
                      permutations = 999)
  
  message(sprintf("[beta] PERMANOVA R² = %.4f, p = %.4f", 
                 permanova$R2[1], permanova$`Pr(>F)`[1]))
  
  # Ordination (PCoA)
  ord <- ordinate(ps, method = "PCoA", distance = distance)
  
  # Save results
  outdir <- cfg$output$directory
  saveRDS(dist_matrix, file.path(outdir, sprintf("beta_%s_distance.rds", distance)))
  write.csv(as.data.frame(permanova), 
           file.path(outdir, sprintf("beta_%s_permanova_primary.csv", distance)))
  
  return(list(
    distance_matrix = dist_matrix,
    permanova = permanova,
    ordination = ord
  ))
}

#' Run Secondary Beta Diversity Analysis
#' 
#' PERMANOVA for each secondary variable
#' 
#' @param ps phyloseq object
#' @param cfg Configuration list
#' @param plan Analysis plan
#' @param distance Distance metric
#' @return List of results for each secondary variable
#' 
#' @export
run_secondary_beta_analysis <- function(ps, cfg, plan, distance = "bray") {
  
  if (length(plan$secondary_vars) == 0 || !plan$run_beta_stratified) {
    message("[beta] No secondary analyses requested")
    return(NULL)
  }
  
  message(sprintf("[beta] Running secondary analyses: %s", 
                 paste(plan$secondary_vars, collapse = ", ")))
  
  # Calculate distance once
  dist_matrix <- phyloseq::distance(ps, method = distance)
  meta <- as(sample_data(ps), "data.frame")
  
  results <- list()
  
  for (sec_var in plan$secondary_vars) {
    message(sprintf("[beta] Testing: %s", sec_var))
    
    # Check for NAs and skip if too many missing values
    var_data <- meta[[sec_var]]
    n_missing <- sum(is.na(var_data))
    n_total <- length(var_data)
    
    if (n_missing > 0) {
      message(sprintf("[beta]   Warning: %d/%d samples have NA for %s", 
                     n_missing, n_total, sec_var))
      
      if (n_missing >= n_total * 0.5) {
        message(sprintf("[beta]   Skipping: too many missing values (>50%%)"))
        next
      }
      
      # Remove NA samples for this test
      valid_samples <- !is.na(var_data)
      ps_subset <- prune_samples(valid_samples, ps)
      dist_subset <- phyloseq::distance(ps_subset, method = distance)
      meta_subset <- as(sample_data(ps_subset), "data.frame")
    } else {
      dist_subset <- dist_matrix
      meta_subset <- meta
    }
    
    # Check if enough groups remain
    n_groups <- length(unique(meta_subset[[sec_var]]))
    if (n_groups < 2) {
      message(sprintf("[beta]   Skipping: only %d group(s) available", n_groups))
      next
    }
    
    formula_str <- sprintf("dist_subset ~ %s", sec_var)
    
    set.seed(cfg$analysis$seed %||% 42)
    permanova <- tryCatch({
      adonis2(as.formula(formula_str), 
              data = meta_subset, 
              permutations = 999)
    }, error = function(e) {
      warning(sprintf("PERMANOVA failed for %s: %s", sec_var, e$message))
      return(NULL)
    })
    
    if (is.null(permanova)) next
    
    message(sprintf("[beta]   R² = %.4f, p = %.4f", 
                   permanova$R2[1], permanova$`Pr(>F)`[1]))
    
    results[[sec_var]] <- list(
      permanova = permanova,
      n_samples = nrow(meta_subset),
      n_missing = n_missing
    )
    
    # Save
    outdir <- cfg$output$directory
    write.csv(as.data.frame(permanova), 
             file.path(outdir, sprintf("beta_%s_permanova_%s.csv", distance, sec_var)))
  }
  
  return(results)
}

#' Run Multifactor PERMANOVA
#' 
#' Test multiple variables simultaneously with interactions
#' Example: Group + Sex + Group:Sex
#' 
#' @param ps phyloseq object
#' @param cfg Configuration list
#' @param plan Analysis plan
#' @param distance Distance metric
#' @return PERMANOVA results with multiple factors
#' 
#' @export
run_multifactor_beta_analysis <- function(ps, cfg, plan, distance = "bray") {
  
  if (!plan$run_beta_multifactor) {
    message("[beta] Multifactor analysis not enabled")
    return(NULL)
  }
  
  group_col <- plan$primary$group_column
  
  # Check if we have secondary variables
  if (length(plan$secondary_vars) == 0) {
    warning("[beta] Multifactor requested but no secondary variables defined")
    return(NULL)
  }
  
  message(sprintf("[beta] Running multifactor PERMANOVA"))
  
  # Calculate distance
  dist_matrix <- phyloseq::distance(ps, method = distance)
  meta <- as(sample_data(ps), "data.frame")
  
  # Build formula with main effects and interactions
  # Example: Group + Sex + AgeCategory + Group:Sex
  formula_terms <- c(group_col, plan$secondary_vars)
  
  # Add interaction terms (primary with each secondary)
  if (length(plan$secondary_vars) > 0) {
    interactions <- paste(group_col, plan$secondary_vars, sep = ":")
    formula_terms <- c(formula_terms, interactions)
  }
  
  formula_str <- sprintf("dist_matrix ~ %s", paste(formula_terms, collapse = " + "))
  
  message(sprintf("[beta] Formula: %s", formula_str))
  
  # Run PERMANOVA
  set.seed(cfg$analysis$seed %||% 42)
  permanova <- tryCatch({
    adonis2(as.formula(formula_str), 
           data = meta, 
           permutations = 999,
           by = "margin")  # Type III (marginal effects)
  }, error = function(e) {
    warning(sprintf("[beta] Multifactor PERMANOVA failed: %s", e$message))
    return(NULL)
  })
  
  if (!is.null(permanova)) {
    # Print results
    message(sprintf("[beta] Multifactor PERMANOVA results:"))
    print(permanova)
    
    # Save
    outdir <- cfg$output$directory
    write.csv(as.data.frame(permanova), 
             file.path(outdir, sprintf("beta_%s_permanova_multifactor.csv", distance)))
  }
  
  return(permanova)
}

#' Run Pairwise PERMANOVA Comparisons
#' 
#' Compare specific group pairs with Bonferroni correction
#' 
#' @param ps phyloseq object
#' @param cfg Configuration list
#' @param plan Analysis plan
#' @param distance Distance metric
#' @return Data frame with pairwise results
#' 
#' @export
run_pairwise_beta_comparisons <- function(ps, cfg, plan, distance = "bray") {
  
  if (!plan$run_pairwise_comparisons || length(plan$pairwise_contrasts) == 0) {
    message("[beta] Pairwise comparisons not requested")
    return(NULL)
  }
  
  group_col <- plan$primary$group_column
  message(sprintf("[beta] Running pairwise PERMANOVA comparisons"))
  
  # Calculate distance
  dist_matrix <- phyloseq::distance(ps, method = distance)
  meta <- as(sample_data(ps), "data.frame")
  
  results <- data.frame()
  
  for (contrast in plan$pairwise_contrasts) {
    group1 <- contrast$group1
    group2 <- contrast$group2
    
    message(sprintf("[beta]   %s vs %s", group1, group2))
    
    # Subset to these two groups
    keep_samples <- meta[[group_col]] %in% c(group1, group2)
    ps_subset <- prune_samples(keep_samples, ps)
    meta_subset <- as(sample_data(ps_subset), "data.frame")
    
    # Subset distance matrix
    dist_subset <- as.dist(as.matrix(dist_matrix)[keep_samples, keep_samples])
    
    # PERMANOVA
    formula_str <- sprintf("dist_subset ~ %s", group_col)
    
    set.seed(cfg$analysis$seed %||% 42)
    permanova <- adonis2(as.formula(formula_str), 
                        data = meta_subset, 
                        permutations = 999)
    
    results <- rbind(results, data.frame(
      comparison = sprintf("%s vs %s", group1, group2),
      group1 = group1,
      group2 = group2,
      n1 = sum(meta_subset[[group_col]] == group1),
      n2 = sum(meta_subset[[group_col]] == group2),
      R2 = permanova$R2[1],
      F_statistic = permanova$F[1],
      p_value = permanova$`Pr(>F)`[1]
    ))
  }
  
  # Bonferroni correction
  results$p_adjusted <- p.adjust(results$p_value, method = "bonferroni")
  results$significance <- ifelse(results$p_adjusted < 0.001, "***",
                                ifelse(results$p_adjusted < 0.01, "**",
                                      ifelse(results$p_adjusted < 0.05, "*", "ns")))
  
  # Save
  outdir <- cfg$output$directory
  write.csv(results, 
           file.path(outdir, sprintf("beta_%s_pairwise_comparisons.csv", distance)),
           row.names = FALSE)
  
  message(sprintf("[beta] Pairwise comparisons complete"))
  print(results)
  
  return(results)
}

#' Plot Ordination by Variable
#' 
#' Create PCoA plot colored by specified variable
#' 
#' @param ps phyloseq object
#' @param ord Ordination object
#' @param color_by Variable name for coloring
#' @param shape_by Optional variable for shape
#' @param title Plot title
#' @return ggplot object
#' 
#' @export
plot_ordination_by_variable <- function(ps, ord, color_by, 
                                       shape_by = NULL, 
                                       title = NULL) {
  
  library(ggplot2)
  
  if (is.null(title)) {
    title <- sprintf("PCoA - %s", color_by)
  }
  
  p <- plot_ordination(ps, ord, color = color_by, shape = shape_by) +
    theme_bw() +
    ggtitle(title) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
  
  return(p)
}

#' Test Dispersion Homogeneity
#' 
#' PERMDISP to check if group dispersions are equal
#' Important assumption for PERMANOVA interpretation
#' 
#' @param ps phyloseq object
#' @param group_col Grouping variable
#' @param distance Distance metric
#' @return PERMDISP results
#' 
#' @export
test_beta_dispersion <- function(ps, group_col, distance = "bray") {
  
  message(sprintf("[beta] Testing dispersion homogeneity: %s", group_col))
  
  # Calculate distance
  dist_matrix <- phyloseq::distance(ps, method = distance)
  
  # Get groups
  meta <- as(sample_data(ps), "data.frame")
  groups <- meta[[group_col]]
  
  # PERMDISP
  disp <- betadisper(dist_matrix, groups)
  disp_test <- permutest(disp, permutations = 999)
  
  message(sprintf("[beta] PERMDISP p = %.4f", disp_test$tab$`Pr(>F)`[1]))
  
  if (disp_test$tab$`Pr(>F)`[1] < 0.05) {
    warning(sprintf("[beta] Significant dispersion differences detected (p=%.4f)", 
                   disp_test$tab$`Pr(>F)`[1]))
    warning("[beta] PERMANOVA results may be confounded by dispersion effects")
  }
  
  return(list(
    betadisper = disp,
    test = disp_test
  ))
}
