# Analysis Plan Module
# Parses configuration and creates execution plan for flexible multi-group analysis
# Author: Taimoor Khan
# Date: 2025-11-01

#' Get Analysis Plan from Configuration
#' 
#' Determines what analyses to run based on config mode and metadata structure
#' 
#' @param cfg Configuration list from YAML
#' @return List with analysis plan details
#' 
#' @export
get_analysis_plan <- function(cfg) {
  mode <- tolower(cfg$analysis$mode %||% "standard")
  
  # Validate mode
  valid_modes <- c("simple", "standard", "comprehensive")
  if (!mode %in% valid_modes) {
    warning(sprintf("Invalid mode '%s'. Using 'standard'. Valid modes: %s", 
                   mode, paste(valid_modes, collapse = ", ")))
    mode <- "standard"
  }
  
  message(sprintf("[plan] Analysis mode: %s", toupper(mode)))
  
  # Initialize plan
  plan <- list(
    mode = mode,
    
    # Primary comparison (always required)
    primary = list(
      group_column = cfg$metadata$primary_comparison$group_column %||% 
                     cfg$metadata$group_column %||% "Group",
      levels = cfg$metadata$primary_comparison$levels %||% NULL,
      reference = cfg$metadata$primary_comparison$reference %||% NULL
    ),
    
    # Analysis flags
    run_alpha_primary = TRUE,
    run_beta_primary = TRUE,
    run_composition = TRUE,
    
    # Secondary analyses (standard and comprehensive)
    run_alpha_stratified = FALSE,
    run_beta_stratified = FALSE,
    run_alpha_continuous = FALSE,
    run_beta_multifactor = FALSE,
    run_differential_abundance = FALSE,
    run_pairwise_comparisons = FALSE,
    
    # Variables
    secondary_vars = c(),
    continuous_vars = c(),
    covariates = cfg$metadata$covariates %||% c(),
    
    # Composite scores
    composite_scores = list(),
    
    # Methods
    alpha_metrics = cfg$analysis$alpha$metrics %||% c("Shannon", "Observed", "Simpson"),
    beta_metrics = cfg$analysis$beta$metrics %||% c("bray"),
    da_methods = cfg$analysis$differential_abundance$methods %||% c("deseq2"),
    
    # Options
    pairwise_contrasts = list()
  )
  
  # Parse secondary comparisons (standard and comprehensive modes)
  if (mode %in% c("standard", "comprehensive")) {
    sec_comps <- cfg$metadata$secondary_comparisons
    if (!is.null(sec_comps) && length(sec_comps) > 0) {
      enabled_vars <- sapply(sec_comps, function(x) {
        # Handle both list and character inputs
        if (is.list(x)) {
          if (isTRUE(x$enabled)) x$variable else NULL
        } else if (is.character(x)) {
          x  # Character strings are treated as enabled by default
        } else {
          NULL
        }
      })
      plan$secondary_vars <- unlist(enabled_vars[!sapply(enabled_vars, is.null)])
      
      if (length(plan$secondary_vars) > 0) {
        message(sprintf("[plan] Secondary comparisons: %s", 
                       paste(plan$secondary_vars, collapse = ", ")))
      }
    }
  }
  
  # Comprehensive mode enables advanced analyses
  if (mode == "comprehensive") {
    # Alpha diversity enhancements
    if (isTRUE(cfg$analysis$alpha$stratify_by_secondary)) {
      plan$run_alpha_stratified <- TRUE
    }
    if (isTRUE(cfg$analysis$alpha$correlate_continuous)) {
      plan$run_alpha_continuous <- TRUE
    }
    
    # Beta diversity enhancements
    if (isTRUE(cfg$analysis$beta$test_multifactor)) {
      plan$run_beta_multifactor <- TRUE
    }
    if (isTRUE(cfg$analysis$beta$plot_secondary)) {
      plan$run_beta_stratified <- TRUE
    }
    
    # Differential abundance
    if (isTRUE(cfg$analysis$differential_abundance$enabled)) {
      plan$run_differential_abundance <- TRUE
    }
    
    # Pairwise comparisons
    if (isTRUE(cfg$analysis$pairwise_comparisons$enabled)) {
      plan$run_pairwise_comparisons <- TRUE
      contrasts_raw <- cfg$analysis$pairwise_comparisons$contrasts %||% list()
      
      # Parse contrasts - handle both "Group1-Group2" strings and list formats
      plan$pairwise_contrasts <- lapply(contrasts_raw, function(x) {
        if (is.character(x)) {
          # Parse "Group1-Group2" format
          parts <- strsplit(x, "-")[[1]]
          if (length(parts) == 2) {
            list(group1 = trimws(parts[1]), group2 = trimws(parts[2]))
          } else {
            warning(sprintf("Invalid contrast format: %s", x))
            NULL
          }
        } else if (is.list(x)) {
          # Already in list format
          x
        } else {
          NULL
        }
      })
      
      # Remove NULLs
      plan$pairwise_contrasts <- Filter(Negate(is.null), plan$pairwise_contrasts)
      
      if (length(plan$pairwise_contrasts) > 0) {
        message(sprintf("[plan] Pairwise contrasts: %d comparisons", length(plan$pairwise_contrasts)))
      }
    }
    
    # Parse continuous variables
    cont_vars <- cfg$metadata$continuous_variables
    if (!is.null(cont_vars) && length(cont_vars) > 0) {
      enabled_cont <- sapply(cont_vars, function(x) {
        # Handle both list and character inputs
        if (is.list(x)) {
          if (isTRUE(x$enabled)) x$name else NULL
        } else if (is.character(x)) {
          x  # Character strings are treated as enabled by default
        } else {
          NULL
        }
      })
      plan$continuous_vars <- unlist(enabled_cont[!sapply(enabled_cont, is.null)])
      
      if (length(plan$continuous_vars) > 0) {
        message(sprintf("[plan] Continuous variables: %s", 
                       paste(plan$continuous_vars, collapse = ", ")))
      }
    }
    
    # Parse composite scores
    plan$composite_scores <- cfg$metadata$composite_scores %||% list()
  }
  
  # Summary
  message(sprintf("[plan] Analyses enabled:"))
  message(sprintf("  - Alpha diversity (primary): %s", plan$run_alpha_primary))
  message(sprintf("  - Beta diversity (primary): %s", plan$run_beta_primary))
  message(sprintf("  - Alpha stratified: %s", plan$run_alpha_stratified))
  message(sprintf("  - Alpha continuous correlations: %s", plan$run_alpha_continuous))
  message(sprintf("  - Beta multifactor: %s", plan$run_beta_multifactor))
  message(sprintf("  - Differential abundance: %s", plan$run_differential_abundance))
  
  return(plan)
}

#' Validate Analysis Plan Against Metadata
#' 
#' Checks that required columns exist in phyloseq object
#' 
#' @param plan Analysis plan from get_analysis_plan()
#' @param ps phyloseq object with sample_data
#' @return Modified plan with validated variables only
#' 
#' @export
validate_analysis_plan <- function(plan, ps) {
  meta <- as(sample_data(ps), "data.frame")
  meta_cols <- colnames(meta)
  
  # Check primary group column
  if (!plan$primary$group_column %in% meta_cols) {
    stop(sprintf("Primary group column '%s' not found in metadata. Available: %s",
                plan$primary$group_column, paste(meta_cols, collapse = ", ")))
  }
  
  # Validate secondary variables
  if (length(plan$secondary_vars) > 0) {
    valid_secondary <- plan$secondary_vars[plan$secondary_vars %in% meta_cols]
    invalid_secondary <- setdiff(plan$secondary_vars, meta_cols)
    
    if (length(invalid_secondary) > 0) {
      warning(sprintf("Removing invalid secondary variables: %s", 
                     paste(invalid_secondary, collapse = ", ")))
    }
    
    plan$secondary_vars <- valid_secondary
  }
  
  # Validate continuous variables
  if (length(plan$continuous_vars) > 0) {
    valid_continuous <- plan$continuous_vars[plan$continuous_vars %in% meta_cols]
    invalid_continuous <- setdiff(plan$continuous_vars, meta_cols)
    
    if (length(invalid_continuous) > 0) {
      warning(sprintf("Removing invalid continuous variables: %s", 
                     paste(invalid_continuous, collapse = ", ")))
    }
    
    plan$continuous_vars <- valid_continuous
  }
  
  # Validate covariates
  if (length(plan$covariates) > 0) {
    valid_covariates <- plan$covariates[plan$covariates %in% meta_cols]
    invalid_covariates <- setdiff(plan$covariates, meta_cols)
    
    if (length(invalid_covariates) > 0) {
      warning(sprintf("Removing invalid covariates: %s", 
                     paste(invalid_covariates, collapse = ", ")))
    }
    
    plan$covariates <- valid_covariates
  }
  
  message("[plan] Validation complete")
  return(plan)
}

#' Print Analysis Plan Summary
#' 
#' @param plan Analysis plan object
#' @export
print_analysis_plan <- function(plan) {
  cat("\n")
  cat("=" , rep("=", 60), "=\n", sep = "")
  cat("  ANALYSIS PLAN SUMMARY\n")
  cat("=" , rep("=", 60), "=\n", sep = "")
  cat(sprintf("Mode: %s\n", toupper(plan$mode)))
  cat(sprintf("Primary comparison: %s\n", plan$primary$group_column))
  
  if (length(plan$secondary_vars) > 0) {
    cat(sprintf("Secondary comparisons: %s\n", paste(plan$secondary_vars, collapse = ", ")))
  }
  
  if (length(plan$continuous_vars) > 0) {
    cat(sprintf("Continuous variables: %s\n", paste(plan$continuous_vars, collapse = ", ")))
  }
  
  if (length(plan$covariates) > 0) {
    cat(sprintf("Covariates: %s\n", paste(plan$covariates, collapse = ", ")))
  }
  
  cat("\nAnalyses to run:\n")
  analyses <- c(
    "Alpha diversity (primary)" = plan$run_alpha_primary,
    "Alpha diversity (stratified)" = plan$run_alpha_stratified,
    "Alpha-continuous correlations" = plan$run_alpha_continuous,
    "Beta diversity (primary)" = plan$run_beta_primary,
    "Beta diversity (multifactor)" = plan$run_beta_multifactor,
    "Differential abundance" = plan$run_differential_abundance,
    "Pairwise comparisons" = plan$run_pairwise_comparisons
  )
  
  for (i in seq_along(analyses)) {
    status <- ifelse(analyses[i], "[x]", "[ ]")
    cat(sprintf("  %s %s\n", status, names(analyses)[i]))
  }
  
  cat("=" , rep("=", 60), "=\n", sep = "")
  cat("\n")
}
