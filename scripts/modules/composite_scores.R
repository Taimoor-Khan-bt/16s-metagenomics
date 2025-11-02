# Composite Scores Module
# Creates derived metadata variables from component indicators
# Author: Taimoor Khan
# Date: 2025-11-01

#' Create Composite Scores from Components
#' 
#' Generates new metadata columns by combining binary/categorical indicators
#' 
#' @param ps phyloseq object with sample_data
#' @param cfg Configuration list with composite_scores definition
#' @return Modified phyloseq object with new metadata columns
#' 
#' @export
create_composite_scores <- function(ps, cfg) {
  composite_defs <- cfg$metadata$composite_scores
  
  if (is.null(composite_defs) || length(composite_defs) == 0) {
    message("[composite] No composite scores defined")
    return(ps)
  }
  
  meta <- as(sample_data(ps), "data.frame")
  meta_cols <- colnames(meta)
  
  for (score_def in composite_defs) {
    # Skip if disabled
    if (!isTRUE(score_def$enabled)) {
      next
    }
    
    score_name <- score_def$name
    components <- score_def$components
    method <- score_def$method %||% "count_positive"
    positive_values <- score_def$positive_values %||% c("yes", "1", "positive", "abnormal", "severe")
    
    message(sprintf("[composite] Creating score '%s' from: %s", 
                   score_name, paste(components, collapse = ", ")))
    
    # Validate components exist
    missing_comps <- setdiff(components, meta_cols)
    if (length(missing_comps) > 0) {
      warning(sprintf("Skipping '%s': missing components %s", 
                     score_name, paste(missing_comps, collapse = ", ")))
      next
    }
    
    # Calculate score based on method
    if (method == "count_positive") {
      score_values <- calculate_count_positive(meta, components, positive_values)
    } else if (method == "sum") {
      score_values <- calculate_sum(meta, components)
    } else if (method == "mean") {
      score_values <- calculate_mean(meta, components)
    } else if (method == "any") {
      score_values <- calculate_any(meta, components, positive_values)
    } else if (method == "all") {
      score_values <- calculate_all(meta, components, positive_values)
    } else {
      warning(sprintf("Unknown method '%s' for score '%s'. Using count_positive.", 
                     method, score_name))
      score_values <- calculate_count_positive(meta, components, positive_values)
    }
    
    # Add to metadata
    meta[[score_name]] <- score_values
    
    # Optionally categorize
    if (!is.null(score_def$categories)) {
      cat_name <- paste0(score_name, "_Category")
      meta[[cat_name]] <- categorize_score(score_values, score_def$categories)
      message(sprintf("[composite] Created categorical version: %s", cat_name))
    }
    
    # Summary
    if (is.numeric(score_values)) {
      message(sprintf("[composite] %s: min=%.2f, max=%.2f, mean=%.2f", 
                     score_name, min(score_values, na.rm = TRUE), 
                     max(score_values, na.rm = TRUE), 
                     mean(score_values, na.rm = TRUE)))
    } else {
      tbl <- table(score_values, useNA = "ifany")
      message(sprintf("[composite] %s distribution: %s", 
                     score_name, paste(names(tbl), "=", tbl, collapse = ", ")))
    }
  }
  
  # Update phyloseq object
  sample_data(ps) <- meta
  
  return(ps)
}

#' Count Positive Indicators
#' 
#' @param meta Data frame
#' @param components Column names
#' @param positive_values Values considered "positive"
#' @return Numeric vector with counts
calculate_count_positive <- function(meta, components, positive_values) {
  # Normalize positive values for comparison
  pos_vals_lower <- tolower(as.character(positive_values))
  
  counts <- apply(meta[, components, drop = FALSE], 1, function(row) {
    # Convert to lowercase strings for comparison
    row_lower <- tolower(as.character(row))
    
    # Count matches
    sum(row_lower %in% pos_vals_lower, na.rm = TRUE)
  })
  
  return(counts)
}

#' Sum Numeric Components
#' 
#' @param meta Data frame
#' @param components Column names (must be numeric)
#' @return Numeric vector with sums
calculate_sum <- function(meta, components) {
  # Coerce to numeric
  comp_data <- as.data.frame(lapply(meta[, components, drop = FALSE], function(x) {
    as.numeric(as.character(x))
  }))
  
  sums <- rowSums(comp_data, na.rm = TRUE)
  return(sums)
}

#' Mean of Numeric Components
#' 
#' @param meta Data frame
#' @param components Column names (must be numeric)
#' @return Numeric vector with means
calculate_mean <- function(meta, components) {
  comp_data <- as.data.frame(lapply(meta[, components, drop = FALSE], function(x) {
    as.numeric(as.character(x))
  }))
  
  means <- rowMeans(comp_data, na.rm = TRUE)
  return(means)
}

#' Check if ANY Component is Positive
#' 
#' @param meta Data frame
#' @param components Column names
#' @param positive_values Values considered "positive"
#' @return Logical vector
calculate_any <- function(meta, components, positive_values) {
  pos_vals_lower <- tolower(as.character(positive_values))
  
  any_positive <- apply(meta[, components, drop = FALSE], 1, function(row) {
    row_lower <- tolower(as.character(row))
    any(row_lower %in% pos_vals_lower, na.rm = TRUE)
  })
  
  return(ifelse(any_positive, "Yes", "No"))
}

#' Check if ALL Components are Positive
#' 
#' @param meta Data frame
#' @param components Column names
#' @param positive_values Values considered "positive"
#' @return Character vector ("Yes" or "No")
calculate_all <- function(meta, components, positive_values) {
  pos_vals_lower <- tolower(as.character(positive_values))
  
  all_positive <- apply(meta[, components, drop = FALSE], 1, function(row) {
    row_lower <- tolower(as.character(row))
    # Must have all components as positive (no NAs allowed)
    all(!is.na(row)) && all(row_lower %in% pos_vals_lower)
  })
  
  return(ifelse(all_positive, "Yes", "No"))
}

#' Categorize Continuous Score
#' 
#' @param scores Numeric vector
#' @param categories List with cutoffs and labels
#' @return Factor with categories
#' 
#' @examples
#' categories <- list(
#'   cutoffs = c(0, 1, 2, 4),
#'   labels = c("None", "Low", "Moderate", "High")
#' )
categorize_score <- function(scores, categories) {
  cutoffs <- categories$cutoffs
  labels <- categories$labels
  
  if (length(cutoffs) != length(labels) + 1) {
    stop("categories must have length(labels) + 1 cutoffs")
  }
  
  cats <- cut(scores, 
             breaks = cutoffs, 
             labels = labels,
             include.lowest = TRUE,
             right = FALSE)
  
  return(as.character(cats))
}

#' Create Nutritional Risk Score
#' 
#' Convenience function for common nutritional risk composite
#' 
#' @param ps phyloseq object
#' @param stunting_col Column name for stunting status
#' @param underweight_col Column name for underweight status
#' @param wasting_col Column name for wasting status
#' @param anemia_col Column name for anemia status
#' @return Modified phyloseq with NutritionalRisk column
#' 
#' @export
add_nutritional_risk_score <- function(ps, 
                                       stunting_col = "Stunting",
                                       underweight_col = "Underweight", 
                                       wasting_col = "Wasting",
                                       anemia_col = "Anemia") {
  
  meta <- as(sample_data(ps), "data.frame")
  
  components <- c(stunting_col, underweight_col, wasting_col, anemia_col)
  missing <- setdiff(components, colnames(meta))
  
  if (length(missing) > 0) {
    stop(sprintf("Missing columns: %s", paste(missing, collapse = ", ")))
  }
  
  # Count number of risk factors (assuming "yes" or 1 = positive)
  positive_vals <- c("yes", "1", "positive", "abnormal")
  
  risk_score <- calculate_count_positive(meta, components, positive_vals)
  
  # Categorize: 0 = None, 1 = Low, 2 = Moderate, 3-4 = High
  risk_category <- cut(risk_score, 
                      breaks = c(-1, 0, 1, 2, 4),
                      labels = c("None", "Low", "Moderate", "High"),
                      include.lowest = TRUE)
  
  meta$NutritionalRisk <- risk_score
  meta$NutritionalRisk_Category <- as.character(risk_category)
  
  sample_data(ps) <- meta
  
  message(sprintf("[composite] Created NutritionalRisk score (0-4) and category"))
  message(sprintf("[composite] Distribution: %s", 
                 paste(table(meta$NutritionalRisk_Category), collapse = ", ")))
  
  return(ps)
}
