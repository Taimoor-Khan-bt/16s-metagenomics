#!/usr/bin/env Rscript
# Input Validation Module
# Provides comprehensive input validation for 16S pipeline
# Author: Pipeline Development Team
# Last Updated: 2024-11-02

# Load required packages
suppressPackageStartupMessages({
  library(dada2)
  library(yaml)
})

#' Validate Sample IDs
#' 
#' Check that sample IDs follow expected format and conventions
#' 
#' @param sample_ids Character vector of sample IDs
#' @param allow_patterns Character vector of allowed regex patterns (optional)
#' @return List with 'valid' (logical) and 'messages' (character vector)
#' @export
validate_sample_ids <- function(sample_ids, allow_patterns = NULL) {
  messages <- character()
  valid <- TRUE
  
  # Check for empty input
  if (length(sample_ids) == 0) {
    messages <- c(messages, "ERROR: No sample IDs provided")
    return(list(valid = FALSE, messages = messages))
  }
  
  # Check for NAs
  if (any(is.na(sample_ids))) {
    n_na <- sum(is.na(sample_ids))
    messages <- c(messages, sprintf("ERROR: %d sample ID(s) are NA", n_na))
    valid <- FALSE
  }
  
  # Check for empty strings
  if (any(sample_ids == "")) {
    n_empty <- sum(sample_ids == "")
    messages <- c(messages, sprintf("ERROR: %d sample ID(s) are empty strings", n_empty))
    valid <- FALSE
  }
  
  # Check for duplicates
  if (any(duplicated(sample_ids))) {
    dupes <- sample_ids[duplicated(sample_ids)]
    messages <- c(messages, sprintf("ERROR: %d duplicate sample ID(s): %s", 
                                   length(dupes), 
                                   paste(head(dupes, 3), collapse = ", ")))
    valid <- FALSE
  }
  
  # Check for whitespace
  has_whitespace <- grepl("\\s", sample_ids)
  if (any(has_whitespace)) {
    n_whitespace <- sum(has_whitespace)
    messages <- c(messages, sprintf("WARNING: %d sample ID(s) contain whitespace", n_whitespace))
  }
  
  # Check for special characters that could cause issues
  has_special <- grepl("[^A-Za-z0-9._-]", sample_ids)
  if (any(has_special)) {
    bad_ids <- sample_ids[has_special]
    messages <- c(messages, sprintf("WARNING: %d sample ID(s) contain special characters: %s",
                                   length(bad_ids),
                                   paste(head(bad_ids, 3), collapse = ", ")))
  }
  
  # Check against allowed patterns if provided
  if (!is.null(allow_patterns)) {
    matches_pattern <- sapply(sample_ids, function(id) {
      any(sapply(allow_patterns, function(pat) grepl(pat, id)))
    })
    
    if (!all(matches_pattern)) {
      bad_ids <- sample_ids[!matches_pattern]
      messages <- c(messages, sprintf("WARNING: %d sample ID(s) don't match allowed patterns: %s",
                                     length(bad_ids),
                                     paste(head(bad_ids, 3), collapse = ", ")))
    }
  }
  
  return(list(valid = valid, messages = messages))
}

#' Validate Metadata Structure
#' 
#' Check that metadata data frame has required columns and proper structure
#' 
#' @param metadata Data frame containing sample metadata
#' @param required_cols Character vector of required column names
#' @param sample_id_col Name of sample ID column (default: "SampleID")
#' @return List with 'valid' (logical) and 'messages' (character vector)
#' @export
validate_metadata_structure <- function(metadata, required_cols = NULL, sample_id_col = "SampleID") {
  messages <- character()
  valid <- TRUE
  
  # Check that input is a data frame
  if (!is.data.frame(metadata)) {
    messages <- c(messages, sprintf("ERROR: Metadata must be a data frame, got: %s", class(metadata)[1]))
    return(list(valid = FALSE, messages = messages))
  }
  
  # Check for empty metadata
  if (nrow(metadata) == 0) {
    messages <- c(messages, "ERROR: Metadata has 0 rows")
    return(list(valid = FALSE, messages = messages))
  }
  
  if (ncol(metadata) == 0) {
    messages <- c(messages, "ERROR: Metadata has 0 columns")
    return(list(valid = FALSE, messages = messages))
  }
  
  # Check for sample ID column
  if (!sample_id_col %in% colnames(metadata)) {
    messages <- c(messages, sprintf("ERROR: Sample ID column '%s' not found in metadata", sample_id_col))
    messages <- c(messages, sprintf("       Available columns: %s", paste(colnames(metadata), collapse = ", ")))
    valid <- FALSE
  } else {
    # Validate sample IDs in metadata
    id_validation <- validate_sample_ids(metadata[[sample_id_col]])
    if (!id_validation$valid) {
      messages <- c(messages, "ERROR: Sample ID validation failed:")
      messages <- c(messages, paste("      ", id_validation$messages))
      valid <- FALSE
    } else if (length(id_validation$messages) > 0) {
      messages <- c(messages, id_validation$messages)
    }
  }
  
  # Check for required columns
  if (!is.null(required_cols)) {
    missing_cols <- setdiff(required_cols, colnames(metadata))
    if (length(missing_cols) > 0) {
      messages <- c(messages, sprintf("ERROR: Missing required column(s): %s", 
                                     paste(missing_cols, collapse = ", ")))
      valid <- FALSE
    }
  }
  
  # Check for columns with all NA values
  all_na_cols <- sapply(metadata, function(x) all(is.na(x)))
  if (any(all_na_cols)) {
    na_cols <- names(metadata)[all_na_cols]
    messages <- c(messages, sprintf("WARNING: %d column(s) have all NA values: %s",
                                   length(na_cols),
                                   paste(na_cols, collapse = ", ")))
  }
  
  # Check for factor columns with only one level
  factor_cols <- sapply(metadata, is.factor)
  if (any(factor_cols)) {
    single_level <- sapply(metadata[, factor_cols, drop = FALSE], function(x) nlevels(x) == 1)
    if (any(single_level)) {
      single_cols <- names(metadata[, factor_cols, drop = FALSE])[single_level]
      messages <- c(messages, sprintf("WARNING: %d factor column(s) have only one level: %s",
                                     length(single_cols),
                                     paste(single_cols, collapse = ", ")))
    }
  }
  
  return(list(valid = valid, messages = messages))
}

#' Validate Numeric Ranges
#' 
#' Check that numeric parameters are within valid ranges
#' 
#' @param value Numeric value to validate
#' @param param_name Name of parameter (for error messages)
#' @param min Minimum allowed value (inclusive, NULL for no minimum)
#' @param max Maximum allowed value (inclusive, NULL for no maximum)
#' @param allow_na Logical, whether NA values are allowed (default: FALSE)
#' @return List with 'valid' (logical) and 'messages' (character vector)
#' @export
validate_numeric_range <- function(value, param_name, min = NULL, max = NULL, allow_na = FALSE) {
  messages <- character()
  valid <- TRUE
  
  # Check for NA
  if (is.na(value)) {
    if (!allow_na) {
      messages <- c(messages, sprintf("ERROR: Parameter '%s' cannot be NA", param_name))
      valid <- FALSE
    }
    return(list(valid = valid, messages = messages))
  }
  
  # Check type
  if (!is.numeric(value)) {
    messages <- c(messages, sprintf("ERROR: Parameter '%s' must be numeric, got: %s", 
                                   param_name, class(value)[1]))
    return(list(valid = FALSE, messages = messages))
  }
  
  # Check minimum
  if (!is.null(min) && value < min) {
    messages <- c(messages, sprintf("ERROR: Parameter '%s' (%s) is below minimum allowed value (%s)",
                                   param_name, value, min))
    valid <- FALSE
  }
  
  # Check maximum
  if (!is.null(max) && value > max) {
    messages <- c(messages, sprintf("ERROR: Parameter '%s' (%s) exceeds maximum allowed value (%s)",
                                   param_name, value, max))
    valid <- FALSE
  }
  
  return(list(valid = valid, messages = messages))
}

#' Validate FASTQ File Pairs
#' 
#' Check that forward and reverse FASTQ files are properly paired
#' 
#' @param fastq_dir Directory containing FASTQ files
#' @param forward_pattern Pattern for forward reads (e.g., "_1.fastq.gz")
#' @param reverse_pattern Pattern for reverse reads (e.g., "_2.fastq.gz")
#' @return List with 'valid' (logical), 'messages' (character vector), and 'pairs' (data frame)
#' @export
validate_fastq_pairs <- function(fastq_dir, 
                                 forward_pattern = "_1\\.(fastq|fq)\\.gz$",
                                 reverse_pattern = "_2\\.(fastq|fq)\\.gz$") {
  messages <- character()
  valid <- TRUE
  
  # Check directory exists
  if (!dir.exists(fastq_dir)) {
    messages <- c(messages, sprintf("ERROR: FASTQ directory does not exist: %s", fastq_dir))
    return(list(valid = FALSE, messages = messages, pairs = NULL))
  }
  
  # Get forward and reverse files
  forward_files <- list.files(fastq_dir, pattern = forward_pattern, full.names = TRUE)
  reverse_files <- list.files(fastq_dir, pattern = reverse_pattern, full.names = TRUE)
  
  # Check if files were found
  if (length(forward_files) == 0) {
    messages <- c(messages, sprintf("ERROR: No forward read files found matching pattern: %s", forward_pattern))
    valid <- FALSE
  }
  
  if (length(reverse_files) == 0) {
    messages <- c(messages, sprintf("ERROR: No reverse read files found matching pattern: %s", reverse_pattern))
    valid <- FALSE
  }
  
  if (!valid) {
    return(list(valid = FALSE, messages = messages, pairs = NULL))
  }
  
  # Extract sample names
  forward_samples <- sub(forward_pattern, "", basename(forward_files))
  reverse_samples <- sub(reverse_pattern, "", basename(reverse_files))
  
  # Check for matching pairs
  unpaired_forward <- setdiff(forward_samples, reverse_samples)
  unpaired_reverse <- setdiff(reverse_samples, forward_samples)
  
  if (length(unpaired_forward) > 0) {
    messages <- c(messages, sprintf("ERROR: %d forward read(s) without matching reverse: %s",
                                   length(unpaired_forward),
                                   paste(head(unpaired_forward, 3), collapse = ", ")))
    valid <- FALSE
  }
  
  if (length(unpaired_reverse) > 0) {
    messages <- c(messages, sprintf("ERROR: %d reverse read(s) without matching forward: %s",
                                   length(unpaired_reverse),
                                   paste(head(unpaired_reverse, 3), collapse = ", ")))
    valid <- FALSE
  }
  
  # Create pairs data frame
  paired_samples <- intersect(forward_samples, reverse_samples)
  if (length(paired_samples) > 0) {
    # Find actual files for each paired sample
    pairs <- data.frame(
      sample = paired_samples,
      forward = sapply(paired_samples, function(s) {
        forward_files[grep(paste0("^", s), basename(forward_files))]
      }),
      reverse = sapply(paired_samples, function(s) {
        reverse_files[grep(paste0("^", s), basename(reverse_files))]
      }),
      stringsAsFactors = FALSE
    )
    
    # Check file sizes
    pairs$forward_size <- file.size(pairs$forward)
    pairs$reverse_size <- file.size(pairs$reverse)
    
    # Warn about small files (< 1KB likely empty or corrupt)
    small_files <- (pairs$forward_size < 1000 | pairs$reverse_size < 1000) & 
                   !is.na(pairs$forward_size) & !is.na(pairs$reverse_size)
    if (any(small_files, na.rm = TRUE)) {
      small_samples <- pairs$sample[small_files]
      messages <- c(messages, sprintf("WARNING: %d sample(s) have very small FASTQ files (< 1KB): %s",
                                     sum(small_files),
                                     paste(head(small_samples, 3), collapse = ", ")))
    }
  } else {
    pairs <- NULL
    messages <- c(messages, "ERROR: No paired samples found")
    valid <- FALSE
  }
  
  return(list(valid = valid, messages = messages, pairs = pairs))
}

#' Check Metadata-Sample Alignment
#' 
#' Verify that metadata samples match FASTQ sample IDs
#' 
#' @param metadata_ids Character vector of sample IDs from metadata
#' @param fastq_ids Character vector of sample IDs from FASTQ files
#' @param require_exact Logical, whether metadata and FASTQ IDs must match exactly (default: FALSE)
#' @return List with 'valid' (logical) and 'messages' (character vector)
#' @export
check_metadata_sample_alignment <- function(metadata_ids, fastq_ids, require_exact = FALSE) {
  messages <- character()
  valid <- TRUE
  
  # Find mismatches
  missing_in_fastq <- setdiff(metadata_ids, fastq_ids)
  missing_in_metadata <- setdiff(fastq_ids, metadata_ids)
  
  if (length(missing_in_fastq) > 0) {
    messages <- c(messages, sprintf("%s: %d sample(s) in metadata but not in FASTQ files: %s",
                                   if (require_exact) "ERROR" else "WARNING",
                                   length(missing_in_fastq),
                                   paste(head(missing_in_fastq, 3), collapse = ", ")))
    if (require_exact) valid <- FALSE
  }
  
  if (length(missing_in_metadata) > 0) {
    messages <- c(messages, sprintf("%s: %d sample(s) in FASTQ files but not in metadata: %s",
                                   if (require_exact) "ERROR" else "WARNING",
                                   length(missing_in_metadata),
                                   paste(head(missing_in_metadata, 3), collapse = ", ")))
    if (require_exact) valid <- FALSE
  }
  
  # Report on matches
  n_matched <- length(intersect(metadata_ids, fastq_ids))
  if (n_matched == 0) {
    messages <- c(messages, "ERROR: No samples match between metadata and FASTQ files")
    valid <- FALSE
  } else {
    messages <- c(messages, sprintf("INFO: %d sample(s) matched between metadata and FASTQ files", n_matched))
  }
  
  return(list(valid = valid, messages = messages))
}

#' Validate Pipeline Configuration
#' 
#' Comprehensive validation of pipeline configuration parameters
#' 
#' @param config List containing pipeline configuration
#' @return List with 'valid' (logical) and 'messages' (character vector)
#' @export
validate_pipeline_config <- function(config) {
  messages <- character()
  valid <- TRUE
  
  # Validate threads
  if (!is.null(config$threads)) {
    thread_check <- validate_numeric_range(config$threads, "threads", min = 1, max = 128)
    if (!thread_check$valid) {
      valid <- FALSE
      messages <- c(messages, thread_check$messages)
    }
  }
  
  # Validate rarefaction depth
  if (!is.null(config$rarefaction$depth)) {
    rare_check <- validate_numeric_range(config$rarefaction$depth, "rarefaction.depth", 
                                         min = 100, max = 1000000)
    if (!rare_check$valid) {
      valid <- FALSE
      messages <- c(messages, rare_check$messages)
    }
  }
  
  # Validate trimming parameters
  if (!is.null(config$trimming$truncLen_forward)) {
    trunc_f_check <- validate_numeric_range(config$trimming$truncLen_forward, 
                                            "trimming.truncLen_forward", min = 50, max = 500)
    if (!trunc_f_check$valid) {
      valid <- FALSE
      messages <- c(messages, trunc_f_check$messages)
    }
  }
  
  if (!is.null(config$trimming$truncLen_reverse)) {
    trunc_r_check <- validate_numeric_range(config$trimming$truncLen_reverse,
                                            "trimming.truncLen_reverse", min = 50, max = 500)
    if (!trunc_r_check$valid) {
      valid <- FALSE
      messages <- c(messages, trunc_r_check$messages)
    }
  }
  
  # Validate maxEE values
  if (!is.null(config$trimming$maxEE_forward)) {
    maxee_f_check <- validate_numeric_range(config$trimming$maxEE_forward,
                                            "trimming.maxEE_forward", min = 0, max = 10)
    if (!maxee_f_check$valid) {
      valid <- FALSE
      messages <- c(messages, maxee_f_check$messages)
    }
  }
  
  if (!is.null(config$trimming$maxEE_reverse)) {
    maxee_r_check <- validate_numeric_range(config$trimming$maxEE_reverse,
                                            "trimming.maxEE_reverse", min = 0, max = 10)
    if (!maxee_r_check$valid) {
      valid <- FALSE
      messages <- c(messages, maxee_r_check$messages)
    }
  }
  
  # Check for required paths
  required_paths <- list(
    input_dir = c("input_dir", "io.input_dir"),
    output_dir = c("output_dir", "io.output_dir"),
    metadata_file = c("metadata_file", "io.metadata_csv")
  )
  
  for (path_name in names(required_paths)) {
    paths_to_check <- required_paths[[path_name]]
    found <- FALSE
    
    for (path_key in paths_to_check) {
      # Handle nested keys with dot notation
      if (grepl("\\.", path_key)) {
        parts <- strsplit(path_key, "\\.")[[1]]
        value <- config
        for (part in parts) {
          value <- value[[part]]
          if (is.null(value)) break
        }
        if (!is.null(value)) {
          found <- TRUE
          break
        }
      } else {
        if (!is.null(config[[path_key]])) {
          found <- TRUE
          break
        }
      }
    }
    
    if (!found) {
      messages <- c(messages, sprintf("WARNING: Configuration missing '%s'", path_name))
    }
  }
  
  return(list(valid = valid, messages = messages))
}

# Print validation helper
print_validation_results <- function(results, prefix = "") {
  if (length(results$messages) > 0) {
    for (msg in results$messages) {
      cat(paste0(prefix, msg, "\n"))
    }
  }
  
  if (results$valid) {
    cat(paste0(prefix, "✓ Validation passed\n"))
  } else {
    cat(paste0(prefix, "✗ Validation failed\n"))
  }
  
  invisible(results$valid)
}
