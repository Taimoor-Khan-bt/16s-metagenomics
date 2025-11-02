#!/usr/bin/env Rscript
# Error Handling Module
# Provides comprehensive error handling and informative error messages
# Author: Pipeline Development Team
# Last Updated: 2024-11-02

#' Create Informative Error Message
#' 
#' Generate a detailed error message with context and actionable guidance
#' 
#' @param error_type Type of error (e.g., "file_not_found", "validation_failed")
#' @param context List containing contextual information about the error
#' @return Character string with formatted error message
#' @export
create_error_message <- function(error_type, context = list()) {
  
  messages <- list(
    
    # File errors
    file_not_found = function(ctx) {
      sprintf(
        "File Not Found Error:\n  File: %s\n  Context: %s\n  → Check that the file path is correct\n  → Verify the file exists: ls -lh %s",
        ctx$file %||% "unknown",
        ctx$context %||% "reading input file",
        dirname(ctx$file %||% ".")
      )
    },
    
    directory_not_found = function(ctx) {
      sprintf(
        "Directory Not Found Error:\n  Directory: %s\n  Context: %s\n  → Create the directory: mkdir -p %s\n  → Check parent directory permissions",
        ctx$dir %||% "unknown",
        ctx$context %||% "accessing directory",
        ctx$dir %||% "."
      )
    },
    
    empty_directory = function(ctx) {
      sprintf(
        "Empty Directory Error:\n  Directory: %s\n  Expected: %s\n  → Check that input files were copied correctly\n  → Verify file naming patterns match expectations\n  → List directory contents: ls -lh %s",
        ctx$dir %||% "unknown",
        ctx$expected %||% "FASTQ files",
        ctx$dir %||% "."
      )
    },
    
    # Data errors
    empty_dataset = function(ctx) {
      sprintf(
        "Empty Dataset Error:\n  Stage: %s\n  Possible causes:\n  → All samples filtered out due to low quality\n  → Incorrect filtering parameters (too strict)\n  → Problems in upstream processing\n  Actions:\n  → Review filtering summary: %s\n  → Check quality profiles\n  → Adjust maxEE or truncLen parameters",
        ctx$stage %||% "unknown",
        ctx$summary_file %||% "output/filtering_summary.csv"
      )
    },
    
    mismatched_samples = function(ctx) {
      sprintf(
        "Sample Mismatch Error:\n  Metadata samples: %d\n  FASTQ samples: %d\n  Missing from FASTQ: %s\n  Missing from metadata: %s\n  Actions:\n  → Verify sample ID naming is consistent\n  → Check for typos in metadata file\n  → Ensure all FASTQ files were transferred\n  → Review sample ID standards: docs/SAMPLE_ID_STANDARDS.md",
        ctx$n_metadata %||% 0,
        ctx$n_fastq %||% 0,
        paste(head(ctx$missing_fastq %||% character(), 3), collapse = ", "),
        paste(head(ctx$missing_metadata %||% character(), 3), collapse = ", ")
      )
    },
    
    duplicate_samples = function(ctx) {
      sprintf(
        "Duplicate Sample IDs Error:\n  Duplicates found: %s\n  Location: %s\n  Actions:\n  → Check for copy-paste errors in metadata\n  → Verify sample IDs are unique\n  → Review sample naming scheme",
        paste(head(ctx$duplicates %||% character(), 5), collapse = ", "),
        ctx$location %||% "metadata file"
      )
    },
    
    # Metadata errors
    missing_column = function(ctx) {
      sprintf(
        "Missing Column Error:\n  Required column: %s\n  Available columns: %s\n  Context: %s\n  Actions:\n  → Check metadata column names (case-sensitive)\n  → Verify configuration file metadata.columns section\n  → Review column mapping in config",
        ctx$column %||% "unknown",
        paste(ctx$available %||% character(), collapse = ", "),
        ctx$context %||% "reading metadata"
      )
    },
    
    invalid_metadata_values = function(ctx) {
      sprintf(
        "Invalid Metadata Values Error:\n  Column: %s\n  Invalid values: %s\n  Expected: %s\n  Actions:\n  → Check for typos or inconsistent naming\n  → Verify all categorical values are spelled correctly\n  → Remove or replace NA values if not intended",
        ctx$column %||% "unknown",
        paste(head(ctx$invalid %||% character(), 5), collapse = ", "),
        ctx$expected %||% "valid non-NA values"
      )
    },
    
    # Processing errors
    insufficient_reads = function(ctx) {
      sprintf(
        "Insufficient Reads Error:\n  Sample: %s\n  Read count: %d\n  Minimum required: %d\n  Stage: %s\n  Actions:\n  → Check sample quality profiles\n  → This sample may need to be excluded\n  → Review filtering parameters\n  → Consider sequencing depth issues",
        ctx$sample %||% "unknown",
        ctx$read_count %||% 0,
        ctx$min_reads %||% 1000,
        ctx$stage %||% "processing"
      )
    },
    
    rarefaction_depth_error = function(ctx) {
      sprintf(
        "Rarefaction Depth Error:\n  Requested depth: %d\n  Maximum sample depth: %d\n  Samples below threshold: %d\n  Actions:\n  → Lower rarefaction depth in config\n  → Review read depth distribution\n  → Consider excluding low-depth samples\n  → Check: output/filtering_summary.csv",
        ctx$depth %||% 0,
        ctx$max_depth %||% 0,
        ctx$n_below %||% 0
      )
    },
    
    taxonomy_assignment_failed = function(ctx) {
      sprintf(
        "Taxonomy Assignment Error:\n  Database: %s\n  Error details: %s\n  Actions:\n  → Verify database file exists and is not corrupted\n  → Check database file size matches expected\n  → Ensure sufficient memory available\n  → Review database path in configuration",
        ctx$database %||% "unknown",
        ctx$error %||% "unknown error"
      )
    },
    
    # Configuration errors
    invalid_config = function(ctx) {
      sprintf(
        "Invalid Configuration Error:\n  Config file: %s\n  Issue: %s\n  Actions:\n  → Check YAML syntax (indentation, colons, quotes)\n  → Validate against base_config.yaml structure\n  → Run: scripts/validate_config.R %s\n  → Review documentation: docs/CONFIG.md",
        ctx$file %||% "unknown",
        ctx$issue %||% "unknown issue",
        ctx$file %||% "config.yaml"
      )
    },
    
    missing_config_param = function(ctx) {
      sprintf(
        "Missing Configuration Parameter:\n  Parameter: %s\n  Required for: %s\n  Actions:\n  → Add parameter to your config file\n  → Check base_config.yaml for parameter structure\n  → Verify config inheritance is working\n  → Example: %s",
        ctx$param %||% "unknown",
        ctx$context %||% "processing",
        ctx$example %||% "param: value"
      )
    },
    
    # Statistical errors
    insufficient_groups = function(ctx) {
      sprintf(
        "Insufficient Groups Error:\n  Analysis: %s\n  Groups found: %d\n  Minimum required: %d\n  Actions:\n  → Check grouping variable in metadata\n  → Verify group column has multiple levels\n  → Review metadata.primary_comparison.group_column in config",
        ctx$analysis %||% "statistical test",
        ctx$n_groups %||% 0,
        ctx$min_groups %||% 2
      )
    },
    
    insufficient_replicates = function(ctx) {
      sprintf(
        "Insufficient Replicates Error:\n  Group: %s\n  Replicates: %d\n  Minimum required: %d\n  Actions:\n  → Statistical tests require adequate replication\n  → Consider grouping samples differently\n  → Review experimental design",
        ctx$group %||% "unknown",
        ctx$n_replicates %||% 0,
        ctx$min_replicates %||% 3
      )
    },
    
    # Resource errors
    insufficient_memory = function(ctx) {
      sprintf(
        "Insufficient Memory Error:\n  Required: ~%s GB\n  Available: %s GB\n  Process: %s\n  Actions:\n  → Close other applications\n  → Reduce number of threads in config\n  → Process samples in batches\n  → Consider running on a machine with more RAM",
        ctx$required_gb %||% "unknown",
        ctx$available_gb %||% "unknown",
        ctx$process %||% "pipeline step"
      )
    },
    
    insufficient_disk = function(ctx) {
      sprintf(
        "Insufficient Disk Space Error:\n  Required: ~%s GB\n  Available: %s GB\n  Location: %s\n  Actions:\n  → Clean up temporary files\n  → Move output to different disk\n  → Remove old analysis results\n  → Check: df -h",
        ctx$required_gb %||% "unknown",
        ctx$available_gb %||% "unknown",
        ctx$location %||% "output directory"
      )
    }
  )
  
  # Get error message generator
  generator <- messages[[error_type]]
  
  if (is.null(generator)) {
    return(sprintf("Unknown Error Type: %s\nContext: %s", 
                   error_type, 
                   paste(names(context), context, sep = "=", collapse = ", ")))
  }
  
  # Generate message
  generator(context)
}

#' Null Coalescing Operator
#' 
#' Return first argument if not NULL, otherwise return second argument
#' 
#' @param a First value
#' @param b Default value if a is NULL
#' @return a if not NULL, else b
#' @export
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}

#' Stop with Informative Error
#' 
#' Stop execution with a detailed, actionable error message
#' 
#' @param error_type Type of error
#' @param context Context information
#' @param call. Whether to include call in error (default: FALSE)
#' @export
stop_informative <- function(error_type, context = list(), call. = FALSE) {
  msg <- create_error_message(error_type, context)
  stop(msg, call. = call.)
}

#' Warn with Informative Message
#' 
#' Issue a warning with detailed, actionable information
#' 
#' @param error_type Type of warning
#' @param context Context information
#' @param call. Whether to include call in warning (default: FALSE)
#' @param immediate. Whether to show warning immediately (default: TRUE)
#' @export
warn_informative <- function(error_type, context = list(), call. = FALSE, immediate. = TRUE) {
  msg <- create_error_message(error_type, context)
  warning(msg, call. = call., immediate. = immediate.)
}

#' Try with Informative Errors
#' 
#' Wrap tryCatch with informative error handling
#' 
#' @param expr Expression to evaluate
#' @param error_type Type of error for informative message
#' @param context Context information
#' @param silent Whether to suppress original error message (default: FALSE)
#' @return Result of expr or NULL on error
#' @export
try_informative <- function(expr, error_type = NULL, context = list(), silent = FALSE) {
  tryCatch(
    expr,
    error = function(e) {
      if (!is.null(error_type)) {
        context$error <- conditionMessage(e)
        msg <- create_error_message(error_type, context)
        if (!silent) cat("\n", msg, "\n\n", sep = "")
        stop(msg, call. = FALSE)
      } else {
        if (!silent) cat("\nOriginal error:", conditionMessage(e), "\n")
        stop(e)
      }
    }
  )
}

#' Check File Exists
#' 
#' Check if file exists and provide informative error if not
#' 
#' @param file Path to file
#' @param context Description of what the file is used for
#' @export
check_file_exists <- function(file, context = "reading file") {
  if (!file.exists(file)) {
    stop_informative("file_not_found", list(file = file, context = context))
  }
  invisible(TRUE)
}

#' Check Directory Exists
#' 
#' Check if directory exists and provide informative error if not
#' 
#' @param dir Path to directory
#' @param context Description of what the directory is used for
#' @export
check_dir_exists <- function(dir, context = "accessing directory") {
  if (!dir.exists(dir)) {
    stop_informative("directory_not_found", list(dir = dir, context = context))
  }
  invisible(TRUE)
}

#' Check Directory Not Empty
#' 
#' Check if directory contains expected files
#' 
#' @param dir Path to directory
#' @param pattern File pattern to look for (optional)
#' @param expected Description of expected contents
#' @export
check_dir_not_empty <- function(dir, pattern = NULL, expected = "files") {
  check_dir_exists(dir)
  
  files <- if (is.null(pattern)) {
    list.files(dir)
  } else {
    list.files(dir, pattern = pattern)
  }
  
  if (length(files) == 0) {
    stop_informative("empty_directory", list(
      dir = dir,
      expected = expected
    ))
  }
  invisible(TRUE)
}

#' Print Error Summary
#' 
#' Print a formatted summary of errors and warnings
#' 
#' @param errors Character vector of error messages
#' @param warnings Character vector of warning messages
#' @export
print_error_summary <- function(errors = character(), warnings = character()) {
  if (length(errors) > 0) {
    cat("\n❌ ERRORS:\n")
    cat(paste("  ", seq_along(errors), ". ", errors, sep = "", collapse = "\n"), "\n")
  }
  
  if (length(warnings) > 0) {
    cat("\n⚠️  WARNINGS:\n")
    cat(paste("  ", seq_along(warnings), ". ", warnings, sep = "", collapse = "\n"), "\n")
  }
  
  if (length(errors) == 0 && length(warnings) == 0) {
    cat("\n✓ No errors or warnings\n")
  }
  
  invisible(list(errors = errors, warnings = warnings))
}
