#!/usr/bin/env Rscript
# Pipeline Input Validation Wrapper
# Comprehensive validation of all pipeline inputs before processing
# Author: Pipeline Development Team
# Last Updated: 2024-11-02

suppressPackageStartupMessages({
  library(yaml)
  library(optparse)
})

# Simple null coalescing helper
`%||%` <- function(a, b) if (is.null(a)) b else a

# Source modules with path resolution
script_dir <- tryCatch({
  dirname(sys.frame(1)$ofile)
}, error = function(e) {
  commandArgs(trailingOnly = FALSE)
  this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(this_file) > 0) {
    dirname(this_file)
  } else {
    getwd()
  }
})

if (is.null(script_dir) || length(script_dir) == 0 || script_dir == "" || is.na(script_dir)) {
  script_dir <- file.path(getwd(), "scripts")
}

module_dir <- file.path(script_dir, "modules")
if (!dir.exists(module_dir)) {
  # Try alternative paths
  if (dir.exists(file.path(getwd(), "scripts", "modules"))) {
    module_dir <- file.path(getwd(), "scripts", "modules")
  } else if (dir.exists(file.path(dirname(script_dir), "scripts", "modules"))) {
    module_dir <- file.path(dirname(script_dir), "scripts", "modules")
  }
}

source(file.path(module_dir, "input_validation.R"))
source(file.path(module_dir, "error_handling.R"))

# Also source config inheritance if available
config_inheritance_script <- file.path(dirname(module_dir), "config_inheritance.R")
if (file.exists(config_inheritance_script)) {
  source(config_inheritance_script)
  use_inheritance <- TRUE
} else {
  use_inheritance <- FALSE
}

# Parse command line arguments
option_list <- list(
  make_option(c("-c", "--config"), type = "character", default = NULL,
              help = "Path to configuration file", metavar = "FILE"),
  make_option(c("-m", "--metadata"), type = "character", default = NULL,
              help = "Path to metadata file", metavar = "FILE"),
  make_option(c("-f", "--fastq-dir"), type = "character", default = NULL,
              help = "Path to FASTQ directory", metavar = "DIR"),
  make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
              help = "Print verbose output")
)

parser <- OptionParser(usage = "%prog [options]", option_list = option_list,
                      description = "\nValidate pipeline inputs before processing")

args <- parse_args(parser)

# Helper function for verbose printing
vcat <- function(...) {
  if (args$verbose) {
    cat(...)
  }
}

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  16S Metagenomics Pipeline - Input Validation\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Track validation results
all_errors <- character()
all_warnings <- character()
validation_passed <- TRUE

# Initialize variables
pipeline_config <- NULL
pipeline_metadata <- NULL
fastq_samples <- NULL

# 1. Validate Configuration
cat("1️⃣  Validating Configuration...\n")
if (is.null(args$config)) {
  all_errors <- c(all_errors, "No configuration file specified (use --config)")
  validation_passed <- FALSE
} else {
  tryCatch({
    check_file_exists(args$config, "configuration validation")
    
    # Load config with inheritance if available
    if (use_inheritance) {
      config <- load_config_with_inheritance(args$config)
    } else {
      config <- yaml::read_yaml(args$config)
    }
    
    # Validate configuration parameters
    config_validation <- validate_pipeline_config(config)
    
    if (!config_validation$valid) {
      all_errors <- c(all_errors, config_validation$messages[grepl("^ERROR:", config_validation$messages)])
      validation_passed <- FALSE
    }
    
    warnings <- config_validation$messages[grepl("^WARNING:", config_validation$messages)]
    if (length(warnings) > 0) {
      all_warnings <- c(all_warnings, warnings)
    }
    
    if (config_validation$valid) {
      cat("   ✓ Configuration structure valid\n")
      vcat(sprintf("   • Project: %s\n", config$project$name %||% "not specified"))
      vcat(sprintf("   • Threads: %d\n", config$project$threads %||% config$threads %||% "not specified"))
      vcat(sprintf("   • Rarefaction depth: %s\n", config$diversity$rarefaction$depth %||% config$rarefaction$depth %||% "not specified"))
    } else {
      cat("   ✗ Configuration validation failed\n")
    }
    
    # Store config for later use
    pipeline_config <- config
    
  }, error = function(e) {
    all_errors <- c(all_errors, sprintf("Configuration error: %s", conditionMessage(e)))
    validation_passed <- FALSE
    cat("   ✗ Configuration validation failed\n")
  })
}

cat("\n")

# 2. Validate Metadata
cat("2️⃣  Validating Metadata...\n")

metadata_file <- args$metadata %||% pipeline_config$metadata_file %||% pipeline_config$io$metadata_csv
sample_id_col <- pipeline_config$metadata$columns$sample_id %||% 
                 pipeline_config$metadata$id_column %||% 
                 "SampleID"

if (is.null(metadata_file)) {
  all_errors <- c(all_errors, "No metadata file specified")
  validation_passed <- FALSE
} else {
  tryCatch({
    check_file_exists(metadata_file, "metadata validation")
    metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
    
    # Get required columns from config
    required_cols <- c(sample_id_col)
    group_col <- pipeline_config$metadata$primary_comparison$group_column %||%
                 pipeline_config$metadata$group_column %||%
                 pipeline_config$analysis$compare_by
    
    if (!is.null(group_col)) {
      required_cols <- c(required_cols, group_col)
    }
    
    # Validate metadata structure
    meta_validation <- validate_metadata_structure(metadata, 
                                                   required_cols = required_cols,
                                                   sample_id_col = sample_id_col)
    
    if (!meta_validation$valid) {
      all_errors <- c(all_errors, meta_validation$messages[grepl("^ERROR:", meta_validation$messages)])
      validation_passed <- FALSE
    }
    
    warnings <- meta_validation$messages[grepl("^WARNING:", meta_validation$messages)]
    if (length(warnings) > 0) {
      all_warnings <- c(all_warnings, warnings)
    }
    
    if (meta_validation$valid) {
      cat(sprintf("   ✓ Metadata structure valid\n"))
      cat(sprintf("   • Samples: %d\n", nrow(metadata)))
      cat(sprintf("   • Columns: %d\n", ncol(metadata)))
      vcat(sprintf("   • Column names: %s\n", paste(colnames(metadata), collapse = ", ")))
      
      # Check for grouping variable
      group_col <- pipeline_config$metadata$primary_comparison$group_column %||%
                   pipeline_config$metadata$group_column %||%
                   pipeline_config$analysis$compare_by
      
      if (!is.null(group_col)) {
        if (group_col %in% colnames(metadata)) {
          groups <- unique(metadata[[group_col]])
          cat(sprintf("   • Groups (%s): %s\n", group_col, paste(groups, collapse = ", ")))
          
          if (length(groups) < 2) {
            all_warnings <- c(all_warnings, sprintf("Only %d group(s) found for statistical comparisons", length(groups)))
          }
        }
      }
      
      # Store metadata for later use
      pipeline_metadata <- metadata
      
    } else {
      cat("   ✗ Metadata validation failed\n")
    }
    
  }, error = function(e) {
    all_errors <- c(all_errors, sprintf("Metadata error: %s", conditionMessage(e)))
    validation_passed <- FALSE
    cat("   ✗ Metadata validation failed\n")
  })
}

cat("\n")

# 3. Validate FASTQ Files
cat("3️⃣  Validating FASTQ Files...\n")

fastq_dir <- args$fastq_dir %||% pipeline_config$input_dir %||% pipeline_config$io$input_dir

if (is.null(fastq_dir)) {
  all_errors <- c(all_errors, "No FASTQ directory specified")
  validation_passed <- FALSE
} else {
  tryCatch({
    check_dir_exists(fastq_dir, "FASTQ file validation")
    
    # Get file patterns from config or use defaults
    forward_pattern <- "_1\\.fastq\\.gz$"
    reverse_pattern <- "_2\\.fastq\\.gz$"
    
    # Validate file pairs
    fastq_validation <- validate_fastq_pairs(fastq_dir, forward_pattern, reverse_pattern)
    
    if (!fastq_validation$valid) {
      all_errors <- c(all_errors, fastq_validation$messages[grepl("^ERROR:", fastq_validation$messages)])
      validation_passed <- FALSE
    }
    
    warnings <- fastq_validation$messages[grepl("^WARNING:", fastq_validation$messages)]
    if (length(warnings) > 0) {
      all_warnings <- c(all_warnings, warnings)
    }
    
    if (fastq_validation$valid) {
      cat("   ✓ FASTQ file pairs valid\n")
      cat(sprintf("   • Sample pairs: %d\n", nrow(fastq_validation$pairs)))
      
      if (args$verbose && !is.null(fastq_validation$pairs)) {
        total_size <- sum(fastq_validation$pairs$forward_size + fastq_validation$pairs$reverse_size)
        cat(sprintf("   • Total size: %.2f GB\n", total_size / 1e9))
        cat(sprintf("   • Average size per sample: %.2f MB\n", 
                   total_size / nrow(fastq_validation$pairs) / 1e6))
      }
      
      # Store FASTQ info for later use
      fastq_samples <- fastq_validation$pairs$sample
      
    } else {
      cat("   ✗ FASTQ validation failed\n")
    }
    
  }, error = function(e) {
    all_errors <- c(all_errors, sprintf("FASTQ error: %s", conditionMessage(e)))
    validation_passed <- FALSE
    cat("   ✗ FASTQ validation failed\n")
  })
}

cat("\n")

# 4. Cross-Validate Metadata and FASTQ
cat("4️⃣  Cross-Validating Metadata and FASTQ Samples...\n")

if (exists("pipeline_metadata") && exists("fastq_samples")) {
  tryCatch({
    metadata_ids <- pipeline_metadata[[sample_id_col]]
    
    # Check alignment
    alignment_check <- check_metadata_sample_alignment(metadata_ids, fastq_samples, 
                                                       require_exact = FALSE)
    
    if (!alignment_check$valid) {
      all_errors <- c(all_errors, alignment_check$messages[grepl("^ERROR:", alignment_check$messages)])
      validation_passed <- FALSE
    }
    
    warnings <- alignment_check$messages[grepl("^(WARNING|INFO):", alignment_check$messages)]
    if (length(warnings) > 0) {
      for (w in warnings) {
        if (grepl("^INFO:", w)) {
          cat(sprintf("   %s\n", w))
        } else {
          all_warnings <- c(all_warnings, w)
        }
      }
    }
    
    if (alignment_check$valid) {
      cat("   ✓ Metadata and FASTQ samples aligned\n")
    } else {
      cat("   ✗ Metadata-FASTQ alignment failed\n")
    }
    
  }, error = function(e) {
    all_errors <- c(all_errors, sprintf("Alignment error: %s", conditionMessage(e)))
    validation_passed <- FALSE
    cat("   ✗ Cross-validation failed\n")
  })
} else {
  cat("   ⊘ Skipped (previous validations failed)\n")
}

cat("\n")

# Print Summary
cat("═══════════════════════════════════════════════════════════════\n")
cat("  Validation Summary\n")
cat("═══════════════════════════════════════════════════════════════\n")

print_error_summary(all_errors, all_warnings)

if (validation_passed) {
  cat("\n✅ All validations passed! Pipeline is ready to run.\n\n")
  quit(status = 0)
} else {
  cat("\n❌ Validation failed. Please address the errors above before running the pipeline.\n\n")
  quit(status = 1)
}
