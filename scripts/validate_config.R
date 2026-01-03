#!/usr/bin/env Rscript
# Configuration Validation Script
# Validates pipeline configuration before execution
# Catches errors early to prevent pipeline failures

suppressPackageStartupMessages({
  library(yaml)
})

# Source config inheritance utilities
# Try to get script directory, fallback to current directory
script_dir <- tryCatch({
  dirname(sys.frame(1)$ofile)
}, error = function(e) {
  if (file.exists("scripts/config_inheritance.R")) {
    "scripts"
  } else {
    getwd()
  }
})

# Source the inheritance utilities
inheritance_file <- if (file.exists(file.path(script_dir, "config_inheritance.R"))) {
  file.path(script_dir, "config_inheritance.R")
} else if (file.exists("scripts/config_inheritance.R")) {
  "scripts/config_inheritance.R"
} else {
  stop("Cannot find config_inheritance.R")
}
source(inheritance_file)

validate_config <- function(config_path) {
  errors <- character()
  warnings <- character()
  
  # Check config file exists
  if (!file.exists(config_path)) {
    stop("Configuration file not found: ", config_path)
  }
  
  # Load config with inheritance support
  cfg <- tryCatch({
    # First check if inheritance is used
    if (has_inheritance(config_path)) {
      message("[validate] Loading config with inheritance...")
      load_config_with_inheritance(config_path, dirname(config_path))
    } else {
      yaml::read_yaml(config_path)
    }
  }, error = function(e) {
    stop("Failed to parse YAML config: ", e$message)
  })
  
  # Store original path
  cfg$config_file_path <- normalizePath(config_path)
  
  # === Required Fields ===
  
  # Project section
  if (is.null(cfg$project)) {
    errors <- c(errors, "Missing required section: project")
  } else {
    if (is.null(cfg$project$name)) {
      errors <- c(errors, "Missing required field: project.name")
    }
    if (is.null(cfg$project$sequencing_type)) {
      warnings <- c(warnings, "Missing project.sequencing_type, assuming '16S'")
      cfg$project$sequencing_type <- "16S"
    }
  }
  
  # IO section
  if (is.null(cfg$io)) {
    errors <- c(errors, "Missing required section: io")
  } else {
    if (is.null(cfg$io$input_dir)) {
      errors <- c(errors, "Missing required field: io.input_dir")
    } else {
      # Check input directory exists
      if (!dir.exists(cfg$io$input_dir)) {
        errors <- c(errors, paste0("Input directory not found: ", cfg$io$input_dir))
      } else {
        # Check for FASTQ files
        fastq_files <- list.files(
          cfg$io$input_dir, 
          pattern = "\\.(fastq|fq)\\.gz$", 
          recursive = FALSE
        )
        if (length(fastq_files) == 0) {
          errors <- c(errors, paste0("No FASTQ files found in input directory: ", cfg$io$input_dir))
        }
      }
    }
    
    if (is.null(cfg$io$output_dir)) {
      warnings <- c(warnings, "Missing io.output_dir, using 'output'")
      cfg$io$output_dir <- "output"
    }
    
    # Check metadata if specified
    if (!is.null(cfg$io$metadata_csv)) {
      if (!file.exists(cfg$io$metadata_csv)) {
        errors <- c(errors, paste0("Metadata file not found: ", cfg$io$metadata_csv))
      }
    }
  }
  
  # Amplicon section
  if (is.null(cfg$amplicon)) {
    warnings <- c(warnings, "Missing amplicon section, using defaults")
    cfg$amplicon <- list()
  }
  
  # Check taxonomy databases
  if (!is.null(cfg$amplicon$taxonomy)) {
    if (!is.null(cfg$amplicon$taxonomy$train_set)) {
      if (!file.exists(cfg$amplicon$taxonomy$train_set)) {
        errors <- c(errors, paste0(
          "Taxonomy training set not found: ", cfg$amplicon$taxonomy$train_set,
          "\nDownload from: https://zenodo.org/record/4587955"
        ))
      } else {
        # Check file is not suspiciously small (should be >100MB)
        size_mb <- file.size(cfg$amplicon$taxonomy$train_set) / 1024^2
        if (size_mb < 100) {
          warnings <- c(warnings, paste0(
            "Training set file is small (", round(size_mb, 1), " MB). May be incomplete."
          ))
        }
      }
    }
    
    if (!is.null(cfg$amplicon$taxonomy$species_db)) {
      if (!file.exists(cfg$amplicon$taxonomy$species_db)) {
        warnings <- c(warnings, paste0(
          "Species database not found: ", cfg$amplicon$taxonomy$species_db,
          "\nSpecies-level assignment will be skipped."
        ))
      }
    }
  }
  
  # === Data Type Validation ===
  
  # Check numeric fields
  if (!is.null(cfg$project$threads)) {
    if (!is.numeric(cfg$project$threads) || cfg$project$threads < 1) {
      errors <- c(errors, "project.threads must be a positive integer")
    }
  }
  
  if (!is.null(cfg$plots$dpi)) {
    if (!is.numeric(cfg$plots$dpi) || cfg$plots$dpi < 72) {
      errors <- c(errors, "plots.dpi must be numeric and >= 72")
    }
  }
  
  if (!is.null(cfg$plots$base_size)) {
    if (!is.numeric(cfg$plots$base_size) || cfg$plots$base_size < 6) {
      errors <- c(errors, "plots.base_size must be numeric and >= 6")
    }
  }
  
  # === Metadata Validation ===
  
  if (!is.null(cfg$io$metadata_csv) && file.exists(cfg$io$metadata_csv)) {
    meta <- tryCatch(
      read.csv(cfg$io$metadata_csv, check.names = FALSE, stringsAsFactors = FALSE),
      error = function(e) {
        errors <<- c(errors, paste0("Failed to read metadata CSV: ", e$message))
        return(NULL)
      }
    )
    
    if (!is.null(meta)) {
      # Check ID column exists
      id_col <- cfg$metadata$id_column
      if (is.null(id_col)) {
        id_col <- "SampleID"
        warnings <- c(warnings, "metadata.id_column not specified, assuming 'SampleID'")
      }
      
      if (!id_col %in% colnames(meta)) {
        errors <- c(errors, paste0(
          "ID column '", id_col, "' not found in metadata. ",
          "Available columns: ", paste(colnames(meta), collapse = ", ")
        ))
      } else {
        # Check for duplicates
        if (any(duplicated(meta[[id_col]]))) {
          dupes <- meta[[id_col]][duplicated(meta[[id_col]])]
          errors <- c(errors, paste0(
            "Duplicate sample IDs in metadata: ", 
            paste(head(dupes, 5), collapse = ", ")
          ))
        }
        
        # Check for NAs in ID column
        if (any(is.na(meta[[id_col]]))) {
          errors <- c(errors, "NA values found in sample ID column")
        }
      }
      
      # Check grouping column
      group_col <- NULL
      if (!is.null(cfg$metadata$primary_comparison$group_column)) {
        group_col <- cfg$metadata$primary_comparison$group_column
      } else if (!is.null(cfg$metadata$group_column)) {
        group_col <- cfg$metadata$group_column
        warnings <- c(warnings, "Using deprecated config: metadata.group_column. Use metadata.primary_comparison.group_column")
      }
      
      if (!is.null(group_col)) {
        if (!group_col %in% colnames(meta)) {
          errors <- c(errors, paste0(
            "Group column '", group_col, "' not found in metadata. ",
            "Available columns: ", paste(colnames(meta), collapse = ", ")
          ))
        } else {
          # Check number of groups
          n_groups <- length(unique(meta[[group_col]][!is.na(meta[[group_col]])]))
          if (n_groups < 2) {
            errors <- c(errors, paste0(
              "Need at least 2 groups for comparison. Found ", n_groups, " group(s)."
            ))
          } else {
            message(sprintf("[validate] Found %d groups: %s", 
              n_groups, 
              paste(unique(meta[[group_col]][!is.na(meta[[group_col]])]), collapse = ", ")
            ))
          }
        }
      }
    }
  }
  
  # === Report Results ===
  
  if (length(errors) > 0) {
    attr(cfg, "errors") <- errors
    cat("\n❌ Configuration Validation FAILED:\n")
    cat(paste0("  ", seq_along(errors), ". ", errors, collapse = "\n"), "\n\n")
    stop("Configuration validation failed. Please fix the errors above.")
  }
  
  if (length(warnings) > 0) {
    attr(cfg, "warnings") <- warnings
    cat("\n⚠️  Configuration Warnings:\n")
    cat(paste0("  ", seq_along(warnings), ". ", warnings, collapse = "\n"), "\n\n")
  }
  
  cat("✓ Configuration validation passed\n")
  return(cfg)
}

# Allow running as standalone script (only when directly executed, not sourced)
# Check if this is the main script being executed
if (!interactive()) {
  # Only run standalone mode if this script was called directly
  # When sourced, sys.frames() will have multiple frames
  if (length(sys.frames()) <= 1) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) == 0) {
      cat("Usage: Rscript validate_config.R <config.yaml>\n")
      quit(status = 1)
    }
    
    config_file <- args[1]
    cfg <- validate_config(config_file)
    cat("\nConfiguration is valid and ready for pipeline execution.\n")
  }
}
