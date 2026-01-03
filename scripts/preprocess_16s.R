# Preprocessing for 16S: organize trimmed reads under output/<cohort>/trimmed and generate filtered reads under output/<cohort>/filtered

suppressPackageStartupMessages({
  library(dada2)
  library(yaml)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Smart sample ID matching function
# Handles various naming patterns and extracts the core sample ID
match_sample_ids <- function(file_samples, meta_samples) {
  # Try exact match first
  if (any(file_samples %in% meta_samples)) {
    return(file_samples)
  }
  
  # Common patterns to try (in order of specificity)
  patterns <- c(
    "-LFM[0-9]+$",           # Remove -LFM##### suffix
    "-[A-Z]+[0-9]+$",        # Remove -ABC### suffix (generic)
    "_S[0-9]+$",             # Remove _S## Illumina sample number
    "_[0-9]+$",              # Remove trailing _###
    "^([^_-]+[-_][0-9]+[A-Za-z]*).*", # Extract prefix up to first ID (e.g., KMUN-001 from KMUN-001-LFM4957)
    "^([A-Z]+-[0-9]+[A-Za-z]*).*"     # Extract pattern like ABC-123
  )
  
  for (pattern in patterns) {
    if (grepl("\\(", pattern)) {
      # Extraction pattern
      test_ids <- sub(pattern, "\\1", file_samples)
    } else {
      # Removal pattern
      test_ids <- gsub(pattern, "", file_samples)
    }
    
    # Check if this pattern produces matches
    matches <- sum(test_ids %in% meta_samples)
    if (matches > 0) {
      message(sprintf("[preprocess-16S] Sample ID pattern matched: %s (%d/%d samples)", 
                      pattern, matches, length(test_ids)))
      return(test_ids)
    }
  }
  
  # No pattern worked, return original
  warning("[preprocess-16S] Could not find matching pattern between file names and metadata IDs")
  return(file_samples)
}

run_preprocess_16s <- function(cfg) {
  set.seed(cfg$project$random_seed %||% 1234)
  base_out <- cfg$project$output_dir %||% "output"
  cohort <- cfg$io$cohort
  if (is.null(cohort) || is.na(cohort) || cohort == "") cohort <- basename(cfg$io$input_dir)
  out_cohort <- file.path(base_out, cohort)
  trimmed_dir <- file.path(out_cohort, "trimmed")
  filtered_dir <- file.path(out_cohort, "filtered")
  dir.create(trimmed_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(filtered_dir, recursive = TRUE, showWarnings = FALSE)

  # Locate trimmed reads: prefer <output_dir>/<cohort>/trimmed, then <input_dir>/trimmed, else raw reads
  input_dir <- cfg$io$input_dir
  trimmed_src_input <- file.path(input_dir, "trimmed") # Check for a pre-existing trimmed dir

  # 1. Check if files are already in the destination (e.g. from run_pipeline.sh)
  if (dir.exists(trimmed_dir) && length(list.files(trimmed_dir, pattern = "\\.(fastq|fq)\\.gz$")) > 0) {
     message("[preprocess-16S] Using trimmed files found in cohort output: ", trimmed_dir)
  } else if (dir.exists(trimmed_src_input)) {
    # 2. Check if they are in input/trimmed
    message("[preprocess-16S] Found trimmed files in input subdirectory: ", trimmed_src_input)
    # Copy trimmed FASTQs into cohort output for provenance
    tfiles <- list.files(trimmed_src_input, pattern = "\\.(fastq|fq)\\.gz$", full.names = TRUE)
    if (length(tfiles)) file.copy(tfiles, trimmed_dir, overwrite = TRUE)
  }

  # Use trimmed files if present, else fall back to pattern in input_dir
  fnFs <- sort(list.files(trimmed_dir, pattern = "_R1_trimmed\\.(fastq|fq)\\.gz$", full.names = TRUE))
  fnRs <- sort(list.files(trimmed_dir, pattern = "_R2_trimmed\\.(fastq|fq)\\.gz$", full.names = TRUE))
  if (length(fnFs) == 0 || length(fnRs) == 0) {
    fnFs <- sort(list.files(input_dir, pattern = "(_1|_R1)\\.(fastq|fq)\\.gz$", full.names = TRUE))
    fnRs <- sort(list.files(input_dir, pattern = "(_2|_R2)\\.(fastq|fq)\\.gz$", full.names = TRUE))
  }

  # Sample names
  sample.names <- gsub("_R1_trimmed\\.(fastq|fq)\\.gz$|(_1|_R1)\\.(fastq|fq)\\.gz$", "", basename(fnFs))

  # Filter by metadata if available
  meta_path <- cfg$io$metadata_csv
  if (!is.null(meta_path) && file.exists(meta_path)) {
    meta <- tryCatch(read.csv(meta_path, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(meta)) {
      id_col <- cfg$metadata$id_column %||% "SampleID"
      if (id_col %in% colnames(meta)) {
        meta_samples <- as.character(meta[[id_col]])
        
        # Use smart matching to handle various naming patterns
        matched_ids <- match_sample_ids(sample.names, meta_samples)
        
        # Only keep samples that are in metadata
        keep_idx <- matched_ids %in% meta_samples
        
        if (sum(keep_idx) > 0) {
          message(sprintf("[preprocess-16S] Metadata filtering: %d/%d samples matched", sum(keep_idx), length(sample.names)))
          fnFs <- fnFs[keep_idx]
          fnRs <- fnRs[keep_idx]
          sample.names <- sample.names[keep_idx]
          
          # Store the mapping for later use in analysis
          attr(sample.names, "matched_ids") <- matched_ids[keep_idx]
        } else {
          warning("[preprocess-16S] No samples matched metadata; processing all samples")
        }
      }
    }
  }
  
  # Final check after metadata filtering
  if (length(fnFs) != length(fnRs) || length(fnFs) == 0) {
    stop("[preprocess-16S] Could not detect paired FASTQ files.")
  }
  
  # Verify R1/R2 pairing
  r1_samples <- gsub("_R1_trimmed\\.(fastq|fq)\\.gz$|(_1|_R1)\\.(fastq|fq)\\.gz$", "", basename(fnFs))
  r2_samples <- gsub("_R2_trimmed\\.(fastq|fq)\\.gz$|(_2|_R2)\\.(fastq|fq)\\.gz$", "", basename(fnRs))
  if (!all(r1_samples == r2_samples)) {
    unpaired <- setdiff(union(r1_samples, r2_samples), intersect(r1_samples, r2_samples))
    stop(sprintf("[preprocess-16S] Unmatched R1/R2 pairs detected: %s", paste(head(unpaired, 5), collapse=", ")))
  }
  message(sprintf("[preprocess-16S] Verified %d properly paired FASTQ files", length(fnFs)))

  # Filter and trim (reads are already primer-trimmed)
  filtFs <- file.path(filtered_dir, paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(filtered_dir, paste0(sample.names, "_R_filt.fastq.gz"))
  names(filtFs) <- sample.names
  names(filtRs) <- sample.names

  dir.create(filtered_dir, recursive = TRUE, showWarnings = FALSE)

  out <- suppressWarnings(filterAndTrim(
    fnFs, filtFs,
    fnRs, filtRs,
    truncLen = cfg$amplicon$dada2$truncLen %||% c(0,0),
    maxEE = cfg$amplicon$dada2$maxEE %||% c(3,5),
    truncQ = cfg$amplicon$dada2$truncQ %||% 2,
    minLen = cfg$amplicon$dada2$minLen %||% 100,
    rm.phix = TRUE,
    compress = TRUE,
    multithread = TRUE
  ))

  # Save filtering summary under analysis directory later as well but keep here for traceability
  dir.create(file.path(out_cohort, "analysis"), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(out, file.path(out_cohort, "analysis", "filtering_summary.csv"), row.names = TRUE)

  invisible(list(
    filtered_F = filtFs, 
    filtered_R = filtRs, 
    samples = sample.names, 
    matched_ids = attr(sample.names, "matched_ids"),  # Pass the matched metadata IDs
    out_cohort = out_cohort
  ))
}
