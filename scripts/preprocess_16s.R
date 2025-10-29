# Preprocessing for 16S: organize trimmed reads under output/<cohort>/trimmed and generate filtered reads under output/<cohort>/filtered

suppressPackageStartupMessages({
  library(dada2)
  library(yaml)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

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

  # Locate trimmed reads: prefer <input_dir>/trimmed, else raw reads
  input_dir <- cfg$io$input_dir
  trimmed_src <- file.path(input_dir, "trimmed")
  if (dir.exists(trimmed_src)) {
    # Copy trimmed FASTQs into cohort output for provenance
    tfiles <- list.files(trimmed_src, pattern = "\\.fastq\\.gz$", full.names = TRUE)
    if (length(tfiles)) file.copy(tfiles, trimmed_dir, overwrite = TRUE)
  }

  # Use trimmed files if present, else fall back to pattern in input_dir
  fnFs <- sort(list.files(trimmed_dir, pattern = "_R1_trimmed\\.fastq\\.gz$", full.names = TRUE))
  fnRs <- sort(list.files(trimmed_dir, pattern = "_R2_trimmed\\.fastq\\.gz$", full.names = TRUE))
  if (length(fnFs) == 0 || length(fnRs) == 0) {
    fnFs <- sort(list.files(input_dir, pattern = "_1\\.fastq\\.gz$|_R1\\.fastq\\.gz$", full.names = TRUE))
    fnRs <- sort(list.files(input_dir, pattern = "_2\\.fastq\\.gz$|_R2\\.fastq\\.gz$", full.names = TRUE))
  }
  if (length(fnFs) != length(fnRs) || length(fnFs) == 0) {
    stop("[preprocess-16S] Could not detect paired FASTQ files.")
  }

  # Sample names
  sample.names <- gsub("_R1_trimmed\\.fastq\\.gz$|_R1\\.fastq\\.gz$|_1\\.fastq\\.gz$", "", basename(fnFs))

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
  utils::write.csv(out, file.path(out_cohort, "analysis", "filtering_summary.csv"), row.names = TRUE)

  invisible(list(filtered_F = filtFs, filtered_R = filtRs, samples = sample.names, out_cohort = out_cohort))
}
