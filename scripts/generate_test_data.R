# Test Dataset Generator for 16S Pipeline
# Creates minimal synthetic test data for integration testing

suppressPackageStartupMessages({
  library(dada2)
  library(Biostrings)
})

#' Generate synthetic test dataset
#' 
#' Creates a minimal dataset with 6 samples for testing pipeline
generate_test_dataset <- function(outdir = "test_data") {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  message("Generating synthetic test dataset...")
  
  # Define test samples
  samples <- data.frame(
    sample_id = c("TEST-001", "TEST-002", "TEST-003", 
                  "TEST-004", "TEST-005", "TEST-006"),
    group = c("Control", "Control", "Control", 
              "Treatment", "Treatment", "Treatment"),
    age = c(25, 30, 28, 26, 32, 29),
    bmi = c(22.5, 24.1, 23.8, 21.9, 25.3, 23.2),
    stringsAsFactors = FALSE
  )
  
  # Save metadata
  write.csv(samples, file.path(outdir, "metadata.csv"), row.names = FALSE)
  message("✓ Created metadata.csv")
  
  # Generate synthetic ASV sequences (16S V3-V4 region, ~420bp)
  set.seed(1234)
  n_asvs <- 50
  
  # Base sequence (realistic 16S)
  base_seq <- paste0(
    "TACGTAGGTGGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGAGCGTAGGCGGACGCGCAAGTCTGATGTGAAAGC",
    "CCGGGCTCAACCTTGGAACTGCATTTGGAACTGTCAGGCTTGAATCTCGGAGAGGTAAGCGGAATTCCTAGTGTAGCG",
    "GTGAAATGCGTAGATATTAGGAGGAACACCAGTGGCGAAGGCGGCTTACTGGACGATAACTGACGCTGAGGCTCGAAA",
    "GCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAATGTTAGCCGTGGGAGAGCACACA",
    "CTTAGCGGGTCTGCGAAGGTTGAGCCGATGAAGTGACCGCCTGGGAAGTACGGTCGCAAGGCTGAAACTCAAAGGAA"
  )
  
  # Create variants with mutations
  asv_seqs <- character(n_asvs)
  for (i in 1:n_asvs) {
    seq <- base_seq
    # Add 1-3 random mutations per ASV
    n_muts <- sample(1:3, 1)
    positions <- sample(1:nchar(seq), n_muts)
    for (pos in positions) {
      bases <- c("A", "T", "G", "C")
      substr(seq, pos, pos) <- sample(bases, 1)
    }
    asv_seqs[i] <- seq
  }
  
  # Generate abundance matrix (samples x ASVs)
  # Control group: ASVs 1-30 dominant
  # Treatment group: ASVs 20-50 dominant (shift in community)
  abund_matrix <- matrix(0, nrow = 6, ncol = n_asvs)
  rownames(abund_matrix) <- samples$sample_id
  colnames(abund_matrix) <- paste0("ASV", 1:n_asvs)
  
  for (i in 1:6) {
    if (samples$group[i] == "Control") {
      # Control: more of ASVs 1-30
      abund_matrix[i, 1:30] <- rpois(30, lambda = 100)
      abund_matrix[i, 31:50] <- rpois(20, lambda = 20)
    } else {
      # Treatment: more of ASVs 20-50
      abund_matrix[i, 1:20] <- rpois(20, lambda = 20)
      abund_matrix[i, 21:50] <- rpois(30, lambda = 100)
    }
  }
  
  # Ensure minimum depth
  abund_matrix <- abund_matrix + rpois(length(abund_matrix), lambda = 5)
  
  message("✓ Generated ", n_asvs, " synthetic ASVs")
  message("✓ Mean read depth: ", round(mean(rowSums(abund_matrix))))
  
  # Generate FASTQ files
  fastq_dir <- file.path(outdir, "fastq")
  dir.create(fastq_dir, showWarnings = FALSE)
  
  for (i in 1:nrow(samples)) {
    samp <- samples$sample_id[i]
    
    # R1 and R2 filenames (use _1 and _2 pattern for compatibility with trimming script)
    r1_file <- file.path(fastq_dir, paste0(samp, "_1.fastq.gz"))
    r2_file <- file.path(fastq_dir, paste0(samp, "_2.fastq.gz"))
    
    # Generate reads for this sample
    reads_r1 <- character()
    reads_r2 <- character()
    quals <- character()
    
    read_id <- 1
    for (j in 1:n_asvs) {
      n_reads <- abund_matrix[i, j]
      if (n_reads > 0) {
        for (k in 1:n_reads) {
          # Forward read (first 250bp)
          fwd_seq <- substr(asv_seqs[j], 1, 250)
          # Reverse read (last 250bp, reverse complement)
          rev_seq <- substr(asv_seqs[j], nchar(asv_seqs[j]) - 249, nchar(asv_seqs[j]))
          
          # Add to reads
          reads_r1 <- c(reads_r1, 
                       paste0("@read_", read_id),
                       fwd_seq,
                       "+",
                       paste(rep("I", nchar(fwd_seq)), collapse = ""))
          reads_r2 <- c(reads_r2,
                       paste0("@read_", read_id),
                       rev_seq,
                       "+",
                       paste(rep("I", nchar(rev_seq)), collapse = ""))
          read_id <- read_id + 1
        }
      }
    }
    
    # Write compressed FASTQ
    writeLines(reads_r1, con <- gzfile(r1_file, "w")); close(con)
    writeLines(reads_r2, con <- gzfile(r2_file, "w")); close(con)
  }
  
  message("✓ Created FASTQ files for ", nrow(samples), " samples")
  message("\nTest dataset ready in: ", outdir)
  message("Samples: ", nrow(samples))
  message("ASVs: ", n_asvs)
  message("Total reads: ", sum(abund_matrix))
  
  return(invisible(list(
    samples = samples,
    asv_seqs = asv_seqs,
    abundance = abund_matrix
  )))
}

# Run if executed directly
if (!interactive()) {
  generate_test_dataset("test_data")
}
