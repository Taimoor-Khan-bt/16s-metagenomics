# Analysis for 16S: DADA2 inference, chimera removal, taxonomy, phyloseq creation, optional rarefaction, and stats

suppressPackageStartupMessages({
  library(dada2)
  library(phyloseq)
  library(Biostrings)
  library(vegan)
  library(yaml)
  library(ggplot2)
  library(tidyverse)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

run_analysis_16s <- function(cfg, pre = NULL) {
  set.seed(cfg$project$random_seed %||% 1234)
  base_out <- cfg$project$output_dir %||% "output"
  cohort <- cfg$io$cohort
  if (is.null(cohort) || is.na(cohort) || cohort == "") cohort <- basename(cfg$io$input_dir)
  out_cohort <- file.path(base_out, cohort)
  analysis_dir <- file.path(out_cohort, "analysis")
  dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)

  # Discover filtered files if not provided
  if (is.null(pre)) {
    filtered_dir <- file.path(out_cohort, "filtered")
    filtFs <- sort(list.files(filtered_dir, pattern = "_F_filt\\.fastq\\.gz$", full.names = TRUE))
    filtRs <- sort(list.files(filtered_dir, pattern = "_R_filt\\.fastq\\.gz$", full.names = TRUE))
    sample.names <- gsub("_F_filt\\.fastq\\.gz$", "", basename(filtFs))
  } else {
    filtFs <- pre$filtered_F; filtRs <- pre$filtered_R; sample.names <- pre$samples
  }

  if (length(filtFs) == 0) stop("[analysis-16S] No filtered files found.")

  # Learn error rates
  errF <- suppressWarnings(learnErrors(filtFs, multithread = TRUE))
  errR <- suppressWarnings(learnErrors(filtRs, multithread = TRUE))

  # Inference and merge
  dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
  dadaRs <- dada(filtRs, err = errR, multithread = TRUE)
  mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

  seqtab <- makeSequenceTable(mergers)
  utils::write.csv(dim(seqtab), file.path(analysis_dir, "seqtab_dimensions.csv"), row.names = FALSE)

  seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)
  utils::write.csv(dim(seqtab.nochim), file.path(analysis_dir, "seqtab_nochim_dimensions.csv"), row.names = FALSE)

  # Length summary
  asv_lengths <- nchar(colnames(seqtab.nochim))
  sink(file.path(analysis_dir, "seq_length_summary.txt"));
  cat("ASV length summary (nchar of sequences)\n\n");
  print(summary(asv_lengths));
  sink()

  # Read tracking (core stages)
  getN <- function(x) sum(getUniques(x))
  track <- cbind(
    denoisedF = sapply(dadaFs, getN),
    denoisedR = sapply(dadaRs, getN),
    merged = sapply(mergers, getN),
    nonchim = rowSums(seqtab.nochim)
  )
  rownames(track) <- sample.names
  utils::write.csv(track, file.path(analysis_dir, "read_tracking.csv"))

  # Taxonomy
  silva_train <- cfg$amplicon$taxonomy$silva_train_set %||% "silva_nr99_v138.1_train_set.fa"
  silva_species <- cfg$amplicon$taxonomy$silva_species %||% "silva_species_assignment_v138.1.fa"
  taxa <- assignTaxonomy(seqtab.nochim, silva_train, multithread = TRUE)
  taxa <- addSpecies(taxa, silva_species)

  # Save ASV sequences
  uniqueSeqs <- colnames(seqtab.nochim)
  names(uniqueSeqs) <- paste0("ASV", seq_along(uniqueSeqs))
  writeXStringSet(DNAStringSet(uniqueSeqs), file.path(analysis_dir, "asv_sequences.fasta"))

  # Create phyloseq object
  ps <- phyloseq(
    otu_table(seqtab.nochim, taxa_are_rows = FALSE),
    tax_table(taxa)
  )

  # Optional metadata
  if (!is.null(cfg$io$metadata_csv) && file.exists(cfg$io$metadata_csv)) {
    meta <- tryCatch(read.csv(cfg$io$metadata_csv, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(meta)) {
      rn <- trimws(rownames(meta)); rownames(meta) <- rn
      sn <- sample_names(ps); sn <- trimws(sn)
      common <- intersect(sn, rownames(meta))
      if (length(common) > 0) {
        meta_sub <- meta[sn, , drop = FALSE]
        sample_data(ps) <- sample_data(meta_sub)
      } else {
        message("[analysis-16S] metadata.csv found but sample IDs do not match; proceeding without metadata.")
      }
    }
  }

  # Phylogeny (optional) will be handled in visualization if enabled; we keep ps as-is here

  # Rename ASVs to stable IDs after taxonomy
  taxa_names(ps) <- paste0("ASV", seq_along(taxa_names(ps)))
  saveRDS(ps, file.path(analysis_dir, "phyloseq_object_raw.rds"))

  # Filtering for downstream analysis
  ps_filtered <- subset_taxa(ps,
    Kingdom == "Bacteria" & !is.na(Phylum) & Phylum != "" & Phylum != "uncharacterized" &
    (is.na(Order) | Order != "Chloroplast") & (is.na(Family) | Family != "Mitochondria")
  )

  # Rarefy if requested
  if (!is.null(cfg$analysis$rarefaction_depth) && cfg$analysis$rarefaction_depth > 0) {
    depth <- cfg$analysis$rarefaction_depth
  } else {
    depth <- min(sample_sums(ps_filtered))
  }
  ps_rarefied <- suppressWarnings(rarefy_even_depth(ps_filtered, sample.size = depth, rngseed = cfg$project$random_seed %||% 1234))
  saveRDS(ps_rarefied, file.path(analysis_dir, "phyloseq_rarefied.rds"))

  # Alpha diversity table (saved for analysis outputs)
  alpha_df <- estimate_richness(ps_rarefied, measures = c("Observed", "Shannon", "InvSimpson"))
  alpha_df$Sample <- rownames(alpha_df)
  if (!is.null(sample_data(ps_rarefied, errorIfNULL = FALSE))) {
    md <- as.data.frame(sample_data(ps_rarefied)); md$Sample <- rownames(md)
    alpha_df <- merge(alpha_df, md, by = "Sample", all.x = TRUE)
  }
  utils::write.csv(alpha_df, file.path(analysis_dir, "alpha_diversity.csv"), row.names = FALSE)

  # Optional PERMANOVA (stored as text if applicable)
  md <- as.data.frame(sample_data(ps_rarefied))
  if (!is.null(md) && nrow(md) > 1 && "Group" %in% colnames(md) && length(na.omit(unique(md$Group))) > 1) {
    bray <- phyloseq::distance(ps_rarefied, method = "bray")
    md_df <- data.frame(Group = md$Group, row.names = rownames(md))
    ad <- vegan::adonis2(bray ~ Group, data = md_df, permutations = 999)
    capture.output(ad, file = file.path(analysis_dir, "permanova_group.txt"))
  }

  invisible(list(ps_raw = file.path(analysis_dir, "phyloseq_object_raw.rds"), ps_rarefied = file.path(analysis_dir, "phyloseq_rarefied.rds")))
}
