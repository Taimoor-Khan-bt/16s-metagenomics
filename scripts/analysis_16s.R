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

# Source analysis modules
modules_dir <- file.path(dirname(sys.frame(1)$ofile), "modules")
if (dir.exists(modules_dir)) {
  source(file.path(modules_dir, "analysis_plan.R"))
  source(file.path(modules_dir, "composite_scores.R"))
  source(file.path(modules_dir, "alpha_diversity.R"))
  source(file.path(modules_dir, "beta_diversity.R"))
  message("[analysis-16S] Loaded modular analysis components")
}


# --- metadata helpers ----------------------------------------------------------
metadata_required <- c("SampleID", "Group")
metadata_optional <- c("SubjectID", "Age", "Sex", "CollectionSite", "SequencingType", "Batch", "LibraryPrep", "dmft")

coerce_metadata_types <- function(df) {
  # Ensure presence of required columns
  for (rc in metadata_required) if (!rc %in% colnames(df)) df[[rc]] <- NA
  # Coerce basic types with safe defaults
  if ("SampleID" %in% colnames(df)) df$SampleID <- as.character(df$SampleID)
  if ("SubjectID" %in% colnames(df)) df$SubjectID <- as.character(df$SubjectID)
  if ("Group" %in% colnames(df)) df$Group <- as.factor(df$Group)
  if ("Age" %in% colnames(df)) df$Age <- suppressWarnings(as.numeric(df$Age))
  if ("Sex" %in% colnames(df)) df$Sex <- factor(df$Sex, levels = c("Male","Female","Other","Unknown"))
  if ("CollectionSite" %in% colnames(df)) df$CollectionSite <- as.factor(df$CollectionSite)
  if ("SequencingType" %in% colnames(df)) df$SequencingType <- factor(df$SequencingType, levels = c("16S","shotgun"))
  if ("Batch" %in% colnames(df)) df$Batch <- as.character(df$Batch)
  if ("LibraryPrep" %in% colnames(df)) df$LibraryPrep <- as.character(df$LibraryPrep)
  if ("dmft" %in% colnames(df)) df$dmft <- suppressWarnings(as.numeric(df$dmft))
  df
}

validate_and_align_metadata <- function(meta_path, ps, analysis_dir, cfg) {
  issues <- c()
  meta <- tryCatch(read.csv(meta_path, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(meta)) {
    issues <- c(issues, sprintf("Could not read metadata: %s", meta_path))
    return(list(meta = NULL, issues = issues))
  }
  
  # Map config-defined id/group columns (support both old and new config structure)
  id_col <- cfg$metadata$id_column %||% "SampleID"
  
  # Support both old config (group_column) and new config (primary_comparison$group_column)
  if (!is.null(cfg$metadata$primary_comparison)) {
    group_col <- cfg$metadata$primary_comparison$group_column %||% "Group"
  } else {
    group_col <- cfg$metadata$group_column %||% "Group"
  }
  
  # If the configured id_column doesn't exist but "Sample" does, rename it
  if (!id_col %in% colnames(meta) && "Sample" %in% colnames(meta)) {
    colnames(meta)[colnames(meta) == "Sample"] <- id_col
    issues <- c(issues, sprintf("Renamed 'Sample' column to '%s'", id_col))
  }
  
  # Ensure SampleID column exists
  if (!id_col %in% colnames(meta)) {
    if (!is.null(rownames(meta)) && all(nzchar(rownames(meta)))) {
      meta[[id_col]] <- rownames(meta)
    } else {
      issues <- c(issues, sprintf("Missing '%s' column; using row numbers", id_col))
      meta[[id_col]] <- seq_len(nrow(meta))
    }
  }
  rownames(meta) <- as.character(meta[[id_col]])
  
  # Type coercion
  meta <- coerce_metadata_types(meta)
  
  # Ensure Group column exists
  if (!group_col %in% colnames(meta)) {
    meta[[group_col]] <- factor("Unknown")
    issues <- c(issues, sprintf("Missing '%s' column; filled with 'Unknown'", group_col))
  }

  sn <- trimws(sample_names(ps))
  rn <- trimws(rownames(meta))
  common <- intersect(sn, rn)
  if (length(common) == 0) {
    issues <- c(issues, "No matching sample IDs between metadata and phyloseq object")
    return(list(meta = NULL, issues = issues))
  }
  
  # Align order to sample_names
  meta_aligned <- meta[sn, , drop = FALSE]
  
  # Save validated metadata and any issues
  utils::write.csv(meta_aligned, file.path(analysis_dir, "metadata_validated.csv"), row.names = TRUE)
  if (length(issues)) writeLines(issues, file.path(analysis_dir, "metadata_issues.txt"))
  list(meta = meta_aligned, issues = issues)
}

# --- phylogeny helper (optional, minimal) -------------------------------------
maybe_build_tree <- function(ps, cfg) {
  build <- isTRUE(cfg$amplicon$phylogeny$build_tree)
  if (!build) return(ps)
  if (!requireNamespace("phangorn", quietly = TRUE) || !requireNamespace("DECIPHER", quietly = TRUE)) return(ps)
  max_tips <- cfg$amplicon$phylogeny$max_tips %||% 500
  keep_taxa <- names(sort(taxa_sums(ps), decreasing = TRUE))[seq_len(min(max_tips, ntaxa(ps)))]
  ps_sub <- prune_taxa(keep_taxa, ps)
  taxa_ids <- taxa_names(ps_sub)
  seqs <- taxa_ids
  # Try to locate ASV sequences
  # Locate ASV fasta in cohort-scoped analysis folder (where run_analysis_16s writes it)
  base_out <- cfg$project$output_dir %||% "output"
  cohort <- cfg$io$cohort
  if (is.null(cohort) || is.na(cohort) || cohort == "") cohort <- basename(cfg$io$input_dir)
  fasta_file <- file.path(base_out, cohort, "analysis", "asv_sequences.fasta")
  if (!all(grepl("^[ACGTNacgtn]+$", seqs)) && file.exists(fasta_file)) {
    fas <- Biostrings::readDNAStringSet(fasta_file)
    names(fas) <- gsub("\t.*$", "", names(fas))
    have <- intersect(names(fas), taxa_ids)
    if (length(have) > 0) {
      fas <- fas[have]; seqs <- as.character(fas); names(seqs) <- names(fas)
    }
  }
  names(seqs) <- taxa_ids
  alignment <- DECIPHER::AlignSeqs(Biostrings::DNAStringSet(seqs), processors = cfg$project$threads %||% 2)
  aln_mat <- as.matrix(alignment)
  phang_align <- phangorn::phyDat(aln_mat, type = "DNA")
  dm <- phangorn::dist.ml(phang_align)
  treeNJ <- phangorn::NJ(dm)
  fit <- phangorn::pml(treeNJ, data = phang_align)
  fitGTR <- phangorn::optim.pml(fit, model = "GTR", optInv = TRUE, optGamma = TRUE, k = 4)
  phyloseq::phy_tree(ps_sub) <- fitGTR$tree
  ps_sub
}

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
    matched_ids <- NULL
  } else {
    filtFs <- pre$filtered_F
    filtRs <- pre$filtered_R
    sample.names <- pre$samples
    matched_ids <- pre$matched_ids  # Get the metadata-matched IDs
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
  # Rename samples to match metadata
  if (!is.null(matched_ids)) {
    rownames(seqtab) <- matched_ids
  }
  utils::write.csv(dim(seqtab), file.path(analysis_dir, "seqtab_dimensions.csv"), row.names = FALSE)

  seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)
  utils::write.csv(dim(seqtab.nochim), file.path(analysis_dir, "seqtab_nochim_dimensions.csv"), row.names = FALSE)

  # Update sample.names to matched IDs for consistency
  if (!is.null(matched_ids)) {
    sample.names <- matched_ids
  }

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

  # Optional metadata with validation
  if (!is.null(cfg$io$metadata_csv) && file.exists(cfg$io$metadata_csv)) {
    va <- validate_and_align_metadata(cfg$io$metadata_csv, ps, analysis_dir, cfg)
    if (!is.null(va$meta)) sample_data(ps) <- sample_data(va$meta)
    if (length(va$issues)) message("[analysis-16S] Metadata notes: ", paste(va$issues, collapse = "; "))
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

  # Beta distances as configured (save to analysis_dir)
  metrics <- cfg$analysis$beta_metrics %||% c("bray")
  for (m in metrics) {
    mm <- tolower(m)
    if (mm == "bray") {
      d <- phyloseq::distance(ps_rarefied, method = "bray")
      saveRDS(d, file.path(analysis_dir, "beta_bray.rds"))
    } else if (grepl("unifrac", mm)) {
      ps_for_uni <- ps_rarefied
      if (is.null(phy_tree(ps_for_uni, errorIfNULL = FALSE))) {
        ps_for_uni <- try(maybe_build_tree(ps_for_uni, cfg), silent = TRUE)
      }
      if (!is.null(phy_tree(ps_for_uni, errorIfNULL = FALSE))) {
        d <- phyloseq::distance(ps_for_uni, method = "unifrac")
        saveRDS(d, file.path(analysis_dir, "beta_unifrac.rds"))
      } else {
        message("[analysis-16S] UniFrac requested but no tree available and build not possible; skipping.")
      }
    } else {
      message("[analysis-16S] Beta metric not implemented: ", m)
    }
  }

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

#' Run Analysis by Plan (NEW - for modular multi-variable analysis)
#' 
#' Executes analyses based on analysis plan from config
#' This function orchestrates: composite scores, alpha diversity, beta diversity
#' 
#' @param ps_path Path to phyloseq object (rarefied)
#' @param cfg Configuration list
#' @return NULL (results saved to files)
#' 
#' @export
run_analysis_by_plan <- function(ps_path, cfg) {
  message("\n[analysis-plan] Starting plan-based analysis")
  
  # Load phyloseq
  ps <- readRDS(ps_path)
  
  # Create analysis plan
  plan <- get_analysis_plan(cfg)
  plan <- validate_analysis_plan(plan, ps)
  print_analysis_plan(plan)
  
  # Add output directory to cfg for module access
  cohort <- cfg$io$cohort %||% basename(cfg$io$input_dir)
  analysis_dir <- file.path(cfg$project$output_dir %||% "output", cohort, "analysis")
  cfg$output <- list(directory = analysis_dir)
  
  # Create composite scores if defined
  if (length(plan$composite_scores) > 0) {
    message("\n[composite] Creating composite scores")
    ps <- create_composite_scores(ps, cfg)
    # Save updated phyloseq
    saveRDS(ps, ps_path)
  }
  
  # Get distance metrics
  distance_metrics <- cfg$analysis$beta$metrics %||% c("bray")
  alpha_metrics <- cfg$analysis$alpha$metrics %||% c("Shannon", "Observed", "Simpson")
  
  # === ALPHA DIVERSITY ===
  message("\n=== ALPHA DIVERSITY ANALYSES ===")
  
  # Primary
  if (plan$run_alpha_primary) {
    run_primary_alpha_analysis(ps, cfg, plan, metrics = alpha_metrics)
  }
  
  # Secondary
  if (length(plan$secondary_vars) > 0) {
    run_secondary_alpha_analysis(ps, cfg, plan, metrics = alpha_metrics)
  }
  
  # Stratified
  if (plan$run_alpha_stratified) {
    run_stratified_alpha_analysis(ps, cfg, plan, metrics = alpha_metrics)
  }
  
  # Continuous correlations
  if (plan$run_alpha_continuous) {
    run_alpha_continuous_correlations(ps, cfg, plan, metrics = alpha_metrics)
  }
  
  # === BETA DIVERSITY ===
  message("\n=== BETA DIVERSITY ANALYSES ===")
  
  for (dist_metric in distance_metrics) {
    message(sprintf("\n[beta] Using distance: %s", dist_metric))
    
    # Primary
    if (plan$run_beta_primary) {
      run_primary_beta_analysis(ps, cfg, plan, distance = dist_metric)
    }
    
    # Secondary
    if (plan$run_beta_stratified) {
      run_secondary_beta_analysis(ps, cfg, plan, distance = dist_metric)
    }
    
    # Multifactor
    if (plan$run_beta_multifactor) {
      run_multifactor_beta_analysis(ps, cfg, plan, distance = dist_metric)
    }
    
    # Pairwise
    if (plan$run_pairwise_comparisons) {
      run_pairwise_beta_comparisons(ps, cfg, plan, distance = dist_metric)
    }
    
    # Dispersion test
    if (isTRUE(cfg$analysis$beta$test_dispersion)) {
      test_beta_dispersion(ps, plan$primary$group_column, distance = dist_metric)
    }
  }
  
  # === DIFFERENTIAL ABUNDANCE ===
  message(sprintf("[DEBUG] plan$run_differential_abundance = %s", plan$run_differential_abundance %||% "NULL"))
  
  if (isTRUE(plan$run_differential_abundance)) {
    message("\n=== DIFFERENTIAL ABUNDANCE ANALYSIS ===")
    
    # Source the DA module
    if (!exists("run_deseq2_analysis")) {
      message("[DA] Sourcing differential_abundance.R module...")
      source(file.path("scripts", "modules", "differential_abundance.R"))
    }
    
    # Run DESeq2
    group_col <- plan$primary$group_column
    rank <- cfg$analysis$differential_abundance$taxonomic_rank %||% "Genus"
    
    message(sprintf("[DA] Running DESeq2 analysis for rank: %s, grouping by: %s", rank, group_col))
    
    da_results <- run_deseq2_analysis(ps, cfg, analysis_dir, group_col, rank)
    
    if (!is.null(da_results)) {
      message("[DA] Differential abundance analysis complete")
    } else {
      message("[DA] WARNING: DA analysis returned NULL")
    }
  }
  
  # === GROUP COMPARISONS ===
  if (isTRUE(cfg$analysis$group_comparisons$enabled)) {
    message("\n=== TAXONOMIC GROUP COMPARISONS ===")
    message("[group-comp] Group comparison plots will be generated during visualization")
  }
  
  message("\n[analysis-plan] Plan-based analysis complete!")
  message(sprintf("[analysis-plan] Results saved to: %s", cfg$output$directory %||% "output"))
  
  invisible(NULL)
}

