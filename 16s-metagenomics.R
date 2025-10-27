# ===============================================================
# 16S rRNA Data Analysis Pipeline (R-only version)
# - ASV Generation, Taxonomy Assignment, Diversity & Phylogenetics
# ===============================================================
# Author: Taimoor Khan
# Email: taimoorkhan007.tk@gmail.com
# ===============================================================

# -------------------------
# 1. Install and Load Required Packages
# -------------------------
# Function to install packages if not available
install_if_missing <- function(pkgs) {
  new_pkgs <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
  if (length(new_pkgs) > 0) {
    if (!require("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(new_pkgs)
  }
}

required_pkgs <- c(
  "dada2", "phyloseq", "ShortRead", "Biostrings",
  "DECIPHER", "vegan", "DESeq2", "tidyverse", "ape"
)

# Install missing packages
install_if_missing(required_pkgs)

# Load packages
invisible(lapply(required_pkgs, library, character.only = TRUE))

# -------------------------
# 2. Set Project Variables
# -------------------------
input_path  <- "HF/trimmed"        # Folder containing trimmed FASTQ files
filt_path   <- "HF/filtered"       # Folder to store filtered FASTQs
output_dir  <- "output"            # Output folder for results
ref_fasta   <- "silva_nr99_v138.1_train_set.fa"  # Reference DB
species_fasta <- "silva_species_assignment_v138.1.fa"
metadata_file <- "metadata.csv"    # Optional metadata (set NULL if not available)

dir.create(filt_path, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# 3. Identify Input Files
# -------------------------
fnFs <- sort(list.files(input_path, pattern = "_R1_trimmed.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(input_path, pattern = "_R2_trimmed.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

stopifnot(length(fnFs) == length(fnRs))

# -------------------------
# 4. Filtering & Trimming
# -------------------------
filtFs <- file.path(filt_path, paste0(sample.names, "_R1_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R2_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen = c(0, 0),          # no truncation; reads already trimmed by cutadapt
  maxEE = c(3, 5),             # relax expected errors to retain more reads
  truncQ = 2,
  minLen = 100,                # drop extremely short reads only
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)

write.csv(out, file.path(output_dir, "filtering_summary.csv"))
cat("Filtering complete.\n")

# Keep only samples that produced filtered files
keep <- file.exists(filtFs) & file.exists(filtRs)
if (!all(keep)) {
  dropped <- sample.names[!keep]
  warning(sprintf("Some samples had no reads pass the filter and were dropped: %s", paste(dropped, collapse = ", ")))
}
filtFs <- filtFs[keep]
filtRs <- filtRs[keep]
sample.names <- sample.names[keep]
if (length(filtFs) == 0) {
  stop("No filtered read files were generated. Consider relaxing filter parameters (increase maxEE, lower truncQ, or adjust truncLen).")
}

# -------------------------
# 5. Learn Error Rates
# -------------------------
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# Optional: visualize error rates
# plotErrors(errF, nominalQ=TRUE)
# plotErrors(errR, nominalQ=TRUE)

# -------------------------
# 6. Denoising and Merging
# -------------------------
names(filtFs) <- names(filtRs) <- sample.names

dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# -------------------------
# 7. Chimera Removal & Length Filtering
# -------------------------
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)
seq_lengths <- nchar(colnames(seqtab.nochim))
# Write sequence length summary to help choose appropriate band
len_summary <- summary(seq_lengths)
write.table(as.data.frame(len_summary), file.path(output_dir, "seq_length_summary.txt"), col.names = NA, quote = FALSE)
# Do not apply a hard length filter by default; keep all non-chimeric ASVs
seqtab.filt <- seqtab.nochim
cat("Retained", ncol(seqtab.filt), "ASVs after chimera removal (no length filter).\n")

# -------------------------
# 8. Taxonomy Assignment
# -------------------------
taxa <- assignTaxonomy(seqtab.filt, refFasta = ref_fasta, minBoot = 50, multithread = TRUE)
taxa <- as.matrix(taxa)
rownames(taxa) <- colnames(seqtab.filt)

# Optional: add species assignment if available
taxa <- addSpecies(taxa, species_fasta)

# -------------------------
# 9. Export ASV Sequences (for MAFFT / FastTree)
# -------------------------
asv_seqs <- colnames(seqtab.filt)
names(asv_seqs) <- paste0("ASV_", seq_along(asv_seqs))
writeXStringSet(DNAStringSet(asv_seqs), file.path(output_dir, "asv_sequences.fasta"))

# -------------------------
# 10. Run MAFFT & FastTree Externally (commented)
# -------------------------
# Run these commands in your terminal (not inside R):
#
# mafft --auto output/asv_sequences.fasta > output/asv_aligned.fasta
# FastTreeMP -nt output/asv_aligned.fasta > output/asv_tree.nwk
#
# After running, continue from next section.

# -------------------------
# 11. Import Tree & Create Phyloseq Object (robust)
# -------------------------

# Ensure taxa and ASV names are consistent
rownames(taxa) <- colnames(seqtab.filt)
asv_tab <- otu_table(as.matrix(seqtab.filt), taxa_are_rows = FALSE)
TAX <- tax_table(as.matrix(taxa))

ps <- phyloseq(asv_tab, TAX)

tree_path <- file.path(output_dir, "asv_tree.nwk")
if (file.exists(tree_path)) {
  message("Found tree: ", tree_path, " — attempting to import and match tips to ASVs.")
  tree <- read.tree(tree_path)
  tips <- tree$tip.label
  asv_names <- taxa_names(ps)                # ASV names in phyloseq (should be sequence strings or ASV IDs)
  # If asv_sequences (from step 9) exist in environment, use mapping seq -> ASV_ID
  # asv_seqs was created earlier as: names(asv_seqs) <- paste0("ASV_", seq_along(asv_seqs))
  if (!exists("asv_seqs")) {
    # Fallback: construct mapping from seqtab column names (they are the sequences)
    seq_vec <- colnames(seqtab.filt)
    id_vec <- paste0("ASV_", seq_along(seq_vec))
    names(seq_vec) <- id_vec
    asv_seqs <- seq_vec
    rm(seq_vec, id_vec)
  }
  seq_to_id <- setNames(names(asv_seqs), asv_seqs) # names(asv_seqs) are ASV IDs; values are sequences
  
  # Helper to try matching in several ways
  match_tips_to_asv <- function(tips, asv_names, seq_to_id) {
    # 1) direct match: tree tip labels are exactly ASV names
    if (all(tips %in% asv_names)) {
      return(list(mode = "tip_is_asv", common = intersect(tips, asv_names)))
    }
    # 2) tree tips are the raw sequences that match ASV sequences
    seq_matches <- intersect(tips, names(seq_to_id))
    if (length(seq_matches) > 0) {
      mapped_ids <- seq_to_id[seq_matches]
      return(list(mode = "tip_is_seq", map = setNames(mapped_ids, seq_matches)))
    }
    # 3) tips have extra header text: take first token before whitespace and try again
    tokens <- sapply(strsplit(tips, "\\s+"), `[`, 1)
    token_matches_seq <- intersect(tokens, names(seq_to_id))
    if (length(token_matches_seq) > 0) {
      mapped_ids <- seq_to_id[token_matches_seq]
      return(list(mode = "tip_token_is_seq", map = setNames(mapped_ids, token_matches_seq), tokens = tokens))
    }
    # 4) maybe the tree tips are ASV IDs but with small differences: try to strip non-alphanumeric
    clean_tips <- gsub("[^A-Za-z0-9_]", "", tips)
    clean_match <- intersect(clean_tips, asv_names)
    if (length(clean_match) > 0) {
      return(list(mode = "tip_cleaned_to_asv", cleaned = clean_tips, common = intersect(clean_tips, asv_names)))
    }
    # nothing matched
    return(NULL)
  }
  
  mt <- match_tips_to_asv(tips, asv_names, seq_to_id)
  
  if (is.null(mt)) {
    warning("Unable to match tree tip labels to ASV names. Tree will NOT be merged. ",
            "First few tree tips:\n", paste(head(tips, 10), collapse = ", "), "\n",
            "First few ASV names:\n", paste(head(asv_names, 10), collapse = ", "))
  } else {
    if (mt$mode == "tip_is_asv") {
      common <- mt$common
      tree_pruned <- ape::keep.tip(tree, common)
      ps <- merge_phyloseq(ps, phy_tree(tree_pruned))
      message("Merged tree: tips matched directly to ASV names. Kept ", length(common), " tips.")
    } else if (mt$mode == "tip_is_seq") {
      # map sequence tips -> ASV IDs and relabel tree tips to ASV IDs
      seqs_found <- names(mt$map)
      mapped_ids <- unname(mt$map)
      # replace tip labels
      idx <- match(seqs_found, tree$tip.label)
      tree$tip.label[idx] <- mapped_ids
      common <- intersect(tree$tip.label, asv_names)
      tree_pruned <- ape::keep.tip(tree, common)
      ps <- merge_phyloseq(ps, phy_tree(tree_pruned))
      message("Merged tree: tree tips were sequences and were remapped to ASV IDs. Kept ", length(common), " tips.")
    } else if (mt$mode == "tip_token_is_seq") {
      # tokens matched: build mapping from full tip -> mapped id using tokens
      tokens <- mt$tokens
      token_matches <- names(mt$map)
      mapped_ids <- unname(mt$map)
      idx <- which(tokens %in% token_matches)
      tree$tip.label[idx] <- mapped_ids[match(tokens[idx], token_matches)]
      common <- intersect(tree$tip.label, asv_names)
      tree_pruned <- ape::keep.tip(tree, common)
      ps <- merge_phyloseq(ps, phy_tree(tree_pruned))
      message("Merged tree: tree tip tokens mapped to ASV sequences and relabeled. Kept ", length(common), " tips.")
    } else if (mt$mode == "tip_cleaned_to_asv") {
      # use cleaned tips
      cleaned <- mt$cleaned
      # set cleaned labels and prune
      tree$tip.label <- cleaned
      common <- intersect(tree$tip.label, asv_names)
      tree_pruned <- ape::keep.tip(tree, common)
      ps <- merge_phyloseq(ps, phy_tree(tree_pruned))
      message("Merged tree after cleaning tip labels. Kept ", length(common), " tips.")
    } else {
      warning("Unexpected mapping mode: ", mt$mode, "; not merging tree.")
    }
  }
} else {
  message("No tree file found at ", tree_path, " — skipping tree import.")
}

# Add metadata if provided
if (!is.null(metadata_file) && file.exists(metadata_file)) {
  meta <- read.csv(metadata_file, row.names = 1)
  ps <- merge_phyloseq(ps, sample_data(meta))
}

# Verify the object
if (!is.null(ps)) {
  validObject(ps)
  saveRDS(ps, file.path(output_dir, "phyloseq_object_raw.rds"))
  message("Phyloseq object created and saved: ", file.path(output_dir, "phyloseq_object_raw.rds"))
} else {
  stop("phyloseq object 'ps' is NULL after tree merge attempt.")
}

# -------------------------
# 12. Normalization (Rarefaction)
# -------------------------
min_reads <- min(sample_sums(ps))
ps_rarefied <- rarefy_even_depth(ps, sample.size = min_reads, rngseed = 42, replace = FALSE)
saveRDS(ps_rarefied, file.path(output_dir, "phyloseq_rarefied.rds"))

# -------------------------
# 13. Alpha & Beta Diversity
# -------------------------
alpha <- estimate_richness(ps_rarefied, measures = c("Observed", "Shannon", "Simpson"))
write.csv(alpha, file.path(output_dir, "alpha_diversity.csv"))

bray_dist <- phyloseq::distance(ps_rarefied, method = "bray")
saveRDS(bray_dist, file.path(output_dir, "beta_bray.rds"))

# -------------------------
# 14. Differential Abundance (if Group column exists)
# -------------------------
if ("sample_data" %in% slotNames(ps) && "Group" %in% colnames(sample_data(ps))) {
  dds <- phyloseq_to_deseq2(ps, ~ Group)
  dds <- DESeq(dds)
  res <- results(dds)
  res_sig <- subset(res, padj < 0.05)
  write.csv(as.data.frame(res_sig), file.path(output_dir, "deseq2_results.csv"))
}

# -------------------------
# 15. Save Outputs & Session Info
# -------------------------
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
cat("✅ Analysis complete. All outputs saved in:", output_dir, "\n")
