# =================================================================
# Complete 16S rRNA Sequence Analysis Pipeline
# From trimmed reads to advanced analysis
# 
# Author: Taimoor Khan
# Email: taimoorkhan007.tk@gmail.com
# Version: 1.0.0
# Last Updated: October 20, 2025
# =================================================================

# Function to safely load a package
safe_library <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    message(sprintf("Package %s is not available. Some functionality may be limited.", pkg))
    return(FALSE)
  }
  return(TRUE)
}

# Load required packages
pkg_status <- list()

# Core packages
pkg_status$core <- c(
  dada2 = safe_library("dada2"),
  phyloseq = safe_library("phyloseq"),
  vegan = safe_library("vegan"),
  tidyverse = safe_library("tidyverse"),
  ape = safe_library("ape"),
  DESeq2 = safe_library("DESeq2"),
  DECIPHER = safe_library("DECIPHER"),
  phangorn = safe_library("phangorn"),
  ggtree = safe_library("ggtree"),
  RColorBrewer = safe_library("RColorBrewer"),
  ComplexHeatmap = safe_library("ComplexHeatmap"),
  indicspecies = safe_library("indicspecies")
)

# Create output directory
output_dir <- "output"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ==================================================================
# 1. Initial Data Processing with DADA2
# ==================================================================
cat("🔍 Starting DADA2 processing pipeline...\n")

# File paths
path <- "HF/trimmed"
# Match only the expected gzipped files to avoid accidental suffix issues
fnFs <- sort(list.files(path, pattern = "_R1_trimmed\\.fastq\\.gz$", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_trimmed\\.fastq\\.gz$", full.names = TRUE))

# Extract sample names robustly (remove exact suffix)
sample.names <- sub("_R1_trimmed\\.fastq\\.gz$", "", basename(fnFs))

# Quality filter and trim
cat("Filtering and trimming reads...\n")
filtFs <- file.path(output_dir, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(output_dir, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

dir.create(file.path(output_dir, "filtered"), recursive = TRUE, showWarnings = FALSE)

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen = c(0, 0),          # reads already primer-trimmed by cutadapt
  maxEE = c(3, 5),             # relax to retain more reads
  truncQ = 2,
  minLen = 100,                # drop extremely short reads only
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)

write.csv(out, file.path(output_dir, "filtering_summary.csv"))

# Keep only samples that produced filtered files (some may have been dropped entirely)
keep <- file.exists(filtFs) & file.exists(filtRs)
if (!all(keep)) {
  dropped <- sample.names[!keep]
  warning(sprintf("Some samples had no reads pass the filter and were dropped: %s", paste(dropped, collapse = ", ")))
  sample.names <- sample.names[keep]
  filtFs <- filtFs[keep]
  filtRs <- filtRs[keep]
}
if (length(filtFs) == 0) {
  stop("No filtered read files were generated. Consider relaxing filter parameters (e.g., increase maxEE, lower truncQ, or adjust minLen).")
}

# Learn error rates
cat("Learning error rates...\n")
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Sample inference
cat("Performing sample inference...\n")
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# Merge paired reads
cat("Merging paired reads...\n")
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
write.csv(dim(seqtab), file.path(output_dir, "seqtab_dimensions.csv"))

# Remove chimeras
cat("Removing chimeras...\n")
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
write.csv(dim(seqtab.nochim), file.path(output_dir, "seqtab_nochim_dimensions.csv"))
# Summarize ASV lengths to guide optional length filtering
asv_lengths <- nchar(colnames(seqtab.nochim))
# Write a simple summary to text to avoid data.frame/dimname issues
sink(file.path(output_dir, "seq_length_summary.txt"))
cat("ASV length summary (nchar of sequences)\n\n")
print(summary(asv_lengths))
sink()

# Track reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, 
               sapply(dadaFs, getN),
               sapply(dadaRs, getN),
               sapply(mergers, getN),
               rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
write.csv(track, file.path(output_dir, "read_tracking.csv"))

if (ncol(seqtab.nochim) == 0) {
  stop("No ASVs after chimera removal. Check trimming/filtering parameters or disable strict filters.")
}

# Assign taxonomy
cat("Assigning taxonomy...\n")
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa", multithread=TRUE)
taxa <- addSpecies(taxa, "silva_species_assignment_v138.1.fa")

# Save ASV sequences
uniqueSeqs <- colnames(seqtab.nochim)
names(uniqueSeqs) <- paste0("ASV", seq_along(uniqueSeqs))
writeXStringSet(DNAStringSet(uniqueSeqs), file.path(output_dir, "asv_sequences.fasta"))

# Create phyloseq object with basic components
cat("Creating phyloseq object...\n")
ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows=FALSE),
  tax_table(taxa)
)

# Attach sample metadata if available
if (file.exists("metadata.csv")) {
  meta <- tryCatch(read.csv("metadata.csv", row.names = 1, check.names = FALSE), error = function(e) NULL)
  if (!is.null(meta)) {
    # Keep only samples present in the phyloseq object and order accordingly
    common_samples <- intersect(sample_names(ps), rownames(meta))
    if (length(common_samples) > 0) {
      meta_sub <- meta[common_samples, , drop = FALSE]
      # Reorder rows to match sample order in ps
      meta_sub <- meta_sub[sample_names(ps), , drop = FALSE]
      ps <- merge_phyloseq(ps, sample_data(meta_sub))
    } else {
      warning("metadata.csv found, but sample IDs do not match trimmed FASTQ sample names. Skipping metadata merge.")
    }
  } else {
    warning("Unable to read metadata.csv. Skipping metadata merge.")
  }
}

# Add phylogenetic tree if phangorn is available
if (requireNamespace("phangorn", quietly = TRUE)) {
  cat("Creating phylogenetic tree...\n")
  seqs <- getSequences(seqtab.nochim)
  names(seqs) <- seqs
  alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
  phang.align <- phangorn::phyDat(as(alignment, "matrix"), type="DNA")
  dm <- phangorn::dist.ml(phang.align)
  treeNJ <- phangorn::NJ(dm)
  fit <- phangorn::pml(treeNJ, data=phang.align)
  fitGTR <- phangorn::update(fit, k=4, inv=0.2)
  tree <- phy_tree(fitGTR$tree)
  ps <- merge_phyloseq(ps, tree)
} else {
  message("Package 'phangorn' not available. Skipping phylogenetic tree creation.")
}

# Rename ASVs
taxa_names(ps) <- paste0("ASV", seq_along(taxa_names(ps)))

# Save raw phyloseq object
saveRDS(ps, file.path(output_dir, "phyloseq_object_raw.rds"))

# ==================================================================
# 2. Data Processing and Normalization
# ==================================================================
cat("Processing and normalizing data...\n")

# Remove non-bacterial sequences, chloroplast and mitochondria; require known Phylum
ps_filtered <- subset_taxa(ps,
  Kingdom == "Bacteria" &
  !is.na(Phylum) & Phylum != "" & Phylum != "uncharacterized" &
  (is.na(Order) | Order != "Chloroplast") &
  (is.na(Family) | Family != "Mitochondria")
)

# Save a brief taxa summary after filtering
tax_after <- as.data.frame(tax_table(ps_filtered))
write.csv(as.data.frame(table(tax_after$Phylum)), file.path(output_dir, "phylum_counts_after_filter.csv"), row.names = FALSE)

# Rarefy to even depth
sample_sums_df <- data.frame(Sum = sample_sums(ps_filtered))
rarefaction_depth <- min(sample_sums(ps_filtered))
ps_rarefied <- rarefy_even_depth(ps_filtered, sample.size = rarefaction_depth, rngseed = 123)

saveRDS(ps_rarefied, file.path(output_dir, "phyloseq_rarefied.rds"))

# ==================================================================
# 3. Basic Analysis
# ==================================================================
cat("Performing basic analysis...\n")

alpha_df <- estimate_richness(ps_rarefied, measures = c("Observed", "Shannon", "InvSimpson"))
alpha_df$Sample <- rownames(alpha_df)
if (!is.null(sample_data(ps_rarefied, errorIfNULL = FALSE))) {
  md_alpha <- as.data.frame(sample_data(ps_rarefied))
  md_alpha$Sample <- rownames(md_alpha)
  alpha_df <- merge(alpha_df, md_alpha, by = "Sample", all.x = TRUE)
}
write.csv(alpha_df, file.path(output_dir, "alpha_diversity.csv"), row.names = FALSE)

# Alpha diversity plots
alpha_long <- tidyr::pivot_longer(alpha_df, cols = c("Observed","Shannon","InvSimpson"), names_to = "Metric", values_to = "Value")
alpha_plot <- ggplot(alpha_long, aes(x = Metric, y = Value)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_minimal() +
  labs(x = "Metric", y = "Value", title = "Alpha Diversity Metrics") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(output_dir, "alpha_diversity_plot.pdf"), alpha_plot, width = 8, height = 6)
ggsave(file.path(output_dir, "alpha_diversity_plot.tiff"), alpha_plot, width = 8, height = 6, dpi = 600)

# Alpha diversity stats by Group (if present)
if ("Group" %in% colnames(alpha_df)) {
  kw_obs <- kruskal.test(Observed ~ Group, data = alpha_df)
  kw_sha <- kruskal.test(Shannon ~ Group, data = alpha_df)
  kw_inv <- kruskal.test(InvSimpson ~ Group, data = alpha_df)
  sink(file.path(output_dir, "alpha_diversity_stats.txt"))
  cat("Kruskal–Wallis tests by Group\n\n")
  print(kw_obs); cat("\n"); print(kw_sha); cat("\n"); print(kw_inv); cat("\n")
  sink()
  # Pairwise Wilcoxon with BH correction
  pw_obs <- pairwise.wilcox.test(alpha_df$Observed, alpha_df$Group, p.adjust.method = "BH")
  pw_sha <- pairwise.wilcox.test(alpha_df$Shannon, alpha_df$Group, p.adjust.method = "BH")
  pw_inv <- pairwise.wilcox.test(alpha_df$InvSimpson, alpha_df$Group, p.adjust.method = "BH")
  capture.output(pw_obs, file = file.path(output_dir, "alpha_pairwise_observed.txt"))
  capture.output(pw_sha, file = file.path(output_dir, "alpha_pairwise_shannon.txt"))
  capture.output(pw_inv, file = file.path(output_dir, "alpha_pairwise_invsimpson.txt"))
}

# Beta diversity (Bray-Curtis)
bray_dist <- phyloseq::distance(ps_rarefied, method = "bray")
pcoa <- ordinate(ps_rarefied, method = "PCoA", distance = bray_dist)

# PCoA plot (color by Group if available)
pcoa_df <- data.frame(pcoa$vectors)
pcoa_df$Sample <- rownames(pcoa_df)
if (!is.null(sample_data(ps_rarefied, errorIfNULL = FALSE))) {
  meta_df <- as.data.frame(sample_data(ps_rarefied))
  meta_df$Sample <- rownames(meta_df)
  pcoa_df <- merge(pcoa_df, meta_df, by = "Sample", all.x = TRUE)
}

if ("Group" %in% colnames(pcoa_df)) {
  pcoa_plot <- ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color = Group)) +
    geom_point(size = 3, alpha = 0.9) +
    theme_minimal() +
    labs(x = paste0("PC1 (", round(pcoa$values$Relative_eig[1] * 100, 1), "%)"),
         y = paste0("PC2 (", round(pcoa$values$Relative_eig[2] * 100, 1), "%)"),
         title = "PCoA of Bray-Curtis Distances")
} else {
  pcoa_plot <- ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2)) +
    geom_point(size = 3, alpha = 0.9) +
    theme_minimal() +
    labs(x = paste0("PC1 (", round(pcoa$values$Relative_eig[1] * 100, 1), "%)"),
         y = paste0("PC2 (", round(pcoa$values$Relative_eig[2] * 100, 1), "%)"),
         title = "PCoA of Bray-Curtis Distances")
}

ggsave(file.path(output_dir, "beta_diversity_pcoa.pdf"), pcoa_plot, width = 8, height = 6)
ggsave(file.path(output_dir, "beta_diversity_pcoa.tiff"), pcoa_plot, width = 8, height = 6, dpi = 600)

# PERMANOVA and dispersion if Group present
if (!is.null(sample_data(ps_rarefied, errorIfNULL = FALSE))) {
  md <- as.data.frame(sample_data(ps_rarefied))
  if ("Group" %in% colnames(md)) {
    ad <- vegan::adonis2(as.matrix(bray_dist) ~ Group, data = md, permutations = 999)
    capture.output(ad, file = file.path(output_dir, "permanova_group.txt"))
    # Test for homogeneity of dispersion
    bd <- vegan::betadisper(bray_dist, md$Group)
    bd_perm <- vegan::permutest(bd, permutations = 999)
    capture.output(list(betadisper = bd, permutest = bd_perm), file = file.path(output_dir, "betadisper_group.txt"))
  }
}

# Taxonomic composition at Phylum level
tax_phylum <- tax_glom(ps_rarefied, taxrank = "Phylum")
phylum_counts <- transform_sample_counts(tax_phylum, function(x) x/sum(x))
phylum_data <- psmelt(phylum_counts)

# Top 10 phyla plot (sample-wise)
top_phyla <- phylum_data %>%
  group_by(Phylum) %>%
  summarise(mean_abundance = mean(Abundance)) %>%
  arrange(desc(mean_abundance)) %>%
  head(10)

phylum_plot <- ggplot(subset(phylum_data, Phylum %in% top_phyla$Phylum),
                     aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sample", y = "Relative Abundance", title = "Top 10 Phyla Composition")

ggsave(file.path(output_dir, "phylum_composition.pdf"), phylum_plot, width = 10, height = 6)
ggsave(file.path(output_dir, "phylum_composition.tiff"), phylum_plot, width = 10, height = 6, dpi = 600)

# Group-aggregated composition (if Group present)
if (!is.null(sample_data(ps_rarefied, errorIfNULL = FALSE))) {
  md <- as.data.frame(sample_data(ps_rarefied))
  if ("Group" %in% colnames(md)) {
    ps_grouped <- ps_rarefied %>%
      merge_samples("Group") %>%
      transform_sample_counts(function(x) x / sum(x))
    phylum_grp <- psmelt(tax_glom(ps_grouped, taxrank = "Phylum"))
    phylum_grp_plot <- ggplot(phylum_grp, aes(x = Sample, y = Abundance, fill = Phylum)) +
      geom_bar(stat = "identity") +
      scale_fill_brewer(palette = "Set3") +
      theme_minimal() +
      labs(x = "Group", y = "Relative Abundance", title = "Phylum Composition by Group")
    ggsave(file.path(output_dir, "phylum_composition_by_group.tiff"), phylum_grp_plot, width = 8, height = 6, dpi = 600)
  }
}

# Rarefaction curves from pre-rarefaction counts
otu_mat_full <- as(otu_table(ps_filtered), "matrix")
if (!taxa_are_rows(ps_filtered)) otu_mat_full <- t(otu_mat_full)
pdf(file.path(output_dir, "rarefaction_curves.pdf"), width = 10, height = 8)
suppressWarnings(vegan::rarecurve(otu_mat_full, step = 100, cex = 0.5, label = TRUE))
dev.off()

# Optional: ANCOM-BC differential abundance by Group
if (requireNamespace("ANCOMBC", quietly = TRUE)) {
  if (!is.null(sample_data(ps_filtered, errorIfNULL = FALSE))) {
    md <- as.data.frame(sample_data(ps_filtered))
    if ("Group" %in% colnames(md) && length(unique(md$Group)) > 1) {
      suppressMessages({
        res <- try(ANCOMBC::ancombc2(data = ps_filtered, formula = "Group", p_adj_method = "BH", lib_cut = 0, prv_cut = 0.1, 
                                     group = "Group", struc_zero = TRUE, neg_lb = TRUE), silent = TRUE)
      })
      if (!inherits(res, "try-error")) {
        # Extract results table
        tt <- res$res
        write.csv(tt, file.path(output_dir, "ancombc_results.csv"), row.names = TRUE)
        # Simple volcano plot if lfc/pvals exist
        if (all(c("lfc","`p_val`","`q_val`") %in% colnames(tt[[1]]))) {
          # Skipping complex plotting due to nested results structure differences across versions
        }
      } else {
        message("ANCOM-BC failed; skipping differential abundance.")
      }
    }
  }
}

# ==================================================================
# 4. Advanced Analysis
# ==================================================================
cat("Performing advanced analysis...\n")

# Create advanced analysis directory
advanced_dir <- file.path(output_dir, "advanced_analysis")
dir.create(advanced_dir, recursive = TRUE, showWarnings = FALSE)

# Core microbiome analysis
ps_rel <- transform_sample_counts(ps_rarefied, function(x) x/sum(x))
prevalence_threshold <- 0.5
abundance_threshold <- 0.001

# Calculate prevalence and abundance
otu_mat <- as(otu_table(ps_rel), "matrix")
if (taxa_are_rows(ps_rel)) {
    prevalence <- rowSums(otu_mat > 0) / ncol(otu_mat)
    abundance <- rowMeans(otu_mat)
} else {
    prevalence <- colSums(otu_mat > 0) / nrow(otu_mat)
    abundance <- colMeans(otu_mat)
}

core_taxa <- names(which(prevalence >= prevalence_threshold & abundance >= abundance_threshold))
if (length(core_taxa) > 0) {
    ps_core <- prune_taxa(core_taxa, ps_rel)
} else {
    message("No core taxa found with the given thresholds. Adjusting thresholds...")
    prevalence_threshold <- 0.3  # Lower threshold
    abundance_threshold <- 0.0005  # Lower threshold
    core_taxa <- names(which(prevalence >= prevalence_threshold & abundance >= abundance_threshold))
    if (length(core_taxa) > 0) {
        ps_core <- prune_taxa(core_taxa, ps_rel)
        message("Using adjusted thresholds: prevalence >= 0.3, abundance >= 0.0005")
    } else {
        stop("No core taxa found even with adjusted thresholds.")
    }
}

# Save core microbiome results
tax_mat <- as(tax_table(ps_core), "matrix")
core_tax <- data.frame(
  ASV = rownames(tax_mat),
  as.data.frame(tax_mat, stringsAsFactors = FALSE),
  Prevalence = prevalence[rownames(tax_mat)],
  MeanAbundance = abundance[rownames(tax_mat)],
  stringsAsFactors = FALSE
)
write.csv(core_tax, file.path(advanced_dir, "core_microbiome.csv"))

# Core microbiome prevalence plot
prevalence_plot <- ggplot(data.frame(
    ASV = names(prevalence),
    Prevalence = prevalence,
    Abundance = abundance
  ), aes(x = Abundance, y = Prevalence)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = prevalence_threshold, linetype = 2, color = "red") +
  geom_vline(xintercept = abundance_threshold, linetype = 2, color = "red") +
  scale_x_log10() +
  theme_minimal() +
  labs(
    title = "Core Microbiome Analysis",
    x = "Mean Relative Abundance",
    y = "Prevalence (Proportion of Samples)"
  )

ggsave(file.path(advanced_dir, "core_microbiome_prevalence.pdf"), prevalence_plot, width = 8, height = 6)

# Enhanced heatmap of core taxa
if (nrow(otu_table(ps_core)) > 0 && ncol(otu_table(ps_core)) > 0) {
    core_abundance <- as.matrix(otu_table(ps_core))
    # Make sure we have the correct orientation
    if (!taxa_are_rows(ps_core)) {
        core_abundance <- t(core_abundance)
    }
    tax_mat <- as(tax_table(ps_core), "matrix")
    genera <- tax_mat[, "Genus"]
    genera[is.na(genera)] <- "Unknown genus"
    rownames(core_abundance) <- genera

    pdf(file.path(advanced_dir, "core_taxa_heatmap.pdf"), width = 10, height = 8)
    Heatmap(core_abundance,
            name = "Abundance",
            column_title = "Core Taxa Abundance Across Samples",
            clustering_distance_rows = "euclidean",
            clustering_distance_columns = "euclidean",
            show_row_names = TRUE,
            show_column_names = TRUE)
    dev.off()
} else {
    message("No core taxa found for heatmap generation.")
}

# Additional Taxa Abundance Analysis
if (nsamples(ps_rarefied) > 1) {
  # Get counts of each taxonomic level
  tax_counts <- lapply(c("Phylum", "Class", "Order", "Family", "Genus"), function(level) {
    ps_level <- tax_glom(ps_rarefied, taxrank = level)
    counts <- transform_sample_counts(ps_level, function(x) x/sum(x))
    return(psmelt(counts))
  })
  names(tax_counts) <- c("Phylum", "Class", "Order", "Family", "Genus")
  
  # Save counts at each taxonomic level
  for (level in names(tax_counts)) {
    write.csv(tax_counts[[level]], 
              file.path(advanced_dir, paste0(tolower(level), "_abundance.csv")))
  }
  
  # Create plots for top taxa at different levels
  for (level in c("Phylum", "Class", "Order", "Family", "Genus")) {
    # Get data for current taxonomic level
    level_data <- tax_counts[[level]]
    
    # Calculate mean abundance for each taxa
    top_taxa <- level_data %>%
      group_by(!!sym(level)) %>%
      summarise(mean_abundance = mean(Abundance)) %>%
      arrange(desc(mean_abundance)) %>%
      head(20)  # Top 20 taxa
    
    # Create abundance plot
    taxa_plot <- level_data %>%
      filter(!!sym(level) %in% top_taxa[[level]]) %>%
      ggplot(aes(x = Sample, y = Abundance, fill = !!sym(level))) +
      geom_bar(stat = "identity") +
      scale_fill_brewer(palette = "Set3") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "right",
        legend.title = element_blank()
      ) +
      labs(
        y = "Relative Abundance",
        title = paste("Top 20", level, "Composition")
      )
    
    # Save plot
    ggsave(
      file.path(advanced_dir, paste0("top_", tolower(level), "_abundance.pdf")),
      taxa_plot,
      width = 10,
      height = 6
    )
    
    # Create heatmap data
    taxa_mat <- level_data %>%
      filter(!!sym(level) %in% top_taxa[[level]]) %>%
      reshape2::dcast(as.formula(paste(level, "~ Sample")), 
                     value.var = "Abundance", 
                     fill = 0)
    rownames(taxa_mat) <- taxa_mat[[level]]
    taxa_mat[[level]] <- NULL
    
    # Create and save heatmap
    pdf(file.path(advanced_dir, paste0("top_", tolower(level), "_heatmap.pdf")),
        width = 12, height = 8)
    draw(Heatmap(as.matrix(taxa_mat),
           name = "Relative\nAbundance",
           column_title = paste("Top 20", level, "Abundance"),
           clustering_distance_rows = "euclidean",
           clustering_distance_columns = "euclidean",
           show_row_names = TRUE,
           show_column_names = TRUE,
           row_names_gp = gpar(fontsize = 8),
           column_names_gp = gpar(fontsize = 8)))
    dev.off()
  }
}

# Save session info
writeLines(capture.output(sessionInfo()), 
           file.path(output_dir, "sessionInfo.txt"))

cat("Analysis complete! Results saved in output directory.\n")