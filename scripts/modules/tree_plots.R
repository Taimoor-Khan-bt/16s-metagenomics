# Phylogenetic Tree Plotting Module
# Functions for building and visualizing phylogenetic trees

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(ggtree)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Build phylogenetic tree from ASV sequences
#'
#' @param ps Phyloseq object
#' @param cfg Configuration list
#' @return Phyloseq object with tree attached
#' @export
build_phylogenetic_tree <- function(ps, cfg) {
  
  message("[tree] Building phylogenetic tree...")
  
  # Check required packages
  if (!requireNamespace("phangorn", quietly = TRUE) || 
      !requireNamespace("DECIPHER", quietly = TRUE)) {
    message("[tree] Missing phangorn or DECIPHER packages; cannot build tree")
    return(ps)
  }
  
  tryCatch({
    # Limit to top taxa to avoid memory issues
    max_tips <- cfg$amplicon$phylogeny$max_tips %||% 500
    keep_taxa <- names(sort(taxa_sums(ps), decreasing = TRUE))[1:min(max_tips, ntaxa(ps))]
    ps_sub <- prune_taxa(keep_taxa, ps)
    
    message("[tree] Building tree for ", ntaxa(ps_sub), " most abundant taxa")
    
    # Look for ASV sequences
    base_out <- cfg$project$output_dir %||% "output"
    cohort <- cfg$io$cohort %||% basename(cfg$io$input_dir)
    fasta_file <- file.path(base_out, cohort, "analysis", "asv_sequences.fasta")
    
    if (!file.exists(fasta_file)) {
      message("[tree] ASV sequences file not found: ", fasta_file)
      return(ps)
    }
    
    # Read sequences
    fas <- Biostrings::readDNAStringSet(fasta_file)
    names(fas) <- gsub("\t.*$", "", names(fas))
    
    # Match to taxa
    taxa_ids <- taxa_names(ps_sub)
    have <- intersect(names(fas), taxa_ids)
    
    if (length(have) == 0) {
      message("[tree] No matching ASV IDs in fasta file")
      return(ps)
    }
    
    if (length(have) < length(taxa_ids)) {
      message("[tree] Only ", length(have), "/", length(taxa_ids), " taxa have sequences")
      ps_sub <- prune_taxa(have, ps_sub)
      taxa_ids <- taxa_names(ps_sub)
    }
    
    # Prepare sequences
    seqs <- as.character(fas[have])
    names(seqs) <- have
    seqs_ds <- Biostrings::DNAStringSet(seqs)
    
    # Align sequences
    message("[tree] Aligning sequences...")
    alignment <- DECIPHER::AlignSeqs(seqs_ds, processors = cfg$project$threads %||% 2)
    
    # Convert to phangorn format
    phang_align <- phangorn::phyDat(as.matrix(alignment), type = "DNA")
    
    # Build neighbor-joining tree
    message("[tree] Building neighbor-joining tree...")
    dm <- phangorn::dist.ml(phang_align)
    treeNJ <- phangorn::NJ(dm)
    
    # Optimize tree with GTR model
    message("[tree] Optimizing with GTR+G+I model...")
    fit <- phangorn::pml(treeNJ, data = phang_align)
    fitGTR <- phangorn::optim.pml(
      fit, 
      model = "GTR", 
      optInv = TRUE, 
      optGamma = TRUE,
      rearrangement = "stochastic", 
      control = phangorn::pml.control(trace = 0)
    )
    
    # Attach tree to phyloseq object
    phy_tree(ps_sub) <- fitGTR$tree
    
    message("[tree] Phylogenetic tree built successfully")
    message("[tree] Tree has ", length(fitGTR$tree$tip.label), " tips")
    
    return(ps_sub)
    
  }, error = function(e) {
    message("[tree] Error building tree: ", e$message)
    return(ps)
  })
}

#' Plot rectangular phylogenetic tree
#'
#' @param ps Phyloseq object with tree
#' @param cfg Configuration list
#' @param outdir Output directory
#' @export
plot_tree_rectangular <- function(ps, cfg, outdir) {
  
  message("[tree] Creating rectangular tree plot...")
  
  tryCatch({
    # Check if tree exists
    tree <- phy_tree(ps, errorIfNULL = FALSE)
    if (is.null(tree)) {
      message("[tree] No phylogenetic tree found")
      return(invisible(NULL))
    }
    
    if (length(tree$tip.label) < 2) {
      message("[tree] Tree too small (< 2 tips)")
      return(invisible(NULL))
    }
    
    # Prepare data
    taxa_abund <- taxa_sums(ps)
    tt <- as(tax_table(ps), "matrix")
    
    df <- data.frame(
      label = names(taxa_abund),
      Abundance = taxa_abund,
      Phylum = tt[names(taxa_abund), "Phylum"],
      stringsAsFactors = FALSE
    )
    
    df$Phylum[is.na(df$Phylum)] <- "Unassigned"
    df$AbundanceScaled <- scales::rescale(log1p(df$Abundance), to = c(1, 5))
    df <- df[df$label %in% tree$tip.label, ]
    
    n_phyla <- length(unique(df$Phylum))
    subtitle_text <- paste(nrow(df), "ASVs from", n_phyla, "phyla; sized by log abundance")
    
    # Get color palette
    colors <- get_tree_palette(n_phyla, cfg)
    
    # Create tree plot
    p <- ggtree(tree, layout = "rectangular", linewidth = 0.4) %<+% df +
      geom_tippoint(aes(size = AbundanceScaled, color = Phylum), alpha = 0.8) +
      scale_color_manual(values = colors) +
      scale_size_continuous(range = c(1, 4)) +
      guides(
        size = "none", 
        color = guide_legend(override.aes = list(size = 3))
      ) +
      theme_tree2() +
      theme(
        legend.position = "right",
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        plot.caption = element_text(size = 10, hjust = 1)
      ) +
      labs(
        title = "Phylogenetic Tree (Rectangular Layout)",
        subtitle = subtitle_text,
        caption = "Maximum likelihood tree (GTR+G+I model); tips colored by Phylum"
      )
    
    # Save
    filename <- file.path(outdir, "phylo_tree_rectangular.tiff")
    ggsave(filename, p, width = 12, height = 10, dpi = 600, compression = "lzw")
    message("[tree] Saved rectangular tree: ", basename(filename))
    
    # Also save PDF
    filename_pdf <- file.path(outdir, "phylo_tree_rectangular.pdf")
    ggsave(filename_pdf, p, width = 12, height = 10)
    
    return(p)
    
  }, error = function(e) {
    message("[tree] Error creating rectangular tree plot: ", e$message)
  })
}

#' Plot circular phylogenetic tree
#'
#' @param ps Phyloseq object with tree
#' @param cfg Configuration list
#' @param outdir Output directory
#' @export
plot_tree_circular <- function(ps, cfg, outdir) {
  
  message("[tree] Creating circular tree plot...")
  
  tryCatch({
    # Check if tree exists
    tree <- phy_tree(ps, errorIfNULL = FALSE)
    if (is.null(tree)) {
      message("[tree] No phylogenetic tree found")
      return(invisible(NULL))
    }
    
    if (length(tree$tip.label) < 2) {
      message("[tree] Tree too small (< 2 tips)")
      return(invisible(NULL))
    }
    
    # Prepare data
    taxa_abund <- taxa_sums(ps)
    tt <- as(tax_table(ps), "matrix")
    
    df <- data.frame(
      label = names(taxa_abund),
      Abundance = taxa_abund,
      Phylum = tt[names(taxa_abund), "Phylum"],
      stringsAsFactors = FALSE
    )
    
    df$Phylum[is.na(df$Phylum)] <- "Unassigned"
    df$AbundanceScaled <- scales::rescale(log1p(df$Abundance), to = c(1, 5))
    df <- df[df$label %in% tree$tip.label, ]
    
    n_phyla <- length(unique(df$Phylum))
    subtitle_text <- paste(nrow(df), "ASVs from", n_phyla, "phyla; sized by log abundance")
    
    # Get color palette
    colors <- get_tree_palette(n_phyla, cfg)
    
    # Create tree plot
    p <- ggtree(tree, layout = "circular", linewidth = 0.3) %<+% df +
      geom_tippoint(aes(size = AbundanceScaled, color = Phylum), alpha = 0.75) +
      scale_color_manual(values = colors) +
      scale_size_continuous(range = c(1, 4)) +
      guides(
        size = "none", 
        color = guide_legend(override.aes = list(size = 3))
      ) +
      theme_tree() +
      theme(
        legend.position = "right",
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        plot.caption = element_text(size = 10, hjust = 1)
      ) +
      labs(
        title = "Phylogenetic Tree (Circular Layout)",
        subtitle = subtitle_text,
        caption = "Maximum likelihood tree (GTR+G+I model); tips colored by Phylum"
      )
    
    # Save
    filename <- file.path(outdir, "phylo_tree_circular.tiff")
    ggsave(filename, p, width = 10, height = 10, dpi = 600, compression = "lzw")
    message("[tree] Saved circular tree: ", basename(filename))
    
    # Also save PDF
    filename_pdf <- file.path(outdir, "phylo_tree_circular.pdf")
    ggsave(filename_pdf, p, width = 10, height = 10)
    
    return(p)
    
  }, error = function(e) {
    message("[tree] Error creating circular tree plot: ", e$message)
  })
}

#' Get color palette for tree plots
#' @keywords internal
get_tree_palette <- function(n, cfg) {
  pal <- tolower(cfg$visualization$color_palette %||% "okabe-ito")
  
  if (pal == "okabe-ito") {
    cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
              "#0072B2", "#D55E00", "#CC79A7", "#999999")
    return(rep(cols, length.out = n))
  } else if (pal == "viridis") {
    return(viridisLite::viridis(n))
  } else if (pal == "tableau10") {
    cols <- c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
              "#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC")
    return(rep(cols, length.out = n))
  } else {
    return(scales::hue_pal()(n))
  }
}
