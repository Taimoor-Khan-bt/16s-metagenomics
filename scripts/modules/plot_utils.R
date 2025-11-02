# Plotting Utilities Module
# Common functions for plot theming, saving, and data preparation

suppressPackageStartupMessages({
  library(ggplot2)
  library(phyloseq)
  library(viridisLite)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Get publication-ready ggplot2 theme
#'
#' @param cfg Configuration list with plots settings
#' @return ggplot2 theme object
#' @export
plot_theme <- function(cfg) {
  theme_sel <- tolower(cfg$plots$theme %||% "classic")
  base_size <- cfg$plots$base_size %||% 12
  font_family <- cfg$plots$font_family %||% "sans"
  
  base <- switch(
    theme_sel,
    minimal = theme_minimal(base_size = base_size, base_family = font_family),
    bw = theme_bw(base_size = base_size, base_family = font_family),
    classic = theme_classic(base_size = base_size, base_family = font_family),
    theme_classic(base_size = base_size, base_family = font_family)
  )
  
  base + theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold", size = base_size * 1.1),
    axis.text = element_text(size = base_size * 0.9),
    legend.title = element_text(face = "bold", size = base_size),
    legend.text = element_text(size = base_size * 0.8),
    plot.title = element_text(face = "bold", size = base_size * 1.2, hjust = 0.5),
    plot.subtitle = element_text(size = base_size * 0.9, hjust = 0.5),
    plot.caption = element_text(size = base_size * 0.7, hjust = 1)
  )
}

#' Get color palette for plots
#'
#' @param n Number of colors needed
#' @param cfg Configuration list with plots settings
#' @return Character vector of colors
#' @export
palette_vals <- function(n, cfg) {
  pal <- tolower(cfg$plots$color_palette %||% "viridis")
  
  if (pal == "okabe-ito") {
    cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
              "#0072B2", "#D55E00", "#CC79A7", "#999999")
    return(rep(cols, length.out = n))
  } else if (pal == "tableau10") {
    cols <- c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
              "#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC")
    return(rep(cols, length.out = n))
  } else if (pal == "viridis") {
    return(viridisLite::viridis(n))
  } else if (pal == "plasma") {
    return(viridisLite::plasma(n))
  } else {
    warning("[viz] Unknown palette '", pal, "'; falling back to viridis")
    return(viridisLite::viridis(n))
  }
}

#' Save plot with consistent settings
#'
#' @param p ggplot object
#' @param path Output file path
#' @param cfg Configuration list
#' @param width Plot width in inches (overrides cfg)
#' @param height Plot height in inches (overrides cfg)
#' @param format Output format: "tiff", "pdf", "png"
#' @export
save_plot <- function(p, path, cfg, width = NULL, height = NULL, format = "tiff") {
  w <- width %||% cfg$plots$width %||% 8
  h <- height %||% cfg$plots$height %||% 6
  dpi <- cfg$plots$dpi %||% 600
  
  # Create directory if needed
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  
  if (format == "tiff") {
    ggsave(filename = path, plot = p, width = w, height = h, 
           dpi = dpi, compression = "lzw")
  } else if (format == "pdf") {
    ggsave(filename = path, plot = p, width = w, height = h, device = "pdf")
  } else {
    ggsave(filename = path, plot = p, width = w, height = h, dpi = dpi)
  }
  
  message("[viz] Saved: ", basename(path))
}

#' Validate and prepare phyloseq metadata
#'
#' @param ps Phyloseq object
#' @param cfg Configuration list
#' @return Phyloseq object with validated metadata
#' @export
prepare_metadata <- function(ps, cfg) {
  meta <- tryCatch(
    as(sample_data(ps), "data.frame"), 
    error = function(e) NULL
  )
  
  if (is.null(meta) || ncol(meta) == 0) {
    meta <- data.frame(
      SampleID = sample_names(ps), 
      Group = "Unknown", 
      row.names = sample_names(ps), 
      check.names = FALSE
    )
  }
  
  id_col <- cfg$metadata$id_column %||% "SampleID"
  group_col <- cfg$metadata$group_column %||% "Group"
  
  # Add missing required columns
  missing_cols <- setdiff(c(id_col, group_col), colnames(meta))
  for (mc in missing_cols) {
    meta[[mc]] <- if (mc == group_col) "Unknown" else rownames(meta)
  }
  
  sample_data(ps) <- sample_data(meta)
  return(ps)
}

#' Get top N taxa by abundance
#'
#' @param ps Phyloseq object
#' @param n Number of top taxa to return
#' @param rank Taxonomic rank
#' @return Character vector of top taxa names
#' @export
get_top_taxa <- function(ps, n = 20, rank = "Genus") {
  # Agglomerate to rank
  ps_agg <- tryCatch(
    tax_glom(ps, taxrank = rank, NArm = FALSE),
    error = function(e) {
      warning("[viz] Could not agglomerate to ", rank, ": ", e$message)
      return(ps)
    }
  )
  
  # Get abundance
  abund <- taxa_sums(ps_agg)
  top_taxa <- names(sort(abund, decreasing = TRUE)[1:min(n, length(abund))])
  
  return(top_taxa)
}

#' Check if statistical package is available
#'
#' @param package Package name
#' @return Logical indicating availability
#' @export
has_package <- function(package) {
  requireNamespace(package, quietly = TRUE)
}
