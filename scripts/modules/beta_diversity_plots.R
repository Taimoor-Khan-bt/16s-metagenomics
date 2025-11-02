# Beta Diversity Plotting Module
# Functions for ordination and beta diversity visualizations

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(vegan)
})

# Source utilities - handle both direct execution and sourcing from other scripts
if (!exists("plot_theme")) {
  script_dir <- tryCatch(
    dirname(sys.frame(1)$ofile),
    error = function(e) "scripts/modules"
  )
  if (file.exists(file.path(script_dir, "plot_utils.R"))) {
    source(file.path(script_dir, "plot_utils.R"))
  } else if (file.exists("scripts/modules/plot_utils.R")) {
    source("scripts/modules/plot_utils.R")
  }
}

#' Plot ordination (PCoA, NMDS, etc.)
#'
#' @param ps Phyloseq object
#' @param cfg Configuration list
#' @param outdir Output directory
#' @param color_by Grouping variable for colors
#' @param shape_by Variable for point shapes (optional)
#' @export
plot_ordination <- function(ps, cfg, outdir, color_by = NULL, shape_by = NULL) {
  tryCatch({
    # Determine grouping variable
    if (is.null(color_by)) {
      if (!is.null(cfg$metadata$primary_comparison)) {
        group_col <- cfg$metadata$primary_comparison$group_column %||% "Group"
      } else {
        group_col <- cfg$metadata$group_column %||% "Group"
      }
    } else {
      group_col <- color_by
    }
    
    # Get metadata
    meta <- as(sample_data(ps), "data.frame")
    if (!group_col %in% colnames(meta)) {
      message("[beta] Skipping: '", group_col, "' not in metadata")
      return(invisible(NULL))
    }
    
    # Distance methods to try
    methods <- cfg$plots$ordination_methods %||% c("PCoA", "NMDS")
    distances <- cfg$plots$beta_distances %||% c("bray", "jaccard")
    
    for (dist_method in distances) {
      for (ord_method in methods) {
        tryCatch({
          # Calculate ordination
          if (tolower(ord_method) == "pcoa") {
            ord <- ordinate(ps, method = "PCoA", distance = dist_method)
            method_label <- "PCoA"
          } else if (tolower(ord_method) == "nmds") {
            ord <- ordinate(ps, method = "NMDS", distance = dist_method)
            method_label <- "NMDS"
          } else {
            next
          }
          
          # Create plot
          p <- plot_ordination(ps, ord, color = group_col, 
                              shape = shape_by, type = "samples") +
            scale_color_manual(values = palette_vals(
              length(unique(meta[[group_col]])), cfg
            )) +
            plot_theme(cfg) +
            labs(
              title = paste(method_label, "Ordination"),
              subtitle = paste("Distance:", dist_method, "| Colored by:", group_col),
              caption = paste("Each point represents one sample")
            ) +
            theme(legend.position = "right")
          
          # Add ellipses if multiple groups
          n_groups <- length(unique(meta[[group_col]][!is.na(meta[[group_col]])]))
          if (n_groups >= 2 && n_groups <= 10) {
            p <- p + stat_ellipse(aes(group = .data[[group_col]]), 
                                 type = "norm", linetype = 2, alpha = 0.5)
          }
          
          # Save
          filename <- paste0("ordination_", tolower(ord_method), "_", 
                            dist_method, "_", 
                            tolower(gsub("[^[:alnum:]]", "_", group_col)), 
                            ".tiff")
          save_plot(p, file.path(outdir, filename), cfg)
          
        }, error = function(e) {
          message("[beta] Failed ", ord_method, " with ", dist_method, ": ", e$message)
        })
      }
    }
    
    message("[beta] Generated ordination plots")
    
  }, error = function(e) {
    warning("[beta] Failed to generate ordination plots: ", e$message)
  })
}

#' Check PERMANOVA assumptions (homogeneity of dispersions)
#'
#' @param ps Phyloseq object
#' @param group_column Column name for grouping variable
#' @param distance Distance method
#' @return List with betadisper result and interpretation
#' @export
check_permanova_assumptions <- function(ps, group_column, distance = "bray") {
  tryCatch({
    # Get distance matrix
    dist_mat <- phyloseq::distance(ps, method = distance)
    
    # Get metadata
    meta <- as(sample_data(ps), "data.frame")
    
    # Check if group column exists
    if (!group_column %in% colnames(meta)) {
      warning("[beta] Group column '", group_column, "' not found in metadata")
      return(NULL)
    }
    
    # Get grouping factor
    group_factor <- factor(meta[[group_column]])
    
    # Check for sufficient groups
    n_groups <- length(levels(group_factor))
    if (n_groups < 2) {
      warning("[beta] Need at least 2 groups for betadisper, found: ", n_groups)
      return(NULL)
    }
    
    # Run betadisper to test homogeneity of dispersions
    cat("\n")
    message("[beta] Checking PERMANOVA assumptions (betadisper)...")
    message("[beta] Testing homogeneity of multivariate dispersions")
    message("[beta] Distance: ", distance, " | Grouping: ", group_column)
    
    bd <- betadisper(dist_mat, group_factor)
    bd_test <- permutest(bd, permutations = 999)
    
    # Interpret results
    p_value <- bd_test$tab$"Pr(>F)"[1]
    assumption_met <- p_value > 0.05
    
    cat("\n")
    cat("══════════════════════════════════════════════════════════\n")
    cat("PERMANOVA Assumption Check: Homogeneity of Dispersions\n")
    cat("══════════════════════════════════════════════════════════\n")
    cat("Distance method:", distance, "\n")
    cat("Grouping variable:", group_column, "\n")
    cat("Groups tested:", paste(levels(group_factor), collapse = ", "), "\n")
    cat("\nBetadisper Test:\n")
    print(bd_test)
    cat("\n")
    
    if (assumption_met) {
      cat("✓ ASSUMPTION MET: No significant difference in dispersions\n")
      cat("  (p = ", sprintf("%.4f", p_value), " > 0.05)\n", sep = "")
      cat("  → PERMANOVA is appropriate for these data\n")
      cat("  → Group differences reflect location, not dispersion\n")
    } else {
      cat("⚠ WARNING: ASSUMPTION VIOLATED\n")
      cat("  (p = ", sprintf("%.4f", p_value), " ≤ 0.05)\n", sep = "")
      cat("  → Dispersions differ significantly among groups\n")
      cat("  → PERMANOVA results may be unreliable\n")
      cat("  → Consider:\n")
      cat("     • Interpreting PERMANOVA with caution\n")
      cat("     • Using alternative tests (ANOSIM, MRPP)\n")
      cat("     • Transforming data to stabilize variance\n")
      cat("     • Reporting both location and dispersion effects\n")
    }
    cat("══════════════════════════════════════════════════════════\n\n")
    
    # Average distances to centroid by group
    cat("Average distance to centroid by group:\n")
    print(bd$group.distances)
    cat("\n")
    
    return(list(
      betadisper = bd,
      permutest = bd_test,
      p_value = p_value,
      assumption_met = assumption_met,
      distance = distance,
      group_column = group_column
    ))
    
  }, error = function(e) {
    warning("[beta] Assumption checking failed: ", e$message)
    return(NULL)
  })
}

#' Perform and visualize PERMANOVA test with assumption checking
#'
#' @param ps Phyloseq object
#' @param cfg Configuration list
#' @param formula_str Formula for PERMANOVA (e.g., "~ Group")
#' @param distance Distance method
#' @param check_assumptions Logical, whether to check assumptions first (default: TRUE)
#' @return List with PERMANOVA result and assumption check
#' @export
run_permanova <- function(ps, cfg, formula_str, distance = "bray", check_assumptions = TRUE) {
  tryCatch({
    # Extract group variable from formula
    group_var <- gsub("^.*~\\s*", "", formula_str)
    group_var <- trimws(strsplit(group_var, "\\+")[[1]][1])
    
    # Check assumptions if requested
    assumption_check <- NULL
    if (check_assumptions) {
      assumption_check <- check_permanova_assumptions(ps, group_var, distance)
      
      if (!is.null(assumption_check) && !assumption_check$assumption_met) {
        message("[beta] ⚠ Proceeding with PERMANOVA despite violated assumptions")
        message("[beta]   Interpret results with caution!")
      }
    }
    
    # Get distance matrix
    dist_mat <- phyloseq::distance(ps, method = distance)
    
    # Get metadata
    meta <- as(sample_data(ps), "data.frame")
    
    # Run PERMANOVA
    formula_obj <- as.formula(formula_str)
    perm <- adonis2(formula_obj, data = meta, permutations = 999)
    
    cat("\n")
    cat("══════════════════════════════════════════════════════════\n")
    cat("PERMANOVA Results\n")
    cat("══════════════════════════════════════════════════════════\n")
    message("[beta] Distance: ", distance, " | Formula: ", formula_str)
    print(perm)
    cat("══════════════════════════════════════════════════════════\n\n")
    
    return(list(
      permanova = perm,
      assumptions = assumption_check,
      distance = distance,
      formula = formula_str
    ))
    
  }, error = function(e) {
    warning("[beta] PERMANOVA failed: ", e$message)
    return(NULL)
  })
}

#' Plot betadisper results
#'
#' @param assumption_check Result from check_permanova_assumptions()
#' @param outdir Output directory
#' @param cfg Configuration list
#' @export
plot_betadisper <- function(assumption_check, outdir, cfg) {
  if (is.null(assumption_check)) {
    return(invisible(NULL))
  }
  
  tryCatch({
    bd <- assumption_check$betadisper
    
    # Create boxplot of distances to centroid
    distances_df <- data.frame(
      Group = bd$group,
      Distance_to_Centroid = bd$distances
    )
    
    p <- ggplot(distances_df, aes(x = Group, y = Distance_to_Centroid, fill = Group)) +
      geom_boxplot(alpha = 0.7, outlier.shape = 16) +
      geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
      scale_fill_manual(values = palette_vals(nlevels(bd$group), cfg)) +
      plot_theme(cfg) +
      labs(
        title = "Homogeneity of Dispersions Test (Betadisper)",
        subtitle = paste0("Distance: ", assumption_check$distance, 
                         " | p = ", sprintf("%.4f", assumption_check$p_value)),
        x = assumption_check$group_column,
        y = "Distance to Group Centroid",
        caption = if (assumption_check$assumption_met) {
          "✓ Assumption met: dispersions are homogeneous (p > 0.05)"
        } else {
          "⚠ Assumption violated: dispersions differ significantly (p ≤ 0.05)"
        }
      ) +
      theme(
        legend.position = "none",
        plot.caption = element_text(hjust = 0, face = "italic")
      )
    
    # Save plot
    filename <- paste0("betadisper_", assumption_check$distance, "_",
                      tolower(gsub("[^[:alnum:]]", "_", assumption_check$group_column)),
                      ".tiff")
    save_plot(p, file.path(outdir, filename), cfg)
    
    message("[beta] Saved betadisper plot: ", filename)
    
  }, error = function(e) {
    warning("[beta] Failed to plot betadisper: ", e$message)
  })
}
