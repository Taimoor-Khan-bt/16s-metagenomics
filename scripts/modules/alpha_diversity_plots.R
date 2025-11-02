# Alpha Diversity Plotting Module
# Functions for alpha diversity visualizations

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(dplyr)
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

#' Plot alpha diversity boxplots with statistical tests
#'
#' @param ps Phyloseq object
#' @param cfg Configuration list
#' @param outdir Output directory
#' @param color_by Grouping variable (defaults to cfg$metadata$group_column)
#' @export
plot_alpha_diversity <- function(ps, cfg, outdir, color_by = NULL) {
  tryCatch({
    # Get diversity measures
    requested_measures <- cfg$plots$alpha_measures %||% c("Shannon", "Observed", "Simpson")
    available_measures <- c("Shannon", "Observed", "Simpson", "InvSimpson", "Chao1", "ACE")
    measures <- intersect(requested_measures, available_measures)
    if (length(measures) == 0) measures <- c("Shannon", "Observed")
    
    # Calculate diversity
    df_alpha <- estimate_richness(ps, measures = measures)
    
    # Fix rowname matching (R converts hyphens to periods)
    alpha_ids <- rownames(df_alpha)
    alpha_ids_fixed <- gsub("\\.", "-", alpha_ids)
    
    # Get metadata
    meta <- as(sample_data(ps), "data.frame")
    meta_ids <- rownames(meta)
    
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
    
    # Add grouping column
    df_alpha[[group_col]] <- meta[[group_col]][match(alpha_ids_fixed, meta_ids)]
    
    # Validate grouping variable
    if (!group_col %in% colnames(df_alpha)) {
      message("[alpha] Skipping: column '", group_col, "' not found")
      return(invisible(NULL))
    }
    if (all(is.na(df_alpha[[group_col]]))) {
      message("[alpha] Skipping: all values NA in '", group_col, "'")
      return(invisible(NULL))
    }
    
    # Check for statistical package
    has_ggpubr <- has_package("ggpubr")
    n_groups <- length(unique(df_alpha[[group_col]][!is.na(df_alpha[[group_col]])]))
    
    # Plot each measure
    for (meas in measures) {
      # Transform if needed
      if (meas == "Observed") {
        df_alpha$Observed_log1p <- log1p(df_alpha$Observed)
        y_col <- "Observed_log1p"
        y_label <- "log1p(Observed)"
      } else {
        y_col <- meas
        y_label <- meas
      }
      
      # Create plot
      p <- ggplot(df_alpha, aes(x = .data[[group_col]], y = .data[[y_col]], 
                                 fill = .data[[group_col]])) +
        geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +
        geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
        scale_fill_manual(values = palette_vals(n_groups, cfg)) +
        plot_theme(cfg) +
        labs(
          title = paste(meas, "Diversity"),
          subtitle = paste("Grouped by", group_col, "| N =", nrow(df_alpha)),
          x = group_col, 
          y = y_label
        ) +
        theme(legend.position = "none")
      
      # Add statistics if available
      if (has_ggpubr && n_groups >= 2) {
        test_method <- if (n_groups == 2) "wilcox.test" else "kruskal.test"
        
        if (n_groups <= 4) {
          p <- p + ggpubr::stat_compare_means(
            method = test_method,
            label = "p.signif",
            hide.ns = FALSE
          )
        } else {
          p <- p + ggpubr::stat_compare_means(method = test_method)
        }
        caption_text <- paste0("Test: ", test_method)
      } else {
        caption_text <- "Boxplots show median and IQR"
      }
      
      p <- p + labs(caption = caption_text)
      
      # Save plot
      filename_suffix <- if (is.null(color_by)) "primary" else 
        tolower(gsub("[^[:alnum:]]", "_", group_col))
      save_plot(
        p, 
        file.path(outdir, paste0("alpha_", tolower(meas), "_", filename_suffix, ".tiff")), 
        cfg
      )
    }
    
    message("[alpha] Generated ", length(measures), " alpha diversity plots")
    
  }, error = function(e) {
    warning("[alpha] Failed to generate plots: ", e$message)
  })
}

#' Plot alpha diversity vs continuous variables
#'
#' @param ps Phyloseq object
#' @param cfg Configuration list
#' @param outdir Output directory
#' @param continuous_vars Character vector of continuous variable names
#' @export
plot_alpha_continuous <- function(ps, cfg, outdir, continuous_vars) {
  tryCatch({
    # Calculate diversity
    requested_measures <- cfg$plots$alpha_measures %||% c("Shannon", "Observed")
    available_measures <- c("Shannon", "Observed", "Simpson", "InvSimpson")
    measures <- intersect(requested_measures, available_measures)
    
    df_alpha <- estimate_richness(ps, measures = measures)
    
    # Get metadata
    meta <- as(sample_data(ps), "data.frame")
    
    # Match and merge
    alpha_ids_fixed <- gsub("\\.", "-", rownames(df_alpha))
    meta_ids <- rownames(meta)
    
    for (cont_var in continuous_vars) {
      if (!cont_var %in% colnames(meta)) {
        message("[alpha] Skipping '", cont_var, "': not found in metadata")
        next
      }
      
      df_alpha[[cont_var]] <- meta[[cont_var]][match(alpha_ids_fixed, meta_ids)]
      
      # Skip if not numeric or all NA
      if (!is.numeric(df_alpha[[cont_var]]) || all(is.na(df_alpha[[cont_var]]))) {
        message("[alpha] Skipping '", cont_var, "': not numeric or all NA")
        next
      }
      
      # Plot each measure
      for (meas in measures) {
        p <- ggplot(df_alpha, aes(x = .data[[cont_var]], y = .data[[meas]])) +
          geom_point(alpha = 0.6, size = 3) +
          geom_smooth(method = "lm", se = TRUE, color = "blue", alpha = 0.3) +
          plot_theme(cfg) +
          labs(
            title = paste(meas, "vs", cont_var),
            subtitle = paste("N =", sum(!is.na(df_alpha[[cont_var]]))),
            x = cont_var,
            y = meas
          )
        
        # Add correlation if enough data
        if (sum(!is.na(df_alpha[[meas]]) & !is.na(df_alpha[[cont_var]])) >= 3) {
          cor_test <- cor.test(df_alpha[[meas]], df_alpha[[cont_var]], 
                                method = "spearman")
          p <- p + labs(
            caption = sprintf("Spearman ρ = %.3f, p = %.3g", 
                              cor_test$estimate, cor_test$p.value)
          )
        }
        
        save_plot(
          p,
          file.path(outdir, paste0("alpha_", tolower(meas), "_vs_", 
                                    tolower(cont_var), ".tiff")),
          cfg
        )
      }
    }
    
    message("[alpha] Generated continuous correlation plots")
    
  }, error = function(e) {
    warning("[alpha] Failed to generate continuous plots: ", e$message)
  })
}
