#!/usr/bin/env Rscript
# =============================================================================
# core_microbiome.R — Core microbiome analysis (Multi-Group Support)
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: core_microbiome.R <table_tsv> <taxonomy_tsv> <metadata> <group_col> <prevalence> <out_dir> [tax_level]")
}

table_tsv  <- args[1]
tax_tsv    <- args[2]
meta_file  <- args[3]
group_col  <- args[4]
prevalence <- as.numeric(args[5])
out_dir    <- args[6]
tax_level  <- if (length(args) >= 7 && nchar(trimws(args[7])) > 0) args[7] else "ASV"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── 1. Load Data ─────────────────────────────────────────────────────────────
tab <- read.table(table_tsv, header = TRUE, sep = "\t", skip = 1, 
                  row.names = 1, check.names = FALSE, comment.char = "", quote = "")

tax_df <- read.table(tax_tsv, header = TRUE, sep = "\t", 
                     row.names = 1, check.names = FALSE, quote = "")

meta <- read.table(meta_file, header = TRUE, sep = "\t", 
                   row.names = 1, check.names = FALSE, comment.char = "", quote = "")

# ── 2. Sync Samples ──────────────────────────────────────────────────────────
common_samples <- intersect(colnames(tab), rownames(meta))
if (length(common_samples) < 2) stop("Fewer than 2 common samples.")

tab  <- tab[, common_samples, drop = FALSE]
meta <- meta[common_samples, , drop = FALSE]
grp  <- as.factor(meta[[group_col]])
grp_levels <- levels(grp)

# ── 3. Identify Core Taxa per Group ──────────────────────────────────────────
core_ids <- list()
for (g in grp_levels) {
  samp_g    <- common_samples[grp == g]
  tab_g     <- tab[, samp_g, drop = FALSE]
  prev_g    <- rowSums(tab_g > 0) / ncol(tab_g)
  core_ids[[g]] <- names(prev_g[prev_g >= prevalence])
  message("Group ", g, ": ", length(core_ids[[g]]), " core features.")
}

all_core_ids <- unique(unlist(core_ids))

# ── 4. Helper: Map Taxonomy ──────────────────────────────────────────────────
get_tax_name <- function(id) {
  if (id %in% rownames(tax_df)) {
    return(gsub("[a-z]__", "", tax_df[id, "Taxon"]))
  }
  return(id)
}

# ── 5. Multi-Group Statistical Test (Omnibus) ────────────────────────────────
if (length(all_core_ids) > 0) {
  stats_rows <- lapply(all_core_ids, function(tx_id) {
    # Create 2xN Contingency Table (Present/Absent vs Groups)
    count_matrix <- matrix(NA, nrow = 2, ncol = length(grp_levels))
    for (i in seq_along(grp_levels)) {
      g <- grp_levels[i]
      count_matrix[1, i] <- sum(tab[tx_id, grp == g] > 0)  # Present
      count_matrix[2, i] <- sum(tab[tx_id, grp == g] == 0) # Absent
    }
    
    # Choose test based on number of groups and expected cell counts:
    #   2 groups  → Fisher's Exact (exact p-value for 2×2 tables)
    #   3+ groups with any expected count < 5 → Fisher's with Monte Carlo simulation
    #   3+ groups, all expected counts ≥ 5     → Chi-square approximation is valid
    test_used <- "chi_square"
    p_val <- tryCatch({
      if (length(grp_levels) == 2) {
        test_used <<- "fisher_exact"
        fisher.test(count_matrix)$p.value
      } else {
        expected <- suppressWarnings(chisq.test(count_matrix)$expected)
        if (any(expected < 5, na.rm = TRUE)) {
          test_used <<- "fisher_exact_simulated"
          fisher.test(count_matrix, simulate.p.value = TRUE, B = 2000)$p.value
        } else {
          chisq.test(count_matrix)$p.value
        }
      }
    }, error = function(e) { message("Test failed for ", tx_id, ": ", e$message); NA })

    res <- data.frame(FeatureID = tx_id, Taxon = get_tax_name(tx_id),
                      test_used = test_used, p_value = p_val)
    # Add Boolean flags for each group
    for (g in grp_levels) { res[[paste0("is_core_", g)]] <- tx_id %in% core_ids[[g]] }
    return(res)
  })
  
  stats_df <- do.call(rbind, stats_rows)
  stats_df$p_adj     <- p.adjust(stats_df$p_value, method = "BH")
  stats_df$tax_level <- tax_level
  write.table(stats_df, file.path(out_dir, "core_stats.tsv"), 
              sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  # Write an empty file to satisfy Snakemake if no core taxa exist
  write.table(data.frame(Note="No core taxa found"), 
              file.path(out_dir, "core_stats.tsv"), sep="\t")
}

# ── 6. Plotting ─────────────────────────────────────────────────────────────

# Helper: extract readable short name from SILVA taxonomy string.
# Returns the deepest non-empty rank; if the result would be duplicated across
# the table it appends the genus + a numeric suffix to keep rows distinct.
short_tax <- function(x) {
  parts <- strsplit(x, ";")[[1]]
  parts <- trimws(gsub("[a-z]__", "", parts))
  parts <- parts[nchar(parts) > 0]
  if (length(parts) == 0) return(x)
  parts[length(parts)]
}

make_unique_labels <- function(ids, taxa) {
  labels <- sapply(taxa, short_tax)
  # For duplicated short names, prepend genus and append an index
  dupes <- names(which(table(labels) > 1))
  if (length(dupes) > 0) {
    for (d in dupes) {
      idx <- which(labels == d)
      labels[idx] <- paste0(d, " ASV", seq_along(idx))
    }
  }
  labels
}

pdf(file.path(out_dir, "core_plots.pdf"), width = 13, height = 8)

# ── Page 1: Prevalence barplot + significance annotation ─────────────────────
if (exists("stats_df") && nrow(stats_df) > 0) {
  # Build per-group prevalence data
  prev_rows <- lapply(all_core_ids, function(tx_id) {
    do.call(rbind, lapply(grp_levels, function(g) {
      samp_g <- common_samples[grp == g]
      prev   <- sum(tab[tx_id, samp_g] > 0) / length(samp_g)
      data.frame(FeatureID = tx_id, Group = g, Prevalence = prev,
                 stringsAsFactors = FALSE)
    }))
  })
  prev_df <- do.call(rbind, prev_rows)

  # Attach p_adj and unique short names (disambiguates same-genus ASVs)
  unique_labels <- make_unique_labels(stats_df$FeatureID, stats_df$Taxon)
  id_to_label   <- setNames(unique_labels, stats_df$FeatureID)

  prev_df <- merge(prev_df,
                   stats_df[, c("FeatureID", "Taxon", "test_used", "p_value", "p_adj")],
                   by = "FeatureID")
  prev_df$ShortName <- id_to_label[prev_df$FeatureID]
  prev_df$Significant <- ifelse(!is.na(prev_df$p_adj) & prev_df$p_adj < 0.05,
                                "p_adj < 0.05", "p_adj ≥ 0.05")

  # Build label: show p_adj once per taxon (on the group with highest prevalence)
  label_df <- prev_df[ave(prev_df$Prevalence, prev_df$FeatureID, FUN = max) ==
                        prev_df$Prevalence, ]
  label_df <- label_df[!duplicated(label_df$FeatureID), ]
  label_df$label <- paste0("p_adj=", signif(label_df$p_adj, 2))

  # Cap at top 30 taxa by max prevalence across groups for readability
  top_ids <- names(sort(tapply(prev_df$Prevalence, prev_df$FeatureID, max),
                        decreasing = TRUE))[seq_len(min(30, length(all_core_ids)))]
  plot_df  <- prev_df[prev_df$FeatureID %in% top_ids, ]
  label_df <- label_df[label_df$FeatureID %in% top_ids, ]

  # Order taxa by descending max prevalence
  tax_order <- names(sort(tapply(plot_df$Prevalence, plot_df$ShortName, max),
                          decreasing = TRUE))
  plot_df$ShortName  <- factor(plot_df$ShortName,  levels = rev(tax_order))
  label_df$ShortName <- factor(label_df$ShortName, levels = rev(tax_order))

  p1 <- ggplot(plot_df, aes(x = Prevalence, y = ShortName, fill = Group)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.65, alpha = 0.85) +
    # One label per taxon, pinned at x = 1.01 (right margin), no dodge
    geom_text(data = label_df,
              aes(x = 1.01, y = ShortName, label = label, color = Significant),
              inherit.aes = FALSE, hjust = 0, size = 3) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40") +
    scale_x_continuous(labels = scales::percent_format(), limits = c(0, 1.30),
                       expand = c(0, 0)) +
    scale_color_manual(values = c("p_adj < 0.05" = "#C62828", "p_adj \u2265 0.05" = "grey50"),
                       name = "Significance") +
    scale_fill_brewer(palette = "Set2", name = group_col) +
    theme_bw(base_size = 11) +
    theme(legend.position   = "right",
          axis.text.y       = element_text(size = 9),
          panel.grid.major.y = element_blank()) +
    labs(x     = paste0("Prevalence (dashed = ", prevalence * 100, "% threshold)"),
         y     = "Taxon",
         title = paste("Core Microbiome Prevalence by", group_col),
         subtitle = paste0("Taxa present in \u2265", prevalence * 100,
                           "% of samples in at least one group  |  ",
                           "Taxonomy level: ", tax_level, "  |  ",
                           sum(!is.na(stats_df$p_adj) & stats_df$p_adj < 0.05),
                           " taxa with p_adj < 0.05"))
  print(p1)
}

# ── Page 2: UpSetR intersection diagram ──────────────────────────────────────
if (requireNamespace("UpSetR", quietly = TRUE) && length(grp_levels) > 1) {
  upset_list <- lapply(core_ids, function(ids) if (length(ids) > 0) ids else character(0))
  if (sum(sapply(upset_list, length)) > 0) {
    print(UpSetR::upset(UpSetR::fromList(upset_list), order.by = "freq",
                        main.bar.color = "steelblue", sets.bar.color = "darkgreen",
                        text.scale = 1.4))
  }
}

dev.off()
if (exists("p1")) {
  ggsave(file.path(out_dir, "core_plots.png"),
         plot = p1, width = 13, height = 8, units = "in", dpi = 600)
  message("Saved: core_plots.png (600 DPI)")
}
message("Saved: core_plots.pdf")