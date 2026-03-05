#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("Usage: lefse_analysis.R <table_tsv> <taxonomy_tsv> <metadata> <group_col> <out_dir> [strategy] [dual_plots]")

table_tsv  <- args[1]
tax_tsv    <- args[2]
meta_file  <- args[3]
group_col  <- args[4]
out_dir    <- args[5]
strategy   <- if (length(args) >= 6) args[6] else "rename"
dual_plots <- isTRUE(tolower(if (length(args) >= 7) args[7] else "false") == "true")

# =============================================================================
# format_taxon_label()
# Converts raw SILVA taxonomy strings to publication-ready labels.
#   In:  "k__Bacteria;p__Firmicutes;c__Clostridia;f__Lachnospiraceae;g__"
#   Out: "Unclassified Lachnospiraceae"
# Logic:
#   1. Normalise separator, split, strip rank prefixes (d__/p__/...).
#   2. Drop empty / 'uncultured' tokens.
#   3. If the LAST raw token was blank (trailing __), wrap last named token
#      in "Unclassified [...]".
#   4. Otherwise return the last non-empty, cleaned token.
# =============================================================================
format_taxon_label <- function(x) {
  if (is.na(x) || nchar(trimws(x)) == 0) return("Unclassified")
  parts <- strsplit(gsub("\\|", ";", x), ";")[[1]]
  parts <- trimws(parts)
  clean <- sub("^[dpcofgskDPCOFGSK]__", "", parts)
  non_empty <- clean[nchar(clean) > 0 & tolower(clean) != "uncultured"]
  if (length(non_empty) == 0) return("Unclassified")
  last_raw_clean <- sub("^[dpcofgskDPCOFGSK]__", "", parts[length(parts)])
  if (nchar(trimws(last_raw_clean)) == 0) {
    return(paste("Unclassified", non_empty[length(non_empty)]))
  }
  return(non_empty[length(non_empty)])
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── 1. Load Data ─────────────────────────────────────────────────────────────
tab <- read.table(table_tsv, header = TRUE, sep = "\t", skip = 1, 
                  row.names = 1, check.names = FALSE, comment.char = "", quote = "")

tax_df <- read.table(tax_tsv, header = TRUE, sep = "\t", 
                     row.names = 1, check.names = FALSE, quote = "")

meta <- read.table(meta_file, header = TRUE, sep = "\t", 
                   row.names = 1, check.names = FALSE, comment.char = "", quote = "")

# ── 2. Collapse Taxonomy (Prevent Duplicates) ────────────────────────────────
# Map hashes to names
tax_strings <- tax_df[rownames(tab), "Taxon"]
tax_strings <- gsub(";", "|", tax_strings) 
tax_strings[is.na(tax_strings)] <- paste0("Unassigned_", rownames(tab)[is.na(tax_strings)])

# Add tax column and aggregate
tab$TaxonomyGroup <- tax_strings
tab_collapsed <- tab %>%
  group_by(TaxonomyGroup) %>%
  summarise(across(everything(), sum)) %>%
  as.data.frame()

rownames(tab_collapsed) <- tab_collapsed$TaxonomyGroup
tab_final <- tab_collapsed[, -1] # Remove the TaxonomyGroup column

# ── 3. Sync Samples (Fix 50 vs 49 error) ─────────────────────────────────────
common_samples <- intersect(colnames(tab_final), rownames(meta))
tab_final <- tab_final[, common_samples, drop = FALSE]
meta <- meta[common_samples, , drop = FALSE]

message("Final check: Table has ", ncol(tab_final), " samples. Metadata has ", nrow(meta), " samples.")

# ── 4. microbiomeMarker LEfSe ────────────────────────────────────────────────
if (requireNamespace("microbiomeMarker", quietly = TRUE)) {
  library(microbiomeMarker)
  library(phyloseq)

  otu <- otu_table(as.matrix(tab_final), taxa_are_rows = TRUE)
  sam <- sample_data(meta)
  
  # Create tax_table to satisfy microbiomeMarker requirements
  tax_mat <- as.matrix(data.frame(Taxon = rownames(tab_final), row.names = rownames(tab_final)))
  tax_tab <- tax_table(tax_mat)
  
  ps <- phyloseq(otu, sam, tax_tab)

  message("Running LEfSe...")
  mm <- tryCatch(
    run_lefse(ps, group = group_col, 
              multigrp_strat = TRUE, 
              lda_cutoff = 2.0, 
              wilcoxon_cutoff = 0.05),
    error = function(e) { message("LEfSe failed: ", e$message); NULL }
  )

  if (!is.null(mm) && !is.null(marker_table(mm)) && nrow(marker_table(mm)) > 0) {
    res_df <- marker_table(mm) %>% as.data.frame()
    write.table(res_df, file.path(out_dir, "lefse_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

    p_raw <- plot_ef_bar(mm) + theme_bw() + labs(title = "LEfSe Analysis (raw labels)")

    # Cleaned plot: apply format_taxon_label to y-axis labels
    # NOTE: ggplot passes the full label vector at once, so sapply() is required.
    p_clean <- p_raw +
      scale_y_discrete(labels = if (strategy != "none") function(x) sapply(x, format_taxon_label) else identity) +
      labs(title = "LEfSe Analysis")

    # Cleaned PDF (primary output)
    pdf(file.path(out_dir, "lefse_plots.pdf"), width = 12, height = 8)
    print(p_clean)
    dev.off()

    # Raw PDF (always written so Snakemake output is satisfied)
    pdf(file.path(out_dir, "lefse_plots_raw.pdf"), width = 12, height = 8)
    print(p_raw)
    dev.off()

    message("LEfSe complete.")
    quit(save = "no", status = 0)
  }
}

# ── 5. Fallback Analysis ─────────────────────────────────────────────────────
message("Using Fallback Kruskal-Wallis...")
results <- lapply(rownames(tab_final), function(taxon) {
  dat <- data.frame(val = as.numeric(tab_final[taxon,]), grp = as.factor(meta[[group_col]]))
  p <- tryCatch(kruskal.test(val ~ grp, data = dat)$p.value, error = function(e) NA)
  data.frame(taxon = taxon, p_value = p, stringsAsFactors = FALSE)
})

res_df <- do.call(rbind, results)
res_df$p_adj <- p.adjust(res_df$p_value, method = "BH")
write.table(res_df, file.path(out_dir, "lefse_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

top_sig <- head(res_df[!is.na(res_df$p_adj) & res_df$p_adj < 0.1, ], 20)
if (nrow(top_sig) == 0) top_sig <- head(res_df, 20)

make_kw_plot <- function(df, clean) {
  # ggplot passes the full label vector at once; sapply() vectorises the function.
  lab_fn <- if (clean && strategy != "none") function(x) sapply(x, format_taxon_label) else identity
  ggplot(df, aes(x = reorder(taxon, p_adj), y = -log10(p_adj + 1e-10))) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    theme_bw(base_size = 12) +
    scale_x_discrete(labels = lab_fn) +
    labs(
      title = if (clean) "Kruskal-Wallis Significance" else "Kruskal-Wallis Significance (raw labels)",
      x = NULL, y = expression(-log[10](p[adj]))
    )
}

# Raw PDF
pdf(file.path(out_dir, "lefse_plots_raw.pdf"), width = 10, height = 6)
if (nrow(top_sig) > 0) print(make_kw_plot(top_sig, clean = FALSE))
dev.off()

# Cleaned PDF
pdf(file.path(out_dir, "lefse_plots.pdf"), width = 10, height = 6)
if (nrow(top_sig) > 0) print(make_kw_plot(top_sig, clean = TRUE))
dev.off()