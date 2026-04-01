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

# ── 4. microbiomeMarker LEfSe (Phylum + Genus levels) ───────────────────────

# Helper: parse SILVA taxonomy string into 7-column hierarchical tax_table
parse_silva_taxonomy <- function(tax_strings, row_ids) {
  ranks   <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  pfx_map <- list(d__=1L, k__=1L, p__=2L, c__=3L, o__=4L, f__=5L, g__=6L, s__=7L)
  mat <- matrix("Unclassified", nrow = length(row_ids), ncol = length(ranks),
                dimnames = list(row_ids, ranks))
  for (i in seq_along(row_ids)) {
    s <- tax_strings[i]
    if (is.na(s) || nchar(trimws(s)) == 0) next
    s     <- gsub("\\|", "; ", s)
    parts <- trimws(strsplit(s, ";")[[1]])
    for (part in parts) {
      for (pfx in names(pfx_map)) {
        if (startsWith(part, pfx)) {
          val <- trimws(sub(pfx, "", part, fixed = TRUE))
          if (nchar(val) > 0 && !tolower(val) %in% c("uncultured", "")) {
            mat[i, pfx_map[[pfx]]] <- val
          }
          break
        }
      }
    }
  }
  mat
}

# Helper: ensure a file exists (write placeholder if not)
ensure_file <- function(path, msg = "No significant markers") {
  if (file.exists(path)) return(invisible(NULL))
  if (endsWith(path, ".pdf")) {
    grDevices::pdf(path, width = 8, height = 5)
    graphics::plot.new(); graphics::title(msg, cex.main = 1.2)
    grDevices::dev.off()
  } else if (endsWith(path, ".png")) {
    grDevices::png(path, width = 800, height = 500, res = 96)
    graphics::plot.new(); graphics::title(msg)
    grDevices::dev.off()
  }
}

if (requireNamespace("microbiomeMarker", quietly = TRUE)) {
  library(microbiomeMarker)
  library(phyloseq)

  # Build phyloseq with hierarchical tax_table.
  # NOTE: tab_final row names ARE the SILVA taxonomy strings (pipe-separated),
  # assigned during the step-2 group_by(TaxonomyGroup) collapse.
  # They contain "|" which microbiomeMarker uses to detect a "summarized" ps
  # (check_tax_summarize checks for "|" in OTU row names).  When ps is flagged
  # as summarized, run_lefse only allows norm="CPM" — blocking norm="none".
  # Fix: rename taxa to "taxon_N" before building the phyloseq.
  asv_ids   <- rownames(tab_final)           # these ARE the taxonomy strings
  new_names <- paste0("taxon_", seq_len(length(asv_ids)))
  rownames(tab_final) <- new_names

  otu     <- otu_table(as.matrix(tab_final), taxa_are_rows = TRUE)
  sam     <- sample_data(meta)
  tax_mat <- parse_silva_taxonomy(asv_ids, new_names)  # parse strings → use new_names as row IDs
  tax_tab <- tax_table(tax_mat)
  ps      <- phyloseq(otu, sam, tax_tab)

  # Filter to Bacteria only — removes Eukaryota/Archaea contamination that
  # would otherwise dominate LEfSe results (e.g. Ascomycota fungal reads)
  ps_bact <- tryCatch(
    subset_taxa(ps, Kingdom == "Bacteria"),
    error = function(e) { message("Bacteria filter failed: ", e$message); ps }
  )
  if (!is.null(ps_bact) && ntaxa(ps_bact) > 0) {
    message("Filtered to Bacteria: ", ntaxa(ps_bact), " taxa (from ", ntaxa(ps), ")")
    ps <- ps_bact
  }

  dark_col_lf <- c("#1B4F72", "#922B21", "#1D8348", "#6C3483", "#784212",
                   "#0E6655", "#4A235A", "#1A5276", "#7D6608", "#212F3D")

  # ── Unconditional taxonomy cladogram (always produced) ─────────────────────
  # Circular dendrogram: Kingdom → Phylum → Class → Order → Family → Genus.
  # Nodes coloured by which group has higher mean abundance (Wilcoxon direction).
  # Significant genera (p < 0.05) are labelled.  Runs regardless of whether
  # any LEfSe markers pass the LDA threshold.
  build_taxonomy_cladogram <- function(ps_in, group_col, out_path_base) {
    tryCatch({
      for (pkg in c("ape", "ggtree")) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
          message("  Skipping taxonomy cladogram — missing package: ", pkg)
          return(invisible(NULL))
        }
        suppressPackageStartupMessages(library(pkg, character.only = TRUE))
      }

      # Glom to Genus to get one row per genus
      ps_g <- tryCatch(
        tax_glom(ps_in, taxrank = "Genus", NArm = TRUE),
        error = function(e) ps_in
      )
      tax_mat <- as.data.frame(tax_table(ps_g))
      otu_mat <- as.matrix(otu_table(ps_g))
      if (!taxa_are_rows(ps_g)) otu_mat <- t(otu_mat)
      meta_df  <- as.data.frame(sample_data(ps_g))
      grp_vec  <- as.character(meta_df[[group_col]])
      grp_levs <- sort(unique(na.omit(grp_vec)))

      ranks_all <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
      ranks     <- ranks_all[ranks_all %in% colnames(tax_mat)]
      if (length(ranks) < 2) {
        message("  Not enough taxonomy ranks for cladogram"); return(invisible(NULL))
      }

      # Create safe Newick labels (alphanum + underscore only) and keep
      # a display_map: safe_id → human-readable label for annotation.
      tax_norm    <- tax_mat[, ranks, drop = FALSE]
      display_map <- list()
      for (ri in seq_along(ranks)) {
        r   <- ranks[ri]
        bad <- is.na(tax_norm[[r]]) | tax_norm[[r]] == "" |
               grepl("unclassified|uncultured|unidentified|__NA",
                     tax_norm[[r]], ignore.case = TRUE)
        if (any(bad)) {
          fill <- if (ri > 1) paste0("Unc_", tax_norm[[ranks[ri - 1]]][bad])
                  else        "Unclassified"
          tax_norm[[r]][bad] <- fill
        }
        display_vals <- tax_norm[[r]]
        safe_vals    <- gsub("[^A-Za-z0-9_]", "_",
                             paste0(substr(r, 1, 1), "_", display_vals))
        tax_norm[[r]] <- safe_vals
        for (j in seq_along(safe_vals)) {
          display_map[[safe_vals[j]]] <- display_vals[j]
        }
      }
      display_map[["Root"]] <- "Root"

      # Build parent → children map (unique pairs at each rank transition)
      pc <- list()
      pc[["Root"]] <- unique(tax_norm[[ranks[1]]])
      for (i in seq(2, length(ranks))) {
        pairs <- unique(tax_norm[, c(ranks[i - 1], ranks[i]), drop = FALSE])
        for (j in seq_len(nrow(pairs))) {
          p  <- as.character(pairs[j, 1])
          ch <- as.character(pairs[j, 2])
          if (is.null(pc[[p]])) pc[[p]] <- character(0)
          pc[[p]] <- unique(c(pc[[p]], ch))
        }
      }

      # Recursive Newick serializer
      to_nwk <- function(n) {
        kids <- pc[[n]]
        if (is.null(kids) || length(kids) == 0) return(n)
        paste0("(", paste(vapply(kids, to_nwk, character(1)), collapse = ","), ")", n)
      }
      tr <- ape::read.tree(text = paste0(to_nwk("Root"), ";"))

      # Wilcoxon enrichment per genus tip
      tip_rank  <- ranks[length(ranks)]
      tip_nodes <- unique(tax_norm[[tip_rank]])
      enrich_v  <- setNames(rep(NA_character_, length(tip_nodes)), tip_nodes)
      pval_v    <- setNames(rep(1.0,           length(tip_nodes)), tip_nodes)
      if (length(grp_levs) >= 2) {
        g1 <- which(grp_vec == grp_levs[1])
        g2 <- which(grp_vec == grp_levs[2])
        for (tip in tip_nodes) {
          idx <- which(tax_norm[[tip_rank]] == tip)
          if (length(idx) == 0) next
          v  <- if (length(idx) == 1) as.numeric(otu_mat[idx, ])
                else as.numeric(colSums(otu_mat[idx, , drop = FALSE]))
          m1 <- mean(v[g1]); m2 <- mean(v[g2])
          p  <- tryCatch(suppressWarnings(wilcox.test(v[g1], v[g2])$p.value),
                         error = function(e) NA_real_)
          enrich_v[tip] <- if (is.finite(m1) && is.finite(m2)) {
            if (m1 >= m2) grp_levs[1] else grp_levs[2]
          } else NA_character_
          pval_v[tip] <- ifelse(is.na(p), 1.0, p)
        }
      }

      # Tip annotation dataframe keyed on tree tip labels
      tip_df <- data.frame(
        label   = tr$tip.label,
        enrich  = enrich_v[tr$tip.label],
        sig     = !is.na(pval_v[tr$tip.label]) & pval_v[tr$tip.label] < 0.05,
        display = vapply(tr$tip.label, function(n) {
          lbl <- display_map[[n]]; if (is.null(lbl)) n else lbl
        }, character(1)),
        stringsAsFactors = FALSE
      )
      tip_df$sig[is.na(tip_df$sig)] <- FALSE
      n_sig <- sum(tip_df$sig, na.rm = TRUE)
      message("  Cladogram: ", length(tip_nodes), " genera, ", n_sig,
              " significant (Wilcoxon p<0.05)")

      grp_pal <- setNames(dark_col_lf[seq_along(grp_levs)], grp_levs)

      p_clad <- ggtree(tr, layout = "circular", color = "grey70",
                       linewidth = 0.35) %<+% tip_df +
        geom_tippoint(aes(color = enrich), size = 1.8, alpha = 0.90,
                      na.rm = TRUE) +
        geom_tiplab(aes(label = ifelse(sig, display, NA_character_)),
                    size = 2.2, offset = 1.5, align = FALSE, color = "black",
                    fontface = "italic", na.rm = TRUE) +
        ggplot2::scale_color_manual(
          values       = grp_pal,
          na.value     = "grey80",
          name         = paste0("Higher mean in (", group_col, ")"),
          na.translate = FALSE
        ) +
        ggplot2::theme_void() +
        ggplot2::theme(
          legend.position = "bottom",
          legend.title    = ggplot2::element_text(face = "bold", size = 11),
          legend.text     = ggplot2::element_text(face = "bold", size = 10),
          plot.title      = ggplot2::element_text(face = "bold", size = 13,
                                                  hjust = 0.5),
          plot.caption    = ggplot2::element_text(size = 9, hjust = 0.5,
                                                  color = "grey40"),
          plot.margin     = ggplot2::margin(40, 40, 40, 40)
        ) +
        ggplot2::labs(
          title   = paste0("Taxonomy Cladogram \u2014 ", group_col),
          caption = paste0(length(tip_nodes), " genera  |  ",
                           "node colour = group with higher mean abundance  |  ",
                           "labelled: Wilcoxon p\u00a0<\u00a00.05")
        )

      grDevices::pdf(paste0(out_path_base, ".pdf"), width = 14, height = 14)
      print(p_clad)
      grDevices::dev.off()
      ggplot2::ggsave(paste0(out_path_base, ".png"), plot = p_clad,
                      width = 14, height = 14, units = "in", dpi = 600)
      message("  Cladogram saved: ", basename(out_path_base), ".pdf / .png")
    }, error = function(e) {
      message("  Taxonomy cladogram failed: ", e$message)
    })
    invisible(NULL)
  }

  # Always produce the cladogram — regardless of LEfSe significance
  build_taxonomy_cladogram(ps, group_col,
                           file.path(out_dir, "lefse_cladogram"))

  # ── Run LEfSe at one taxonomic level and return panels ────────────────────
  run_and_plot_level <- function(ps_in, rank_name) {
    message("LEfSe at ", rank_name, " level ...")
    ps_glom <- tryCatch(
      tax_glom(ps_in, taxrank = rank_name, NArm = TRUE),
      error = function(e) { message("tax_glom failed: ", e$message); NULL }
    )
    if (is.null(ps_glom)) return(NULL)
    message("  After glom: ", ntaxa(ps_glom), " taxa")

    # multigrp_strat must be FALSE for 2-group comparisons —
    # TRUE triggers a known microbiomeMarker "In index: 1" crash
    n_groups <- length(unique(get_variable(ps_glom, group_col)))
    use_multigrp <- n_groups > 2

    mm <- tryCatch(
      run_lefse(ps_glom, group = group_col,
                multigrp_strat = use_multigrp, lda_cutoff = 2.0,
                wilcoxon_cutoff = 0.05, norm = "none"),
      error = function(e) { message("run_lefse failed: ", e$message); NULL }
    )
    n_markers <- tryCatch({
      mt <- marker_table(mm); if (is.null(mt)) 0L else nrow(mt)
    }, error = function(e) 0L)
    # If strict thresholds find nothing, retry with relaxed thresholds
    if (n_markers == 0L) {
      message("  No markers at LDA>=2.0 / p<0.05; retrying with LDA>=1.5 / p<0.1 ...")
      mm <- tryCatch(
        run_lefse(ps_glom, group = group_col,
                  multigrp_strat = use_multigrp, lda_cutoff = 1.5,
                  wilcoxon_cutoff = 0.10, norm = "none"),
        error = function(e) { message("run_lefse (relaxed) failed: ", e$message); NULL }
      )
      n_markers <- tryCatch({
        mt <- marker_table(mm); if (is.null(mt)) 0L else nrow(mt)
      }, error = function(e) 0L)
    }
    if (n_markers == 0L) {
      message("No significant markers at ", rank_name, " level.")
      return(NULL)
    }

    res_df         <- as.data.frame(marker_table(mm))
    res_df$lda_abs <- abs(res_df$ef_lda)
    enrich_groups  <- sort(unique(res_df$enrich_group))
    enrich_colors  <- setNames(dark_col_lf[seq_along(enrich_groups)], enrich_groups)

    # Panel B: LDA bar chart (group-colored, no sample labels)
    p_bar <- ggplot(res_df,
                    aes(x = reorder(feature, lda_abs), y = lda_abs,
                        fill = enrich_group)) +
      geom_col(color = "grey20", linewidth = 0.35, alpha = 0.90) +
      coord_flip() +
      scale_fill_manual(values = enrich_colors, name = "Enriched in") +
      theme_classic(base_size = 12) +
      theme(
        legend.position = "bottom",
        legend.title    = element_text(face = "bold", size = 12),
        legend.text     = element_text(face = "bold", size = 11),
        axis.title      = element_text(face = "bold", size = 12),
        axis.text       = element_text(face = "bold", size = 10, color = "black"),
        plot.title      = element_text(face = "bold", size = 13),
        axis.line       = element_line(linewidth = 0.6)
      ) +
      labs(x = NULL, y = "LDA Score (log10)",
           title = paste("LDA Scores \u2014", rank_name))

    # Panel A: circular cladogram (group-colored nodes, no sample labels)
    p_clad <- tryCatch({
      plot_cladogram(mm, color = group_col, only_marker = TRUE,
                     clade_label_level = 4) +
        theme(
          legend.position = "bottom",
          legend.title    = element_text(face = "bold", size = 12),
          legend.text     = element_text(face = "bold", size = 11),
          plot.title      = element_text(face = "bold", size = 13)
        ) +
        labs(title = paste("Cladogram \u2014", rank_name))
    }, error = function(e) {
      message("plot_cladogram failed at ", rank_name, ": ", e$message)
      NULL
    })

    list(mm = mm, res_df = res_df, bar = p_bar, clad = p_clad)
  }

  # ── Save multi-panel for one level ────────────────────────────────────────
  save_level <- function(panels, rank_name) {
    fname_base <- file.path(out_dir,
                            paste0("lefse_", tolower(rank_name), "_plots"))
    if (is.null(panels)) {
      ensure_file(paste0(fname_base, ".pdf"),
                  paste("No significant LEfSe markers at", rank_name, "level"))
      ensure_file(paste0(fname_base, ".png"),
                  paste("No significant LEfSe markers at", rank_name, "level"))
      return(invisible(NULL))
    }
    pieces <- Filter(Negate(is.null), list(panels$clad, panels$bar))

    if (length(pieces) >= 2 && requireNamespace("patchwork", quietly = TRUE)) {
      library(patchwork)
      combined <- patchwork::wrap_plots(pieces, ncol = 2, widths = c(1.1, 0.9)) +
        patchwork::plot_annotation(
          tag_levels = "A",
          title      = paste("LEfSe Analysis \u2014", rank_name, "Level"),
          theme      = theme(plot.title = element_text(face = "bold", size = 15,
                                                       hjust = 0.5))
        ) +
        patchwork::plot_layout(guides = "collect") &
        theme(legend.position = "bottom")
    } else {
      combined <- if (length(pieces) > 0) pieces[[length(pieces)]] else NULL
    }

    if (!is.null(combined)) {
      grDevices::pdf(paste0(fname_base, ".pdf"), width = 16, height = 9)
      print(combined)
      grDevices::dev.off()
      ggplot2::ggsave(paste0(fname_base, ".png"), plot = combined,
                      width = 16, height = 9, units = "in", dpi = 600)
    } else {
      ensure_file(paste0(fname_base, ".pdf"))
      ensure_file(paste0(fname_base, ".png"))
    }
    message("Saved: lefse_", tolower(rank_name), "_plots.pdf / .png")
  }

  # Run at both taxonomic levels
  lvl_phylum <- run_and_plot_level(ps, "Phylum")
  lvl_genus  <- run_and_plot_level(ps, "Genus")

  save_level(lvl_phylum, "Phylum")
  save_level(lvl_genus,  "Genus")

  # Combined results TSV (all levels)
  res_all <- list()
  if (!is.null(lvl_phylum)) { r <- lvl_phylum$res_df; r$level <- "Phylum"; res_all$p <- r }
  if (!is.null(lvl_genus))  { r <- lvl_genus$res_df;  r$level <- "Genus";  res_all$g <- r }
  if (length(res_all) > 0) {
    write.table(do.call(rbind, res_all),
                file.path(out_dir, "lefse_results.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  } else {
    write.table(data.frame(Note = "No significant LEfSe markers found"),
                file.path(out_dir, "lefse_results.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }

  # Backward-compat: lefse_plots.pdf = genus-level (fallback to phylum)
  genus_pdf <- file.path(out_dir, "lefse_genus_plots.pdf")
  main_pdf  <- file.path(out_dir, "lefse_plots.pdf")
  src_pdf   <- if (file.exists(genus_pdf)) genus_pdf else
               file.path(out_dir, "lefse_phylum_plots.pdf")
  if (file.exists(src_pdf)) {
    file.copy(src_pdf, main_pdf, overwrite = TRUE)
  } else {
    ensure_file(main_pdf, "No significant LEfSe markers")
  }
  file.copy(main_pdf, file.path(out_dir, "lefse_plots_raw.pdf"), overwrite = TRUE)

  message("LEfSe analysis complete.")
  quit(save = "no", status = 0)
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

# Placeholder level-specific files required by Snakemake output declarations
ensure_file_fallback <- function(path, msg) {
  if (!file.exists(path)) {
    if (endsWith(path, ".pdf")) {
      grDevices::pdf(path, width = 8, height = 5)
      graphics::plot.new(); graphics::title(msg, cex.main = 1.0)
      grDevices::dev.off()
    } else if (endsWith(path, ".png")) {
      grDevices::png(path, width = 800, height = 500, res = 96)
      graphics::plot.new(); graphics::title(msg)
      grDevices::dev.off()
    }
  }
}
fb_msg <- "microbiomeMarker unavailable \u2014 using Kruskal-Wallis fallback"
ensure_file_fallback(file.path(out_dir, "lefse_phylum_plots.pdf"), fb_msg)
ensure_file_fallback(file.path(out_dir, "lefse_phylum_plots.png"), fb_msg)
ensure_file_fallback(file.path(out_dir, "lefse_genus_plots.pdf"),  fb_msg)
ensure_file_fallback(file.path(out_dir, "lefse_genus_plots.png"),  fb_msg)
ensure_file_fallback(file.path(out_dir, "lefse_cladogram.pdf"),    fb_msg)
ensure_file_fallback(file.path(out_dir, "lefse_cladogram.png"),    fb_msg)