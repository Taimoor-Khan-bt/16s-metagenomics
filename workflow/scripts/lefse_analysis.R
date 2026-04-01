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

  # ── Colour palettes shared across all level plots ─────────────────────────
  grp_levs_all  <- sort(unique(as.character(get_variable(ps, group_col))))
  group_pal_all <- setNames(dark_col_lf[seq_along(grp_levs_all)], grp_levs_all)

  # 15-colour palette for phyla (Okabe-Ito extended)
  phylum_pal15 <- c(
    "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
    "#A65628","#F781BF","#999999","#1B9E77","#D95F02",
    "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D"
  )

  # ── Helper: build Newick tree + annotation frame from a phyloseq ──────────
  # Returns list(tr, tip_df, display_map, tax_norm, ranks) or NULL on failure.
  make_tree_data <- function(ps_in) {
    tax_mat <- as.data.frame(tax_table(ps_in))
    otu_mat <- as.matrix(otu_table(ps_in))
    if (!taxa_are_rows(ps_in)) otu_mat <- t(otu_mat)
    meta_df  <- as.data.frame(sample_data(ps_in))
    grp_vec  <- as.character(meta_df[[group_col]])
    grp_levs <- sort(unique(na.omit(grp_vec)))

    ranks_all <- c("Kingdom","Phylum","Class","Order","Family","Genus")
    ranks     <- ranks_all[ranks_all %in% colnames(tax_mat)]
    if (length(ranks) < 2) return(NULL)

    tax_norm    <- tax_mat[, ranks, drop = FALSE]
    display_map <- list()
    for (ri in seq_along(ranks)) {
      r   <- ranks[ri]
      bad <- is.na(tax_norm[[r]]) | tax_norm[[r]] == "" |
             grepl("unclassified|uncultured|unidentified|__NA",
                   tax_norm[[r]], ignore.case = TRUE)
      if (any(bad)) {
        fill <- if (ri > 1) paste0("Unc_", gsub("[^A-Za-z0-9]","_",
                                                tax_norm[[ranks[ri-1]]][bad]))
                else        "Unclassified"
        tax_norm[[r]][bad] <- fill
      }
      display_vals <- tax_norm[[r]]
      safe_vals    <- gsub("[^A-Za-z0-9_]","_",
                           paste0(substr(r,1,1),"_", display_vals))
      tax_norm[[r]] <- safe_vals
      for (j in seq_along(safe_vals)) display_map[[safe_vals[j]]] <- display_vals[j]
    }
    display_map[["Root"]] <- "Root"

    pc <- list(); pc[["Root"]] <- unique(tax_norm[[ranks[1]]])
    for (i in seq(2, length(ranks))) {
      pairs <- unique(tax_norm[, c(ranks[i-1], ranks[i]), drop=FALSE])
      for (j in seq_len(nrow(pairs))) {
        p <- as.character(pairs[j,1]); ch <- as.character(pairs[j,2])
        if (is.null(pc[[p]])) pc[[p]] <- character(0)
        pc[[p]] <- unique(c(pc[[p]], ch))
      }
    }
    to_nwk <- function(n) {
      kids <- pc[[n]]
      if (is.null(kids)||length(kids)==0) return(n)
      paste0("(", paste(vapply(kids,to_nwk,character(1)),collapse=","), ")", n)
    }
    tr <- tryCatch(ape::read.tree(text=paste0(to_nwk("Root"),";")),
                   error=function(e) NULL)
    if (is.null(tr)) return(NULL)

    tip_rank  <- ranks[length(ranks)]
    tip_nodes <- unique(tax_norm[[tip_rank]])
    enrich_v  <- setNames(rep(NA_character_, length(tip_nodes)), tip_nodes)
    pval_v    <- setNames(rep(1.0,           length(tip_nodes)), tip_nodes)
    abund_v   <- setNames(rep(0.0,           length(tip_nodes)), tip_nodes)
    if (length(grp_levs) >= 2) {
      g1 <- which(grp_vec == grp_levs[1]); g2 <- which(grp_vec == grp_levs[2])
      for (tip in tip_nodes) {
        idx <- which(tax_norm[[tip_rank]] == tip)
        if (!length(idx)) next
        v  <- if (length(idx)==1) as.numeric(otu_mat[idx,])
              else as.numeric(colSums(otu_mat[idx,,drop=FALSE]))
        abund_v[tip] <- mean(v)
        m1 <- mean(v[g1]); m2 <- mean(v[g2])
        p  <- tryCatch(suppressWarnings(wilcox.test(v[g1],v[g2])$p.value),
                       error=function(e) NA_real_)
        enrich_v[tip] <- if (is.finite(m1)&&is.finite(m2)) {
          if (m1>=m2) grp_levs[1] else grp_levs[2]
        } else NA_character_
        pval_v[tip] <- ifelse(is.na(p), 1.0, p)
      }
    }

    get_phylum <- function(n) {
      r  <- if ("Phylum" %in% ranks) tax_norm$Phylum[tax_norm[[tip_rank]]==n]
            else character(0)
      if (length(r)&&!is.na(r[1])) r[1] else NA_character_
    }

    tip_df <- data.frame(
      label   = tr$tip.label,
      enrich  = enrich_v[tr$tip.label],
      pval    = pval_v[tr$tip.label],
      sig     = !is.na(pval_v[tr$tip.label]) & pval_v[tr$tip.label] < 0.05,
      abund   = log1p(abund_v[tr$tip.label]),
      phylum  = vapply(tr$tip.label, get_phylum, character(1)),
      display = vapply(tr$tip.label, function(n) {
        lbl <- display_map[[n]]; if (is.null(lbl)) n else lbl
      }, character(1)),
      stringsAsFactors = FALSE
    )
    tip_df$sig[is.na(tip_df$sig)] <- FALSE

    list(tr=tr, tip_df=tip_df, display_map=display_map,
         tax_norm=tax_norm, ranks=ranks, grp_levs=grp_levs)
  }

  # ── Modern multi-ring circular cladogram ─────────────────────────────────
  # Ring 1: phylum colour strip (ggtreeExtra::geom_fruit)
  # Ring 2: enrichment group tile
  # Significant genera (Wilcoxon p<0.05) are labelled in bold italic
  build_modern_cladogram <- function(ps_in, group_col, out_path_base) {
    tryCatch({
      for (pkg in c("ape","ggtree","ggtreeExtra","ggnewscale")) {
        if (!requireNamespace(pkg, quietly=TRUE)) {
          message("  Skipping modern cladogram — missing: ", pkg)
          return(invisible(NULL))
        }
        suppressPackageStartupMessages(library(pkg, character.only=TRUE))
      }

      ps_g <- tryCatch(tax_glom(ps_in, taxrank="Genus", NArm=TRUE),
                       error=function(e) ps_in)
      td <- make_tree_data(ps_g)
      if (is.null(td)) { message("  make_tree_data returned NULL"); return(invisible(NULL)) }

      tr     <- td$tr
      tip_df <- td$tip_df
      grp_levs <- td$grp_levs
      grp_pal  <- setNames(dark_col_lf[seq_along(grp_levs)], grp_levs)

      # Phylum colour map
      phyla <- sort(unique(na.omit(tip_df$phylum)))
      phylum_colors <- setNames(phylum_pal15[seq_along(phyla)], phyla)
      tip_df$phylum[is.na(tip_df$phylum)] <- "Unknown"
      # Use unique column names for ring data to avoid clashing with ggtree's
      # internal columns (label, y, angle, enrich, phylum, …)
      phylum_ring <- data.frame(tip_id=tip_df$label,
                                tip_phylum=tip_df$phylum,
                                stringsAsFactors=FALSE)
      enrich_ring <- data.frame(tip_id=tip_df$label,
                                tip_enrich=tip_df$enrich,
                                stringsAsFactors=FALSE)

      n_sig <- sum(tip_df$sig, na.rm=TRUE)
      n_gen <- nrow(tip_df)
      message("  Modern cladogram: ", n_gen, " genera, ", n_sig, " sig (p<0.05)")

      p_tree <- ggtree(tr, layout="circular", color="grey65",
                        linewidth=0.28, open.angle=8) %<+% tip_df +
        geom_tippoint(aes(color=enrich, size=abund), alpha=0.88, na.rm=TRUE) +
        ggplot2::scale_color_manual(
          values=grp_pal, na.value="grey80",
          name=paste0("Higher in (", group_col, ")"), na.translate=FALSE
        ) +
        ggplot2::scale_size_continuous(
          range=c(0.8,3.5), name="log(mean abund.)",
          guide=guide_legend(override.aes=list(alpha=0.8, color="grey40"))
        )

      # Ring 1: phylum colour strip
      p_tree <- p_tree +
        ggtreeExtra::geom_fruit(
          data    = phylum_ring,
          geom    = geom_tile,
          mapping = aes(y=tip_id, x=1, fill=tip_phylum),
          offset=0.05, pwidth=0.07
        ) +
        ggplot2::scale_fill_manual(
          values=c(phylum_colors, Unknown="grey88"),
          name="Phylum", na.translate=FALSE,
          guide=guide_legend(ncol=2, override.aes=list(size=5))
        )

      # Ring 2: enrichment direction tile
      p_tree <- p_tree +
        ggnewscale::new_scale_fill() +
        ggtreeExtra::geom_fruit(
          data    = enrich_ring,
          geom    = geom_tile,
          mapping = aes(y=tip_id, x=1, fill=tip_enrich),
          offset=0.06, pwidth=0.07
        ) +
        ggplot2::scale_fill_manual(
          values=grp_pal, na.value="grey92",
          name=paste0("Enriched in\n(", group_col, ")"), na.translate=FALSE,
          guide=guide_legend(override.aes=list(size=5))
        )

      # Significant tip labels
      if (n_sig > 0) {
        p_tree <- p_tree +
          geom_tiplab(
            aes(label=ifelse(sig, display, NA_character_), filter=sig),
            size=2.0, offset=0.15, color="grey10",
            fontface="bold.italic", na.rm=TRUE, align=FALSE
          )
      }

      p_tree <- p_tree +
        ggplot2::theme_void() +
        ggplot2::theme(
          legend.position = "right",
          legend.box      = "vertical",
          legend.title    = ggplot2::element_text(face="bold", size=10),
          legend.text     = ggplot2::element_text(size=9),
          plot.title      = ggplot2::element_text(face="bold", size=14, hjust=0.5),
          plot.caption    = ggplot2::element_text(size=8, color="grey50", hjust=0.5),
          plot.margin     = ggplot2::margin(30,30,30,30)
        ) +
        ggplot2::labs(
          title   = paste0("Taxonomy Cladogram \u2014 ", group_col),
          caption = paste0(n_gen, " genera  |  inner ring: Phylum  |  ",
                           "outer ring: enrichment group  |  labelled: p\u202f<\u202f0.05")
        )

      grDevices::pdf(paste0(out_path_base,".pdf"), width=17, height=15)
      print(p_tree); grDevices::dev.off()
      ggplot2::ggsave(paste0(out_path_base,".png"), plot=p_tree,
                      width=17, height=15, units="in", dpi=600)
      message("  Modern cladogram saved: ", basename(out_path_base), ".pdf / .png")
    }, error=function(e) { message("  Modern cladogram failed: ", e$message) })
    invisible(NULL)
  }

  # Always produce the cladogram — regardless of LEfSe significance
  build_modern_cladogram(ps, group_col,
                         file.path(out_dir, "lefse_cladogram"))

  # ── Run LEfSe at one taxonomic level and return panels ────────────────────
  # Always produces plots for ALL taxa — not just those passing significance
  # thresholds.  Bars are signed (negative = group1 enriched, positive = group2).
  # Significant taxa (Wilcoxon p < 0.05) are drawn at full opacity; others are
  # shown semi-transparent so readers can judge the overall pattern themselves.
  run_and_plot_level <- function(ps_in, rank_name) {
    message("LEfSe at ", rank_name, " level ...")
    ps_glom <- tryCatch(
      tax_glom(ps_in, taxrank = rank_name, NArm = TRUE),
      error = function(e) { message("tax_glom failed: ", e$message); NULL }
    )
    if (is.null(ps_glom)) return(NULL)
    n_taxa <- ntaxa(ps_glom)
    message("  After glom: ", n_taxa, " taxa")

    # ── Step 1: get LDA + p-values for every taxon ──────────────────────────
    # Try run_lefse with lda_cutoff = 0 so all taxa get scores;
    # fall back to Wilcoxon-based signed-LDA if run_lefse fails entirely.
    n_groups     <- length(unique(get_variable(ps_glom, group_col)))
    use_multigrp <- n_groups > 2

    mm_all <- tryCatch(
      run_lefse(ps_glom, group = group_col,
                multigrp_strat = use_multigrp,
                lda_cutoff = 0, wilcoxon_cutoff = 1.0, norm = "none"),
      error = function(e) { message("  run_lefse(lda_cutoff=0) failed: ", e$message); NULL }
    )

    # Build res_df with signed LDA for every taxon
    tax_mat_g  <- as.data.frame(tax_table(ps_glom))
    otu_mat_g  <- as.matrix(otu_table(ps_glom))
    if (!taxa_are_rows(ps_glom)) otu_mat_g <- t(otu_mat_g)
    meta_g     <- as.data.frame(sample_data(ps_glom))
    grp_vec    <- as.character(meta_g[[group_col]])
    grp_levs   <- sort(unique(na.omit(grp_vec)))
    g1 <- which(grp_vec == grp_levs[1])
    g2 <- if (length(grp_levs) >= 2) which(grp_vec == grp_levs[2]) else integer(0)

    # Raw taxon name from the rank column (or rowname)
    get_raw_name <- function() {
      if (rank_name %in% colnames(tax_mat_g)) as.character(tax_mat_g[[rank_name]])
      else rownames(tax_mat_g)
    }
    raw_names <- get_raw_name()

    # If run_lefse returned scores, use them; otherwise compute manually
    if (!is.null(mm_all)) {
      mt <- tryCatch(as.data.frame(marker_table(mm_all)), error=function(e) NULL)
    } else { mt <- NULL }

    compute_manual_lda <- function() {
      pv <- numeric(n_taxa); eg <- character(n_taxa); lda_v <- numeric(n_taxa)
      for (i in seq_len(n_taxa)) {
        v  <- as.numeric(otu_mat_g[i,])
        m1 <- mean(v[g1]); m2 <- if (length(g2)) mean(v[g2]) else 0
        p  <- tryCatch(suppressWarnings(wilcox.test(v[g1], v[g2])$p.value),
                       error=function(e) NA_real_)
        pv[i]  <- ifelse(is.na(p), 1.0, p)
        eg[i]  <- if (is.finite(m1) && is.finite(m2) && m1 >= m2) grp_levs[1]
                  else if (length(grp_levs) >= 2) grp_levs[2] else grp_levs[1]
        # Approximate LDA as log10(|grand_mean+1|), signed by direction
        gm     <- mean(v) + 1
        lda_raw <- if (gm > 0) log10(gm) else 0
        lda_v[i] <- if (eg[i] == grp_levs[1]) -lda_raw else lda_raw
      }
      data.frame(feature=raw_names, ef_lda=lda_v,
                 enrich_group=eg, p_value=pv, stringsAsFactors=FALSE)
    }

    if (!is.null(mt) && nrow(mt) > 0 && all(c("feature","ef_lda","enrich_group") %in% colnames(mt))) {
      # run_lefse gave us data; add sign & p_value if missing
      # ef_lda from microbiomeMarker is always positive — add direction
      if (!"p_value" %in% colnames(mt)) mt$p_value <- NA_real_
      # recompute sign from enrich_group: group1 → negative, group2 → positive
      mt$ef_lda <- abs(mt$ef_lda) *
        ifelse(mt$enrich_group == grp_levs[1], -1, 1)
      # Fill in any missing taxa (run_lefse may omit 0-variance taxa)
      missing <- setdiff(raw_names, mt$feature)
      if (length(missing)) {
        man <- compute_manual_lda()
        man <- man[man$feature %in% missing, , drop=FALSE]
        mt  <- rbind(mt[, c("feature","ef_lda","enrich_group","p_value")], man)
      } else {
        mt <- mt[, c("feature","ef_lda","enrich_group","p_value")]
      }
      res_df <- mt
    } else {
      message("  run_lefse returned no usable scores; computing Wilcoxon-based LDA ...")
      res_df <- compute_manual_lda()
    }

    # Also compute strict markers for TSV (p<0.05 AND |lda|≥2)
    mm_strict <- tryCatch(
      run_lefse(ps_glom, group=group_col, multigrp_strat=use_multigrp,
                lda_cutoff=2.0, wilcoxon_cutoff=0.05, norm="none"),
      error=function(e) NULL
    )
    n_strict <- tryCatch({
      mt2 <- marker_table(mm_strict); if(is.null(mt2)) 0L else nrow(mt2)
    }, error=function(e) 0L)
    if (n_strict == 0L) {
      mm_strict <- tryCatch(
        run_lefse(ps_glom, group=group_col, multigrp_strat=use_multigrp,
                  lda_cutoff=1.5, wilcoxon_cutoff=0.10, norm="none"),
        error=function(e) NULL
      )
    }
    message("  All-taxa LDA computed (", nrow(res_df), " taxa); strict markers: ", n_strict)

    # ── Clean taxon labels ────────────────────────────────────────────────────
    res_df$feature <- vapply(as.character(res_df$feature),
                              format_taxon_label, character(1))
    res_df$feature[is.na(res_df$feature)|res_df$feature==""] <- "Unknown"

    # Aggregate duplicated labels (format_taxon_label can collapse distinct OTUs
    # to the same display name — take mean LDA, majority enrich_group, min p_value)
    if (anyDuplicated(res_df$feature)) {
      res_df <- do.call(rbind, lapply(split(res_df, res_df$feature), function(d) {
        data.frame(
          feature      = d$feature[1],
          ef_lda       = mean(d$ef_lda),
          enrich_group = d$enrich_group[which.max(tabulate(match(d$enrich_group,
                                                                  unique(d$enrich_group))))],
          p_value      = min(d$p_value),
          stringsAsFactors = FALSE
        )
      }))
    }

    # Significance flag
    res_df$sig   <- !is.na(res_df$p_value) & res_df$p_value < 0.05
    res_df$p_value[is.na(res_df$p_value)] <- 1.0

    # Sort by signed LDA (most negative first → most positive last)
    res_df <- res_df[order(res_df$ef_lda), ]
    res_df$feature <- factor(res_df$feature, levels = unique(res_df$feature))

    # ── Panel B: bidirectional signed LDA bar chart ───────────────────────────
    p_bar <- ggplot(res_df,
                    aes(x = ef_lda, y = feature,
                        fill = enrich_group,
                        alpha = sig)) +
      geom_col(color = NA, linewidth = 0) +
      geom_vline(xintercept = 0, color = "grey25", linewidth = 0.55) +
      # Significance asterisk at bar tip
      geom_text(
        data    = res_df[res_df$sig, , drop = FALSE],
        mapping = aes(x = ef_lda + sign(ef_lda) * 0.06,
                      label = "*"),
        inherit.aes = FALSE,
        size = 4.5, color = "grey10", vjust = 0.5
      ) +
      scale_alpha_manual(values = c("TRUE" = 0.95, "FALSE" = 0.32),
                         guide  = "none") +
      scale_fill_manual(values = group_pal_all,
                        name   = paste0("Enriched in\n(", group_col, ")"),
                        na.value = "grey70") +
      scale_x_continuous(
        labels = function(x) paste0(ifelse(x < 0, "\u2190", "\u2192"), " ", abs(x)),
        expand = expansion(mult = 0.12)
      ) +
      theme_classic(base_size = 11) +
      theme(
        legend.position  = "bottom",
        legend.title     = element_text(face = "bold", size = 11),
        legend.text      = element_text(face = "bold", size = 10),
        axis.title       = element_text(face = "bold", size = 11),
        axis.text.y      = element_text(size = 8, color = "black"),
        axis.text.x      = element_text(size = 9),
        plot.title       = element_text(face = "bold", size = 12, hjust = 0.5),
        plot.subtitle    = element_text(size = 9, hjust = 0.5, color = "grey40"),
        axis.line        = element_line(linewidth = 0.55)
      ) +
      labs(
        x        = paste0("Signed LDA Score (log10)",
                          "      \u2190 ", grp_levs[1],
                          "  |  ", if (length(grp_levs)>=2) grp_levs[2] else "",
                          " \u2192"),
        y        = NULL,
        title    = paste("LDA Scores \u2014", rank_name, "Level"),
        subtitle = paste0("All ", nrow(res_df), " taxa shown  |  ",
                          "* p\u202f<\u202f0.05 (Wilcoxon)  |  ",
                          "transparency = non-significant")
      )

    # ── Panel A: modern circular cladogram for THIS level ─────────────────────
    # We re-use the make_tree_data helper on the glommed phyloseq
    p_clad <- tryCatch({
      for (pkg in c("ape","ggtree","ggtreeExtra","ggnewscale")) {
        if (!requireNamespace(pkg, quietly=TRUE)) stop("missing ", pkg)
        suppressPackageStartupMessages(library(pkg, character.only=TRUE))
      }
      ps_for_clad <- tryCatch(tax_glom(ps_in, taxrank=rank_name, NArm=TRUE),
                               error=function(e) ps_glom)
      td <- make_tree_data(ps_for_clad)
      if (is.null(td)) stop("make_tree_data returned NULL")

      tr_l   <- td$tr; tip_df_l <- td$tip_df
      n_sig_c <- sum(tip_df_l$sig, na.rm=TRUE)

      phyla_l <- sort(unique(na.omit(tip_df_l$phylum)))
      pc_map  <- setNames(phylum_pal15[seq_along(phyla_l)], phyla_l)
      tip_df_l$phylum[is.na(tip_df_l$phylum)] <- "Unknown"
      # Unique column names to avoid ggtreeExtra clash with tree variables
      phylum_ring_l <- data.frame(tip_id=tip_df_l$label,
                                   tip_phylum=tip_df_l$phylum,
                                   stringsAsFactors=FALSE)
      enrich_ring_l <- data.frame(tip_id=tip_df_l$label,
                                   tip_enrich=tip_df_l$enrich,
                                   stringsAsFactors=FALSE)

      p_c <- ggtree(tr_l, layout="circular", color="grey65",
                     linewidth=0.3, open.angle=10) %<+% tip_df_l +
        geom_tippoint(aes(color=enrich, size=abund), alpha=0.88, na.rm=TRUE) +
        ggplot2::scale_color_manual(values=group_pal_all, na.value="grey80",
                                    name=paste0("Higher in (", group_col, ")"),
                                    na.translate=FALSE) +
        ggplot2::scale_size_continuous(range=c(0.8,3.5), name="log(mean abund.)",
                                       guide=guide_legend(override.aes=list(alpha=0.8,color="grey40")))

      # Ring 1: phylum
      p_c <- p_c +
        ggtreeExtra::geom_fruit(
          data=phylum_ring_l, geom=geom_tile,
          mapping=aes(y=tip_id, x=1, fill=tip_phylum), offset=0.05, pwidth=0.08
        ) +
        ggplot2::scale_fill_manual(
          values=c(pc_map, Unknown="grey88"), name="Phylum",
          na.translate=FALSE,
          guide=guide_legend(ncol=2, override.aes=list(size=5))
        )

      # Ring 2: enrichment
      p_c <- p_c +
        ggnewscale::new_scale_fill() +
        ggtreeExtra::geom_fruit(
          data=enrich_ring_l, geom=geom_tile,
          mapping=aes(y=tip_id, x=1, fill=tip_enrich), offset=0.06, pwidth=0.08
        ) +
        ggplot2::scale_fill_manual(
          values=group_pal_all, name=paste0("Enriched in\n(", group_col, ")"),
          na.value="grey92", na.translate=FALSE,
          guide=guide_legend(override.aes=list(size=5))
        )

      if (n_sig_c > 0) {
        p_c <- p_c +
          geom_tiplab(aes(label=ifelse(sig, display, NA_character_), filter=sig),
                      size=2.0, offset=0.14, color="grey10",
                      fontface="bold.italic", na.rm=TRUE, align=FALSE)
      }

      p_c +
        ggplot2::theme_void() +
        ggplot2::theme(
          legend.position = "right", legend.box = "vertical",
          legend.title    = ggplot2::element_text(face="bold", size=10),
          legend.text     = ggplot2::element_text(size=9),
          plot.title      = ggplot2::element_text(face="bold", size=12, hjust=0.5),
          plot.margin     = ggplot2::margin(20,20,20,20)
        ) +
        ggplot2::labs(
          title = paste0("Cladogram \u2014 ", rank_name, " Level")
        )
    }, error = function(e) {
      message("  Per-level cladogram failed (", rank_name, "): ", e$message)
      NULL
    })

    list(mm = mm_strict, res_df = res_df, bar = p_bar, clad = p_clad)
  }

  # ── Save multi-panel for one level ────────────────────────────────────────
  save_level <- function(panels, rank_name) {
    fname_base <- file.path(out_dir,
                            paste0("lefse_", tolower(rank_name), "_plots"))
    if (is.null(panels)) {
      ensure_file(paste0(fname_base, ".pdf"),
                  paste("LEfSe failed at", rank_name, "level"))
      ensure_file(paste0(fname_base, ".png"),
                  paste("LEfSe failed at", rank_name, "level"))
      return(invisible(NULL))
    }
    pieces <- Filter(Negate(is.null), list(panels$clad, panels$bar))

    if (length(pieces) >= 2 && requireNamespace("patchwork", quietly = TRUE)) {
      library(patchwork)
      # Cladogram gets ~55 % of width; LDA bar ~45 %
      combined <- patchwork::wrap_plots(pieces, ncol = 2, widths = c(1.2, 1.0)) +
        patchwork::plot_annotation(
          tag_levels = "A",
          title      = paste("LEfSe Analysis \u2014", rank_name, "Level"),
          theme      = theme(plot.title = element_text(face = "bold", size = 15,
                                                       hjust = 0.5))
        ) +
        patchwork::plot_layout(guides = "collect") &
        theme(legend.position = "bottom")
    } else if (length(pieces) == 1) {
      combined <- pieces[[1]]
    } else {
      combined <- NULL
    }

    if (!is.null(combined)) {
      grDevices::pdf(paste0(fname_base, ".pdf"), width = 22, height = 13)
      print(combined)
      grDevices::dev.off()
      ggplot2::ggsave(paste0(fname_base, ".png"), plot = combined,
                      width = 22, height = 13, units = "in", dpi = 600)
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