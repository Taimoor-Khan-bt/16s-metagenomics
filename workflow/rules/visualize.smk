# =============================================================================
# rules/visualize.smk — Step 6: Summary QZV visualizations
# =============================================================================

# ── Core pipeline PNG plots (added to rule all) ────────────────────────────────
rule core_plots:
    """
    Generate all core QC/diversity/taxonomy PNGs from already-exported flat files.
    No QIIME 2 Docker needed — runs pure R after the export step.
    Outputs (600 DPI):
      qc/dada2_stats.png             read-retention waterfall
      qc/read_depth.png              per-sample total read count
      taxonomy/taxonomy_barplot.png    phylum stacked bars per sample
      taxonomy/top_taxa_abundance.png  top-20 phyla lollipop
    """
    input:
        dada2_stats   = f"{OUT}/exported/dada2_stats/stats.tsv",
        metadata      = config["metadata_file"],
        feature_table = f"{OUT}/exported/feature_table/feature-table.tsv",
        taxonomy      = f"{OUT}/exported/taxonomy/taxonomy.tsv",
        # Ensure all per-metric alpha TSVs exist before this rule starts
        alpha_tsvs    = expand(
            f"{OUT}/exported/alpha_diversity/{{metric}}/alpha-diversity.tsv",
            metric=config["diversity"]["alpha_metrics"],
        ),
        # Ensure all PCoA ordinations are exported before this rule starts
        beta_ordinations = expand(
            f"{OUT}/exported/beta_diversity/{{beta_metric}}/ordination.txt",
            beta_metric=config["diversity"]["beta_metrics"],
        ),
        # Bray-Curtis distance matrix for within/between dissimilarity boxplot
        bray_dm = f"{OUT}/exported/beta_diversity/bray_curtis/distance_matrix.tsv",
    output:
        dada2_png    = f"{OUT_VIZ}/qc/dada2_stats.png",
        depth_png    = f"{OUT_VIZ}/qc/read_depth.png",
        tax_bar_png  = f"{OUT_VIZ}/taxonomy/taxonomy_barplot.png",
        top_taxa_png = f"{OUT_VIZ}/taxonomy/top_taxa_abundance.png",
    params:
        rscript   = RSCRIPT,
        group_col = config["analysis"]["group_column"],
        viz_root  = OUT_VIZ,
        alpha_dir = f"{OUT}/exported/alpha_diversity",
        beta_dir  = f"{OUT}/exported/beta_diversity",
    log:
        f"{OUT}/logs/core_plots.log",
    shell:
        """
        mkdir -p "{params.viz_root}/qc" \
                 "{params.viz_root}/diversity" \
                 "{params.viz_root}/taxonomy" \
                 $(dirname {log})
        {params.rscript} workflow/scripts/core_plots.R \
            '{input.dada2_stats}' \
            '{params.alpha_dir}' \
            '{input.metadata}' \
            '{params.group_col}' \
            '{params.viz_root}' \
            '{input.feature_table}' \
            '{input.taxonomy}' \
            '{params.beta_dir}' \
            '{input.bray_dm}' \
            2>&1 | tee {log}
        """


rule viz_dada2_stats:
    """Tabulate DADA2 per-sample denoising stats (reads input/filtered/merged)."""
    input:
        stats = f"{OUT}/dada2_stats.qza",
    output:
        viz = f"{OUT_VIZ}/qc/dada2_stats.qzv",
    params:
        docker = DOCKER,
    log:
        f"{OUT}/logs/viz_dada2_stats.log",
    shell:
        """
        mkdir -p $(dirname {log}) $(dirname {output.viz})
        {params.docker} qiime metadata tabulate \
            --m-input-file    '{input.stats}' \
            --o-visualization '{output.viz}' \
            2>&1 | tee {log}
        """


rule viz_table_summary:
    """
    Feature table summary: per-sample read counts histogram
    and interactive depth selection for rarefaction.
    """
    input:
        table    = f"{OUT}/table.qza",
        metadata = config["metadata_file"],
    output:
        viz = f"{OUT_VIZ}/qc/table_summary.qzv",
    params:
        docker = DOCKER,
    log:
        f"{OUT}/logs/viz_table_summary.log",
    shell:
        """
        mkdir -p $(dirname {log}) $(dirname {output.viz})
        {params.docker} qiime feature-table summarize \
            --i-table                '{input.table}' \
            --m-sample-metadata-file '{input.metadata}' \
            --o-visualization        '{output.viz}' \
            2>&1 | tee {log}
        """


rule viz_rep_seqs:
    """
    Representative sequences visualization — BLAST-searchable in QIIME 2 View.
    """
    input:
        rep_seqs = f"{OUT}/rep_seqs.qza",
    output:
        viz = f"{OUT_VIZ}/qc/rep_seqs.qzv",
    params:
        docker = DOCKER,
    log:
        f"{OUT}/logs/viz_rep_seqs.log",
    shell:
        """
        mkdir -p $(dirname {log}) $(dirname {output.viz})
        {params.docker} qiime feature-table tabulate-seqs \
            --i-data          '{input.rep_seqs}' \
            --o-visualization '{output.viz}' \
            2>&1 | tee {log}
        """
