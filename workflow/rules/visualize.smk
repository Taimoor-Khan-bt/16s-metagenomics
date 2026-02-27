# =============================================================================
# rules/visualize.smk — Step 6: Summary QZV visualizations
# =============================================================================

rule viz_dada2_stats:
    """Tabulate DADA2 per-sample denoising stats (reads input/filtered/merged)."""
    input:
        stats = f"{OUT}/dada2_stats.qza",
    output:
        viz = f"{OUT}/dada2_stats.qzv",
    params:
        docker = DOCKER,
    log:
        f"{OUT}/logs/viz_dada2_stats.log",
    shell:
        """
        mkdir -p $(dirname {log})
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
        viz = f"{OUT}/table_summary.qzv",
    params:
        docker = DOCKER,
    log:
        f"{OUT}/logs/viz_table_summary.log",
    shell:
        """
        mkdir -p $(dirname {log})
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
        viz = f"{OUT}/rep_seqs.qzv",
    params:
        docker = DOCKER,
    log:
        f"{OUT}/logs/viz_rep_seqs.log",
    shell:
        """
        mkdir -p $(dirname {log})
        {params.docker} qiime feature-table tabulate-seqs \
            --i-data          '{input.rep_seqs}' \
            --o-visualization '{output.viz}' \
            2>&1 | tee {log}
        """
