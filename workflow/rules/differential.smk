# =============================================================================
# rules/differential.smk — D: Differential abundance (ANCOM-BC + LEfSe)
# =============================================================================

_DIFF = f"{OUT}/differential"
_COMP = f"{OUT}/composition"

# ── ANCOM-BC via QIIME 2 ─────────────────────────────────────────────────────

rule ancombc:
    """
    ANCOM-BC differential abundance analysis for a collapsed taxonomy level.
    Adjusts for covariates via formula (e.g. "Stunting + AgeMonths + Sex").
    FDR correction (Benjamini-Hochberg) is applied internally by ANCOM-BC.

    Wildcard 'level': 2=Phylum, 6=Genus
    """
    input:
        table    = f"{_COMP}/{{level}}_table.qza",
        metadata = config["metadata_file"],
    output:
        differentials = f"{_DIFF}/ancombc_level{{level}}.qza",
    params:
        docker  = DOCKER,
        formula = config["analysis"]["ancombc_formula"],
    log:
        f"{OUT}/logs/ancombc_level{{level}}.log",
    shell:
        """
        mkdir -p $(dirname {log})
        {params.docker} qiime composition ancombc \
            --i-table        '{input.table}' \
            --m-metadata-file '{input.metadata}' \
            --p-formula      '{params.formula}' \
            --o-differentials '{output.differentials}' \
            --verbose \
            2>&1 | tee {log}
        """


rule viz_ancombc:
    """Tabulate ANCOM-BC differentials as an interactive QZV."""
    input:
        differentials = f"{_DIFF}/ancombc_level{{level}}.qza",
    output:
        viz = f"{_DIFF}/ancombc_level{{level}}.qzv",
    params:
        docker = DOCKER,
    log:
        f"{OUT}/logs/viz_ancombc_level{{level}}.log",
    shell:
        """
        mkdir -p $(dirname {log})
        {params.docker} qiime composition tabulate \
            --i-data          '{input.differentials}' \
            --o-visualization '{output.viz}' \
            2>&1 | tee {log}
        """


rule export_ancombc:
    """Export ANCOM-BC differentials to TSV."""
    input:
        differentials = f"{_DIFF}/ancombc_level{{level}}.qza",
        metadata      = config["metadata_file"],
    output:
        tsv = f"{_DIFF}/ancombc_level{{level}}.tsv",
    params:
        docker  = DOCKER,
        out_dir = f"{_DIFF}/ancombc_level{{level}}_export",
    log:
        f"{OUT}/logs/export_ancombc_level{{level}}.log",
    shell:
        """
        mkdir -p $(dirname {log})
        {params.docker} qiime tools export \
            --input-path  '{input.differentials}' \
            --output-path '{params.out_dir}' \
            2>&1 | tee {log}
        python3 workflow/scripts/merge_ancombc_csvs.py \
            '{params.out_dir}' \
            '{output.tsv}' \
            2>&1 | tee -a {log}
        """



# ── LEfSe via R (microbiomeMarker) ───────────────────────────────────────────

rule r_lefse:
    """
    LEfSe (Linear Discriminant Analysis Effect Size) differential analysis.
    Uses R package microbiomeMarker as the backend.
    Outputs: LDA score barplot + significant taxa table.

    R packages: microbiomeMarker, ggplot2, phyloseq, dplyr
    Install: mamba install -n qiime2 -c conda-forge -c bioconda r-microbiomemarker
             (or via BiocManager::install('microbiomeMarker') in R)
    """
    input:
        genus_tsv = f"{_COMP}/6_relfreq.tsv",
        taxonomy  = f"{OUT}/exported/taxonomy/taxonomy.tsv",
        metadata  = config["metadata_file"],
    output:
        results = f"{_DIFF}/lefse_results.tsv",
        plots   = f"{_DIFF}/lefse_plots.pdf",
    params:
        group_col = config["analysis"]["group_column"],
        out_dir   = _DIFF,
    log:
        f"{OUT}/logs/r_lefse.log",
    shell:
        """
        mkdir -p $(dirname {log})
        Rscript workflow/scripts/lefse_analysis.R \
            '{input.genus_tsv}' \
            '{input.taxonomy}' \
            '{input.metadata}' \
            '{params.group_col}' \
            '{params.out_dir}' \
            2>&1 | tee {log}
        """
