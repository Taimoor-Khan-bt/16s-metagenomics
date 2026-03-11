# =============================================================================
# rules/differential.smk — D: Differential abundance (ANCOM-BC + LEfSe)
# =============================================================================

_DIFF = f"{OUT}/differential"
_COMP = f"{OUT}/composition"

# Constrain 'level' to digits only so viz_ancombc and da_barplot_ancombc
# don't ambiguously match each other's output filenames (e.g. "2_barplot").
wildcard_constraints:
    level = r"\d+"

# ── ANCOM-BC via QIIME 2 ─────────────────────────────────────────────────────

rule ancombc:
    """
    ANCOM-BC differential abundance analysis for a collapsed taxonomy level.
    Adjusts for covariates via formula (e.g. "Stunting + AgeMonths + Sex").
    FDR correction (Benjamini-Hochberg) is applied internally by ANCOM-BC.

    Wildcard 'level': 2=Phylum, 6=Genus
    """
    input:
        table    = f"{_COMP}/{{level}}_table.qza",  # collapsed taxonomy table → readable feature IDs
        metadata = config["metadata_file"],
    output:
        differentials = f"{_DIFF}/ancombc_level{{level}}.qza",
    params:
        docker  = DOCKER,
        formula = config["analysis"]["ancombc_formula"],
        reference_levels = config["analysis"]["ancombc_reference_level"],
    log:
        f"{OUT}/logs/ancombc_level{{level}}.log",
    shell:
        """
        mkdir -p $(dirname {log})
        {params.docker} qiime composition ancombc \
            --i-table        '{input.table}' \
            --m-metadata-file '{input.metadata}' \
            --p-formula      '{params.formula}' \
            --p-reference-levels  {params.reference_levels} \
            --o-differentials '{output.differentials}' \
            --verbose \
            2>&1 | tee {log}
        """


rule viz_ancombc:
    """Tabulate ANCOM-BC differentials as an interactive QZV."""
    input:
        differentials = f"{_DIFF}/ancombc_level{{level}}.qza",
    output:
        viz = f"{OUT_VIZ}/differential/ancombc_level{{level}}.qzv",
    params:
        docker = DOCKER,
    log:
        f"{OUT}/logs/viz_ancombc_level{{level}}.log",
    shell:
        """
        mkdir -p $(dirname {log}) $(dirname {output.viz})
        {params.docker} qiime composition tabulate \
            --i-data          '{input.differentials}' \
            --o-visualization '{output.viz}' \
            2>&1 | tee {log}
        """


rule da_barplot_ancombc:
    """
    ANCOM-BC differential abundance bar chart (log fold change per feature).
    Input is the collapsed taxonomy table so feature IDs are already taxonomy
    strings. --p-level-delimiter ';' trims them to the lowest resolved rank.
    """
    input:
        differentials = f"{_DIFF}/ancombc_level{{level}}.qza",
    output:
        viz = f"{OUT_VIZ}/differential/ancombc_level{{level}}_barplot.qzv",
    params:
        docker = DOCKER,
    log:
        f"{OUT}/logs/da_barplot_ancombc_level{{level}}.log",
    shell:
        """
        mkdir -p $(dirname {log}) $(dirname {output.viz})
        {params.docker} qiime composition da-barplot \
            --i-data                '{input.differentials}' \
            --p-level-delimiter     ';' \
            --o-visualization       '{output.viz}' \
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
    """
    input:
        # Use filtered table (mito/chloro/unassigned removed) for LEfSe
        table_tsv = f"{OUT}/exported/feature_table_filtered/feature-table.tsv",
        taxonomy  = f"{OUT}/exported/taxonomy/taxonomy.tsv",
        metadata  = config["metadata_file"],
    output:
        results   = f"{_DIFF}/lefse_results.tsv",
        plots     = f"{OUT_VIZ}/differential/lefse_plots.pdf",
        plots_raw = f"{OUT_VIZ}/differential/lefse_plots_raw.pdf",
    params:
        group_col  = config["analysis"]["group_column"],
        strategy   = config.get("taxa_processing", {}).get("strategy", "rename"),
        dual_plots = "true" if config.get("taxa_processing", {}).get("generate_dual_plots", False) else "false",
        out_dir    = _DIFF,
        viz_dir    = f"{OUT_VIZ}/differential",
    log:
        f"{OUT}/logs/r_lefse.log",
    shell:
        """
        mkdir -p $(dirname {log}) '{params.viz_dir}'
        {RSCRIPT} workflow/scripts/lefse_analysis.R \
            '{input.table_tsv}' \
            '{input.taxonomy}' \
            '{input.metadata}' \
            '{params.group_col}' \
            '{params.out_dir}' \
            '{params.strategy}' \
            '{params.dual_plots}' \
            2>&1 | tee {log}
        mv '{params.out_dir}/lefse_plots.pdf'     '{output.plots}'
        mv '{params.out_dir}/lefse_plots_raw.pdf' '{output.plots_raw}'
        """
