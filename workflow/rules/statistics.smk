# =============================================================================
# rules/statistics.smk — A: Alpha + B: Beta diversity statistics + E: Core
# =============================================================================

_STATS = f"{OUT}/stats"
_CORE  = f"{OUT}/core_microbiome"
_CDIV  = f"{OUT}/core_diversity"

# ── A: Alpha diversity QIIME 2 group significance ─────────────────────────────
# Works when group has ≥2 samples per category. If it fails, the R Wilcoxon
# fallback in r_alpha_stats covers the same test without sample-count limits.

rule alpha_group_significance:
    """
    Kruskal-Wallis group significance for each alpha metric (QIIME 2).
    Requires ≥2 samples per group — will error gracefully if not met.
    """
    input:
        vector   = f"{_CDIV}/{{metric}}_vector.qza",
        metadata = config["metadata_file"],
    output:
        viz = f"{OUT_VIZ}/diversity/{{metric}}_group_significance.qzv",
    params:
        docker    = DOCKER,
        group_col = config["analysis"]["group_column"],
    log:
        f"{OUT}/logs/alpha_group_sig_{{metric}}.log",
    shell:
        """
        mkdir -p $(dirname {log}) $(dirname {output.viz})
        {params.docker} qiime diversity alpha-group-significance \
            --i-alpha-diversity '{input.vector}' \
            --m-metadata-file   '{input.metadata}' \
            --o-visualization   '{output.viz}' \
            2>&1 | tee {log} || \
        echo "WARNING: alpha-group-significance failed (likely < 2 samples/group). See {log}" | tee -a {log}
        # Create empty output if QIIME 2 failed so Snakemake doesn't retry
        [ -f '{output.viz}' ] || touch '{output.viz}'
        """


# ── B: Beta diversity QIIME 2 PERMANOVA ───────────────────────────────────────

rule beta_group_significance:
    """
    PERMANOVA via QIIME 2 for each beta diversity distance matrix.
    Note: pairwise PERMANOVA does not adjust for covariates.
    Use r_beta_stats for covariate-adjusted vegan::adonis2.
    """
    input:
        distance = f"{_CDIV}/{{beta_metric}}_distance_matrix.qza",
        metadata = config["metadata_file"],
    output:
        viz = f"{OUT_VIZ}/diversity/{{beta_metric}}_permanova.qzv",
    params:
        docker    = DOCKER,
        group_col = config["analysis"]["group_column"],
    log:
        f"{OUT}/logs/beta_permanova_{{beta_metric}}.log",
    shell:
        """
        mkdir -p $(dirname {log}) $(dirname {output.viz})
        {params.docker} qiime diversity beta-group-significance \
            --i-distance-matrix     '{input.distance}' \
            --m-metadata-file       '{input.metadata}' \
            --m-metadata-column     '{params.group_col}' \
            --p-method              permanova \
            --p-pairwise \
            --o-visualization       '{output.viz}' \
            2>&1 | tee {log} || \
        echo "WARNING: beta-group-significance failed. See {log}" | tee -a {log}
        [ -f '{output.viz}' ] || touch '{output.viz}'
        """


# ── B: Export distance matrices to TSV for R analysis ────────────────────────

rule export_distance_matrix:
    """Export a QIIME 2 distance matrix QZA to a flat TSV for R."""
    input:
        matrix = f"{_CDIV}/{{beta_metric}}_distance_matrix.qza",
    output:
        tsv = f"{_STATS}/beta/{{beta_metric}}_distance_matrix.tsv",
    params:
        docker  = DOCKER,
        out_dir = f"{_STATS}/beta/{{beta_metric}}_dm_export",
    log:
        f"{OUT}/logs/export_dm_{{beta_metric}}.log",
    shell:
        """
        mkdir -p $(dirname {log})
        {params.docker} qiime tools export \
            --input-path  '{input.matrix}' \
            --output-path '{params.out_dir}' \
            2>&1 | tee {log}
        mv '{params.out_dir}/distance-matrix.tsv' '{output.tsv}'
        """


# ── A + B: R comprehensive statistical analysis ───────────────────────────────

rule r_alpha_stats:
    """
    R: Wilcoxon rank-sum test + GLM (covariate adjustment) for alpha diversity.
    Outputs: TSV statistics + multi-panel bar chart PDF + 600 DPI PNG.

    R packages: ggplot2, dplyr, tidyr, broom, patchwork
    Install: mamba install -n qiime2 -c conda-forge r-ggplot2 r-dplyr r-tidyr r-broom r-patchwork
    """
    input:
        alpha_tsvs = expand(
            f"{OUT}/exported/alpha_diversity/{{metric}}/alpha-diversity.tsv",
            metric=config["diversity"]["alpha_metrics"]
        ),
        metadata  = config["metadata_file"],
    output:
        stats         = f"{_STATS}/alpha/alpha_statistics.tsv",
        dunn_tsv      = f"{_STATS}/alpha/dunn_posthoc.tsv",
        plots         = f"{OUT_VIZ}/diversity/alpha_plots.pdf",
        plots_png     = f"{OUT_VIZ}/diversity/alpha_plots.png",
        plots_glm_png = f"{OUT_VIZ}/diversity/alpha_plots_glm.png",
    params:
        group_col  = config["analysis"]["group_column"],
        covariates = ",".join(config["analysis"]["covariates"]),
        alpha_dir  = f"{OUT}/exported/alpha_diversity",
        out_dir    = f"{_STATS}/alpha",
        viz_dir    = f"{OUT_VIZ}/diversity",
    log:
        f"{OUT}/logs/r_alpha_stats.log",
    shell:
        """
        mkdir -p $(dirname {log}) '{params.viz_dir}'
        {RSCRIPT} workflow/scripts/alpha_stats.R \
            '{input.metadata}' \
            '{params.alpha_dir}' \
            '{params.group_col}' \
            '{params.covariates}' \
            '{params.out_dir}' \
            2>&1 | tee {log}
        mv '{params.out_dir}/alpha_plots.pdf' '{output.plots}'
        mv '{params.out_dir}/alpha_plots.png' '{output.plots_png}'
        [ -f '{params.out_dir}/alpha_plots_glm.png' ] && mv '{params.out_dir}/alpha_plots_glm.png' '{output.plots_glm_png}' || touch '{output.plots_glm_png}'
        [ -f '{params.out_dir}/dunn_posthoc.tsv' ] && mv '{params.out_dir}/dunn_posthoc.tsv' '{output.dunn_tsv}' || touch '{output.dunn_tsv}'
        """



rule r_beta_stats:
    """
    R: CLR PCoA + PERMANOVA + Bray-Curtis dissimilarity boxplot.
    Outputs: PERMANOVA results TSV + multi-panel PCoA PDF + 600 DPI PNG.

    R packages: vegan, ggplot2, ape, dplyr, patchwork
    Install: mamba install -n qiime2 -c conda-forge r-vegan r-ggplot2 r-ape r-dplyr r-patchwork
    """
    input:
        table_tsv   = f"{OUT}/exported/feature_table_filtered/feature-table.tsv",
        metadata    = config["metadata_file"],
        bray_dm     = f"{_STATS}/beta/bray_curtis_distance_matrix.tsv",
        wunifrac_dm = f"{_STATS}/beta/weighted_unifrac_distance_matrix.tsv",
    output:
        permanova    = f"{_STATS}/beta/permanova_results.tsv",
        pairwise_tsv = f"{_STATS}/beta/permanova_pairwise.tsv",
        permdisp_png = f"{OUT_VIZ}/diversity/permdisp_plots.png",
        plots        = f"{OUT_VIZ}/diversity/pcoa_plots.pdf",
        plots_png    = f"{OUT_VIZ}/diversity/pcoa_plots.png",
    params:
        group_col  = config["analysis"]["group_column"],
        covariates = ",".join(config["analysis"]["covariates"]),
        out_dir    = f"{_STATS}/beta",
        viz_dir    = f"{OUT_VIZ}/diversity",
    log:
        f"{OUT}/logs/r_beta_stats.log",
    shell:
        """
        mkdir -p $(dirname {log}) '{params.viz_dir}'
        {RSCRIPT} workflow/scripts/beta_stats.R \
            '{input.table_tsv}' \
            '{input.metadata}' \
            '{input.bray_dm}' \
            '{input.wunifrac_dm}' \
            '{params.group_col}' \
            '{params.covariates}' \
            '{params.out_dir}' \
            2>&1 | tee {log}
        mv '{params.out_dir}/pcoa_plots.pdf' '{output.plots}'
        mv '{params.out_dir}/pcoa_plots.png' '{output.plots_png}'
        [ -f '{params.out_dir}/permdisp_plots.png' ] && mv '{params.out_dir}/permdisp_plots.png' '{output.permdisp_png}' || touch '{output.permdisp_png}'
        [ -f '{params.out_dir}/permanova_pairwise.tsv' ] && mv '{params.out_dir}/permanova_pairwise.tsv' '{output.pairwise_tsv}' || touch '{output.pairwise_tsv}'
        """


# ── E: Core microbiome ────────────────────────────────────────────────────────

rule core_features:
    """
    Identify core microbiome: ASVs present in ≥X% of samples.
    Prevalence threshold set in config (analysis.core_prevalence).
    """
    input:
        table = f"{OUT}/table.qza",
    output:
        viz = f"{OUT_VIZ}/core_microbiome/core_features.qzv",
    params:
        docker     = DOCKER,
        prevalence = config["analysis"]["core_prevalence"],
    log:
        f"{OUT}/logs/core_features.log",
    shell:
        """
        mkdir -p $(dirname {log}) $(dirname {output.viz})
        {params.docker} qiime feature-table core-features \
            --i-table                        '{input.table}' \
            --p-min-fraction                 {params.prevalence} \
            --o-visualization                '{output.viz}' \
            2>&1 | tee {log}
        """


rule r_core_microbiome:
    """
    R: Chi-square / Fisher exact test comparing core taxa between groups.
    Outputs: TSV stats table (with tax_level column) + prevalence barplot + UpSetR PDF.
    """
    input:
        table_tsv = f"{OUT}/exported/feature_table_filtered/feature-table.tsv",
        taxonomy  = f"{OUT}/exported/taxonomy/taxonomy.tsv",
        metadata  = config["metadata_file"],
    output:
        stats     = f"{_CORE}/core_stats.tsv",
        plots     = f"{OUT_VIZ}/core_microbiome/core_plots.pdf",
        plots_png = f"{OUT_VIZ}/core_microbiome/core_plots.png",
        upset_png = f"{OUT_VIZ}/core_microbiome/core_upset.png",
    params:
        group_col     = config["analysis"]["group_column"],
        prevalence    = config["analysis"]["core_prevalence"],
        tax_level     = config.get("analysis", {}).get("core_taxonomy_level", "ASV"),
        min_abundance = config.get("analysis", {}).get("core_min_abundance", 0.001),
        out_dir       = _CORE,
        viz_dir       = f"{OUT_VIZ}/core_microbiome",
    log:
        f"{OUT}/logs/r_core_microbiome.log",
    shell:
        """
        mkdir -p $(dirname {log}) '{params.viz_dir}'
        {RSCRIPT} workflow/scripts/core_microbiome.R \
            '{input.table_tsv}' \
            '{input.taxonomy}' \
            '{input.metadata}' \
            '{params.group_col}' \
            '{params.prevalence}' \
            '{params.out_dir}' \
            '{params.tax_level}' \
            '{params.min_abundance}' \
            2>&1 | tee {log}
        mv '{params.out_dir}/core_plots.pdf' '{output.plots}'
        [ -f '{params.out_dir}/core_plots.png' ] && mv '{params.out_dir}/core_plots.png' '{output.plots_png}' || touch '{output.plots_png}'
        [ -f '{params.out_dir}/core_upset.png' ] && mv '{params.out_dir}/core_upset.png' '{output.upset_png}' || touch '{output.upset_png}'
        """
