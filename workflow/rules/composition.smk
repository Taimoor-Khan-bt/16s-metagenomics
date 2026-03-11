# =============================================================================
# rules/composition.smk — C: Taxonomy collapse + relative abundance + F: KRONA
# =============================================================================
# Normalization: relative abundance (no rarefaction)

_COMP = f"{OUT}/composition"

# ── Krona directory + level config ──────────────────────────────────────────
_KRONA_VIZ         = f"{OUT_VIZ}/composition/krona"
_KRONA_LEVEL_NAMES = {"2": "phylum", "3": "class", "4": "order",
                      "5": "family", "6": "genus", "7": "species"}
_KRONA_LEVELS      = config.get("krona", {}).get("levels", ["2", "3", "6"])

# ── Taxonomy collapse at phylum / class / genus ───────────────────────────────

rule collapse_taxonomy:
    """
    Collapse feature table to a given taxonomic level.
    SILVA levels: 2=Phylum, 3=Class, 4=Order, 5=Family, 6=Genus, 7=Species
    """
    input:
        table    = f"{OUT}/table_filtered.qza",  # mito/chloro/unassigned removed
        taxonomy = f"{OUT}/taxonomy.qza",
    output:
        collapsed = f"{_COMP}/{{level}}_table.qza",
    params:
        docker = DOCKER,
        level  = lambda wc: wc.level,
    log:
        f"{OUT}/logs/collapse_taxonomy_{{level}}.log",
    shell:
        """
        mkdir -p $(dirname {log})
        {params.docker} qiime taxa collapse \
            --i-table           '{input.table}' \
            --i-taxonomy        '{input.taxonomy}' \
            --p-level           {params.level} \
            --o-collapsed-table '{output.collapsed}' \
            2>&1 | tee {log}
        """


rule relative_frequency:
    """Convert a collapsed table to relative abundance (row sums to 1.0)."""
    input:
        table = f"{_COMP}/{{level}}_table.qza",
    output:
        relfreq = f"{_COMP}/{{level}}_relfreq.qza",
    params:
        docker = DOCKER,
    log:
        f"{OUT}/logs/relfreq_{{level}}.log",
    shell:
        """
        mkdir -p $(dirname {log})
        {params.docker} qiime feature-table relative-frequency \
            --i-table             '{input.table}' \
            --o-relative-frequency-table '{output.relfreq}' \
            2>&1 | tee {log}
        """


rule export_relfreq:
    """Export relative frequency table to BIOM then convert to TSV."""
    input:
        relfreq = f"{_COMP}/{{level}}_relfreq.qza",
    output:
        biom = f"{_COMP}/{{level}}_relfreq.biom",
        tsv  = f"{_COMP}/{{level}}_relfreq.tsv",
    params:
        docker  = DOCKER,
        out_dir = lambda wc: f"{_COMP}/{wc.level}_relfreq_export",
    log:
        f"{OUT}/logs/export_relfreq_{{level}}.log",
    shell:
        """
        mkdir -p $(dirname {log})
        {params.docker} qiime tools export \
            --input-path  '{input.relfreq}' \
            --output-path '{params.out_dir}' \
            2>&1 | tee {log}

        mv '{params.out_dir}/feature-table.biom' '{output.biom}'

        {params.docker} biom convert \
            -i '{output.biom}' \
            -o '{output.tsv}' \
            --to-tsv \
            2>&1 | tee -a {log}
        """


rule viz_taxa_barplot_level:
    """QIIME 2 taxa barplot forced to a specific level (for QZV download)."""
    input:
        table    = f"{OUT}/table_filtered.qza",  # mito/chloro/unassigned removed
        taxonomy = f"{OUT}/taxonomy.qza",
        metadata = config["metadata_file"],
    output:
        viz = f"{OUT_VIZ}/composition/{{level}}_barplot.qzv",
    params:
        docker = DOCKER,
    log:
        f"{OUT}/logs/taxa_barplot_{{level}}.log",
    shell:
        """
        mkdir -p $(dirname {log}) $(dirname {output.viz})
        {params.docker} qiime taxa barplot \
            --i-table         '{input.table}' \
            --i-taxonomy      '{input.taxonomy}' \
            --m-metadata-file '{input.metadata}' \
            --o-visualization '{output.viz}' \
            2>&1 | tee {log}
        """


# ── F: KRONA interactive plots ────────────────────────────────────────────────

rule krona_multisample:
    """
    Multi-level, per-sample annotated Krona HTML plots (via build_krona.py).
    Produces one merged HTML (all samples) + one per collapsed taxonomy level.
    LEfSe LDA scores and alpha diversity are overlaid when available.
    Requires: mamba install -n qiime2 -c bioconda krona && ktUpdateTaxonomy.sh
    """
    input:
        taxonomy_tsv   = f"{OUT}/exported/taxonomy/taxonomy.tsv",
        table_tsv      = f"{OUT}/exported/feature_table_filtered/feature-table.tsv",
        collapse_tsvs  = expand(f"{_COMP}/{{level}}_relfreq.tsv", level=_KRONA_LEVELS),
        alpha_shannon  = f"{OUT}/exported/alpha_diversity/shannon/alpha-diversity.tsv",
        alpha_obs_feat = f"{OUT}/exported/alpha_diversity/observed_features/alpha-diversity.tsv",
    output:
        all_samples = f"{_KRONA_VIZ}/krona_all_samples.html",
        level_htmls = [f"{_KRONA_VIZ}/krona_level{l}_{_KRONA_LEVEL_NAMES.get(l, l)}.html"
                       for l in _KRONA_LEVELS],
    params:
        out_dir       = _KRONA_VIZ,
        alpha_dir     = f"{OUT}/exported/alpha_diversity",
        lefse_tsv     = f"{OUT}/differential/lefse_results.tsv",
        krona_env     = config.get("krona", {}).get("krona_env", ""),
        per_sample    = "--per-sample" if config.get("krona", {}).get("per_sample", True) else "",
        metadata_file = config["metadata_file"],
        collapse_args = " ".join(
            [f"--collapse {l}:{_KRONA_LEVEL_NAMES.get(l, l)}:{_COMP}/{l}_relfreq.tsv"
             for l in _KRONA_LEVELS]
        ),
    log:
        f"{OUT}/logs/krona_multisample.log",
    shell:
        """
        mkdir -p $(dirname {log}) '{params.out_dir}'
        python3 workflow/scripts/build_krona.py \
            --taxonomy   '{input.taxonomy_tsv}' \
            --table      '{input.table_tsv}' \
            --metadata   '{params.metadata_file}' \
            --alpha-dir  '{params.alpha_dir}' \
            --lefse      '{params.lefse_tsv}' \
            --out-dir    '{params.out_dir}' \
            --krona-env  '{params.krona_env}' \
            {params.per_sample} \
            {params.collapse_args} \
            2>&1 | tee {log}
        """


# ── R composition barplots ────────────────────────────────────────────────────

rule r_composition_plots:
    """
    Stacked and grouped barplots at phylum, class, and genus level.
    Outputs: PDF barplots per level.
    R packages required: ggplot2, dplyr, tidyr, scales, RColorBrewer
    """
    input:
        phylum_tsv = f"{_COMP}/2_relfreq.tsv",
        class_tsv  = f"{_COMP}/3_relfreq.tsv",
        genus_tsv  = f"{_COMP}/6_relfreq.tsv",
        metadata   = config["metadata_file"],
    output:
        plots     = f"{OUT_VIZ}/composition/composition_plots.pdf",
        plots_raw = f"{OUT_VIZ}/composition/composition_plots_raw.pdf",
    params:
        group_col  = config["analysis"]["group_column"],
        strategy   = config.get("taxa_processing", {}).get("strategy", "rename"),
        dual_plots = "true" if config.get("taxa_processing", {}).get("generate_dual_plots", False) else "false",
        out_dir    = _COMP,
        viz_dir    = f"{OUT_VIZ}/composition",
    log:
        f"{OUT}/logs/r_composition_plots.log",
    shell:
        """
        mkdir -p $(dirname {log}) '{params.viz_dir}'
        {RSCRIPT} workflow/scripts/composition_plots.R \
            '{input.phylum_tsv}' \
            '{input.class_tsv}' \
            '{input.genus_tsv}' \
            '{input.metadata}' \
            '{params.group_col}' \
            '{params.out_dir}' \
            '{params.strategy}' \
            '{params.dual_plots}' \
            2>&1 | tee {log}
        mv '{params.out_dir}/composition_plots.pdf'     '{output.plots}'
        mv '{params.out_dir}/composition_plots_raw.pdf' '{output.plots_raw}'
        """
