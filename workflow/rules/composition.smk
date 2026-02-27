# =============================================================================
# rules/composition.smk — C: Taxonomy collapse + relative abundance + F: KRONA
# =============================================================================
# Normalization: relative abundance (no rarefaction)

_COMP = f"{OUT}/composition"

# ── Taxonomy collapse at phylum / class / genus ───────────────────────────────

rule collapse_taxonomy:
    """
    Collapse feature table to a given taxonomic level.
    SILVA levels: 2=Phylum, 3=Class, 4=Order, 5=Family, 6=Genus, 7=Species
    """
    input:
        table    = f"{OUT}/table.qza",
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
        table    = f"{_COMP}/{{level}}_table.qza",
        taxonomy = f"{OUT}/taxonomy.qza",
        metadata = config["metadata_file"],
    output:
        viz = f"{_COMP}/{{level}}_barplot.qzv",
    params:
        docker = DOCKER,
    log:
        f"{OUT}/logs/taxa_barplot_{{level}}.log",
    shell:
        """
        mkdir -p $(dirname {log})
        {params.docker} qiime taxa barplot \
            --i-table         '{input.table}' \
            --i-taxonomy      '{input.taxonomy}' \
            --m-metadata-file '{input.metadata}' \
            --o-visualization '{output.viz}' \
            2>&1 | tee {log}
        """


# ── F: KRONA interactive plot ─────────────────────────────────────────────────

rule krona_plot:
    """
    Generate a KRONA interactive HTML plot from taxonomy + feature table.
    Requires: mamba install -n qiime2 -c bioconda krona
    """
    input:
        taxonomy_tsv = f"{OUT}/exported/taxonomy/taxonomy.tsv",
        table_tsv    = f"{OUT}/exported/feature_table/feature-table.tsv",
    output:
        html = f"{_COMP}/krona.html",
    params:
        out_dir = _COMP,
    log:
        f"{OUT}/logs/krona.log",
    shell:
        """
        mkdir -p $(dirname {log}) {params.out_dir}
        python3 workflow/scripts/make_krona_input.py \
            '{input.taxonomy_tsv}' \
            '{input.table_tsv}' \
            '{params.out_dir}/krona_input.txt' \
            2>&1 | tee {log}
        ktImportText \
            -o '{output.html}' \
            '{params.out_dir}/krona_input.txt' \
            2>&1 | tee -a {log}
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
        plots = f"{_COMP}/composition_plots.pdf",
    params:
        group_col  = config["analysis"]["group_column"],
        out_dir    = _COMP,
    log:
        f"{OUT}/logs/r_composition_plots.log",
    shell:
        """
        mkdir -p $(dirname {log})
        Rscript workflow/scripts/composition_plots.R \
            '{input.phylum_tsv}' \
            '{input.class_tsv}' \
            '{input.genus_tsv}' \
            '{input.metadata}' \
            '{params.group_col}' \
            '{output.plots}' \
            2>&1 | tee {log}
        """
