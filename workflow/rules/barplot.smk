# =============================================================================
# rules/barplot.smk — Step 9: Taxonomy barplot
# =============================================================================

rule taxa_barplot:
    """
    Interactive taxonomy barplot showing community composition at all 7
    taxonomic levels (Domain → Species). Upload the QZV to view.qiime2.org.
    """
    input:
        table    = f"{OUT}/table.qza",
        taxonomy = f"{OUT}/taxonomy.qza",
        metadata = config["metadata_file"],
    output:
        viz = f"{OUT}/taxa_barplot.qzv",
    params:
        docker = DOCKER,
    log:
        f"{OUT}/logs/taxa_barplot.log",
    shell:
        """
        mkdir -p $(dirname {log})
        {params.docker} qiime taxa barplot \
            --i-table             '{input.table}' \
            --i-taxonomy          '{input.taxonomy}' \
            --m-metadata-file     '{input.metadata}' \
            --o-visualization     '{output.viz}' \
            2>&1 | tee {log}
        """
