# =============================================================================
# rules/phylogeny.smk — Step 5: Phylogenetic tree (MAFFT + FastTree)
# =============================================================================

_TREE_EXP = f"{OUT}/exported/tree"

rule build_tree:
    """
    Build a phylogenetic tree using MAFFT (alignment) and FastTree (tree).
    Both the unrooted and midpoint-rooted trees are produced.
    The rooted tree is required for phylogenetic diversity metrics (Faith's PD,
    UniFrac).
    """
    input:
        rep_seqs = f"{OUT}/rep_seqs.qza",
    output:
        alignment        = f"{OUT}/aligned_rep_seqs.qza",
        masked_alignment = f"{OUT}/masked_aligned_rep_seqs.qza",
        unrooted_tree    = f"{OUT}/unrooted_tree.qza",
        rooted_tree      = f"{OUT}/rooted_tree.qza",
    params:
        docker  = DOCKER,
        threads = config["dada2"]["threads"],
    log:
        f"{OUT}/logs/build_tree.log",
    shell:
        """
        mkdir -p $(dirname {log})
        {params.docker} qiime phylogeny align-to-tree-mafft-fasttree \
            --i-sequences        '{input.rep_seqs}' \
            --p-n-threads        {params.threads} \
            --o-alignment        '{output.alignment}' \
            --o-masked-alignment '{output.masked_alignment}' \
            --o-tree             '{output.unrooted_tree}' \
            --o-rooted-tree      '{output.rooted_tree}' \
            --verbose \
            2>&1 | tee {log}
        """


rule filter_rep_seqs:
    """
    Prune rep-seqs to only ASVs present in the filtered feature table
    (removes mitochondria / chloroplast / unassigned sequences).
    Used to export a clean Newick tree for ggtree visualisation.
    """
    input:
        rep_seqs = f"{OUT}/rep_seqs.qza",
        table    = f"{OUT}/table_filtered.qza",
    output:
        rep_seqs_filtered = f"{OUT}/rep_seqs_filtered.qza",
    params:
        docker = DOCKER,
    log:
        f"{OUT}/logs/filter_rep_seqs.log",
    shell:
        """
        mkdir -p $(dirname {log})
        {params.docker} qiime feature-table filter-seqs \
            --i-data         '{input.rep_seqs}' \
            --i-table        '{input.table}' \
            --o-filtered-data '{output.rep_seqs_filtered}' \
            2>&1 | tee {log}
        """


rule export_filtered_tree:
    """
    Export the rooted tree pruned to only filtered ASVs as a Newick file.
    This is the tree used for ggtree visualisation — no contaminant tips.
    """
    input:
        table_filtered = f"{OUT}/table_filtered.qza",
        rooted_tree    = f"{OUT}/rooted_tree.qza",
    output:
        nwk = f"{_TREE_EXP}/tree_filtered.nwk",
    params:
        docker  = DOCKER,
        out_dir = f"{_TREE_EXP}/filtered_export",
    log:
        f"{OUT}/logs/export_filtered_tree.log",
    shell:
        """
        mkdir -p $(dirname {log}) '{params.out_dir}'
        # Filter the rooted tree to only ASVs present in the filtered feature table
        {params.docker} qiime phylogeny filter-tree \
            --i-tree    '{input.rooted_tree}' \
            --i-table   '{input.table_filtered}' \
            --o-filtered-tree '{params.out_dir}/filtered_tree.qza' \
            2>&1 | tee {log}
        {params.docker} qiime tools export \
            --input-path  '{params.out_dir}/filtered_tree.qza' \
            --output-path '{params.out_dir}/nwk' \
            2>&1 | tee -a {log}
        cp '{params.out_dir}/nwk/tree.nwk' '{output.nwk}'
        """


rule r_tree_plots:
    """
    R/ggtree: 5 publication-quality phylogenetic tree visualisations.
      01_tree_basic.pdf              — rectangular cladogram with genus labels
      02_tree_abundance_heatmap.pdf  — per-sample relative abundance heatmap ring
      03_tree_phylum_colorstrip.pdf  — phylum colour strip
      04_tree_differential.pdf       — LEfSe LDA bar chart on tree
      05_tree_circular.pdf           — circular layout combining all annotations

    R packages: ape, ggtree, ggtreeExtra, ggplot2, dplyr, tidyr, scales,
                RColorBrewer, aplot, treeio  (auto-installed on first run)
    Install:  mamba install -n qiime2 -c conda-forge -c bioconda \\
                r-ape bioconductor-ggtree bioconductor-ggtreeextra r-ggplot2 \\
                r-dplyr r-tidyr r-scales r-rcolorbrewer r-aplot \\
                bioconductor-treeio
    """
    input:
        tree_nwk     = f"{_TREE_EXP}/tree_filtered.nwk",
        taxonomy_tsv = f"{OUT}/exported/taxonomy/taxonomy.tsv",
        table_tsv    = f"{OUT}/exported/feature_table_filtered/feature-table.tsv",
        metadata_tsv = config["metadata_file"],
        # LEfSe results are optional — script handles missing file gracefully
    output:
        plot01 = f"{OUT_VIZ}/phylogeny/01_tree_basic.pdf",
        plot02 = f"{OUT_VIZ}/phylogeny/02_tree_abundance_heatmap.pdf",
        plot03 = f"{OUT_VIZ}/phylogeny/03_tree_phylum_colorstrip.pdf",
        plot04 = f"{OUT_VIZ}/phylogeny/04_tree_differential.pdf",
        plot05 = f"{OUT_VIZ}/phylogeny/05_tree_circular.pdf",
    params:
        group_col = config["analysis"]["group_column"],
        lefse_tsv = f"{OUT}/differential/lefse_results.tsv",
        out_dir   = f"{OUT_VIZ}/phylogeny",
    log:
        f"{OUT}/logs/r_tree_plots.log",
    shell:
        """
        mkdir -p $(dirname {log}) '{params.out_dir}'
        {RSCRIPT} workflow/scripts/tree_plots.R \
            '{input.tree_nwk}' \
            '{input.taxonomy_tsv}' \
            '{input.table_tsv}' \
            '{input.metadata_tsv}' \
            '{params.group_col}' \
            '{params.lefse_tsv}' \
            '{params.out_dir}' \
            2>&1 | tee {log}
        """
