# =============================================================================
# rules/phylogeny.smk — Step 5: Phylogenetic tree (MAFFT + FastTree)
# =============================================================================

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
