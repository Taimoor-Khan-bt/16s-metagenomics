# =============================================================================
# rules/export.smk — Steps 7, 11–12: Export artifacts to flat files
# =============================================================================

_CORE_DIV = f"{OUT}/core_diversity"
_EXP      = f"{OUT}/exported"

# ── Feature table: BIOM → TSV ─────────────────────────────────────────────────

rule export_feature_table:
    """
    Export the ASV feature table from QIIME 2 artifact to BIOM and convert
    to a human-readable tab-separated file.
    """
    input:
        table = f"{OUT}/table.qza",
    output:
        biom = f"{_EXP}/feature_table/feature-table.biom",
        tsv  = f"{_EXP}/feature_table/feature-table.tsv",
    params:
        docker  = DOCKER,
        out_dir = f"{_EXP}/feature_table",
    log:
        f"{OUT}/logs/export_feature_table.log",
    shell:
        """
        mkdir -p $(dirname {log})
        # Export BIOM artifact
        {params.docker} qiime tools export \
            --input-path  '{input.table}' \
            --output-path '{params.out_dir}' \
            2>&1 | tee {log}
        # Convert BIOM → TSV
        {params.docker} biom convert \
            -i '{output.biom}' \
            -o '{output.tsv}' \
            --to-tsv \
            2>&1 | tee -a {log}
        """


# ── Taxonomy TSV ──────────────────────────────────────────────────────────────

rule export_taxonomy:
    """Export taxonomy assignments to a TSV (Feature ID | Taxon | Confidence)."""
    input:
        taxonomy = f"{OUT}/taxonomy.qza",
    output:
        tsv = f"{_EXP}/taxonomy/taxonomy.tsv",
    params:
        docker  = DOCKER,
        out_dir = f"{_EXP}/taxonomy",
    log:
        f"{OUT}/logs/export_taxonomy.log",
    shell:
        """
        mkdir -p $(dirname {log})
        {params.docker} qiime tools export \
            --input-path  '{input.taxonomy}' \
            --output-path '{params.out_dir}' \
            2>&1 | tee {log}
        """


# ── Phylogenetic tree (Newick) ────────────────────────────────────────────────

rule export_tree:
    """Export the rooted phylogenetic tree in Newick format."""
    input:
        tree = f"{OUT}/rooted_tree.qza",
    output:
        nwk = f"{_EXP}/tree/tree.nwk",
    params:
        docker  = DOCKER,
        out_dir = f"{_EXP}/tree",
    log:
        f"{OUT}/logs/export_tree.log",
    shell:
        """
        mkdir -p $(dirname {log})
        {params.docker} qiime tools export \
            --input-path  '{input.tree}' \
            --output-path '{params.out_dir}' \
            2>&1 | tee {log}
        """


# ── Representative sequences (FASTA) ─────────────────────────────────────────

rule export_rep_seqs:
    """Export ASV representative sequences as FASTA."""
    input:
        rep_seqs = f"{OUT}/rep_seqs.qza",
    output:
        fasta = f"{_EXP}/rep_seqs/dna-sequences.fasta",
    params:
        docker  = DOCKER,
        out_dir = f"{_EXP}/rep_seqs",
    log:
        f"{OUT}/logs/export_rep_seqs.log",
    shell:
        """
        mkdir -p $(dirname {log})
        {params.docker} qiime tools export \
            --input-path  '{input.rep_seqs}' \
            --output-path '{params.out_dir}' \
            2>&1 | tee {log}
        """


# ── DADA2 statistics (TSV) ────────────────────────────────────────────────────

rule export_dada2_stats:
    """Export per-sample DADA2 denoising statistics as TSV."""
    input:
        stats = f"{OUT}/dada2_stats.qza",
    output:
        tsv = f"{_EXP}/dada2_stats/stats.tsv",
    params:
        docker  = DOCKER,
        out_dir = f"{_EXP}/dada2_stats",
    log:
        f"{OUT}/logs/export_dada2_stats.log",
    shell:
        """
        mkdir -p $(dirname {log})
        {params.docker} qiime tools export \
            --input-path  '{input.stats}' \
            --output-path '{params.out_dir}' \
            2>&1 | tee {log}
        """


# ── Alpha diversity vectors (TSV) ─────────────────────────────────────────────

rule export_alpha_diversity:
    """Export each alpha diversity vector as a TSV file (wildcard over metrics)."""
    input:
        vector = f"{_CORE_DIV}/{{metric}}_vector.qza",
    output:
        tsv = f"{_EXP}/alpha_diversity/{{metric}}/alpha-diversity.tsv",
    params:
        docker  = DOCKER,
        out_dir = f"{_EXP}/alpha_diversity/{{metric}}",
    log:
        f"{OUT}/logs/export_alpha_{{metric}}.log",
    shell:
        """
        mkdir -p $(dirname {log})
        {params.docker} qiime tools export \
            --input-path  '{input.vector}' \
            --output-path '{params.out_dir}' \
            2>&1 | tee {log}
        """
