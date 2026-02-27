# =============================================================================
# rules/diversity.smk — Steps 8–10: Core diversity + alpha correlation
# =============================================================================

# ── Step 8: Core diversity metrics (phylogenetic) ──────────────────────────────
# core-metrics-phylogenetic writes to an --output-dir rather than individual files.
# Snakemake needs explicit output declarations; we declare the key outputs we
# will use downstream. The rule removes the output dir before running (QIIME 2
# will not overwrite an existing directory).

_CORE_DIV = f"{OUT}/core_diversity"
_EXP      = f"{OUT}/exported"

rule core_diversity:
    """
    Compute phylogenetic alpha and beta diversity metrics and Emperor PCoA plots.
    Rarefies the feature table to sampling_depth reads per sample (samples with
    fewer reads are dropped).

    Key outputs:
      *_vector.qza     — alpha diversity per sample
      *_emperor.qzv    — 3D PCoA visualization (requires ≥3 samples for Emperor)
      *_distance_matrix.qza
    """
    input:
        phylogeny = f"{OUT}/rooted_tree.qza",
        table     = f"{OUT}/table.qza",
        metadata  = config["metadata_file"],
        tsv       = f"{_EXP}/feature_table/feature-table.tsv",
    output:
        # Declare the outputs we use downstream (Snakemake watches these)
        rarefied_table      = f"{_CORE_DIV}/rarefied_table.qza",
        faith_pd            = f"{_CORE_DIV}/faith_pd_vector.qza",
        shannon             = f"{_CORE_DIV}/shannon_vector.qza",
        evenness            = f"{_CORE_DIV}/evenness_vector.qza",
        observed_features   = f"{_CORE_DIV}/observed_features_vector.qza",
        bray_curtis_matrix  = f"{_CORE_DIV}/bray_curtis_distance_matrix.qza",
        jaccard_matrix      = f"{_CORE_DIV}/jaccard_distance_matrix.qza",
        unweighted_matrix   = f"{_CORE_DIV}/unweighted_unifrac_distance_matrix.qza",
        weighted_matrix     = f"{_CORE_DIV}/weighted_unifrac_distance_matrix.qza",
        bray_curtis_emperor = f"{_CORE_DIV}/bray_curtis_emperor.qzv",
        jaccard_emperor     = f"{_CORE_DIV}/jaccard_emperor.qzv",
        unweighted_emperor  = f"{_CORE_DIV}/unweighted_unifrac_emperor.qzv",
        weighted_emperor    = f"{_CORE_DIV}/weighted_unifrac_emperor.qzv",
    params:
        docker         = DOCKER,
        sampling_depth = config["diversity"]["sampling_depth"],
        threads        = config["dada2"]["threads"],
        out_dir        = _CORE_DIV,
    log:
        f"{OUT}/logs/core_diversity.log",
    shell:
        """
        mkdir -p $(dirname {log})
        
        DEPTH="{params.sampling_depth}"
        if [ "$DEPTH" = "auto" ] || [ "$DEPTH" = "0" ]; then
            echo "Calculating minimum sequencing depth to avoid dropping any samples..."
            DEPTH=$(python3 -c "import pandas as pd; print(int(pd.read_csv('{input.tsv}', sep='\t', skiprows=1, index_col=0).sum(axis=0).min()))")
            echo "Auto sampling depth dynamic value: $DEPTH"
        fi

        # QIIME 2 refuses to write to an existing directory
        rm -rf '{params.out_dir}'
        {params.docker} qiime diversity core-metrics-phylogenetic \
            --i-phylogeny           '{input.phylogeny}' \
            --i-table               '{input.table}' \
            --p-sampling-depth      $DEPTH \
            --m-metadata-file       '{input.metadata}' \
            --p-n-jobs-or-threads   {params.threads} \
            --output-dir            '{params.out_dir}' \
            2>&1 | tee {log}
        """


# ── Step 10: Alpha diversity correlation (works with any sample size) ──────────

rule alpha_correlation:
    """
    Compute Spearman correlation between alpha diversity and numeric metadata
    variables. Works with any number of samples (unlike group-significance
    which needs ≥2 samples per group).
    """
    input:
        vector   = f"{_CORE_DIV}/{{metric}}_vector.qza",
        metadata = config["metadata_file"],
    output:
        viz = f"{_CORE_DIV}/{{metric}}_correlation.qzv",
    params:
        docker = DOCKER,
    log:
        f"{OUT}/logs/alpha_correlation_{{metric}}.log",
    shell:
        """
        mkdir -p $(dirname {log})
        {params.docker} qiime diversity alpha-correlation \
            --i-alpha-diversity '{input.vector}' \
            --m-metadata-file   '{input.metadata}' \
            --p-method          spearman \
            --o-visualization   '{output.viz}' \
            2>&1 | tee {log}
        """
