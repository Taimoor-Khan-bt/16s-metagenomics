# =============================================================================
# rules/denoise.smk — Step 3: DADA2 denoising (error correction + ASV calling)
# =============================================================================

rule denoise_dada2:
    """
    Run DADA2 denoise-paired to correct errors, remove chimeras, merge PE reads, 
    and produce an ASV feature table, representative sequences, and denoising statistics.

    Key parameters (edit in config/config.yaml):
      dada2.trunc_len_f
      dada2.trunc_len_r
      dada2.trim_left_f
      dada2.trim_left_r
      dada2.threads
    """
    input:
        seqs = f"{OUT}/sequences.qza",
    output:
        table       = f"{OUT}/table.qza",
        rep_seqs    = f"{OUT}/rep_seqs.qza",
        stats       = f"{OUT}/dada2_stats.qza",
    params:
        docker      = DOCKER,
        trunc_len_f = config["dada2"]["trunc_len_f"],
        trunc_len_r = config["dada2"]["trunc_len_r"],
        trim_left_f = config["dada2"]["trim_left_f"],
        trim_left_r = config["dada2"]["trim_left_r"],
        threads     = config["dada2"]["threads"],
    log:
        f"{OUT}/logs/denoise_dada2.log",
    shell:
        """
        mkdir -p $(dirname {log})
        {params.docker} qiime dada2 denoise-paired \
            --i-demultiplexed-seqs '{input.seqs}' \
            --p-trunc-len-f {params.trunc_len_f} \
            --p-trunc-len-r {params.trunc_len_r} \
            --p-trim-left-f {params.trim_left_f} \
            --p-trim-left-r {params.trim_left_r} \
            --p-n-threads   {params.threads} \
            --o-table                    '{output.table}' \
            --o-representative-sequences '{output.rep_seqs}' \
            --o-denoising-stats          '{output.stats}' \
            --verbose \
            2>&1 | tee {log}
        """
