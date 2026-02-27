# =============================================================================
# rules/picrust.smk — G: Functional prediction using PICRUSt2
# =============================================================================
# PICRUSt2 requires its own conda environment (incompatible with most envs).
#
# Install PICRUSt2:
#   mamba create -n picrust2 -c conda-forge -c bioconda picrust2 -y
#
# PICRUSt2 predicts gene families (EC numbers) and MetaCyc pathways from
# 16S rRNA amplicon data using a reference phylogenetic tree placement.
# Input: rep_seqs FASTA + feature-table BIOM (raw counts, NOT rarefied)

_PICRUST = f"{OUT}/picrust2"

rule picrust2_pipeline:
    """
    PICRUSt2 full pipeline:
      1. Place ASV sequences onto reference tree
      2. Hidden state prediction for EC + KO gene families
      3. MetaCyc pathway abundance prediction

    Output files:
      EC_metagenome/ec_pred_metagenome_unstrat.tsv.gz
      KO_metagenome/ko_pred_metagenome_unstrat.tsv.gz
      pathways_out/path_abun_unstrat.tsv.gz
    """
    input:
        fasta = f"{OUT}/exported/rep_seqs/dna-sequences.fasta",
        biom  = f"{OUT}/exported/feature_table/feature-table.biom",
    output:
        pathways = f"{_PICRUST}/pathways_out/path_abun_unstrat.tsv.gz",
        ec       = f"{_PICRUST}/EC_metagenome/ec_pred_metagenome_unstrat.tsv.gz",
    params:
        out_dir  = _PICRUST,
        picrust_env = config["picrust"]["env"],
        threads  = config["dada2"]["threads"],
    log:
        f"{OUT}/logs/picrust2_pipeline.log",
    shell:
        """
        mkdir -p $(dirname {log})
        # Remove output dir if it exists (PICRUSt2 won't overwrite)
        rm -rf '{params.out_dir}'
        conda run --no-capture-output -n {params.picrust_env} \
            picrust2_pipeline.py \
                --study_fasta   '{input.fasta}' \
                --input_table   '{input.biom}' \
                --output        '{params.out_dir}' \
                --processes     {params.threads} \
                --verbose \
            2>&1 | tee {log}
        """


rule picrust2_decompress:
    """Decompress PICRUSt2 output TSVs for downstream R analysis."""
    input:
        pathways = f"{_PICRUST}/pathways_out/path_abun_unstrat.tsv.gz",
        ec       = f"{_PICRUST}/EC_metagenome/ec_pred_metagenome_unstrat.tsv.gz",
    output:
        pathways_tsv = f"{_PICRUST}/pathways_out/path_abun_unstrat.tsv",
        ec_tsv       = f"{_PICRUST}/EC_metagenome/ec_pred_metagenome_unstrat.tsv",
    shell:
        """
        gunzip -kf '{input.pathways}'
        gunzip -kf '{input.ec}'
        """


rule r_picrust2_differential:
    """
    R: Wilcoxon rank-sum test on predicted pathway abundances between groups.
    Outputs: TSV of differential pathways + barplot PDF.

    R packages: ggplot2, dplyr, tidyr
    """
    input:
        pathways = f"{_PICRUST}/pathways_out/path_abun_unstrat.tsv",
        metadata = config["metadata_file"],
    output:
        results = f"{_PICRUST}/pathway_differential.tsv",
        plots   = f"{_PICRUST}/pathway_plots.pdf",
    params:
        group_col = config["analysis"]["group_column"],
        out_dir   = _PICRUST,
    log:
        f"{OUT}/logs/r_picrust2_differential.log",
    shell:
        """
        mkdir -p $(dirname {log})
        Rscript workflow/scripts/picrust2_stats.R \
            '{input.pathways}' \
            '{input.metadata}' \
            '{params.group_col}' \
            '{params.out_dir}' \
            2>&1 | tee {log}
        """
