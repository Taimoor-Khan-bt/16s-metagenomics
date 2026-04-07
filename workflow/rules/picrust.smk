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
    PICRUSt2 full pipeline. 
    Note: Output paths updated to match PICRUSt2's native naming convention.
    """
    input:
        fasta = f"{OUT}/exported/rep_seqs/dna-sequences.fasta",
        biom  = f"{OUT}/exported/feature_table/feature-table.biom",
    output:
        # PATH UPDATE: PICRUSt2 creates 'EC_metagenome_out' and 'pathways_out'
        pathways = f"{_PICRUST}/pathways_out/path_abun_unstrat.tsv.gz",
        ec       = f"{_PICRUST}/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz",
        nsti_gz  = f"{_PICRUST}/marker_predicted_and_nsti.tsv.gz",
    params:
        out_dir  = _PICRUST,
        picrust_env = config["picrust"]["env"],
        # FIX: Ensure threads is at least 1 to avoid joblib error
        threads  = config["dada2"]["threads"] if config["dada2"]["threads"] > 0 else 1
    log:
        f"{OUT}/logs/picrust2_pipeline.log",
    shell:
        """
        mkdir -p $(dirname {log})
        rm -rf '{params.out_dir}'
        
        mamba run -n {params.picrust_env} bash -c "picrust2_pipeline.py \
                --study_fasta {input.fasta} \
                --input {input.biom} \
                --output {params.out_dir} \
                --processes {params.threads} \
                --verbose" 2>&1 | tee {log}
        """

rule picrust2_decompress:
    """Decompress PICRUSt2 output TSVs for downstream R analysis."""
    input:
        pathways = f"{_PICRUST}/pathways_out/path_abun_unstrat.tsv.gz",
        ec       = f"{_PICRUST}/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz",
        nsti     = f"{_PICRUST}/marker_predicted_and_nsti.tsv.gz",
    output:
        pathways_tsv = f"{_PICRUST}/pathways_out/path_abun_unstrat.tsv",
        ec_tsv       = f"{_PICRUST}/EC_metagenome_out/pred_metagenome_unstrat.tsv",
        nsti_tsv     = f"{_PICRUST}/marker_predicted_and_nsti.tsv",
    shell:
        """
        gunzip -kf '{input.pathways}'
        gunzip -kf '{input.ec}'
        gunzip -kf '{input.nsti}'
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
        nsti_tsv = f"{_PICRUST}/marker_predicted_and_nsti.tsv",
    output:
        results               = f"{_PICRUST}/pathway_differential.tsv",
        nsti_summary          = f"{_PICRUST}/nsti_summary.tsv",
        plots                 = f"{OUT_VIZ}/picrust2/pathway_plots.pdf",
        plots_raw             = f"{OUT_VIZ}/picrust2/pathway_plots_raw.pdf",
        plots_png             = f"{OUT_VIZ}/picrust2/pathway_plots.png",
        plots_raw_png         = f"{OUT_VIZ}/picrust2/pathway_plots_raw.png",
        plots_heatmap_png     = f"{OUT_VIZ}/picrust2/pathway_plots_heatmap.png",
        plots_raw_heatmap_png = f"{OUT_VIZ}/picrust2/pathway_plots_raw_heatmap.png",
    params:
        group_col  = config["analysis"]["group_column"],
        strategy   = config.get("taxa_processing", {}).get("strategy", "rename"),
        dual_plots = "true" if config.get("taxa_processing", {}).get("generate_dual_plots", False) else "false",
        nsti_max   = config.get("picrust", {}).get("nsti_max", 2.0),
        out_dir    = _PICRUST,
        viz_dir    = f"{OUT_VIZ}/picrust2",
    log:
        f"{OUT}/logs/r_picrust2_differential.log",
    shell:
        """
        mkdir -p $(dirname {log}) '{params.viz_dir}'
        {RSCRIPT} workflow/scripts/picrust2_stats.R \
            '{input.pathways}' \
            '{input.metadata}' \
            '{params.group_col}' \
            '{params.out_dir}' \
            '{params.strategy}' \
            '{params.dual_plots}' \
            '{input.nsti_tsv}' \
            2>&1 | tee {log}
        mv '{params.out_dir}/pathway_plots.pdf'             '{output.plots}'
        mv '{params.out_dir}/pathway_plots_raw.pdf'         '{output.plots_raw}'
        mv '{params.out_dir}/pathway_plots.png'             '{output.plots_png}'
        mv '{params.out_dir}/pathway_plots_raw.png'         '{output.plots_raw_png}'
        mv '{params.out_dir}/pathway_plots_heatmap.png'     '{output.plots_heatmap_png}'
        mv '{params.out_dir}/pathway_plots_raw_heatmap.png' '{output.plots_raw_heatmap_png}'
        [ -f '{params.out_dir}/nsti_summary.tsv' ] && mv '{params.out_dir}/nsti_summary.tsv' '{output.nsti_summary}' || touch '{output.nsti_summary}'
        """
