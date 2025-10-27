# 16s-metagenomics – Release Notes

## v1.0.0 (2025-10-28)

First public release of the 16S rRNA metagenomics analysis pipeline.

Highlights
- Single entrypoint script: `16s-metagenomics-complete.R`
- Robust trimming with `cutadapt-16s-trim.sh` (paired-end, primer-aware) and detailed trimming stats
- DADA2-based denoising with relaxed, cutadapt-aware filtering (no hard truncLen; maxEE tuned; minLen=100)
- SILVA taxonomy assignment + optional species-level assignment
- Phyloseq objects saved (`phyloseq_object_raw.rds`, `phyloseq_rarefied.rds`)
- Contaminant filtering (non-Bacteria, Chloroplast, Mitochondria)
- Alpha/beta diversity with statistical testing (Kruskal–Wallis, pairwise Wilcoxon; PERMANOVA + dispersion)
- Publication-ready figures (TIFF, 600 dpi)
- Rarefaction curves and sequence-length summary to guide filtering

Inputs
- FASTQ files: `HF/*_1.fastq.gz`, `HF/*_2.fastq.gz`
- Metadata: `metadata.csv` with columns: Sample, Insect, Host, Extraction, Group

Key outputs
- Trimming: `HF/trimmed/trimming_statistics.tsv`, `HF/trimmed/cutadapt_summary.log`
- DADA2/QC: `output/filtering_summary.csv`, `output/read_tracking.csv`, `output/seq_length_summary.txt`
- Diversity: `output/alpha_diversity.csv`, `output/alpha_diversity_plot.tiff`, `output/beta_diversity_pcoa.tiff`
- Stats: `output/alpha_diversity_stats.txt`, `output/permanova_group.txt`, `output/betadisper_group.txt`
- Composition: `output/phylum_composition.tiff`, `output/phylum_composition_by_group.tiff`
- Rarefaction: `output/rarefaction_curves.pdf`
- Objects: `output/phyloseq_object_raw.rds`, `output/phyloseq_rarefied.rds`

Notes
- If you install `phangorn`, the pipeline can construct a phylogenetic tree and enable UniFrac.
- Optional differential abundance via ANCOM-BC will run automatically if the package is installed.
