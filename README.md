# 16S rRNA Metagenomics — KMU 2024

QIIME 2-based 16S rRNA Snakemake pipeline. QIIME 2 runs via Docker — no conda QIIME 2 installation needed.

---

## Requirements

| Tool | Version | Install |
|------|---------|---------|
| Docker | any | [docs.docker.com](https://docs.docker.com/get-docker/) |
| Snakemake | ≥9 | `mamba create -n qiime2 snakemake -y` |
| QIIME 2 image | 2024.10 | `docker pull quay.io/qiime2/amplicon:2024.10` |

---

## Quick Start (3 steps)

```bash
# 1. Setup — checks Docker, pulls image, installs Snakemake
bash setup.sh

# 2. Configure
#    a) Edit config/samples.tsv   — add your samples (sample_id + r1 path)
#    b) Edit config/metadata.tsv  — QIIME 2 metadata (#SampleID header)
#    c) Edit config/config.yaml   — set trunc_len, sampling_depth, classifier

# 3. Run
mamba activate qiime2
snakemake --cores all --snakefile workflow/Snakefile
```

---

## Pipeline Steps (24 jobs)

```
samples.tsv ──► build_manifest ──► import_reads ──► denoise_dada2
                                                         │
                              ┌──────────────────────────┤
                              │                          │
                        classify_taxonomy          build_tree
                              │                          │
                        taxa_barplot           core_diversity
                              │                     (17 outputs)
                         export_*           alpha_correlation×4
```

## Key Outputs (`output/qiime2_run/`)

| File | Description |
|------|-------------|
| `taxa_barplot.qzv` ⭐ | Interactive taxonomy barplot |
| `dada2_stats.qzv` | Denoising QC summary |
| `table_summary.qzv` | Per-sample read depths |
| `core_diversity/*_emperor.qzv` | 3D PCoA plots |
| `core_diversity/*_correlation.qzv` | Alpha diversity correlations |
| `exported/feature_table/feature-table.tsv` | ASV table (flat) |
| `exported/taxonomy/taxonomy.tsv` | Taxonomic assignments |
| `exported/tree/tree.nwk` | Phylogenetic tree (Newick) |

**View QZVs** → upload to [https://view.qiime2.org](https://view.qiime2.org)

---

## Useful Commands

```bash
# Dry run (no execution — validate config/rules)
snakemake -n --snakefile workflow/Snakefile

# Run only a specific step
snakemake output/qiime2_run/table.qza --snakefile workflow/Snakefile --cores all

# Force re-run a step
snakemake output/qiime2_run/taxonomy.qza --snakefile workflow/Snakefile --cores all --forcerun classify_taxonomy

# Visualize the DAG
snakemake --dag --snakefile workflow/Snakefile | dot -Tpng > dag.png

# Run with a different config
snakemake --cores all --snakefile workflow/Snakefile --configfile config/my_study_config.yaml
```

---

## Project Structure

```
workflow/
  Snakefile               ← entry point
  rules/
    common.smk            ← Docker helper
    import.smk            ← manifest + import
    denoise.smk           ← DADA2
    taxonomy.smk          ← SILVA classification
    phylogeny.smk         ← MAFFT + FastTree
    visualize.smk         ← summary QZVs
    diversity.smk         ← core metrics + alpha correlation  
    barplot.smk           ← taxonomy barplot
    export.smk            ← flat file exports
config/
  config.yaml             ← all parameters
  samples.tsv             ← sample_id | r1_path
  metadata.tsv            ← QIIME 2 metadata
data/
  raw/                    ← KMU Shahzad 2024 (30 samples)
  example/                ← 2-sample test set
  HF/                     ← HF cohort
refs/
  silva-138-99-nb-classifier.qza
docs/
  qiime2_manual_run.md    ← manual step-by-step reference
```

---

## Classifier Download

```bash
wget https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza \
     -O refs/silva-138-99-nb-classifier.qza
```