# 16S rRNA Metagenomics Pipeline

A reproducible, configuration-driven Snakemake pipeline for 16S rRNA amplicon
sequencing analysis. All QIIME 2 steps run inside the official Docker image —
no QIIME 2 conda installation required. Downstream statistical analysis and
visualizations are produced by modular R scripts.

> **Configuration first.** Every analysis parameter lives in `config/config.yaml`.
> Read `config/config_template.yaml` for a fully-annotated reference of every
> available option before editing your project config.

---

## Table of Contents

- [16S rRNA Metagenomics Pipeline](#16s-rrna-metagenomics-pipeline)
  - [Table of Contents](#table-of-contents)
  - [Requirements](#requirements)
  - [Quick Start](#quick-start)
  - [What the Pipeline Does](#what-the-pipeline-does)
  - [Pipeline Steps](#pipeline-steps)
    - [Step 1 — Import \& Denoise (DADA2)](#step-1--import--denoise-dada2)
    - [Step 2 — Taxonomic Classification](#step-2--taxonomic-classification)
    - [Step 3 — Phylogenetic Tree](#step-3--phylogenetic-tree)
    - [Step 4 — Alpha Diversity](#step-4--alpha-diversity)
    - [Step 5 — Beta Diversity](#step-5--beta-diversity)
    - [Step 6 — Taxonomic Composition](#step-6--taxonomic-composition)
    - [Step 7 — Differential Abundance](#step-7--differential-abundance)
    - [Step 8 — Core Microbiome](#step-8--core-microbiome)
    - [Step 9 — Functional Prediction (PICRUSt2)](#step-9--functional-prediction-picrust2)
    - [Step 10 — Tree Visualization (ggtree)](#step-10--tree-visualization-ggtree)
  - [Execution Commands](#execution-commands)
  - [Outputs](#outputs)
    - [Core Pipeline (`rule all`)](#core-pipeline-rule-all)
    - [Extended Analysis (`rule all_extended`)](#extended-analysis-rule-all_extended)
  - [Project Structure](#project-structure)
  - [Classifier Download](#classifier-download)
  - [Software Versions](#software-versions)

---

## Requirements

| Tool | Version | Install |
|------|---------|---------|
| Docker | any | [docs.docker.com](https://docs.docker.com/get-docker/) |
| Snakemake | ≥9 | `mamba create -n snakemake snakemake -y` |
| QIIME 2 image | 2024.10 | `docker pull quay.io/qiime2/amplicon:2024.10` |
| R | ≥4.3 | included in the QIIME 2 conda env or system R |
| PICRUSt2 | ≥2.5 | `mamba create -n picrust2 -c bioconda picrust2 -y` |

> **Preflight check:** Run `bash preflight.sh` to verify Docker, Snakemake, R
> packages, and disk space before starting.

---

## Quick Start

```bash
# 1. Clone and set up
git clone https://github.com/your-org/16s-metagenomics.git
cd 16s-metagenomics
bash setup.sh                  # checks Docker, pulls QIIME 2 image

# 2. Prepare inputs
#    a) config/samples.tsv    — one row per sample: sample_id | r1 | r2
#    b) config/metadata.tsv   — QIIME 2 metadata (#SampleID header)
#    c) config/config.yaml    — analysis parameters (see config_template.yaml)

# 3. Run
mamba activate snakemake
snakemake --cores all --snakefile workflow/Snakefile            # core pipeline
snakemake all_extended --cores all --snakefile workflow/Snakefile  # + full stats
```

---

## What the Pipeline Does

This pipeline processes raw paired-end 16S rRNA amplicon FASTQ files through
eleven analytical stages and produces publication-ready statistics and
visualizations.

```
Raw FASTQs
    │
    ▼
[1] Import + DADA2 denoising     → ASV table + representative sequences
    │
    ▼
[2] SILVA taxonomic classification → taxonomy assignments per ASV
    │
    ▼
[3] MAFFT alignment + FastTree   → rooted phylogenetic tree
    │
    ├──► [4] Alpha diversity       → richness/evenness + GLM-adjusted tests
    │
    ├──► [5] Beta diversity        → PCoA ordination + PERMANOVA
    │
    ├──► [6] Taxonomic composition → relative abundance at Phylum/Class/Genus
    │
    ├──► [7] Differential abundance → ANCOM-BC (primary) + LEfSe (exploratory)
    │
    ├──► [8] Core microbiome       → prevalence-based core + Fisher's Exact
    │
    ├──► [9] Functional prediction → PICRUSt2 MetaCyc pathways + KW test
    │
    └──► [10] Tree visualization   → annotated ggtree plots (5 styles)
```

Primary group comparison (e.g., case vs. control) is set via
`analysis.group_column` in `config/config.yaml`. Covariates listed under
`analysis.covariates` are included in all multivariable models (GLM, PERMANOVA,
ANCOM-BC), providing covariate-adjusted results.

---

## Pipeline Steps

### Step 1 — Import & Denoise (DADA2)

Raw paired-end FASTQ files are imported into QIIME 2 and denoised using DADA2.
DADA2 performs per-base quality filtering, chimera removal, read merging, and
produces an Amplicon Sequence Variant (ASV) table. Forward and reverse read
truncation lengths are set in `config.yaml` under `dada2.trunc_len_f` and
`dada2.trunc_len_r` — these should be chosen based on the quality drop-off in
your FastQC reports. Sequences classified as mitochondria, chloroplast, or
unassigned are removed from the filtered feature table used in all downstream
analyses.

**Key config options:** `dada2.trunc_len_f`, `dada2.trunc_len_r`,
`dada2.trim_left_f/r`, `dada2.threads`

**Run only denoising:**
```bash
snakemake {output_dir}/table.qza {output_dir}/rep_seqs.qza \
  --cores all --snakefile workflow/Snakefile
```

---

### Step 2 — Taxonomic Classification

ASVs are classified against the SILVA 138.1 reference database (99% identity)
using a pre-trained naive Bayes classifier (`classify-sklearn`). The classifier
file path is set in `config.yaml` under `classifier`. A QIIME 2 interactive
taxonomy barplot is generated for each taxonomy collapse level
(Phylum / Class / Genus by default).

**Key config options:** `classifier`, `analysis.collapse_levels`

**Run only classification:**
```bash
snakemake {output_dir}/taxonomy.qza --cores all --snakefile workflow/Snakefile
```

---

### Step 3 — Phylogenetic Tree

A multiple sequence alignment of representative ASV sequences is generated
using MAFFT. Hypervariable and uninformative columns are masked, and a
maximum-likelihood phylogenetic tree is inferred with FastTree 2. The tree is
midpoint-rooted to produce `rooted_tree.qza`, which is required for
phylogeny-aware diversity metrics (Faith's PD, UniFrac distances).

**Key config options:** `phylogeny.export_filtered_tree`

**Run only phylogeny:**
```bash
snakemake {output_dir}/rooted_tree.qza --cores all --snakefile workflow/Snakefile
```

---

### Step 4 — Alpha Diversity

Alpha diversity (within-sample richness and evenness) is computed from a
rarefied feature table. Four metrics are calculated by default:

| Metric | What it measures |
|--------|-----------------|
| Faith's PD | Phylogenetic diversity (tree branch length covered) |
| Shannon index | Species richness weighted by relative abundance |
| Pielou's evenness | How evenly reads are distributed across taxa |
| Observed ASVs | Raw count of unique ASVs per sample |

Between-group differences are tested with the Wilcoxon rank-sum test (2 groups)
or Kruskal-Wallis (3+ groups). A multivariable General Linear Model (GLM) of
the form `metric ~ group + covariate1 + covariate2 ...` is also fitted,
providing covariate-adjusted effect size estimates (β ± 95% CI), visualized
as a forest plot on page 2 of `alpha_plots.pdf`. All p-values are corrected
with Benjamini-Hochberg FDR (threshold: q < 0.05).

Rarefaction depth is auto-detected to the minimum library size that retains all
samples (set `diversity.sampling_depth: auto`), or specified manually.

**Key config options:** `diversity.sampling_depth`, `diversity.alpha_metrics`,
`analysis.group_column`, `analysis.covariates`

**Run alpha diversity + statistics:**
```bash
snakemake {output_dir}/stats/alpha/alpha_statistics.tsv \
  {viz_dir}/diversity/alpha_plots.pdf \
  --cores all --snakefile workflow/Snakefile
```

---

### Step 5 — Beta Diversity

Beta diversity (between-sample dissimilarity) is assessed using five distance
metrics: Bray-Curtis, Jaccard, unweighted UniFrac, weighted UniFrac, and
Aitchison (CLR-transformed). Principal Coordinates Analysis (PCoA) ordination
plots are produced for each metric and assembled into `pcoa_plots.pdf`.

Between-group differences are tested using PERMANOVA (vegan `adonis2`,
999 permutations) with the covariate-adjusted formula
`distance ~ group + covariate1 + covariate2 ...`. Results are written to
`stats/beta/permanova_results.tsv`.

**Key config options:** `diversity.beta_metrics`, `analysis.group_column`,
`analysis.covariates`

**Run beta diversity + PERMANOVA:**
```bash
snakemake {output_dir}/stats/beta/permanova_results.tsv \
  {viz_dir}/diversity/pcoa_plots.pdf \
  --cores all --snakefile workflow/Snakefile
```

---

### Step 6 — Taxonomic Composition

The feature table is collapsed to each taxonomy level in
`analysis.collapse_levels` (2 = Phylum, 3 = Class, 6 = Genus in SILVA notation).
Relative frequencies are computed per sample and exported as TSV files.
Stacked bar plots are generated per level showing group-stratified relative
abundance. Interactive Krona charts are also produced when `krona` is enabled
in the config.

**Key config options:** `analysis.collapse_levels`, `krona`

**Run composition:**
```bash
snakemake {viz_dir}/composition/composition_plots.pdf \
  --cores all --snakefile workflow/Snakefile
```

---

### Step 7 — Differential Abundance

Two complementary methods are used:

**ANCOM-BC** (primary, covariate-adjusted): Tests for differentially abundant
taxa at each level in `analysis.ancombc_levels`, correcting for sampling
fraction bias and including covariates via the formula in
`analysis.ancombc_formula`. The reference category for the primary group and
each covariate is set in `analysis.ancombc_reference_level`. Results are
exported as TSV files and visualized as bar plots.

**LEfSe** (exploratory): Linear discriminant analysis effect size, applied at
all taxonomy levels. An LDA score cutoff of 2.0 is used. LEfSe is
complementary — it does not control for covariates. Plots are saved in
`visualizations/differential/`.

**Key config options:** `analysis.ancombc_formula`,
`analysis.ancombc_reference_level`, `analysis.ancombc_levels`

**Run differential abundance:**
```bash
snakemake {viz_dir}/differential/lefse_plots.pdf \
  --cores all --snakefile workflow/Snakefile
```

---

### Step 8 — Core Microbiome

The core microbiome is defined as the set of ASVs present in ≥ N% of samples
within each group (threshold set by `analysis.core_prevalence`, default 50%).
Between-group differences in per-ASV core membership are tested using:
- Fisher's Exact Test for 2-group comparisons
- Fisher's Exact Test with Monte Carlo simulation (B = 2000) when expected cell
  counts are < 5 in 3+ group comparisons
- Chi-square where expected counts are sufficient

All p-values are BH-FDR corrected. Results are written to
`core_microbiome/core_stats.tsv` (includes a `test_used` column).
Two plots are produced: a prevalence bar plot (page 1) and an UpSetR
intersection diagram (page 2), saved in `core_plots.pdf`.

**Key config options:** `analysis.core_prevalence`

**Run core microbiome:**
```bash
snakemake {viz_dir}/core_microbiome/core_plots.pdf \
  {output_dir}/core_microbiome/core_stats.tsv \
  --cores all --snakefile workflow/Snakefile
```

---

### Step 9 — Functional Prediction (PICRUSt2)

PICRUSt2 predicts MetaCyc pathway abundances from the ASV table using
phylogenetic placement. Differentially abundant pathways between groups are
identified with the Kruskal-Wallis test followed by BH-FDR correction.
Results are written to `picrust2/pathway_differential.tsv` and visualized
in `pathway_plots.pdf`.

PICRUSt2 runs inside its own conda environment specified by `picrust.env`.

**Key config options:** `picrust.env`

**Run PICRUSt2:**
```bash
snakemake {viz_dir}/picrust2/pathway_plots.pdf \
  --cores all --snakefile workflow/Snakefile
```

---

### Step 10 — Tree Visualization (ggtree)

Five annotated phylogenetic tree plots are generated using ggtree in R:

| File | Description |
|------|-------------|
| `01_tree_basic.pdf` | Rectangular cladogram with genus labels, node size = abundance |
| `02_tree_abundance_heatmap.pdf` | Per-sample relative abundance ring heatmap |
| `03_tree_phylum_colorstrip.pdf` | Outer phylum colour strip |
| `04_tree_differential.pdf` | LEfSe LDA score bars overlaid on tree |
| `05_tree_circular.pdf` | Circular layout with all annotations combined |

**Key config options:** `phylogeny.export_filtered_tree`

**Run tree visualization:**
```bash
snakemake {viz_dir}/phylogeny/05_tree_circular.pdf \
  --cores all --snakefile workflow/Snakefile
```

---

## Execution Commands

Replace `{output_dir}` and `{viz_dir}` with the values set in your
`config/config.yaml` (e.g., `output/my_study` and `output/my_study/visualizations`).

```bash
# ── Activate environment ───────────────────────────────────────────────────────
mamba activate snakemake

# ── Validate config and DAG without running anything ─────────────────────────
snakemake -n --snakefile workflow/Snakefile
snakemake -n --snakefile workflow/Snakefile all_extended

# ── Full core pipeline (import → export) ─────────────────────────────────────
snakemake --cores all --snakefile workflow/Snakefile

# ── Full pipeline including all statistics and visualizations ─────────────────
snakemake all_extended --cores all --snakefile workflow/Snakefile

# ── Run with a custom config file ─────────────────────────────────────────────
snakemake all_extended --cores all --snakefile workflow/Snakefile \
  --configfile config/my_study_config.yaml

# ── Step-by-step (target individual outputs) ──────────────────────────────────

# Denoising only
snakemake output/my_study/table.qza output/my_study/rep_seqs.qza \
  --cores all --snakefile workflow/Snakefile

# Taxonomy only
snakemake output/my_study/taxonomy.qza \
  --cores all --snakefile workflow/Snakefile

# Diversity (alpha + beta) + statistics
snakemake output/my_study/stats/alpha/alpha_statistics.tsv \
  output/my_study/stats/beta/permanova_results.tsv \
  --cores all --snakefile workflow/Snakefile

# Differential abundance (ANCOM-BC + LEfSe)
snakemake output/my_study/visualizations/differential/lefse_plots.pdf \
  --cores all --snakefile workflow/Snakefile

# Core microbiome
snakemake output/my_study/core_microbiome/core_stats.tsv \
  --cores all --snakefile workflow/Snakefile

# PICRUSt2 functional prediction
snakemake output/my_study/visualizations/picrust2/pathway_plots.pdf \
  --cores all --snakefile workflow/Snakefile

# Tree visualization
snakemake output/my_study/visualizations/phylogeny/05_tree_circular.pdf \
  --cores all --snakefile workflow/Snakefile

# ── Utilities ─────────────────────────────────────────────────────────────────

# Force re-run a specific rule (even if output exists)
snakemake --forcerun r_alpha_stats \
  --cores all --snakefile workflow/Snakefile

# Visualize the full DAG as PNG
snakemake --dag --snakefile workflow/Snakefile | dot -Tpng > dag.png

# List all jobs that would run (dry run verbose)
snakemake -n --snakefile workflow/Snakefile --reason
```

---

## Outputs

All outputs are written to `output_dir` and `viz_dir` as set in `config.yaml`.

### Core Pipeline (`rule all`)

| Path | Description |
|------|-------------|
| `table.qza` | Raw ASV feature table |
| `table_filtered.qza` | Filtered table (no mito/chloro/unassigned) |
| `rep_seqs.qza` | Representative ASV sequences |
| `taxonomy.qza` | Taxonomic assignments (SILVA 138.1) |
| `rooted_tree.qza` | Rooted maximum-likelihood phylogenetic tree |
| `dada2_stats.qza` | Per-sample denoising statistics |
| `exported/feature_table/feature-table.tsv` | Flat ASV table |
| `exported/taxonomy/taxonomy.tsv` | Flat taxonomy table |
| `exported/tree/tree.nwk` | Newick tree |
| `visualizations/qc/dada2_stats.qzv` | Interactive denoising QC |
| `visualizations/taxonomy/taxa_barplot.qzv` | Interactive taxonomy barplot |

**View QZV files** → upload to [https://view.qiime2.org](https://view.qiime2.org)

### Extended Analysis (`rule all_extended`)

| Path | Description |
|------|-------------|
| `stats/alpha/alpha_statistics.tsv` | Wilcoxon/KW p-values + GLM coefficients |
| `visualizations/diversity/alpha_plots.pdf` | Violin plots + GLM forest plot |
| `stats/beta/permanova_results.tsv` | PERMANOVA results (all beta metrics) |
| `visualizations/diversity/pcoa_plots.pdf` | PCoA ordination plots |
| `composition/{level}_relfreq.tsv` | Relative abundances per taxonomy level |
| `visualizations/composition/composition_plots.pdf` | Stacked bar plots |
| `differential/ancombc_level{N}.tsv` | ANCOM-BC results per taxonomy level |
| `differential/lefse_results.tsv` | LEfSe LDA scores |
| `visualizations/differential/lefse_plots.pdf` | LEfSe cladogram + bar plot |
| `core_microbiome/core_stats.tsv` | Core ASV Fisher's Exact results |
| `visualizations/core_microbiome/core_plots.pdf` | Prevalence barplot + UpSetR |
| `picrust2/pathway_differential.tsv` | KW-tested MetaCyc pathways |
| `visualizations/picrust2/pathway_plots.pdf` | Top differential pathway plots |
| `visualizations/phylogeny/0[1-5]_tree_*.pdf` | Annotated ggtree plots |

---

## Project Structure

```
workflow/
  Snakefile                  ← pipeline entry point
  rules/
    common.smk               ← shared Docker/Rscript helpers
    import.smk               ← manifest building + QIIME 2 import
    denoise.smk              ← DADA2 denoising
    taxonomy.smk             ← SILVA classification
    phylogeny.smk            ← MAFFT alignment + FastTree + ggtree
    visualize.smk            ← QIIME 2 summary QZVs
    diversity.smk            ← core_diversity plugin + alpha correlations
    barplot.smk              ← QIIME 2 taxa barplot
    export.smk               ← flat file exports
    composition.smk          ← collapse + relative frequency + Krona
    statistics.smk           ← alpha/beta stats + core microbiome (R)
    differential.smk         ← ANCOM-BC + LEfSe (R)
    picrust.smk              ← PICRUSt2 functional prediction
  scripts/
    alpha_stats.R            ← Wilcoxon/KW tests + GLM forest plot
    beta_stats.R             ← PERMANOVA (adonis2) + PCoA plots
    composition_plots.R      ← stacked bar plots
    core_microbiome.R        ← core detection + Fisher's Exact + plots
    lefse_analysis.R         ← LEfSe wrapper + visualization
    picrust2_stats.R         ← KW test on MetaCyc pathways
    tree_plots.R             ← ggtree annotated tree plots
config/
  config.yaml                ← your active project config
  config_template.yaml       ← fully annotated reference template
  samples.tsv                ← sample_id | r1 | r2  (one sample per row)
  metadata.tsv               ← QIIME 2 metadata (#SampleID header)
refs/
  silva-138-99-nb-classifier.qza   ← pre-trained SILVA 138.1 classifier
docs/
  methodology_thesis.txt     ← detailed thesis-style Methods section
  methodology_publication.txt ← concise journal-style Methods section
  qiime2_manual_run.md       ← manual step-by-step QIIME 2 reference
```

---

## Classifier Download

The SILVA 138.1 naive Bayes classifier for QIIME 2 2024.10 must be downloaded
separately:

```bash
wget "https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza" \
     -O refs/silva-138-99-nb-classifier.qza
```

Place the file at the path specified in `config.yaml` under `classifier`.
Alternative classifiers (e.g., region-specific 515F/806R V4 classifiers) are
available at [https://docs.qiime2.org/2024.10/data-resources/](https://docs.qiime2.org/2024.10/data-resources/).

---

## Software Versions

| Software | Version | Citation |
|----------|---------|---------|
| QIIME 2 | 2024.10 | Bolyen et al., Nat Biotechnol (2019) |
| DADA2 | embedded | Callahan et al., Nat Methods (2016) |
| SILVA | 138.1 (99%) | Quast et al., Nucleic Acids Res (2013) |
| MAFFT | 7.x | Katoh & Standley, Mol Biol Evol (2013) |
| FastTree | 2.x | Price et al., PLoS ONE (2010) |
| ANCOM-BC | 2.x | Lin & Peddada, Nat Commun (2020) |
| LEfSe | 1.x | Segata et al., Genome Biol (2011) |
| PICRUSt2 | ≥2.5 | Douglas et al., Nat Biotechnol (2020) |
| vegan | ≥2.6 | Oksanen et al. (2022) |
| ggtree | ≥3.x | Yu et al., Mol Biol Evol (2017) |
| Snakemake | ≥9 | Mölder et al., F1000Research (2021) |