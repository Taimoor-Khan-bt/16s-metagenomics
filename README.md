# Modular 16S Metagenomics Pipeline (R)

A cohesive, config-driven pipeline for 16S rRNA amplicon data. It performs DADA2 preprocessing, taxonomy, optional phylogeny, core diversity analyses, and publication-ready figures. Outputs are organized by cohort for clean multi-study work.

Author: Taimoor Khan  
Email: taimoorkhan007.tk@gmail.com

## Quick start

1) Place input FASTQs and (optionally) `metadata.csv` in a folder (e.g., `HF/`). Metadata is optional; if missing, defaults are used.

2) Edit the config at `config/config.yaml`:
- Set `project.sequencing_type: "16S"`
- Point `io.input_dir` to your folder (e.g., `HF`), optionally `io.cohort` (defaults to the basename of input_dir)
- For plots, toggle features under `plots.enable` (see below)

3) Run the orchestrator:
```bash
Rscript scripts/runner.R --config config/config.yaml
```

Optional: On a fresh machine, set `project.first_run: true` to let `scripts/setup.R` install R/Bioc packages and scaffold directories.

## What’s included

- Preprocessing (DADA2): filter/trim, learnErrors, denoise, merge, chimera removal
- Taxonomy (SILVA): assignTaxonomy, addSpecies
- Optional phylogeny: DECIPHER alignments → phangorn NJ/ML (GTR+G+I)
- Diversity: alpha (Observed, Shannon), beta (Bray PCoA + PERMANOVA)
- Composition: stacked barplots across ranks (Phylum … Genus)
- Heatmap: log-scaled top genera per sample
- Cohort-scoped outputs: `output/<cohort>/{trimmed,filtered,analysis,visualizations}`

## Configuration essentials (config/config.yaml)

project:
- sequencing_type: "16S"
- threads: number of CPU threads
- output_dir: base output folder (default `output`)
- first_run: true|false (install packages + scaffold on first run)

io:
- input_dir: folder holding FASTQs (e.g., `HF`)
- metadata_csv: path to sample metadata (optional)
- cohort: custom name for outputs; if null → basename(input_dir)

amplicon:
- taxonomy.silva_train_set / silva_species: paths to SILVA reference FASTAs
- phylogeny.build_tree: true|false; phylogeny.max_tips: limit tips for readable trees

plots:
- theme: classic|minimal|bw; color_palette: viridis|okabe-ito|tableau10
- ranks: which taxonomic ranks to plot (e.g., [Phylum, Class, Order, Family, Genus])
- top_taxa: N to show in barplots/heatmap ordering
- enable: turn plots on/off (see below)

Example toggles:
```yaml
plots:
  enable:
    alpha: true
    composition: true
    heatmap: true
    ordination: true
    tree_rectangular: true
    tree_circular: true
    per_sample_panels: false
```

## Outputs

```
output/
  <cohort>/
    trimmed/
    filtered/
    analysis/
      phyloseq_object_raw.rds
      phyloseq_rarefied.rds
      alpha_diversity.csv
      read_tracking.csv
      filtering_summary.csv
    visualizations/
      alpha_*.tiff
      composition_*.tiff
      heatmap_top_genera.tiff
      beta_pcoa_bray.tiff
      phylo_tree_{rectangular,circular}.{tiff,png,pdf}
      per_sample_panels.pdf (if enabled)
```

## Biological interpretation guide

- Alpha diversity (Observed, Shannon)
  - What: Within-sample richness (Observed) and richness + evenness (Shannon)
  - Readout: Higher values suggest more diverse communities; group comparisons via non-parametric tests
  - Caveats: Sensitive to depth; prefer consistent library sizes or rarefaction/offsets

- Beta diversity (Bray-Curtis PCoA + PERMANOVA)
  - What: Between-sample compositional dissimilarity and ordination for visualization
  - Readout: Group separation implies compositional differences; PERMANOVA tests explainable variance
  - Caveats: PERMANOVA assumes similar dispersion; check betadisper; compositionality can bias distances

- Composition barplots (ranks: Phylum→Genus)
  - What: Relative abundance of top taxa across samples/groups
  - Readout: Dominant clades and group-level shifts; use for descriptive summaries
  - Caveats: Relative data are compositional—avoid direct percent comparisons as “absolute” increases

- Heatmap (top genera)
  - What: Sample-by-taxon intensity map (log-transformed RA) for pattern discovery
  - Readout: Co-varying taxa, sample clusters, potential outliers
  - Caveats: Choice of top taxa can influence patterns; scale and transformation matter

- Phylogenetic trees (rectangular/circular)
  - What: Approximate evolutionary relationships among abundant ASVs/taxa
  - Readout: Clade-level structure; annotate with abundance for prominence
  - Caveats: Trees are inferred from short regions; topology should be treated as heuristic

- Rarefaction (internally for even-depth objects)
  - What: Subsamples reads to a common depth to compare alpha/beta fairly
  - Readout: Ensures fair within/between-sample comparisons when depths vary
  - Caveats: Discards data; consider depth-aware alternatives for inference

## Requirements

- R (>= 4.1 recommended)
- R/Bioc packages: dada2, phyloseq, DECIPHER, phangorn, vegan, ggplot2, dplyr, tidyr, ggtree (optional for trees)
- SILVA FASTAs: `silva_nr99_v138.1_train_set.fa`, `silva_species_assignment_v138.1.fa`

Tip: Set `project.first_run: true` to let `scripts/setup.R` install packages and scaffold the output structure.

## Running notes

- Metadata is optional. If missing, default Group="Unknown" is used for plot coloring. Matching is done by sample names derived from FASTQ basenames.
- Seeds: `project.random_seed` fixes ordination/rarefaction randomness for reproducibility.
- Performance: Increase `project.threads` for faster DECIPHER alignments and DADA2 steps.

## Project structure (current)

```
.
├── config/
│   └── config.yaml
├── scripts/
│   ├── runner.R              # orchestrator (entrypoint)
│   ├── setup.R               # optional first_run setup
│   ├── preprocess_16s.R      # DADA2 filtering/trim & tracking
│   ├── analysis_16s.R        # inference, taxonomy, phyloseq
│   └── visualization.R       # all figures; toggled by config
├── cutadapt-16s-trim.sh      # optional primer trimming script
├── HF/                       # example input folder
└── output/                   # cohort-scoped outputs
```

## License

MIT License (see `LICENSE`).

## Citation

If this pipeline aids your work, please cite the repository and include the software version/commit hash in Methods for reproducibility.