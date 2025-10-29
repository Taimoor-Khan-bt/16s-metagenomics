# 16s Metagenomics pipeline

Author: Taimoor Khan  
Email: taimoorkhan007.tk@gmail.com

## Requirements

- R (>= 4.1)
- Packages: dada2, phyloseq, Biostrings, DECIPHER, phangorn, vegan, ggplot2, dplyr, tidyr, ggtree (optional for trees)
- References: silva_nr99_v138.1_train_set.fa, silva_species_assignment_v138.1.fa

## Inputs

- FASTQ files under a folder (e.g., `HF/`)
- Optional `metadata.csv` (see Metadata schema)

## Configure

Edit `config/config.yaml`:
- `project.sequencing_type: "16S"`
- `io.input_dir: "HF"` (or your folder)
- `io.metadata_csv: "metadata.csv"` (optional)
- `amplicon.taxonomy.*` paths to SILVA files
- `amplicon.phylogeny.build_tree: true|false`
- plot toggles under `plots.enable`

## Run

```bash
Rscript scripts/runner.R --config config/config.yaml
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
      beta_bray.rds
      beta_unifrac.rds (if tree available)
      metadata_validated.csv
    visualizations/
      alpha_*.tiff
      composition_*.tiff
      heatmap_top_genera.tiff
      beta_pcoa_bray.tiff
      phylo_tree_{rectangular,circular}.{tiff,png,pdf}
```

## Metadata schema

Required columns:
- SampleID
- Group

Recommended optional columns:
- SubjectID, Age, Sex, CollectionSite, SequencingType, Batch, LibraryPrep, dmft

See `docs/metadata_schema.csv` for types and descriptions.

Notes:
- `metadata$id_column` and `metadata$group_column` in the config can map custom column names.
- The pipeline validates and aligns metadata to sample names; a validated copy is saved to `analysis/metadata_validated.csv`.

## Analyses and interpretation

- Alpha diversity (Observed, Shannon): within-sample richness and evenness; compare groups with non-parametric tests.
- Beta diversity (Bray; UniFrac if tree): between-sample dissimilarity; group separation assessed with PERMANOVA.
- Composition barplots (ranks): relative abundance across Phylum–Genus.
- Heatmap (top genera): log-transformed relative abundance by sample and taxon.
- Phylogeny (optional): tree from DECIPHER + phangorn for context and UniFrac.

## Project structure

```
.
├── config/
│   └── config.yaml
├── docs/
│   └── metadata_schema.csv
├── scripts/
│   ├── runner.R
│   ├── setup.R
│   ├── preprocess_16s.R
│   ├── analysis_16s.R
│   └── visualization.R
├── cutadapt-16s-trim.sh
└── output/
```

## License

MIT (see `LICENSE`).