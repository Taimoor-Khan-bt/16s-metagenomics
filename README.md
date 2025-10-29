# 16S rRNA Metagenomics Analysis Pipeline

A production-ready pipeline for 16S rRNA amplicon sequencing analysis. Implements DADA2 workflow with comprehensive quality control, taxonomic classification, diversity analysis, and publication-ready visualizations.

**Author**: Taimoor Khan  
**Contact**: taimoorkhan007.tk@gmail.com  
**Version**: 2.1.0

---

## Purpose

This pipeline processes Illumina-sequenced 16S rRNA amplicons from raw FASTQ files through quality filtering, ASV inference, taxonomic assignment, diversity metrics, and statistical visualization. Outputs include phyloseq objects, diversity tables, and high-resolution publication-ready figures (600 DPI).

**Key Features**:
- DADA2-based ASV inference with quality filtering
- SILVA database taxonomic classification (genus and species level)
- Alpha diversity metrics with automated statistical testing
- Beta diversity analysis with PERMANOVA
- Phylogenetic tree construction for UniFrac distances
- Modular configuration via YAML
- Automated metadata validation

---

## Installation

### System Requirements

- **Operating System**: Linux or macOS
- **R**: Version 4.1 or higher
- **Memory**: Minimum 8 GB RAM (16 GB recommended for large datasets)
- **Disk Space**: 10 GB for reference databases plus dataset storage

### Step 1: Install R Dependencies

```r
# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c(
    "dada2",
    "phyloseq",
    "Biostrings",
    "DECIPHER",
    "ggtree"
))

# Install CRAN packages
install.packages(c(
    "phangorn",
    "vegan",
    "ggplot2",
    "dplyr",
    "tidyr",
    "cowplot",
    "viridisLite",
    "scales",
    "ape",
    "yaml",
    "ggpubr"  # Optional: for statistical annotations
))
```

### Step 2: Download SILVA Reference Databases

Download SILVA v138.1 training sets:

```bash
# Create reference directory
mkdir -p references

# Download taxonomy training set
wget -P references/ https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz
gunzip references/silva_nr99_v138.1_train_set.fa.gz

# Download species assignment database
wget -P references/ https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz
gunzip references/silva_species_assignment_v138.1.fa.gz
```

### Step 3: Clone Repository

```bash
git clone https://github.com/Taimoor-Khan-bt/16s-metagenomics.git
cd 16s-metagenomics
```

### Step 4: Verify Installation

```bash
Rscript -e "library(dada2); library(phyloseq); cat('Installation successful\n')"
```

---

## Setup and Configuration

### Input Data Structure

Place FASTQ files in a dedicated directory:

```
project_dir/
├── HF/                          # Input directory
│   ├── Sample1_R1.fastq.gz
│   ├── Sample1_R2.fastq.gz
│   ├── Sample2_R1.fastq.gz
│   └── Sample2_R2.fastq.gz
├── metadata.csv                 # Sample metadata
└── config/
    └── config.yaml              # Pipeline configuration
```

### Metadata File

Create `metadata.csv` with required columns:

| SampleID | Group     | Age | Sex |
|----------|-----------|-----|-----|
| Sample1  | Control   | 25  | M   |
| Sample2  | Treatment | 30  | F   |

**Required columns**:
- `SampleID`: Must match FASTQ filename prefixes
- `Group`: Experimental groups for comparison

**Optional columns**: SubjectID, Age, Sex, CollectionSite, SequencingType, Batch, LibraryPrep, dmft

See `docs/metadata_schema.csv` for complete specifications.

### Configuration File

Edit `config/config.yaml` to specify paths and parameters:

```yaml
project:
  name: "MyProject"
  cohort_name: "HF"
  sequencing_type: "16S"

io:
  input_dir: "HF"                    # Path to FASTQ directory
  output_dir: "output"
  metadata_csv: "metadata.csv"       # Optional

amplicon:
  taxonomy:
    train_set: "references/silva_nr99_v138.1_train_set.fa"
    species_db: "references/silva_species_assignment_v138.1.fa"
  
  phylogeny:
    build_tree: true                 # Set false to skip tree construction
    max_tips: 100                    # Maximum tips for tree visualization

plots:
  enable:
    alpha_diversity: true
    composition: true
    heatmap: true
    ordination: true
    phylo: true
  dpi: 600                           # Publication quality
  base_size: 14
```

---

## Usage

### Basic Execution

Run the complete pipeline:

```bash
Rscript scripts/runner.R --config config/config.yaml
```

### Pipeline Steps

The pipeline executes the following stages automatically:

1. **Quality Control**: Filter and trim reads based on quality scores
2. **Denoising**: DADA2 error learning and ASV inference
3. **Merging**: Pair-end read merging
4. **Chimera Removal**: Identify and remove chimeric sequences
5. **Taxonomic Assignment**: SILVA-based classification
6. **Phylogenetic Analysis**: Multiple sequence alignment and tree construction
7. **Diversity Metrics**: Calculate alpha and beta diversity
8. **Statistical Testing**: Automated group comparisons
9. **Visualization**: Generate publication-ready figures

### Execution Time

Expected runtime (varies with dataset size):
- Small dataset (10 samples, 50K reads each): ~15 minutes
- Medium dataset (50 samples, 100K reads each): ~1 hour
- Large dataset (100+ samples, 200K reads each): ~3-4 hours

---

## Expected Outputs

### Directory Structure

```
output/
└── HF/                              # Cohort-specific output
    ├── trimmed/                     # Quality-trimmed FASTQ
    ├── filtered/                    # DADA2 filtered FASTQ
    ├── analysis/                    # Core results
    │   ├── phyloseq_object_raw.rds
    │   ├── phyloseq_rarefied.rds
    │   ├── alpha_diversity.csv
    │   ├── read_tracking.csv
    │   ├── filtering_summary.csv
    │   ├── beta_bray.rds
    │   ├── beta_unifrac.rds
    │   ├── beta_wunifrac.rds
    │   ├── metadata_validated.csv
    │   └── sessionInfo.txt
    └── visualizations/              # Publication figures
        ├── alpha_shannon_boxplot.tiff
        ├── alpha_observed_boxplot.tiff
        ├── alpha_simpson_boxplot.tiff
        ├── composition_phylum.tiff
        ├── composition_class.tiff
        ├── composition_order.tiff
        ├── composition_family.tiff
        ├── composition_genus.tiff
        ├── heatmap_top_genus.tiff
        ├── beta_pcoa_bray.tiff
        ├── beta_pcoa_unifrac.tiff
        ├── phylo_tree_rectangular.tiff
        ├── phylo_tree_rectangular.pdf
        ├── phylo_tree_circular.tiff
        └── phylo_tree_circular.pdf
```

### Key Output Files

#### Analysis Files

| File | Description |
|------|-------------|
| `phyloseq_object_raw.rds` | Complete phyloseq object with ASV table, taxonomy, and metadata |
| `phyloseq_rarefied.rds` | Rarefied phyloseq object for diversity analyses |
| `alpha_diversity.csv` | Alpha diversity metrics (Observed, Shannon, Simpson, Chao1) |
| `read_tracking.csv` | Read counts through each processing step |
| `filtering_summary.csv` | Quality filtering statistics per sample |
| `beta_*.rds` | Distance matrices for ordination analyses |

#### Visualization Files

All figures generated at 600 DPI in TIFF format (some also in PDF):

- **Alpha Diversity**: Boxplots with statistical annotations (Wilcoxon/Kruskal-Wallis tests)
- **Composition**: Stacked barplots at taxonomic ranks (Phylum through Genus)
- **Heatmap**: Top 25 genera by abundance across samples
- **Beta Diversity**: PCoA ordination plots with group ellipses
- **Phylogenetic Trees**: Rectangular and circular layouts with abundance-scaled tips

---

## Statistical Analyses

### Alpha Diversity

Within-sample diversity metrics:
- **Observed ASVs**: Species richness
- **Shannon Index**: Richness and evenness combined
- **Simpson Index**: Dominance measure

Statistical comparisons:
- 2 groups: Wilcoxon rank-sum test
- 3+ groups: Kruskal-Wallis test
- P-values displayed on plots with significance symbols

### Beta Diversity

Between-sample dissimilarity:
- **Bray-Curtis**: Abundance-based dissimilarity
- **Jaccard**: Presence/absence-based dissimilarity
- **UniFrac**: Phylogenetic distance (requires tree)
- **Weighted UniFrac**: Abundance-weighted phylogenetic distance

Group separation assessed with PERMANOVA (999 permutations).

### Taxonomic Composition

Relative abundance analysis:
- Aggregated at multiple taxonomic ranks
- Top 10 most abundant taxa displayed
- Low-abundance taxa grouped as "Other"

---

## Troubleshooting

### Common Issues

**Error: "Cannot find metadata column 'Sample'"**
- Solution: Ensure metadata has `SampleID` column or set `metadata.id_column: "Sample"` in config

**Error: "Reference database not found"**
- Solution: Verify SILVA file paths in `config.yaml` and ensure files are decompressed

**Warning: "Low read counts after filtering"**
- Solution: Review quality profiles and adjust filtering parameters in config

**Memory error during tree construction**
- Solution: Reduce `max_tips` in config or set `build_tree: false`

### Getting Help

- Check `output/*/analysis/sessionInfo.txt` for R package versions
- Review `output/*/analysis/read_tracking.csv` for read loss at each step
- Open an issue: https://github.com/Taimoor-Khan-bt/16s-metagenomics/issues

---

## Citation

If you use this pipeline, please cite:

```
Khan T. (2025). 16S rRNA Metagenomics Analysis Pipeline (v2.1.0). 
GitHub repository: https://github.com/Taimoor-Khan-bt/16s-metagenomics
```

Please also cite the core tools:

- **DADA2**: Callahan et al. (2016) Nature Methods 13:581-583
- **phyloseq**: McMurdie & Holmes (2013) PLoS ONE 8(4):e61217
- **SILVA**: Quast et al. (2013) Nucleic Acids Research 41:D590-D596

---

## License

This project is licensed under the MIT License. See `LICENSE` file for details.

---

## Project Structure

```
16s-metagenomics/
├── config/
│   └── config.yaml              # Pipeline configuration
├── docs/
│   └── metadata_schema.csv      # Metadata specification
├── scripts/
│   ├── runner.R                 # Main execution script
│   ├── setup.R                  # Environment setup
│   ├── preprocess_16s.R         # Quality filtering
│   ├── analysis_16s.R           # DADA2 workflow and analysis
│   └── visualization.R          # Figure generation
├── references/                  # SILVA databases (user-added)
├── output/                      # Analysis results (generated)
├── README.md
├── LICENSE
└── RELEASE_v2.1.0.md           # Release notes
```