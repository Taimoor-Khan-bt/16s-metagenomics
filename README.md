# 16S rRNA Metagenomics Pipeline

An automated pipeline for analyzing 16S amplicon sequencing data. Processes raw FASTQ files through quality control, ASV inference, taxonomic classification, diversity analysis, and generates publication-ready figures.

**Author**: Taimoor Khan  
**Contact**: taimoorkhan007.tk@gmail.com  
**Version**: 2.2.0

---

## What This Pipeline Does

- Quality filtering and ASV inference using DADA2
- Taxonomic classification with SILVA database
- Alpha and beta diversity analysis with statistical testing
- Differential abundance testing (DESeq2)
- Phylogenetic tree construction
- Group comparisons across metadata variables
- Correlation analysis with continuous variables
- Publication-ready figures (600 DPI TIFF/PDF)

---

## Installation

### 1. Install R and Required Packages

Ensure R ≥4.3.0 is installed. The pipeline will check and install required packages automatically on first run.

**Required R packages**: dada2, phyloseq, DESeq2, ggplot2, ggtree, vegan, ape, DECIPHER, phangorn

### 2. Download SILVA Reference Files

Download SILVA v138.1 taxonomic databases:

```bash
# Download training set
wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz
gunzip silva_nr99_v138.1_train_set.fa.gz

# Download species assignment database
wget https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz
gunzip silva_species_assignment_v138.1.fa.gz
```

---

## Setup

### Input Data

Place paired-end FASTQ files in a directory:

```
your_project/
├── raw_reads/
│   ├── Sample1_R1.fastq.gz
│   ├── Sample1_R2.fastq.gz
│   ├── Sample2_R1.fastq.gz
│   └── Sample2_R2.fastq.gz
└── metadata.csv
```

### Metadata File

Create `metadata.csv` with sample information:

```csv
SampleID,Group,Age,Sex
Sample1,Control,25,M
Sample2,Treatment,30,F
```

**Required**: `SampleID` column matching FASTQ filenames  
**Required**: `Group` column for primary comparisons  
**Optional**: Additional columns for secondary analyses

### Configuration File

Copy and edit a template from `config/`:

```bash
cp config/kmu_config.yaml config/my_project.yaml
```

Edit key parameters:

```yaml
project:
  name: "MyProject"
  cohort_name: "MyAnalysis"

io:
  input_dir: "raw_reads"           # FASTQ directory
  output_dir: "output"
  metadata_csv: "metadata.csv"

```yaml
amplicon:
  taxonomy:
    train_set: "silva_nr99_v138.1_train_set.fa"
    species_db: "silva_species_assignment_v138.1.fa"
```

metadata:
  primary_comparison:
    group_column: "Group"
  secondary_comparisons:           # Optional
    - "Sex"
    - "AgeCategory"

analysis:
  mode: "comprehensive"            # simple | standard | comprehensive
```

---

## Running the Pipeline

Execute the pipeline:

```bash
./run_pipeline.sh config/my_project.yaml
```

The pipeline runs these steps automatically:
1. Initial quality control (FastQC on raw reads)
2. Adapter/primer trimming (Cutadapt)
3. Post-trim quality control (FastQC)
4. Quality filtering and ASV inference (DADA2)
5. Taxonomic classification (SILVA)
6. Phylogenetic tree construction
7. Diversity calculations
8. Statistical testing
9. Visualization generation
10. Aggregated QC report (MultiQC)

**Runtime**: 15 minutes (10 samples) to 3 hours (100+ samples)

---

## Analysis Modes

### Simple Mode
Basic two-group comparison. Outputs alpha/beta diversity and composition plots.

```yaml
analysis:
  mode: "simple"
```

### Standard Mode
Multiple categorical variables. Generates separate analyses for each grouping variable.

```yaml
analysis:
  mode: "standard"
metadata:
  secondary_comparisons:
    - "Sex"
    - "AgeCategory"
```

### Comprehensive Mode
Full analysis including differential abundance, correlations with continuous variables, group comparisons at multiple taxonomic ranks, and phylogenetic trees.

```yaml
analysis:
  mode: "comprehensive"
metadata:
  secondary_comparisons:
    - "Sex"
    - "AgeCategory"
  continuous_variables:
    - "Age"
    - "BMI"
analysis:
  differential_abundance:
    enabled: true
  group_comparisons:
    enabled: true
```

---

## Outputs

```
output/
├── fastqc_raw/                  # FastQC reports for raw reads
├── trimmed/                     # Adapter-trimmed FASTQ files
├── fastqc_trimmed/              # FastQC reports for trimmed reads
├── multiqc/                     # Aggregated QC report (HTML)
└── MyAnalysis/                  # Cohort-specific results
    ├── analysis/                # Data files
    │   ├── phyloseq_rarefied.rds    # Main phyloseq object
    │   ├── alpha_diversity.csv      # Diversity metrics
    │   ├── deseq2_*.csv             # Differential abundance results
    │   ├── correlation_*.csv        # Correlation results
    │   └── read_tracking.csv        # QC tracking
    └── visualizations/              # Figures (TIFF + PDF)
        ├── alpha_*.tiff             # Alpha diversity boxplots
        ├── beta_*.tiff              # PCoA ordinations
        ├── composition_*.tiff       # Taxonomic barplots
        ├── heatmap_*.tiff           # Clustered heatmaps
        ├── phylo_tree_*.tiff        # Phylogenetic trees
        ├── volcano_*.tiff           # Differential abundance
        └── correlation_*.tiff       # Taxa-variable correlations
```

**Key files**:
- `multiqc/multiqc_report.html`: Start here - comprehensive QC summary
- `phyloseq_rarefied.rds`: Load into R for custom analyses
- `alpha_diversity.csv`: Diversity metrics with statistics
- `deseq2_*.csv`: Differentially abundant taxa
- `correlation_*.csv`: Taxa-variable correlations

---

## Statistical Methods

**Alpha diversity**: Wilcoxon (2 groups) or Kruskal-Wallis (3+ groups)  
**Beta diversity**: PERMANOVA with 999 permutations  
**Differential abundance**: DESeq2 with FDR correction  
**Correlations**: Spearman correlation with FDR correction  
**Group comparisons**: Non-parametric tests with p-value annotations

All p-values are FDR-corrected for multiple testing.

---

## Troubleshooting

**Error: "Cannot find metadata column 'SampleID'"**  
→ Ensure metadata has `SampleID` column matching FASTQ filenames

**Error: "Reference database not found"**  
→ Check SILVA file paths in config are correct

**Warning: "Low read counts after filtering"**  
→ Review `read_tracking.csv` and adjust filtering parameters in config

**Memory error during tree construction**  
→ Set `build_tree: false` in config or reduce `max_tips`

**No differential abundance results**  
→ Check you have at least 3 samples per group and enabled `differential_abundance` in config

---

## Citation

```
Khan T. (2025). 16S rRNA Metagenomics Pipeline (v2.2.0). 
GitHub: https://github.com/Taimoor-Khan-bt/16s-metagenomics
```

**Core tools**:
- DADA2: Callahan et al. (2016) Nat Methods 13:581-583
- phyloseq: McMurdie & Holmes (2013) PLoS ONE 8(4):e61217
- DESeq2: Love et al. (2014) Genome Biol 15:550
- SILVA: Quast et al. (2013) Nucleic Acids Res 41:D590-D596

---

## License

MIT License - see `LICENSE` file