# 16S rRNA Gene Sequence Analysis Pipeline

A comprehensive pipeline for analyzing 16S rRNA gene sequences from metagenomic samples. This pipeline performs quality control, taxonomic assignment, diversity analysis, and visualization of microbial community data.

**Author:** Taimoor Khan  
**Email:** taimoorkhan007.tk@gmail.com

## Features

- Quality filtering and trimming of paired-end reads
- Error rate learning and sample inference using DADA2
- Chimera removal
- Taxonomic assignment using SILVA database
- Alpha diversity analysis (Observed, Shannon, InvSimpson)
- Beta diversity analysis with Bray-Curtis dissimilarity
- Core microbiome analysis
- Advanced taxonomic abundance analysis and visualization
- Publication-ready plots

## Prerequisites

- R (>= 4.0.0)
- Required R packages:
  - dada2
  - phyloseq
  - vegan
  - tidyverse
  - ape
  - DECIPHER
  - ggtree
  - RColorBrewer
  - ComplexHeatmap
  - indicspecies

## Step-by-Step Tutorial

### 1. Setup Environment

```bash
# Clone the repository
git clone [repository-url]
cd 16s-metagenomics

# Create necessary directories
mkdir -p rawFastq/trimmed output
```

### 2. Prepare Input Data

1. Place your raw paired-end fastq files in the `rawFastq` directory
2. Update primer sequences in `cutadapt-16s-trim.sh`:
   ```bash
   # Open the file
   nano cutadapt-16s-trim.sh

   # Update these lines with your primer sequences
   FWD="GTGYCAGCMGCCGCGGTAA"  # 515F primer (update as needed)
   REV="GGACTACNVGGGTWTCTAAT"  # 806R primer (update as needed)
   ```
   **Important:** Using the correct primer sequences is crucial for proper trimming. Make sure to:
   - Use the exact primers used in your sequencing
   - Include any adapters if present
   - Consider primer degeneracy
   - Check primer orientation

3. Run the trimming script:
   ```bash
   bash cutadapt-16s-trim.sh
   ```

### 3. Install Required R Packages

```R
# In R console
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install required packages
BiocManager::install(c(
  "dada2",
  "phyloseq",
  "DECIPHER",
  "phangorn"
))

install.packages(c(
  "tidyverse",
  "vegan",
  "RColorBrewer",
  "ComplexHeatmap",
  "indicspecies"
))
```

### 4. Prepare Reference Databases

1. Download SILVA taxonomy files:
   ```bash
   # Download SILVA files (v138.1)
   wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz
   wget https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz
   
   # Decompress files
   gunzip silva_nr99_v138.1_train_set.fa.gz
   gunzip silva_species_assignment_v138.1.fa.gz
   ```

2. Verify files are in the correct location:
   - `silva_nr99_v138.1_train_set.fa`
   - `silva_species_assignment_v138.1.fa`

### 5. Run the Analysis Pipeline

1. Review and modify parameters if needed:
   - Open `16s-metagenomics-complete.R`
   - Check filtering parameters in the DADA2 section
   - Adjust rarefaction depth if needed
   - Modify visualization parameters if desired

2. Run the pipeline:
   ```bash
   Rscript 16s-metagenomics-complete.R
   ```

3. Monitor the progress:
   - The script will print progress messages
   - Check for any warnings or errors
   - Final message will indicate completion

### 6. Review Output

The pipeline generates various outputs in the `output` directory:
```
output/
├── filtering_summary.csv       # Quality filtering statistics
├── read_tracking.csv          # Read counts through pipeline
├── alpha_diversity.csv        # Alpha diversity metrics
├── advanced_analysis/
│   ├── core_microbiome.csv    # Core microbiome data
│   ├── *_abundance.csv        # Taxonomic abundance data
│   └── *.pdf                  # Visualization plots
└── sessionInfo.txt            # R session information
```

### 7. Best Practices

1. **Quality Control:**
   - Always check the quality reports before proceeding
   - Monitor read loss at each step
   - Verify taxonomy assignment success rate

2. **Parameter Optimization:**
   - Consider adjusting filtering parameters based on your data
   - Test different abundance/prevalence thresholds
   - Validate rarefaction depth choice

3. **Data Validation:**
   - Compare results with expected community composition
   - Verify control samples if available
   - Check for batch effects

## Output

The pipeline generates various output files in the `output` directory:

### Basic Analysis
- Filtering and tracking summaries
- Alpha diversity metrics
- Beta diversity plots
- Taxonomic composition plots

### Advanced Analysis
- Core microbiome analysis
- Taxonomic abundance analysis at multiple levels
- Heatmaps of abundant taxa
- Detailed abundance tables

## File Structure

```
.
├── 16s-metagenomics-complete.R  # Main analysis script
├── rawFastq/
│   └── trimmed/                 # Input fastq files
├── output/                      # Results directory
└── silva_*.fa                   # SILVA taxonomy databases
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

A comprehensive pipeline for analyzing 16S rRNA sequencing data with advanced statistical analysis and publication-ready visualizations.

## Features

### Core Analysis
- ASV generation and taxonomy assignment
- Alpha and beta diversity analyses
- Phylogenetic tree construction
- Differential abundance testing

### Advanced Analysis
- Core microbiome identification
- Phylogenetic diversity metrics
- Microbial co-occurrence networks
- Indicator species analysis
- ANCOM-BC differential abundance

### Visualization
- Publication-ready static plots
- Interactive visualizations
- Advanced heatmaps
- Taxonomic composition plots
- Network visualizations

## Prerequisites

### System Requirements

Before installing the R packages, you need to install system dependencies:

```bash
# For Ubuntu/Debian
sudo apt-get update
sudo apt-get install -y \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    libmpfr-dev \
    libgmp-dev \
    libcairo2-dev \
    libxt-dev \
    libx11-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev

# For CentOS/RHEL
sudo yum install -y \
    libcurl-devel \
    libxml2-devel \
    openssl-devel \
    mpfr-devel \
    gmp-devel \
    cairo-devel \
    libXt-devel \
    libX11-devel \
    harfbuzz-devel \
    fribidi-devel \
    freetype-devel \
    libpng-devel \
    libtiff-devel \
    libjpeg-devel
```

### R Package Installation

First, install BiocManager:
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
```

#### 1. Core Packages (Required)
```R
# Install Bioconductor packages
BiocManager::install(c(
    "dada2",
    "phyloseq",
    "ShortRead",
    "Biostrings",
    "DECIPHER",
    "DESeq2"
))

# Install CRAN packages
install.packages(c(
    "vegan",
    "tidyverse",
    "ape"
))
```

#### 2. Advanced Analysis Packages (Optional)
```R
# Install Bioconductor packages
BiocManager::install(c(
    "microbiome",
    "ANCOMBC"
))

# Install CRAN packages
install.packages(c(
    "plotly",
    "DT",
    "picante",
    "ComplexHeatmap",
    "circlize"
))

# Install GitHub packages
if (!require("devtools")) install.packages("devtools")

# Install SpiecEasi
devtools::install_github("zdk123/SpiecEasi")

# Install NetCoMi (network analysis)
devtools::install_github("stefpeschel/NetCoMi", 
                        dependencies = TRUE,
                        repos = c("https://cloud.r-project.org",
                                BiocManager::repositories()))

# Install metagMisc (additional utilities)
devtools::install_github("vmikk/metagMisc")
```

### Troubleshooting Package Installation

If you encounter issues:

1. **ANCOMBC Installation**: Try installing dependencies first:
   ```R
   BiocManager::install(c("scater", "mia"))
   install.packages("CVXR")
   BiocManager::install("ANCOMBC")
   ```

2. **Cairo/Graphics Issues**: Ensure X11 and Cairo are properly installed:
   ```bash
   # Ubuntu/Debian
   sudo apt-get install libcairo2-dev libxt-dev libx11-dev
   
   # CentOS/RHEL
   sudo yum install cairo-devel libXt-devel libX11-devel
   ```

3. **Network Analysis Packages**: If NetCoMi or SpiecEasi fail:
   ```R
   # Update devtools first
   install.packages("devtools")
   
   # Try installing with specific dependencies
   devtools::install_github("zdk123/SpiecEasi", ref="master")
   devtools::install_github("stefpeschel/NetCoMi", dependencies = TRUE)
   ```

## 🚀 Usage

1. **Basic Analysis Pipeline**
   ```bash
   Rscript 16s-metagenomics-complete.R
   ```

2. **Advanced Analysis**
   ```R
   source("16s-metagenomics-advanced.R")
   ```

3. **Interactive Visualizations**
   ```R
   source("16s-metagenomics-interactive.R")
   ```

## 📁 Project Structure

```
.
├── 16s_pipeline.yaml          # Pipeline configuration
├── 16s-metagenomics-complete.R  # Core analysis script (single entrypoint)
├── 16s-metagenomics-vis.R    # Static visualization script
├── 16s-metagenomics-advanced.R # Advanced analysis script
├── 16s-metagenomics-interactive.R # Interactive visualization script
├── cutadapt-16s-trim.sh      # Trimming script
├── rawFastq/                        # Raw and processed data
│   ├── filtered/
│   └── trimmed/
└── output/                    # Analysis outputs
    ├── advanced_analysis/     # Advanced analysis results and visualizations
```

## 🔍 Analysis Workflow

1. **Data Preprocessing**
   - Quality filtering
   - Trimming
   - ASV generation

2. **Core Analysis**
   - Taxonomy assignment
   - Diversity analyses
   - Basic visualizations

3. **Advanced Analysis**
   - Core microbiome
   - Network analysis
   - Differential abundance
   - Indicator species

4. **Visualization**
   - Static plots
   - Interactive visualizations
   - Publication-ready figures

## 📊 Output Files

- `output/alpha_diversity.csv`: Alpha diversity metrics
- `output/beta_bray.rds`: Bray-Curtis dissimilarity matrix
- `output/advanced_analysis/`: Advanced analysis results and visualizations

## 📝 Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{16S_metagenomics_pipeline,
  author = {Taimoor Khan},
  title = {Advanced 16S rRNA Metagenomic Analysis Pipeline},
  year = {2025},
  url = {https://github.com/taimoorkhan007/16s-metagenomics}
}
```

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

## Publication-ready figures (saved at 600 dpi)

The pipeline exports high-resolution figures suitable for publication:

- output/alpha_diversity_plot.tiff
   - Box+jitter plots for Observed, Shannon, and InvSimpson diversity; can be stratified by Group (Host) when metadata is present.
- output/beta_diversity_pcoa.tiff
   - PCoA of Bray–Curtis distances; points colored by Group (Host) if available.
- output/phylum_composition.tiff
   - Stacked bar plot of relative abundance (Top phyla) per sample.
- output/phylum_composition_by_group.tiff
   - Group-aggregated composition (samples merged by Group, then normalized).
- output/rarefaction_curves.pdf
   - Sample-wise rarefaction curves from pre-rarefaction counts.

Statistical summaries:
- output/alpha_diversity_stats.txt – Kruskal–Wallis tests by Group (+ pairwise Wilcoxon in separate files)
- output/permanova_group.txt – PERMANOVA on Bray–Curtis by Group
- output/betadisper_group.txt – Test for homogeneity of dispersion across Groups

Core data objects and QC:
- output/phyloseq_object_raw.rds, output/phyloseq_rarefied.rds – phyloseq objects for downstream work
- output/filtering_summary.csv, output/read_tracking.csv – DADA2 filtering and read flow
- output/seq_length_summary.txt – ASV sequence-length distribution to guide optional length filtering

Tip: Ensure `metadata.csv` contains Sample IDs matching FASTQ basenames (e.g., HF1) and a Group column (e.g., Host: Buffalo/Cow/Goat) for group-wise plots and tests.