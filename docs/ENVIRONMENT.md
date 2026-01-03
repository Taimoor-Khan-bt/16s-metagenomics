# Environment Management

## Overview
This document describes how to manage the computational environment for the 16S metagenomics pipeline to ensure reproducibility.

## Environment File

### `environment.yaml`
This file contains the complete list of dependencies with **pinned versions** for exact reproducibility. It was generated from a working environment and should reproduce the exact same software stack.

**Generated:** November 2, 2025  
**Conda version:** Based on conda-forge, bioconda channels

### Key Dependencies

**Analysis Tools:**
- FastQC v0.12.1 - Quality control
- Cutadapt v4.6 - Adapter trimming
- MultiQC v1.31 - Report aggregation

**R Environment:**
- R v4.3.3
- DADA2 v1.28.0 - ASV inference
- phyloseq v1.44.0 - Microbiome analysis
- Biostrings v2.68.1 - Sequence manipulation
- ggplot2 v3.5.1 - Visualization
- vegan v2.6.8 - Ecological statistics

**Python Libraries:**
- numpy, pandas, scipy - Data processing
- matplotlib, seaborn - Visualization

## Creating a New Environment

### From Scratch (Recommended for New Users)

1. **Create environment from the locked file:**
   ```bash
   conda env create -f environment.yaml
   ```

2. **Activate the environment:**
   ```bash
   conda activate 16s_pipeline
   ```

3. **Verify installation:**
   ```bash
   ./scripts/preflight_check.sh --input-dir <your_data>
   ```

### Alternative: Manual Installation

If you need to customize or troubleshoot:

```bash
# Create base environment
conda create -n 16s_pipeline python=3.11

# Activate
conda activate 16s_pipeline

# Install core tools
conda install -c bioconda -c conda-forge \
    fastqc=0.12.1 \
    cutadapt=4.6 \
    multiqc=1.31

# Install R and packages
conda install -c conda-forge r-base=4.3.3
conda install -c bioconda \
    bioconductor-dada2=1.28.0 \
    bioconductor-phyloseq=1.44.0 \
    bioconductor-biostrings=2.68.1

conda install -c conda-forge \
    r-ggplot2=3.5.1 \
    r-dplyr=1.1.4 \
    r-tidyr=1.3.1 \
    r-vegan=2.6.8 \
    r-yaml
```

## Updating the Environment

### When to Update
- Adding new features requiring additional packages
- Security patches for existing packages
- Bug fixes in core dependencies

### How to Update

1. **Activate environment:**
   ```bash
   conda activate 16s_pipeline
   ```

2. **Install/update packages:**
   ```bash
   conda install <package>=<version>
   ```

3. **Test thoroughly:**
   ```bash
   # Run preflight checks
   ./scripts/preflight_check.sh --input-dir test_data
   
   # Run pipeline on test dataset
   ./run_pipeline.sh test_config.yaml
   ```

4. **Re-export if successful:**
   ```bash
   conda env export --no-builds > environment.yaml
   ```

5. **Document changes:**
   - Update this file with new versions
   - Note any breaking changes
   - Update pipeline documentation if needed

## Troubleshooting

### Environment Creation Fails

**Issue:** Package conflict during environment creation
```
Solving environment: failed with initial frozen solve. Retrying with flexible solve.
```

**Solutions:**
1. Try creating with a more flexible spec:
   ```bash
   conda env create -f environment.yaml --force
   ```

2. Use mamba (faster solver):
   ```bash
   conda install mamba -n base -c conda-forge
   mamba env create -f environment.yaml
   ```

3. Check channel priorities:
   ```bash
   conda config --show channels
   # Should be: conda-forge, bioconda, defaults
   ```

### Missing R Packages

**Issue:** R packages not found by Bioconductor

**Solution:**
```r
# In R console
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("dada2", "phyloseq", "Biostrings"))
```

### Tool Version Mismatches

**Issue:** Different versions installed than expected

**Solution:**
```bash
# Check installed versions
fastqc --version
cutadapt --version
multiqc --version
Rscript --version

# Compare with environment.yaml
conda list | grep -E "(fastqc|cutadapt|multiqc)"
```

## Environment Isolation Best Practices

### 1. Never Install in Base
```bash
# BAD - installs in base environment
conda install fastqc

# GOOD - always activate target environment first
conda activate 16s_pipeline
conda install fastqc
```

### 2. Use Environment Per Project
```bash
# If you need different versions for different projects
conda create -n project1_16s -f environment_v1.yaml
conda create -n project2_16s -f environment_v2.yaml
```

### 3. Export Regularly
```bash
# After any changes, re-export
conda activate 16s_pipeline
conda env export --no-builds > environment_$(date +%Y%m%d).yaml
```

## Reproducibility Checklist

When sharing your analysis or publishing, ensure:

- [ ] `environment.yaml` is included in repository
- [ ] Environment file is dated or versioned
- [ ] All tool versions documented in methods
- [ ] Pipeline version/commit hash recorded
- [ ] Operating system documented (Linux, macOS, etc.)
- [ ] Hardware specifications noted (RAM, CPU cores)

## System Requirements

### Minimum
- **OS:** Linux or macOS
- **RAM:** 8 GB
- **Disk:** 20 GB free space
- **CPU:** 2 cores

### Recommended
- **OS:** Linux (Ubuntu 20.04+, CentOS 7+)
- **RAM:** 16+ GB
- **Disk:** 100+ GB (for large datasets)
- **CPU:** 8+ cores (for parallel processing)

## Platform-Specific Notes

### Linux (Ubuntu/Debian)
Generally works out of the box. May need:
```bash
sudo apt-get install build-essential
```

### macOS
Some Bioconductor packages may require Xcode:
```bash
xcode-select --install
```

### Windows (WSL2)
Run through Windows Subsystem for Linux:
```bash
# In WSL2 Ubuntu
conda env create -f environment.yaml
```

## Version History

| Date | Version | Changes |
|------|---------|---------|
| 2025-11-02 | 1.0 | Initial locked environment with pinned versions |

## Support

If you encounter environment issues:
1. Check this documentation
2. Run `./scripts/preflight_check.sh` for diagnostics
3. Review `environment.yaml` for version conflicts
4. See GitHub Issues for known problems
