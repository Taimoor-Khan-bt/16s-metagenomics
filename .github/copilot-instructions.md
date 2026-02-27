# Copilot Instructions for 16S Metagenomics Pipeline

## Project Overview
This is an automated 16S rRNA amplicon sequencing analysis pipeline that processes raw FASTQ files through to publication-ready statistical results and visualizations.

**Entry point:** `./run_pipeline.sh config/project.yaml` orchestrates all steps  
**Language stack:** Bash (orchestration) + R (analysis) + Conda (environment)  
**Core principle:** Configuration-driven modular pipeline with comprehensive error handling and provenance tracking

## Architecture & Data Flow

### Pipeline Stages (Bash → R)
1. **QC Raw** (FastQC): `run_pipeline.sh` → FastQC on raw FASTQ
2. **Trimming** (Cutadapt): `scripts/cutadapt-16s-trim.sh` → removes adapters/primers
3. **QC Trimmed** (FastQC): Validates trimming quality
4. **Preprocessing** (R): `scripts/runner.R --step preprocess` → DADA2 ASV inference, taxonomy (SILVA)
5. **Analysis** (R): `scripts/runner.R --step analysis` → diversity metrics, statistical tests
6. **Visualization** (R): `scripts/runner.R --step viz` → publication-quality plots (600 DPI TIFF/PDF)
7. **MultiQC**: Aggregates all QC reports into single HTML

### R Analysis Architecture
- **Orchestrator:** `scripts/runner.R` validates config, manages execution flow, calls modules
- **Modular design:** `scripts/modules/module_*.R` (9 modules: alpha/beta diversity, differential abundance, taxonomy, correlations, core microbiome, QC, utils, PICRUSt)
- **Phyloseq-centric:** All data flows through phyloseq objects (`phyloseq_object_raw.rds`, `phyloseq_rarefied.rds`)
- **Module pattern:** Each module exports standalone functions, sources `module_utils.R` for shared plotting themes

### Configuration System
- **Inheritance:** Configs support `inherit_from: base_config.yaml` for shared defaults (see `docs/CONFIG_INHERITANCE.md`)
- **Three analysis modes:** `simple` (2-group), `standard` (multi-group), `comprehensive` (multi-group + continuous vars + composites)
- **Critical config sections:** `metadata.primary_comparison.group_column` (required), `metadata.secondary_comparisons` (optional), `io.input_dir/metadata_csv`
- **Validation:** `scripts/validate_config.R` runs before execution, checks all required fields

## Critical Developer Workflows

### Running the Pipeline
```bash
# Full pipeline
./run_pipeline.sh config/my_project.yaml

# Specific step only (skip earlier steps)
./run_pipeline.sh config/my_project.yaml analysis

# Available steps: qc, trim, preprocess, analysis, viz, multiqc, all
```

### Testing & Validation
```bash
# Generate synthetic test data → run full pipeline → validate outputs
./scripts/test_pipeline.sh

# Validate environment setup (conda, tools, R packages, disk space)
./scripts/preflight_check.sh --input-dir data/raw --output-dir output

# Validate config structure before running
Rscript scripts/validate_config.R config/my_project.yaml
```

### Debugging
- **Bash errors:** Trapped by custom handler in `run_pipeline.sh` line 26-48, writes to provenance file
- **R errors:** Custom `options(error=...)` in `scripts/runner.R` line 11-22, prints traceback + common issues
- **Check provenance:** `output/<cohort>/analysis/provenance.yaml` logs config, environment, execution time, errors
- **Module-level logs:** Each R module prefixes messages with `[ModuleName]`, e.g., `[DA] Starting DESeq2...`

### Adding New Analysis Modules
1. Create `scripts/modules/module_newanalysis.R`
2. Export functions following pattern: `run_newanalysis <- function(ps, cfg, outdir) { ... }`
3. Source `module_utils.R` for consistent plotting themes: `plot_theme(cfg)`, `get_color_palette(cfg)`
4. Add call in `scripts/runner.R` under appropriate step (preprocess/analysis/viz)
5. Update config schema in `config/config_template.yaml` if adding new parameters

## Project-Specific Conventions

### Metadata Requirements
- **SampleID column:** Must exactly match FASTQ filename prefixes (e.g., `Sample1` → `Sample1_R1.fastq.gz`)
- **Group column:** Required for all analyses, specified in `metadata.primary_comparison.group_column`
- **Sample naming:** Normalize with `scripts/normalize_sample_names.sh` if needed (see `docs/SAMPLE_ID_STANDARDS.md`)

### Statistical Methods (Non-Parametric Default)
- **Alpha diversity:** Wilcoxon (2-group) or Kruskal-Wallis (3+ groups), never t-test/ANOVA
- **Beta diversity:** PERMANOVA with Bray-Curtis distance (999 permutations)
- **Differential abundance:** DESeq2 with cascading fallbacks (parametric → local → mean dispersion fits) for sparse data
- **Multiple testing:** FDR correction (Benjamini-Hochberg) on all p-values
- See `docs/STATISTICAL_METHODS.md` for assumptions, citations, interpretation

### DESeq2 Sparse Data Handling
The pipeline handles zero-inflated microbiome data with three-level fallback strategy (`module_differential_abundance.R` lines 64-84):
1. **Primary:** `estimateSizeFactors(type="poscounts")` (handles zeros natively)
2. **Fallback 1:** Gene-wise dispersion estimates if parametric fit fails
3. **Fallback 2:** Local fit → mean fit if dispersion estimation fails

### Visualization System
- **Modular viz:** Enable with `export USE_MODULAR_VIZ=true` before running
- **Shared themes:** `plot_theme(cfg)` in `module_utils.R` ensures consistent styling across all plots
- **Color palettes:** Okabe-Ito (colorblind-friendly, default), Viridis, Tableau10 (set in `plots.color_palette`)
- **Output formats:** Dual output (TIFF + PDF) at 600 DPI for publication submission

### Output Structure
```
output/<cohort>/
├── analysis/          # CSVs, RDS objects, statistics
│   ├── phyloseq_object_raw.rds
│   ├── phyloseq_rarefied.rds
│   ├── alpha_diversity.csv
│   ├── deseq2_*.csv
│   └── provenance.yaml
├── visualizations/    # TIFF + PDF plots
└── qc/               # FastQC, MultiQC reports
```

## Integration Points

### External Tools (Conda Environment)
- **Conda env name:** `16s_pipeline` (check with `conda activate 16s_pipeline`)
- **Required binaries:** fastqc, cutadapt, multiqc (validated by `preflight_check.sh`)
- **R version:** ≥4.3.0, packages auto-installed on first run

### Reference Databases (SILVA)
- **Required files:** `silva_nr99_v138.1_train_set.fa`, `silva_species_assignment_v138.1.fa`
- **Auto-download:** Run `./setup_pipeline.sh` or set `setup.taxonomy_urls` in config
- **Config paths:** `amplicon.taxonomy.train_set` and `amplicon.taxonomy.species_db`

### Resource Management
```yaml
resources:
  cores: 8                # Sets options(mc.cores)
  ram: "16GB"             # Sets future.globals.maxSize and R_MAX_VSIZE
```

## Key Files Reference
- `run_pipeline.sh`: Main orchestrator, parses YAML, calls bash/R steps
- `scripts/runner.R`: R orchestrator, validates config, executes analysis modules
- `scripts/modules/module_*.R`: Nine modular analysis components
- `config/config_template.yaml`: Comprehensive config schema with comments
- `docs/`: Technical documentation (CONFIG_INHERITANCE, STATISTICAL_METHODS, VALIDATION_GUIDE)

## Common Patterns

### Config Parsing in R
```r
cfg <- validate_config("config/project.yaml")  # Returns nested list
group_col <- cfg$metadata$primary_comparison$group_column
cores <- cfg$resources$cores %||% 4  # Use %||% for defaults
```

### Error Handling Pattern
```r
tryCatch({
  result <- main_analysis()
}, error = function(e) {
  message("[ModuleName] ERROR: ", conditionMessage(e))
  return(NULL)  # Graceful degradation
})
```

### Creating New Configs
```bash
# Start from template, inherit shared settings
cp config/config_template.yaml config/new_project.yaml
# Add inheritance line at top: inherit_from: "base_config.yaml"
```

---
**Documentation:** See `README.md` for user guide, `docs/` for technical details  
**Testing:** Always run `scripts/test_pipeline.sh` after modifying core logic  
**Contact:** taimoorkhan007.tk@gmail.com
