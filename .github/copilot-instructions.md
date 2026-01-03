# Copilot Instructions for 16S Metagenomics Pipeline

## Project Overview
- This is an automated pipeline for 16S rRNA amplicon sequencing analysis, from raw FASTQ to publication-ready figures.
- Main entrypoint: `run_pipeline.sh` (bash). This orchestrates all steps, including R scripts and QC tools.
- Core analysis logic is in R scripts under `scripts/` and `scripts/modules/`.
- Output is organized by cohort in the `output/` directory, with key results and visualizations.

## Architecture & Data Flow
- **Input:** Paired-end FASTQ files + `metadata.csv` (sample info) + YAML config (see `config/`).
- **Pipeline steps:**
  1. QC (FastQC)
  2. Adapter/primer trimming (Cutadapt)
  3. Post-trim QC (FastQC)
  4. ASV inference (DADA2, R)
  5. Taxonomic classification (SILVA, R)
  6. Diversity, stats, and plots (R)
  7. Aggregated QC (MultiQC)
- **R orchestrator:** `scripts/runner.R` loads config, validates, and calls modular analysis/plotting scripts.
- **Modular R code:** Each analysis step is a separate R script in `scripts/modules/` (e.g., `alpha_diversity.R`, `differential_abundance.R`).

## Developer Workflows
- **Run full pipeline:** `./run_pipeline.sh config/my_project.yaml`
- **Test pipeline:** `scripts/test_pipeline.sh` (uses test data)
- **Preflight checks:** `scripts/preflight_check.sh` (validates environment and inputs)
- **Install R packages:** `scripts/install_packages.sh` (auto-run on first use)
- **Generate test data:** `scripts/generate_test_data.R`
- **Validate config:** `scripts/validate_config.R <config.yaml>`

## Project Conventions
- **Config:** All runs are controlled by YAML config files in `config/`. Use `cp config/kmu_config.yaml config/my_project.yaml` to start.
- **Metadata:** `SampleID` in metadata must match FASTQ filenames. `Group` is required for comparisons.
- **Output:** Results are always written to `output/<cohort>/`.
- **Error handling:** Both bash and R scripts have custom error handlers with clear messages and provenance logging.
- **Statistical tests:** Non-parametric by default, FDR correction for all p-values.
- **Visualization:** Modular plotting system in `scripts/modules/` (set `USE_MODULAR_VIZ` env var to enable).

## Integration & Dependencies
- **R ≥4.3.0** required. All R dependencies are auto-installed if missing.
- **External tools:** FastQC, Cutadapt, MultiQC (run via conda env `16s_metagenomics`).
- **Reference DBs:** SILVA files must be downloaded and paths set in config.

## Key Files & Directories
- `run_pipeline.sh`: Main entrypoint, orchestrates all steps
- `scripts/runner.R`: R orchestrator, loads config, runs analysis modules
- `scripts/modules/`: Modular R code for each analysis/plotting step
- `config/`: Config templates and project configs
- `output/`: All results, organized by cohort
- `test_data/`: Example data for testing

## Example: Minimal Run
```bash
cp config/kmu_config.yaml config/my_project.yaml
./run_pipeline.sh config/my_project.yaml
```

## Troubleshooting
- See README and error messages for common issues (e.g., missing files, config errors, memory limits).
- Provenance files are written to output for reproducibility and debugging.

---
For more, see `README.md` and `PIPELINE_FINAL_SUMMARY.md`.
