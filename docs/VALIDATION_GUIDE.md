# Input Validation Guide

## Quick Start

Validate your pipeline inputs before running analysis:

```bash
./scripts/validate_inputs.R --config config/your_config.yaml --verbose
```

This will check:
- ✓ Configuration structure and parameters
- ✓ Metadata file structure and sample IDs
- ✓ FASTQ file pairs and integrity
- ✓ Metadata-FASTQ sample alignment

## Command-Line Options

```bash
./scripts/validate_inputs.R [options]

Options:
  -c, --config FILE       Path to configuration file (required)
  -m, --metadata FILE     Path to metadata file (optional, can be in config)
  -f, --fastq-dir DIR     Path to FASTQ directory (optional, can be in config)
  -v, --verbose           Print detailed output
  -h, --help             Show help message
```

## Usage Examples

### Basic validation
```bash
./scripts/validate_inputs.R --config config/project_config.yaml
```

### Verbose output with explicit paths
```bash
./scripts/validate_inputs.R \
  --config config/project_config.yaml \
  --metadata data/samples/metadata.csv \
  --fastq-dir data/samples/fastq \
  --verbose
```

### In a pipeline script
```bash
# Validate before running
if ./scripts/validate_inputs.R --config config/project_config.yaml; then
    echo "Validation passed, starting pipeline..."
    snakemake --cores 8
else
    echo "Validation failed, fix errors before continuing"
    exit 1
fi
```

## What Gets Validated

### 1. Configuration (config.yaml)
- **Structure:** YAML syntax, required sections
- **Parameters:** Numeric ranges (threads, maxEE, truncLen, rarefaction depth)
- **Paths:** Presence of required path configurations
- **Inheritance:** Config inheritance from base_config.yaml

**Common Issues:**
- Missing required parameters
- Out-of-range values (e.g., threads < 1, maxEE > 10)
- YAML syntax errors
- Inheritance chain broken

### 2. Metadata (metadata.csv)
- **Structure:** Data frame format, non-empty
- **Sample IDs:** Format, duplicates, special characters
- **Columns:** Required columns present (sample_id, group_column)
- **Data Quality:** NA values, single-level factors
- **Groups:** Sufficient groups for statistical comparisons

**Common Issues:**
- Missing sample_id column
- Duplicate sample IDs
- Sample IDs with whitespace or special characters
- Missing grouping variable for statistics
- All NA values in columns

### 3. FASTQ Files
- **Pairing:** Forward and reverse reads matched
- **File Size:** Files not empty (> 1KB)
- **Existence:** Files actually present on disk
- **Naming:** Files follow expected pattern (_1.fastq.gz, _2.fastq.gz)

**Common Issues:**
- Unpaired reads (forward without reverse)
- Empty or corrupted FASTQ files
- Incorrect file naming pattern
- Missing files

### 4. Cross-Validation
- **Alignment:** Metadata samples match FASTQ samples
- **Completeness:** All samples accounted for
- **Discrepancies:** Clear reporting of mismatches

**Common Issues:**
- Samples in metadata but no FASTQ files
- FASTQ files but no metadata entries
- Sample ID naming inconsistencies (TEST-001 vs TEST_001)

## Validation Output

### Success Example
```
═══════════════════════════════════════════════════════════════
  16S Metagenomics Pipeline - Input Validation
═══════════════════════════════════════════════════════════════

1️⃣  Validating Configuration...
   ✓ Configuration structure valid
   • Project: MyProject-2024
   • Threads: 8
   • Rarefaction depth: 5000

2️⃣  Validating Metadata...
   ✓ Metadata structure valid
   • Samples: 24
   • Columns: 6
   • Groups (treatment): Control, DrugA, DrugB

3️⃣  Validating FASTQ Files...
   ✓ FASTQ file pairs valid
   • Sample pairs: 24
   • Total size: 12.5 GB

4️⃣  Cross-Validating Metadata and FASTQ Samples...
   ✓ Metadata and FASTQ samples aligned

═══════════════════════════════════════════════════════════════
✅ All validations passed! Pipeline is ready to run.
```

### Failure Example
```
═══════════════════════════════════════════════════════════════
  Validation Summary
═══════════════════════════════════════════════════════════════

❌ ERRORS:
  1. ERROR: Sample ID column 'SampleID' not found in metadata
  2. ERROR: 2 forward read(s) without matching reverse: S25, S26
  3. ERROR: 3 sample(s) in metadata but not in FASTQ files

⚠️  WARNINGS:
  1. WARNING: 2 sample ID(s) contain whitespace
  2. WARNING: Only 1 group(s) found for statistical comparisons

❌ Validation failed. Please address the errors above before running.
```

## Programmatic Usage

You can also use validation functions in your own R scripts:

```r
# Load modules
source('scripts/modules/input_validation.R')
source('scripts/modules/error_handling.R')

# Validate sample IDs
sample_ids <- c("TEST-001", "TEST-002", "TEST-003")
result <- validate_sample_ids(sample_ids)
if (!result$valid) {
  print_validation_results(result)
  stop("Invalid sample IDs")
}

# Validate metadata
metadata <- read.csv("metadata.csv")
result <- validate_metadata_structure(
  metadata, 
  required_cols = c("SampleID", "Group"),
  sample_id_col = "SampleID"
)

# Validate FASTQ pairs
result <- validate_fastq_pairs("data/fastq")
if (result$valid) {
  cat("Found", nrow(result$pairs), "sample pairs\n")
}

# Validate numeric parameters
result <- validate_numeric_range(
  value = 250,
  param_name = "truncLen_forward",
  min = 50,
  max = 500
)
```

## Integration with Pipeline

### Pre-flight Check Integration

Add to your pipeline runner:

```bash
#!/bin/bash
set -e

# Validate inputs first
echo "Validating pipeline inputs..."
./scripts/validate_inputs.R --config "$CONFIG_FILE" --verbose

# Run preflight checks
echo "Running preflight checks..."
./scripts/preflight_check.sh "$CONFIG_FILE"

# Start pipeline
echo "Starting pipeline..."
./scripts/runner.R "$CONFIG_FILE"
```

### Snakemake Integration

Add as a rule:

```python
rule validate_inputs:
    input:
        config="config/project_config.yaml",
        metadata="data/metadata.csv"
    output:
        touch("results/.validation_passed")
    shell:
        """
        ./scripts/validate_inputs.R \
            --config {input.config} \
            --metadata {input.metadata} \
            --verbose
        """

rule all:
    input:
        "results/.validation_passed",
        "results/final_outputs.rds"
```

## Error Types and Solutions

### Configuration Errors

**Error:** `Parameter 'threads' (128) exceeds maximum allowed value (64)`
- **Solution:** Reduce threads in config to ≤ 64

**Error:** `Configuration missing 'input_dir'`
- **Solution:** Add `input_dir: "path/to/fastq"` to config

### Metadata Errors

**Error:** `Sample ID column 'SampleID' not found in metadata`
- **Solution:** Check column names in CSV (case-sensitive), update config metadata.columns.sample_id

**Error:** `2 duplicate sample ID(s): S001, S002`
- **Solution:** Make sample IDs unique in metadata.csv

### FASTQ Errors

**Error:** `3 forward read(s) without matching reverse`
- **Solution:** Check that each *_1.fastq.gz has matching *_2.fastq.gz

**Error:** `2 sample(s) have very small FASTQ files (< 1KB)`
- **Solution:** Investigate these samples - files may be empty or corrupted

### Cross-Validation Errors

**Error:** `5 sample(s) in metadata but not in FASTQ files`
- **Solution:** Either add missing FASTQ files or remove samples from metadata

**Error:** `No samples match between metadata and FASTQ files`
- **Solution:** Check sample ID naming consistency (underscores vs hyphens, case)

## Tips and Best Practices

1. **Run validation early and often** - Before starting any analysis
2. **Use verbose mode** - Helps debug configuration issues
3. **Check warnings too** - They may indicate data quality issues
4. **Integrate with CI/CD** - Validation exit codes (0/1) work with automation
5. **Validate before compute** - Saves time by catching errors before expensive processing
6. **Keep validation logs** - Redirect output to log file for records

```bash
./scripts/validate_inputs.R \
  --config config/project.yaml \
  --verbose 2>&1 | tee validation_$(date +%Y%m%d).log
```

## Validation Checklist

Before running the pipeline, ensure:

- [ ] Configuration file loads without errors
- [ ] All required config parameters present and in valid ranges
- [ ] Metadata CSV has required columns (sample_id, group_column)
- [ ] No duplicate sample IDs
- [ ] All FASTQ files paired (forward + reverse)
- [ ] Sample IDs match between metadata and FASTQ files
- [ ] At least 2 groups for statistical comparisons
- [ ] Sufficient disk space for outputs
- [ ] All database files exist (SILVA, etc.)

Run validation to check all items:
```bash
./scripts/validate_inputs.R --config your_config.yaml --verbose
```
