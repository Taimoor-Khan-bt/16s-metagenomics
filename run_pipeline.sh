#!/bin/bash
set -euo pipefail

# --- run_pipeline.sh ---
#
# This script is the master entry point for the 16S metagenomics analysis pipeline.
# It automates all steps from raw FASTQ files to final analysis and reports.
#
# Usage:
# ./run_pipeline.sh <config_filename.yaml>
#
# Example:
# ./run_pipeline.sh 16s_pipeline.yaml
#
# Steps:
# 1. Parse the YAML config file to get project parameters.
# 2. Run FastQC on raw reads for initial quality assessment.
# 3. Trim adapters and low-quality bases using Cutadapt.
# 4. Run FastQC on trimmed reads to verify trimming.
# 5. Execute the main R analysis pipeline (DADA2, taxonomy, etc.).
# 6. Aggregate all QC reports into a single MultiQC HTML file.
#
# -------------------------

# Error handler
trap 'error_handler $? $LINENO' ERR

error_handler() {
    local exit_code=$1
    local line_number=$2
    
    # Record failure in provenance if file exists
    if [ -n "${PROVENANCE_FILE:-}" ] && [ -f "$PROVENANCE_FILE" ]; then
        echo "" >> "$PROVENANCE_FILE"
        echo "execution:" >> "$PROVENANCE_FILE"
        echo "  timestamp_end: \"$(date -u +"%Y-%m-%dT%H:%M:%SZ")\"" >> "$PROVENANCE_FILE"
        echo "  status: failed" >> "$PROVENANCE_FILE"
        echo "  error_code: $exit_code" >> "$PROVENANCE_FILE"
        echo "  error_line: $line_number" >> "$PROVENANCE_FILE"
        echo "  duration_seconds: $SECONDS" >> "$PROVENANCE_FILE"
    fi
    
    echo ""
    echo "=========================================="
    echo "ERROR: Pipeline failed at line $line_number with exit code $exit_code"
    echo "=========================================="
    echo "Check the output above for error details."
    echo "Common issues:"
    echo "  - Missing input files"
    echo "  - Conda environment not activated"
    echo "  - Insufficient disk space or memory"
    echo "  - Invalid configuration parameters"
    echo ""
    if [ -n "${PROVENANCE_FILE:-}" ]; then
        echo "Provenance with error details: $PROVENANCE_FILE"
    fi
    exit "$exit_code"
}

# --- Environment and Setup ---
CONDA_ENV="16s_metagenomics"
# Get the absolute path of the directory where this script is located
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

echo "--- 16S Metagenomics Pipeline running in Conda env: $CONDA_ENV ---"
echo "Script directory: $SCRIPT_DIR"

# Check for config file argument
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <path_to_config.yaml>"
    echo "Example: $0 config/kmu_comprehensive.yaml"
    exit 1
fi

CONFIG_FILE=$1
echo "Using config file: $CONFIG_FILE"

# Check if config file exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Config file not found at $CONFIG_FILE"
    exit 1
fi

# Function to parse 2-level nested YAML keys like 'project.output_dir'
get_config_value() {
    local key_path=$1
    local yaml_file=$2
    
    # This parser is intentionally simple and only handles up to 2 levels.
    # It uses sed to find the block for the parent key and then greps the child key.
    if [[ "$key_path" == *.* ]]; then
        local parent_key=$(echo "$key_path" | cut -d. -f1)
        local child_key=$(echo "$key_path" | cut -d. -f2)
        
        # Extract the block under the parent key, then find the child key's value.
        # The sed command finds the parent block, stopping at the next non-indented line.
        # The grep isolates the line with the child key.
        # The final sed cleans up the key and surrounding whitespace/quotes.
        sed -n "/^${parent_key}:/,/^[a-zA-Z]/p" "$yaml_file" | \
            grep -E "^\s*${child_key}:" | \
            sed -E "s/^\s*${child_key}:\s*//" | \
            tr -d '"' | tr -d "'" | head -n 1
    else
        # Top-level key (original behavior)
        grep "^${key_path}:" "$yaml_file" | sed "s/^${key_path}:\s*//" | tr -d '"' | tr -d "'" | head -n 1
    fi
}

# --- Parse Parameters ---
echo "Parsing configuration..."
INPUT_DIR=$(get_config_value "io.input_dir" "$CONFIG_FILE")
OUTPUT_DIR=$(get_config_value "project.output_dir" "$CONFIG_FILE")
COHORT=$(get_config_value "io.cohort" "$CONFIG_FILE")

# Fallback for cohort name if not specified
if [ -z "$COHORT" ] || [ "$COHORT" == "null" ]; then
    COHORT=$(basename "$INPUT_DIR")
fi

# --- Validate Parameters ---
if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: 'input_dir' or 'output_dir' not found in config file."
    exit 1
fi

echo "Input Directory: $INPUT_DIR"
echo "Output Directory: $OUTPUT_DIR"
echo "Cohort: $COHORT"

# Run pre-flight checks before starting pipeline
echo ""
echo "=========================================="
echo "Running pre-flight checks..."
echo "=========================================="
"$SCRIPT_DIR/scripts/preflight_check.sh" --input-dir "$INPUT_DIR" --output-dir "$OUTPUT_DIR"

if [ $? -ne 0 ]; then
    echo ""
    echo "Pre-flight checks failed. Please fix issues before running pipeline."
    exit 1
fi

echo ""
echo "Pre-flight checks passed. Starting pipeline..."
echo ""

# Define Cohort Output Directory
COHORT_OUTPUT_DIR="${OUTPUT_DIR}/${COHORT}"
mkdir -p "$COHORT_OUTPUT_DIR"

# Generate provenance information at start
echo "=========================================="
echo "Recording provenance information..."
echo "=========================================="
PROVENANCE_FILE="$COHORT_OUTPUT_DIR/provenance_$(date +%Y%m%d_%H%M%S).yaml"
"$SCRIPT_DIR/scripts/generate_provenance.sh" --output "$PROVENANCE_FILE" --config "$CONFIG_FILE"
echo "Provenance recorded in: $PROVENANCE_FILE"
echo ""
echo "---------------------------------"
TRIMMED_DIR="${COHORT_OUTPUT_DIR}/trimmed"
FASTQC_RAW_DIR="${COHORT_OUTPUT_DIR}/fastqc_raw"
FASTQC_TRIMMED_DIR="${COHORT_OUTPUT_DIR}/fastqc_trimmed"
MULTIQC_DIR="${COHORT_OUTPUT_DIR}/multiqc"

# Create output directories
mkdir -p "$TRIMMED_DIR"
mkdir -p "$FASTQC_RAW_DIR"
mkdir -p "$FASTQC_TRIMMED_DIR"
mkdir -p "$MULTIQC_DIR"

echo "Pipeline starting..."
echo "---------------------------------"

# --- Step 1: Initial FastQC on Raw Reads ---
echo "STEP 1: Running FastQC on raw reads..."
# Find both .fastq.gz and .fq.gz files
find "${INPUT_DIR}" -maxdepth 1 \( -name "*.fastq.gz" -o -name "*.fq.gz" \) | xargs -r conda run -n "$CONDA_ENV" fastqc -o "$FASTQC_RAW_DIR" --threads "$(nproc)"
echo "FastQC on raw reads complete."
echo "---------------------------------"

# --- Step 2: Trim Primers with Cutadapt ---
echo "STEP 2: Trimming primers with Cutadapt..."
TRIM_SCRIPT="${SCRIPT_DIR}/cutadapt-16s-trim.sh"

if [ ! -f "$TRIM_SCRIPT" ]; then
    echo "Error: Trimming script not found at $TRIM_SCRIPT"
    exit 1
fi

conda run -n "$CONDA_ENV" bash "$TRIM_SCRIPT" "$INPUT_DIR" "$TRIMMED_DIR"

echo "Cutadapt trimming complete."
echo "---------------------------------"

# --- Step 3: FastQC on Trimmed Reads ---
echo "STEP 3: Running FastQC on trimmed reads..."
conda run -n "$CONDA_ENV" fastqc -o "$FASTQC_TRIMMED_DIR" --threads "$(nproc)" "${TRIMMED_DIR}"/*_trimmed.fastq.gz
echo "FastQC on trimmed reads complete."
echo "---------------------------------"

# --- Step 4: Run DADA2 and R Analysis ---
echo "STEP 4: Executing R analysis pipeline..."
R_SCRIPT="${SCRIPT_DIR}/scripts/runner.R"
if [ ! -f "$R_SCRIPT" ]; then
    echo "Error: R script not found at $R_SCRIPT"
    exit 1
fi
Rscript "$R_SCRIPT" --config "$CONFIG_FILE"
echo "R analysis complete."
echo "---------------------------------"

# --- Step 5: Aggregate QC with MultiQC ---
echo "STEP 5: Generating MultiQC report..."
conda run -n "$CONDA_ENV" multiqc "$COHORT_OUTPUT_DIR" -o "$MULTIQC_DIR" --force
echo "MultiQC report generated."
echo "---------------------------------"

echo "Pipeline finished successfully!"

# Record completion time in provenance file
echo "" >> "$PROVENANCE_FILE"
echo "execution:" >> "$PROVENANCE_FILE"
echo "  timestamp_end: \"$(date -u +"%Y-%m-%dT%H:%M:%SZ")\"" >> "$PROVENANCE_FILE"
echo "  status: success" >> "$PROVENANCE_FILE"
echo "  duration_seconds: $SECONDS" >> "$PROVENANCE_FILE"

echo ""
echo "=========================================="
echo "Pipeline execution complete!"
echo "=========================================="
echo "Provenance: $PROVENANCE_FILE"
echo "Results: $OUTPUT_DIR"
echo "=========================================="
