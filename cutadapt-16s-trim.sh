#!/usr/bin/env bash
# ============================================================
# 16S rRNA Primer Trimming Script using Cutadapt
# Author: Taimoor Khan
# Email: taimoorkhan007.tk@gmail.com
# Date: October 20, 2025
# ============================================================

# Exit immediately if a command exits with a non-zero status
set -euo pipefail

# ============================================================
# Configuration Section - MODIFY AS NEEDED
# ============================================================

# Input and output directories
INPUT_DIR="HF"
OUTPUT_DIR="${INPUT_DIR}/trimmed"

# Default primer sequences (Illumina V4 region)
# IMPORTANT: Replace these with your actual primer sequences!
FWD_PRIMER="GTGYCAGCMGCCGCGGTAA"     # 515F primer
REV_PRIMER="GGACTACNVGGGTWTCTAAT"    # 806R primer

# Quality thresholds
MIN_LENGTH=200         # Minimum length after trimming
MIN_QUALITY=20        # Minimum quality score
ERROR_RATE=0.1        # Maximum error rate for matching primers

# Performance settings
THREADS=$(nproc)      # Use all available CPU cores
MAX_MEMORY="4G"       # Maximum memory usage

# ============================================================
# Input Validation
# ============================================================

# Check if input directory exists and contains files
if [[ ! -d "$INPUT_DIR" ]]; then
    echo "Error: Input directory '$INPUT_DIR' does not exist!"
    exit 1
fi

if [[ -z "$(ls -A $INPUT_DIR/*.fastq.gz 2>/dev/null)" ]]; then
    echo "Error: No .fastq.gz files found in '$INPUT_DIR'!"
    echo "Please ensure your input files are gzipped FASTQ files."
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Setup logging
LOGFILE="${OUTPUT_DIR}/cutadapt_summary.log"
STATS_FILE="${OUTPUT_DIR}/trimming_statistics.tsv"

# Initialize log files
{
    echo "=== 16S rRNA Trimming Pipeline ==="
    echo "Started at: $(date)"
    echo "Input directory: ${INPUT_DIR}"
    echo "Output directory: ${OUTPUT_DIR}"
    echo "Forward primer: ${FWD_PRIMER}"
    echo "Reverse primer: ${REV_PRIMER}"
    echo "Minimum length: ${MIN_LENGTH}"
    echo "Quality threshold: ${MIN_QUALITY}"
    echo "Using ${THREADS} threads"
    echo "----------------------------------------"
} > "$LOGFILE"

# Initialize statistics file
echo -e "Sample\tInput_Reads\tPrimer_Found\tPassed_Filters\tPercent_Kept" > "$STATS_FILE"

# Process each pair of FASTQ files
for R1 in ${INPUT_DIR}/*_1.fastq.gz; do
    # Extract sample name
    SAMPLE=$(basename "$R1" _1.fastq.gz)
    R2="${INPUT_DIR}/${SAMPLE}_2.fastq.gz"
    
    # Validate input files
    if [[ ! -f "$R2" ]]; then
        echo "Warning: Missing R2 file for sample ${SAMPLE}, skipping." | tee -a "$LOGFILE"
        continue
    fi
    
    echo "Processing sample: ${SAMPLE}" | tee -a "$LOGFILE"
    
    # Run cutadapt with comprehensive parameters
    CUTADAPT_OUTPUT=$(cutadapt \
        -g "${FWD_PRIMER}" \
        -G "${REV_PRIMER}" \
        -o "${OUTPUT_DIR}/${SAMPLE}_R1_trimmed.fastq.gz" \
        -p "${OUTPUT_DIR}/${SAMPLE}_R2_trimmed.fastq.gz" \
        --minimum-length "${MIN_LENGTH}" \
        --quality-cutoff "${MIN_QUALITY}" \
        --error-rate "${ERROR_RATE}" \
        --cores="${THREADS}" \
        --memory="${MAX_MEMORY}" \
        --discard-untrimmed \
        --trim-n \
        --max-n 0 \
        "$R1" "$R2" 2>&1)
    
    # Extract statistics
    INPUT_READS=$(echo "$CUTADAPT_OUTPUT" | grep "Total read pairs processed:" | grep -o "[0-9,]*" | tr -d ',')
    PRIMER_FOUND=$(echo "$CUTADAPT_OUTPUT" | grep "Reads with adapters:" | grep -o "[0-9,]*" | tr -d ',')
    PASSED_FILTERS=$(echo "$CUTADAPT_OUTPUT" | grep "Pairs written (passing filters):" | grep -o "[0-9,]*" | tr -d ',')
    
    # Calculate percentage
    PERCENT_KEPT=$(awk "BEGIN {print ($PASSED_FILTERS/$INPUT_READS)*100}")
    
    # Log statistics
    echo -e "${SAMPLE}\t${INPUT_READS}\t${PRIMER_FOUND}\t${PASSED_FILTERS}\t${PERCENT_KEPT}" >> "$STATS_FILE"
    
    # Add to log
    {
        echo "Sample ${SAMPLE} completed:"
        echo "  Input reads: ${INPUT_READS}"
        echo "  Primer found: ${PRIMER_FOUND}"
        echo "  Passed filters: ${PASSED_FILTERS}"
        echo "  Percent kept: ${PERCENT_KEPT}%"
        echo "----------------------------------------"
    } >> "$LOGFILE"
    
    echo "✅ Finished ${SAMPLE}" | tee -a "$LOGFILE"
done

# ============================================================
# Final Summary
# ============================================================

# Generate summary statistics
TOTAL_SAMPLES=$(wc -l < "$STATS_FILE")
TOTAL_SAMPLES=$((TOTAL_SAMPLES - 1))  # Subtract header line

echo -e "\nProcessing Complete! 🎉" | tee -a "$LOGFILE"
echo "Processed ${TOTAL_SAMPLES} samples" | tee -a "$LOGFILE"
echo "Summary statistics saved to: ${STATS_FILE}" | tee -a "$LOGFILE"
echo "Detailed log saved to: ${LOGFILE}" | tee -a "$LOGFILE"
echo "Finished at: $(date)" | tee -a "$LOGFILE"

# Validate output
echo -e "\nValidating output files..."
for SAMPLE in $(cut -f1 "$STATS_FILE" | tail -n +2); do
    R1_OUT="${OUTPUT_DIR}/${SAMPLE}_R1_trimmed.fastq.gz"
    R2_OUT="${OUTPUT_DIR}/${SAMPLE}_R2_trimmed.fastq.gz"
    
    if [[ ! -f "$R1_OUT" ]] || [[ ! -f "$R2_OUT" ]]; then
        echo "Warning: Missing output files for sample ${SAMPLE}!"
    fi
done

echo "Done! 🧬"

