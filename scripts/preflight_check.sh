#!/bin/bash
# Pre-flight Checks for 16S Metagenomics Pipeline
# Validates environment setup before running pipeline
# Exit codes: 0 = success, 1 = failure

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

error() {
    echo -e "${RED}✗ ERROR:${NC} $1" >&2
    exit 1
}

warn() {
    echo -e "${YELLOW}⚠ WARNING:${NC} $1" >&2
}

success() {
    echo -e "${GREEN}✓${NC} $1"
}

# === Check 1: Conda Environment ===
check_conda_env() {
    echo "Checking Conda environment..."
    
    if [ -z "$CONDA_DEFAULT_ENV" ]; then
        error "Not in a Conda environment. Run: conda activate 16s_pipeline"
    fi
    
    if [ "$CONDA_DEFAULT_ENV" != "16s_pipeline" ]; then
        warn "Current environment is '$CONDA_DEFAULT_ENV', expected '16s_pipeline'"
        echo "  To use correct environment: conda activate 16s_pipeline"
    fi
    
    success "Conda environment: $CONDA_DEFAULT_ENV"
}

# === Check 2: Required Tools ===
check_tools() {
    echo ""
    echo "Checking required tools..."
    
    local missing_tools=()
    local tools=("fastqc" "cutadapt" "multiqc" "Rscript")
    
    for tool in "${tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            missing_tools+=("$tool")
        else
            local version=""
            case "$tool" in
                fastqc)
                    version=$($tool --version 2>&1 | head -n1)
                    ;;
                cutadapt)
                    version="v$($tool --version 2>&1)"
                    ;;
                multiqc)
                    version="v$($tool --version 2>&1 | grep -oP 'version \K[0-9.]+')"
                    ;;
                Rscript)
                    version=$($tool --version 2>&1 | head -n1)
                    ;;
            esac
            success "$tool found ($version)"
        fi
    done
    
    if [ ${#missing_tools[@]} -gt 0 ]; then
        error "Missing tools: ${missing_tools[*]}"
    fi
}

# === Check 3: R Packages ===
check_r_packages() {
    echo ""
    echo "Checking R packages..."
    
    Rscript --vanilla - <<'EOF'
required_packages <- c(
    "yaml", "dada2", "phyloseq", "Biostrings", 
    "ggplot2", "dplyr", "tidyr", "vegan"
)

missing <- character()
for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        missing <- c(missing, pkg)
    } else {
        version <- packageVersion(pkg)
        cat(sprintf("\033[0;32m✓\033[0m %s (v%s)\n", pkg, version))
    }
}

if (length(missing) > 0) {
    cat(sprintf("\033[0;31m✗ ERROR:\033[0m Missing R packages: %s\n", 
        paste(missing, collapse = ", ")))
    quit(status = 1)
}
EOF
    
    if [ $? -ne 0 ]; then
        error "R package check failed"
    fi
}

# === Check 4: Input Files ===
check_input_files() {
    local input_dir="$1"
    
    echo ""
    echo "Checking input files..."
    
    if [ -z "$input_dir" ]; then
        warn "No input directory specified, skipping file check"
        return 0
    fi
    
    if [ ! -d "$input_dir" ]; then
        error "Input directory not found: $input_dir"
    fi
    
    # Count FASTQ files (support multiple naming patterns)
    local r1_count=$(find "$input_dir" -maxdepth 1 \( -name "*_R1*.fastq.gz" -o -name "*_R1*.fq.gz" -o -name "*_1.fastq.gz" -o -name "*_1.fq.gz" \) | wc -l)
    local r2_count=$(find "$input_dir" -maxdepth 1 \( -name "*_R2*.fastq.gz" -o -name "*_R2*.fq.gz" -o -name "*_2.fastq.gz" -o -name "*_2.fq.gz" \) | wc -l)
    
    if [ "$r1_count" -eq 0 ] || [ "$r2_count" -eq 0 ]; then
        error "No paired-end FASTQ files found in $input_dir (looked for *_R1/_R2 or *_1/_2 patterns)"
    fi
    
    if [ "$r1_count" -ne "$r2_count" ]; then
        warn "Unequal R1 ($r1_count) and R2 ($r2_count) files. Check for missing pairs."
    fi
    
    success "Found $r1_count paired-end samples"
    
    # Validate sample ID format
    echo "Validating sample IDs..."
    local invalid_ids=0
    local checked_ids=()
    
    # Check R1 files for proper naming
    while IFS= read -r file; do
        local basename=$(basename "$file")
        # Extract sample ID (remove _R1/_R2/_1/_2 and extensions)
        local sample_id=$(echo "$basename" | sed -E 's/_(R)?[12](_L[0-9]+)?\.f(ast)?q\.gz$//')
        
        # Check if ID follows recommended pattern: PROJECT-DIGITS or PROJECT-DIGITS-REPLICATE
        if ! echo "$sample_id" | grep -qE '^[A-Z0-9]+-[0-9]{1,}(-[A-Z0-9]+)?$'; then
            if [ "$invalid_ids" -eq 0 ]; then
                warn "Some sample IDs don't follow recommended naming convention:"
                warn "  Expected: PROJECT-DIGITS[-REPLICATE] (e.g., KMUN-001, STUDY-042-A)"
                warn "  See docs/SAMPLE_ID_STANDARDS.md for details"
                echo "  Non-standard samples:"
            fi
            echo "    - $sample_id"
            invalid_ids=$((invalid_ids + 1))
        fi
        
        checked_ids+=("$sample_id")
    done < <(find "$input_dir" -maxdepth 1 \( -name "*_R1*.fastq.gz" -o -name "*_R1*.fq.gz" -o -name "*_1.fastq.gz" -o -name "*_1.fq.gz" \) | head -20)
    
    if [ "$invalid_ids" -gt 0 ]; then
        warn "Found $invalid_ids samples with non-standard naming (showing up to 20)"
        echo "  Pipeline will continue, but consider renaming for consistency."
    else
        success "Sample ID format validation passed"
    fi
    
    # Check for empty files
    local empty_files=$(find "$input_dir" -maxdepth 1 \( -name "*.fastq.gz" -o -name "*.fq.gz" \) -size 0)
    if [ -n "$empty_files" ]; then
        warn "Empty FASTQ files detected:"
        echo "$empty_files" | sed 's/^/  /'
    fi
}

# === Check 5: Disk Space ===
check_disk_space() {
    local output_dir="$1"
    
    echo ""
    echo "Checking disk space..."
    
    if [ -z "$output_dir" ]; then
        output_dir="."
    fi
    
    local available_gb=$(df -BG "$output_dir" | tail -1 | awk '{print $4}' | sed 's/G//')
    
    if [ "$available_gb" -lt 10 ]; then
        warn "Low disk space: ${available_gb}GB available. Recommend at least 10GB."
    else
        success "Disk space: ${available_gb}GB available"
    fi
}

# === Check 6: Memory ===
check_memory() {
    echo ""
    echo "Checking system memory..."
    
    if command -v free &> /dev/null; then
        local total_mem_gb=$(free -g | awk '/^Mem:/{print $2}')
        local avail_mem_gb=$(free -g | awk '/^Mem:/{print $7}')
        
        if [ "$total_mem_gb" -lt 8 ]; then
            warn "System has ${total_mem_gb}GB RAM. Recommend at least 8GB for small datasets."
        else
            success "System memory: ${total_mem_gb}GB total, ${avail_mem_gb}GB available"
        fi
    else
        warn "Cannot check memory (free command not available)"
    fi
}

# === Main Execution ===
main() {
    echo "=========================================="
    echo "16S Metagenomics Pipeline Pre-flight Check"
    echo "=========================================="
    echo ""
    
    local input_dir=""
    local output_dir="output"
    
    # Parse arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            --input-dir)
                input_dir="$2"
                shift 2
                ;;
            --output-dir)
                output_dir="$2"
                shift 2
                ;;
            *)
                shift
                ;;
        esac
    done
    
    # Run all checks
    check_conda_env
    check_tools
    check_r_packages
    check_input_files "$input_dir"
    check_disk_space "$output_dir"
    check_memory
    
    echo ""
    echo "=========================================="
    success "All pre-flight checks passed!"
    echo "=========================================="
    echo ""
}

# Run main function with all arguments
main "$@"
