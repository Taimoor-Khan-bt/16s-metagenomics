#!/bin/bash
# Generate provenance information for pipeline run
# Records git commit, timestamps, versions, and parameters

set -euo pipefail

# Parse arguments
OUTPUT_FILE=""
CONFIG_FILE=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --output)
            OUTPUT_FILE="$2"
            shift 2
            ;;
        --config)
            CONFIG_FILE="$2"
            shift 2
            ;;
        *)
            shift
            ;;
    esac
done

if [ -z "$OUTPUT_FILE" ]; then
    echo "Error: --output required"
    exit 1
fi

# Create provenance file
cat > "$OUTPUT_FILE" << EOF
# Provenance Information
# Generated: $(date -u +"%Y-%m-%dT%H:%M:%SZ")
# This file records metadata about this pipeline execution for reproducibility

execution:
  timestamp_start: "$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  user: "$USER"
  hostname: "$(hostname)"
  working_directory: "$(pwd)"
  conda_environment: "${CONDA_DEFAULT_ENV:-none}"

system:
  os: "$(uname -s)"
  os_version: "$(uname -r)"
  architecture: "$(uname -m)"
  cpu_cores: $(nproc 2>/dev/null || echo "unknown")
  total_memory_gb: $(free -g 2>/dev/null | awk '/^Mem:/{print $2}' || echo "unknown")

git:
  repository: "$(git remote get-url origin 2>/dev/null || echo 'not in git repo')"
  branch: "$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo 'unknown')"
  commit_hash: "$(git rev-parse HEAD 2>/dev/null || echo 'unknown')"
  commit_date: "$(git show -s --format=%ci HEAD 2>/dev/null || echo 'unknown')"
  uncommitted_changes: $(git diff --quiet 2>/dev/null && echo "false" || echo "true")
  status: |
$(git status --short 2>/dev/null | sed 's/^/    /' || echo '    unknown')

software_versions:
  fastqc: "$(fastqc --version 2>&1 | head -n1 || echo 'not found')"
  cutadapt: "$(cutadapt --version 2>&1 || echo 'not found')"
  multiqc: "$(multiqc --version 2>&1 | grep -oP 'version \K[0-9.]+' || echo 'not found')"
  r_version: "$(Rscript --version 2>&1 | head -n1 || echo 'not found')"
EOF

# Add R package versions
cat >> "$OUTPUT_FILE" << 'EOF'
  r_packages: |
EOF

Rscript --vanilla - >> "$OUTPUT_FILE" 2>/dev/null << 'REOF' || echo "    Error getting R package versions" >> "$OUTPUT_FILE"
packages <- c("yaml", "dada2", "phyloseq", "Biostrings", "ggplot2", "dplyr", "tidyr", "vegan")
for (pkg in packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        version <- packageVersion(pkg)
        cat(sprintf("    %s: %s\n", pkg, version))
    }
}
REOF

# Add config file info if provided
if [ -n "$CONFIG_FILE" ] && [ -f "$CONFIG_FILE" ]; then
    cat >> "$OUTPUT_FILE" << EOF

configuration:
  config_file: "$CONFIG_FILE"
  config_md5: "$(md5sum "$CONFIG_FILE" | awk '{print $1}')"
  config_content: |
$(sed 's/^/    /' "$CONFIG_FILE")
EOF
fi

# Add pipeline script checksums for integrity
cat >> "$OUTPUT_FILE" << EOF

pipeline_checksums:
EOF

for script in run_pipeline.sh scripts/*.sh scripts/*.R; do
    if [ -f "$script" ]; then
        echo "  $(basename "$script"): $(md5sum "$script" | awk '{print $1}')" >> "$OUTPUT_FILE"
    fi
done

echo ""
echo "Provenance file generated: $OUTPUT_FILE"
