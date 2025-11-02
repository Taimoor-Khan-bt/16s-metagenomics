#!/bin/bash
# Integration Test Script
# Tests the entire pipeline end-to-end with synthetic data

set -euo pipefail

echo "=========================================="
echo "16S Pipeline Integration Test"
echo "=========================================="
echo ""

# Colors
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

# Test counters
TESTS_RUN=0
TESTS_PASSED=0
TESTS_FAILED=0

# Helper functions
pass() {
    echo -e "${GREEN}✓${NC} $1"
    TESTS_PASSED=$((TESTS_PASSED + 1))
    TESTS_RUN=$((TESTS_RUN + 1))
}

fail() {
    echo -e "${RED}✗${NC} $1"
    TESTS_FAILED=$((TESTS_FAILED + 1))
    TESTS_RUN=$((TESTS_RUN + 1))
}

check_file() {
    if [ -f "$1" ]; then
        pass "Found: $1"
        return 0
    else
        fail "Missing: $1"
        return 1
    fi
}

# Step 1: Generate test data
echo "Step 1: Generating test dataset..."
Rscript scripts/generate_test_data.R
if [ $? -eq 0 ]; then
    pass "Test data generated"
else
    fail "Failed to generate test data"
    exit 1
fi
echo ""

# Step 2: Validate test dataset
echo "Step 2: Validating test dataset..."
check_file "test_data/metadata.csv"
check_file "test_data/fastq/TEST-001_1.fastq.gz"
check_file "test_data/fastq/TEST-001_2.fastq.gz"
echo ""

# Step 3: Run preflight checks
echo "Step 3: Running preflight checks..."
./scripts/preflight_check.sh --input-dir "test_data/fastq" --output-dir "output"
if [ $? -eq 0 ]; then
    pass "Preflight checks passed"
else
    fail "Preflight checks failed"
fi
echo ""

# Step 4: Validate configuration
echo "Step 4: Validating configuration..."
Rscript scripts/validate_config.R config/test_config.yaml
if [ $? -eq 0 ]; then
    pass "Configuration validation passed"
else
    fail "Configuration validation failed"
fi
echo ""

# Step 5: Run trimming
echo "Step 5: Running adapter trimming..."
./cutadapt-16s-trim.sh "test_data/fastq" "test_data/trimmed"
if [ $? -eq 0 ]; then
    pass "Trimming completed"
    check_file "test_data/trimmed/TEST-001_1_trimmed.fastq.gz"
else
    fail "Trimming failed"
fi
echo ""

# Step 6: Run R analysis pipeline
echo "Step 6: Running R analysis pipeline..."
export USE_MODULAR_VIZ=true
Rscript scripts/runner.R --config config/test_config.yaml 2>&1 | tee test_pipeline.log
if [ $? -eq 0 ]; then
    pass "R analysis completed"
else
    fail "R analysis failed"
    echo "Check test_pipeline.log for details"
fi
echo ""

# Step 7: Validate outputs
echo "Step 7: Validating pipeline outputs..."
check_file "output/phyloseq_object_raw.rds"
check_file "output/phyloseq_rarefied.rds"
check_file "output/asv_sequences.fasta"
check_file "output/alpha_diversity.csv"
check_file "output/test/visualizations/alpha_shannon_primary.tiff"
check_file "output/test/visualizations/ordination_pcoa_bray_primary.tiff"
echo ""

# Step 8: Check provenance
echo "Step 8: Checking provenance tracking..."
PROVENANCE_FILE=$(ls -t output/provenance_*.yaml 2>/dev/null | head -1)
if [ -n "$PROVENANCE_FILE" ]; then
    pass "Provenance file created: $PROVENANCE_FILE"
    
    # Check key fields
    if grep -q "timestamp_start" "$PROVENANCE_FILE"; then
        pass "Provenance contains timestamp"
    fi
    if grep -q "git:" "$PROVENANCE_FILE"; then
        pass "Provenance contains git info"
    fi
    if grep -q "software_versions:" "$PROVENANCE_FILE"; then
        pass "Provenance contains software versions"
    fi
else
    fail "No provenance file found"
fi
echo ""

# Summary
echo "=========================================="
echo "Integration Test Summary"
echo "=========================================="
echo "Tests run: $TESTS_RUN"
echo -e "Tests passed: ${GREEN}$TESTS_PASSED${NC}"
if [ $TESTS_FAILED -gt 0 ]; then
    echo -e "Tests failed: ${RED}$TESTS_FAILED${NC}"
else
    echo -e "Tests failed: $TESTS_FAILED"
fi
echo ""

if [ $TESTS_FAILED -eq 0 ]; then
    echo -e "${GREEN}✓ All tests passed!${NC}"
    echo ""
    exit 0
else
    echo -e "${RED}✗ Some tests failed${NC}"
    echo "Check test_pipeline.log for details"
    echo ""
    exit 1
fi
