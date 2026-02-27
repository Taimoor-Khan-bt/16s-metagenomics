#!/usr/bin/env bash
# =============================================================================
# preflight.sh — Read-only preflight check for the QIIME 2 Snakemake pipeline
# =============================================================================
# Checks every tool, package, and file WITHOUT installing anything.
# Reports PASS / WARN / FAIL for each item and prints a fix summary at the end.
#
# Usage: bash preflight.sh
# =============================================================================

set -uo pipefail

GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'
BLUE='\033[0;34m'; BOLD='\033[1m'; NC='\033[0m'

PASS=0; WARN=0; FAIL=0
FIXES=()

pass() { echo -e "  ${GREEN}[PASS]${NC} $*";  ((PASS++));  }
warn() { echo -e "  ${YELLOW}[WARN]${NC} $*"; ((WARN++));  }
fail() { echo -e "  ${RED}[FAIL]${NC} $*";   ((FAIL++));  }
section() { echo -e "\n${BLUE}${BOLD}── $* ──${NC}"; }
fix()  { FIXES+=("$*"); }

QIIME2_ENV="qiime2"
PICRUST_ENV=$(python3 -c "
import yaml
try:
    with open('config/config.yaml') as f:
        c = yaml.safe_load(f)
    print(c.get('picrust',{}).get('env','picrust2'))
except: print('picrust2')
" 2>/dev/null || echo "picrust2")

echo ""
echo -e "${BOLD}============================================================${NC}"
echo -e "${BOLD}  QIIME 2 Pipeline — Preflight Check${NC}"
echo -e "  $(date '+%Y-%m-%d %H:%M')"
echo -e "${BOLD}============================================================${NC}"

# ─────────────────────────────────────────────────────────────────────────────
section "System tools"
# ─────────────────────────────────────────────────────────────────────────────

# Docker
if command -v docker &>/dev/null; then
  pass "docker: $(docker --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+')"
  if docker info &>/dev/null 2>&1; then
    pass "docker daemon: running"
  else
    fail "docker daemon: NOT running"
    fix "Start Docker Desktop or: sudo systemctl start docker"
  fi
else
  fail "docker: not installed"
  fix "Install Docker: https://docs.docker.com/get-docker/"
fi

# mamba / conda
if command -v mamba &>/dev/null; then
  pass "mamba: $(mamba --version 2>/dev/null | head -1)"
elif command -v conda &>/dev/null; then
  warn "mamba not found — conda available (mamba is faster)"
  fix "Install mamba: conda install -n base -c conda-forge mamba -y"
else
  fail "conda/mamba: not installed"
  fix "Install miniforge: https://github.com/conda-forge/miniforge"
fi

# python3
if command -v python3 &>/dev/null; then
  pass "python3: $(python3 --version 2>/dev/null)"
else
  fail "python3: not found"
fi

# ─────────────────────────────────────────────────────────────────────────────
section "QIIME 2 Docker image"
# ─────────────────────────────────────────────────────────────────────────────

IMG="quay.io/qiime2/amplicon:2024.10"
if docker image inspect "$IMG" &>/dev/null 2>&1; then
  SIZE=$(docker image inspect "$IMG" --format='{{.Size}}' 2>/dev/null | \
         awk '{printf "%.1f GB", $1/1073741824}')
  pass "image: $IMG ($SIZE)"
  # Quick functional test
  VER=$(docker run --rm \
    --user "$(id -u):$(id -g)" \
    --volume "${HOME}:${HOME}" \
    "$IMG" qiime --version 2>/dev/null | head -1 || echo "")
  if [[ -n "$VER" ]]; then
    pass "functional test: $VER"
  else
    warn "image present but functional test failed (permissions?)"
    fix "Try: docker run --rm $IMG qiime --version"
  fi
else
  fail "QIIME 2 image: not pulled"
  fix "docker pull $IMG   (~6.5 GB)"
fi

# ─────────────────────────────────────────────────────────────────────────────
section "Snakemake (env: $QIIME2_ENV)"
# ─────────────────────────────────────────────────────────────────────────────

if conda env list 2>/dev/null | grep -q "^${QIIME2_ENV}[[:space:]]"; then
  pass "conda env: $QIIME2_ENV exists"
  if conda run -n "$QIIME2_ENV" snakemake --version &>/dev/null 2>&1; then
    SM_VER=$(conda run -n "$QIIME2_ENV" snakemake --version 2>/dev/null)
    pass "snakemake: $SM_VER"
    SM_MAJOR=$(echo "$SM_VER" | grep -oE '^[0-9]+')
    [[ "${SM_MAJOR:-0}" -ge 8 ]] || warn "snakemake ≥8.0 recommended (found $SM_VER)"
  else
    fail "snakemake: not installed in env $QIIME2_ENV"
    fix "mamba install -n $QIIME2_ENV -c conda-forge snakemake -y"
  fi
else
  fail "conda env '$QIIME2_ENV': does not exist"
  fix "mamba create -n $QIIME2_ENV -c conda-forge snakemake python=3.12 -y"
fi

# ─────────────────────────────────────────────────────────────────────────────
section "R packages (env: $QIIME2_ENV)"
# ─────────────────────────────────────────────────────────────────────────────

check_r_pkg() {
  local pkg="$1"
  local result
  # Use installed.packages() — requireNamespace() can return FALSE for installed
  # packages that have broken optional dependencies.
  result=$(conda run -n "$QIIME2_ENV" Rscript -e \
    "cat('$pkg' %in% rownames(installed.packages()))" 2>/dev/null || echo "FALSE")
  if [[ "$result" == "TRUE" ]]; then
    VER=$(conda run -n "$QIIME2_ENV" Rscript -e \
      "tryCatch(cat(as.character(packageVersion('$pkg'))), error=function(e) cat('installed'))" 2>/dev/null || echo "installed")
    pass "R/$pkg: $VER"
  else
    fail "R/$pkg: not installed"
    fix "mamba install -n $QIIME2_ENV -c conda-forge r-$(echo "$pkg" | tr '[:upper:]' '[:lower:]') -y"
  fi
}


if conda run -n "$QIIME2_ENV" Rscript --version &>/dev/null 2>&1; then
  R_VER=$(conda run -n "$QIIME2_ENV" Rscript --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)
  pass "R: $R_VER"
  
  CRAN_PKGS=(ggplot2 dplyr tidyr vegan ape compositions RColorBrewer scales UpSetR ggpubr broom ggrepel)
  for pkg in "${CRAN_PKGS[@]}"; do
    check_r_pkg "$pkg"
  done
  
  # microbiomeMarker (Bioconductor)
  # Use installed.packages() — requireNamespace() can return FALSE if a
  # dependency is broken even when the package itself is installed.
  MM_CHECK=$(conda run -n "$QIIME2_ENV" Rscript -e \
    "cat('microbiomeMarker' %in% rownames(installed.packages()))" 2>/dev/null || echo "FALSE")
  if [[ "$MM_CHECK" == "TRUE" ]]; then
    MM_VER=$(conda run -n "$QIIME2_ENV" Rscript -e \
      "tryCatch(cat(as.character(packageVersion('microbiomeMarker'))), error=function(e) cat('installed'))" 2>/dev/null || echo "installed")
    pass "R/microbiomeMarker: $MM_VER"
  else
    warn "R/microbiomeMarker: not installed (LEfSe will use Wilcoxon fallback)"
    fix "In R: BiocManager::install('microbiomeMarker')"
    fix "  OR: mamba install -n $QIIME2_ENV -c bioconda bioconductor-microbiomemarker -y"
  fi
else
  fail "R: not installed in env $QIIME2_ENV"
  fix "mamba install -n $QIIME2_ENV -c conda-forge r-base -y"
fi

# ─────────────────────────────────────────────────────────────────────────────
section "KRONA (env: $QIIME2_ENV)"
# ─────────────────────────────────────────────────────────────────────────────

if conda run -n "$QIIME2_ENV" which ktImportText &>/dev/null 2>&1; then
  pass "KRONA (ktImportText): installed"
else
  warn "KRONA: not installed (F: KRONA plots will fail)"
  fix "mamba install -n $QIIME2_ENV -c bioconda krona -y && conda run -n $QIIME2_ENV ktUpdateTaxonomy.sh"
fi

# ─────────────────────────────────────────────────────────────────────────────
section "PICRUSt2 (env: $PICRUST_ENV)"
# ─────────────────────────────────────────────────────────────────────────────

if conda env list 2>/dev/null | grep -q "^${PICRUST_ENV}[[:space:]]"; then
  pass "conda env: $PICRUST_ENV exists"
  if conda run -n "$PICRUST_ENV" picrust2_pipeline.py --version &>/dev/null 2>&1; then
    PC_VER=$(conda run -n "$PICRUST_ENV" picrust2_pipeline.py --version 2>/dev/null || echo "?")
    pass "PICRUSt2: $PC_VER"
  else
    fail "PICRUSt2: env exists but not functional"
    fix "mamba install -n $PICRUST_ENV -c conda-forge -c bioconda picrust2 -y"
  fi
else
  warn "PICRUSt2 env '$PICRUST_ENV': not found (G: functional prediction will fail)"
  fix "mamba create -n $PICRUST_ENV -c conda-forge -c bioconda picrust2 python=3.8 -y"
fi

# ─────────────────────────────────────────────────────────────────────────────
section "Config & input files"
# ─────────────────────────────────────────────────────────────────────────────

for f in config/config.yaml config/samples.tsv config/metadata.tsv; do
  if [[ -f "$f" ]]; then
    LINES=$(wc -l < "$f")
    pass "$f ($LINES lines)"
  else
    fail "$f: missing"
    fix "Copy or create $f (see README.md)"
  fi
done

# Classifier
CLF=$(python3 -c "
import yaml
try:
    c = yaml.safe_load(open('config/config.yaml'))
    print(c.get('classifier','refs/silva-138-99-nb-classifier.qza'))
except: print('refs/silva-138-99-nb-classifier.qza')
" 2>/dev/null || echo "refs/silva-138-99-nb-classifier.qza")

if [[ -f "$CLF" ]]; then
  pass "classifier: $CLF ($(du -sh "$CLF" 2>/dev/null | cut -f1))"
else
  fail "classifier: $CLF not found"
  fix "mkdir -p refs && wget https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza -O $CLF"
fi

# Sample FASTQ files
if [[ -f "config/samples.tsv" ]]; then
  MISSING_FQ=0
  while IFS=$'\t' read -r sid r1; do
    [[ "$sid" == "sample_id" || "$sid" =~ ^# || -z "$sid" ]] && continue
    if [[ -f "$r1" ]]; then
      SIZE=$(du -sh "$r1" 2>/dev/null | cut -f1)
      pass "FASTQ $sid: $r1 ($SIZE)"
    else
      fail "FASTQ $sid: $r1 — NOT FOUND"
      fix "Check path for sample $sid in config/samples.tsv"
      ((MISSING_FQ++))
    fi
  done < config/samples.tsv
fi

# Metadata columns
if [[ -f "config/metadata.tsv" && -f "config/config.yaml" ]]; then
  python3 -c "
import yaml, csv, sys

with open('config/config.yaml') as f:
    cfg = yaml.safe_load(f)

group_col = cfg.get('analysis', {}).get('group_column') or \
            cfg.get('metadata', {}).get('group_column', '')
covariates = cfg.get('analysis', {}).get('covariates', [])

with open('config/metadata.tsv') as f:
    reader = csv.DictReader(f, delimiter='\t')
    header = reader.fieldnames or []
    rows = list(reader)
    sample_col = header[0] if header else '#SampleID'

checks = []
if group_col and group_col not in header:
    checks.append(f'MISSING group_column \"{group_col}\" in metadata (found: {header})')
elif group_col:
    vals = [r[group_col] for r in rows if r.get(group_col)]
    n_groups = len(set(vals))
    checks.append(f'OK: group_column \"{group_col}\" has {n_groups} unique values: {set(vals)}')
for cov in covariates:
    if cov not in header:
        checks.append(f'MISSING covariate \"{cov}\" in metadata')
    else:
        checks.append(f'OK: covariate \"{cov}\" present')
for c in checks:
    print(c)
" 2>/dev/null | while read -r line; do
    if [[ "$line" == OK:* ]]; then
      pass "${line#OK: }"
    elif [[ "$line" == MISSING* ]]; then
      fail "$line"
      fix "Add column '$(echo "$line" | grep -oP '\".*?\"' | head -1 | tr -d '\"')' to config/metadata.tsv"
    fi
  done
fi

# ─────────────────────────────────────────────────────────────────────────────
section "Pipeline dry-run"
# ─────────────────────────────────────────────────────────────────────────────

if conda run -n "$QIIME2_ENV" snakemake --version &>/dev/null 2>&1; then
  echo "  Running base pipeline dry-run…"
  conda run -n "$QIIME2_ENV" snakemake -n --snakefile workflow/Snakefile >/dev/null 2>&1
  if [[ $? -eq 0 ]]; then
    pass "Base pipeline dry-run: OK"
  else
    ERR=$(conda run -n "$QIIME2_ENV" snakemake -n --snakefile workflow/Snakefile 2>&1 | \
          grep -E "Error|Exception|Missing" | head -3)
    fail "Base pipeline dry-run: FAILED"
    [[ -n "$ERR" ]] && fix "Snakefile error: $ERR"
  fi

  echo "  Running extended pipeline dry-run…"
  conda run -n "$QIIME2_ENV" snakemake all_extended -n --snakefile workflow/Snakefile >/dev/null 2>&1
  if [[ $? -eq 0 ]]; then
    pass "Extended pipeline dry-run: OK"
  else
    ERR=$(conda run -n "$QIIME2_ENV" snakemake all_extended -n --snakefile workflow/Snakefile 2>&1 | \
          grep -E "Error|Exception|Missing" | head -3)
    fail "Extended pipeline dry-run: FAILED"
    [[ -n "$ERR" ]] && fix "Snakefile error: $ERR"
  fi

else
  warn "Skipping dry-run (Snakemake not available)"
fi

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────
TOTAL=$((PASS + WARN + FAIL))

echo ""
echo -e "${BOLD}============================================================${NC}"
echo -e "  Results: ${GREEN}${PASS} passed${NC}  |  ${YELLOW}${WARN} warnings${NC}  |  ${RED}${FAIL} failed${NC}  (${TOTAL} checks)"
echo -e "${BOLD}============================================================${NC}"

if [[ ${#FIXES[@]} -gt 0 ]]; then
  echo ""
  echo -e "${BOLD}  Suggested fixes:${NC}"
  for i in "${!FIXES[@]}"; do
    echo -e "  $((i+1)). ${FIXES[$i]}"
  done
  echo ""
fi

if [[ $FAIL -eq 0 && $WARN -eq 0 ]]; then
  echo -e "  ${GREEN}${BOLD}All checks passed — pipeline is ready to run!${NC}"
elif [[ $FAIL -eq 0 ]]; then
  echo -e "  ${YELLOW}No failures — warnings are optional. Pipeline should run.${NC}"
else
  echo -e "  ${RED}${FAIL} failure(s) detected — fix the items above before running.${NC}"
  echo -e "  Run ${BOLD}bash setup.sh${NC} to automatically fix most issues."
fi
echo ""
exit $([[ $FAIL -eq 0 ]] && echo 0 || echo 1)
