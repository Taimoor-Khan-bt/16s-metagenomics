#!/usr/bin/env bash
# =============================================================================
# setup.sh — Complete one-shot setup for the QIIME 2 Snakemake pipeline
# =============================================================================
# Checks and installs all prerequisites:
#   1. Docker + QIIME 2 image
#   2. Snakemake (in 'qiime2' conda env)
#   3. R packages for extended analysis (A–G)
#   4. KRONA (ktImportText) in qiime2 env
#   5. PICRUSt2 conda environment
#   6. Config file validation
#   7. Pipeline dry-run
#
# Usage: bash setup.sh [--skip-picrust] [--skip-r] [--skip-krona]
# =============================================================================

set -euo pipefail

GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; BLUE='\033[0;34m'; NC='\033[0m'
ok()      { echo -e "${GREEN}[✓]${NC} $*"; }
warn()    { echo -e "${YELLOW}[!]${NC} $*"; }
fail()    { echo -e "${RED}[✗]${NC} $*"; exit 1; }
section() { echo -e "\n${BLUE}══${NC} $* ${BLUE}══${NC}"; }

# ── Parse flags ───────────────────────────────────────────────────────────────
SKIP_PICRUST=false
SKIP_R=false
SKIP_KRONA=false
for arg in "$@"; do
  case $arg in
    --skip-picrust) SKIP_PICRUST=true ;;
    --skip-r)       SKIP_R=true ;;
    --skip-krona)   SKIP_KRONA=true ;;
  esac
done

echo "============================================================"
echo "  QIIME 2 Snakemake Pipeline — Full Setup"
echo "  $(date '+%Y-%m-%d %H:%M')"
echo "============================================================"

# ─────────────────────────────────────────────────────────────────────────────
# 1. DOCKER
# ─────────────────────────────────────────────────────────────────────────────
section "1 / 7  Docker"

if ! command -v docker &>/dev/null; then
  fail "Docker not found. Install: https://docs.docker.com/get-docker/"
fi
ok "Docker: $(docker --version)"

if ! docker info &>/dev/null 2>&1; then
  fail "Docker daemon not running. Please start Docker and re-run."
fi
ok "Docker daemon: running"

IMG="quay.io/qiime2/amplicon:2024.10"
if docker image inspect "$IMG" &>/dev/null; then
  ok "QIIME 2 image: $IMG (cached)"
else
  warn "Pulling QIIME 2 image (~6.5 GB)…"
  docker pull "$IMG"
  ok "QIIME 2 image: pulled"
fi

# Smoke test
VER=$(docker run --rm \
  --user "$(id -u):$(id -g)" \
  --volume "${HOME}:${HOME}" \
  --volume "/tmp:/tmp" \
  --workdir "$(pwd)" \
  "$IMG" qiime --version 2>&1 | head -1)
ok "QIIME 2 test: $VER"

# ─────────────────────────────────────────────────────────────────────────────
# 2. SNAKEMAKE
# ─────────────────────────────────────────────────────────────────────────────
section "2 / 7  Snakemake"

QIIME2_ENV="qiime2"

if conda run -n "$QIIME2_ENV" snakemake --version &>/dev/null 2>&1; then
  SM_VER=$(conda run -n "$QIIME2_ENV" snakemake --version 2>/dev/null)
  ok "Snakemake $SM_VER (env: $QIIME2_ENV)"
else
  warn "Snakemake not found in env '$QIIME2_ENV'. Creating / installing…"
  if conda env list | grep -q "^${QIIME2_ENV} "; then
    mamba install -n "$QIIME2_ENV" -c conda-forge -c bioconda snakemake -y
  else
    mamba create -n "$QIIME2_ENV" -c conda-forge -c bioconda snakemake python=3.12 -y
  fi
  ok "Snakemake: $(conda run -n "$QIIME2_ENV" snakemake --version)"
fi

# ─────────────────────────────────────────────────────────────────────────────
# 3. R PACKAGES
# ─────────────────────────────────────────────────────────────────────────────
section "3 / 7  R packages"

if $SKIP_R; then
  warn "Skipping R packages (--skip-r)"
else
  # Check if R is available in the env
  if ! conda run -n "$QIIME2_ENV" Rscript --version &>/dev/null 2>&1; then
    warn "R not found in env '$QIIME2_ENV'. Installing R…"
    mamba install -n "$QIIME2_ENV" -c conda-forge r-base -y
  fi

  R_VERSION=$(conda run -n "$QIIME2_ENV" Rscript --version 2>&1 | head -1)
  ok "R: $R_VERSION"

  CRAN_PKGS=(
    r-ggplot2 r-dplyr r-tidyr r-vegan r-ape
    r-rcolorbrewer r-scales r-upsetr r-ggpubr r-broom
    r-compositions r-ggrepel
  )

  warn "Installing CRAN packages: ${CRAN_PKGS[*]}"
  mamba install -n "$QIIME2_ENV" -c conda-forge "${CRAN_PKGS[@]}" -y
  ok "CRAN R packages installed"

  # microbiomeMarker via BiocManager (bioconductor) — try conda first
  if mamba install -n "$QIIME2_ENV" -c conda-forge -c bioconda bioconductor-microbiomemarker -y 2>/dev/null; then
    ok "microbiomeMarker: installed via conda"
  else
    warn "microbiomeMarker not in conda channels — installing via BiocManager in R…"
    conda run -n "$QIIME2_ENV" Rscript -e "
      if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org')
      if (!requireNamespace('microbiomeMarker', quietly=TRUE)) BiocManager::install('microbiomeMarker', ask=FALSE)
      cat('microbiomeMarker:', as.character(packageVersion('microbiomeMarker')), '\n')
    " 2>&1 || warn "microbiomeMarker install failed — LEfSe will use Wilcoxon fallback"
  fi
fi

# ─────────────────────────────────────────────────────────────────────────────
# 4. KRONA
# ─────────────────────────────────────────────────────────────────────────────
section "4 / 7  KRONA"

if $SKIP_KRONA; then
  warn "Skipping KRONA (--skip-krona)"
else
  if conda run -n "$QIIME2_ENV" ktImportText --version &>/dev/null 2>&1 || \
     conda run -n "$QIIME2_ENV" which ktImportText &>/dev/null 2>&1; then
    ok "KRONA: already installed in env '$QIIME2_ENV'"
  else
    warn "Installing KRONA in env '$QIIME2_ENV'…"
    mamba install -n "$QIIME2_ENV" -c bioconda krona -y
    # Initialize taxonomy database
    conda run -n "$QIIME2_ENV" ktUpdateTaxonomy.sh 2>/dev/null || true
    ok "KRONA: installed"
  fi
fi

# ─────────────────────────────────────────────────────────────────────────────
# 5. PICRUST2
# ─────────────────────────────────────────────────────────────────────────────
section "5 / 7  PICRUSt2"

if $SKIP_PICRUST; then
  warn "Skipping PICRUSt2 (--skip-picrust)"
else
  PICRUST_ENV=$(python3 -c "
import yaml, sys
try:
    with open('config/config.yaml') as f:
        c = yaml.safe_load(f)
    print(c.get('picrust', {}).get('env', 'picrust2'))
except: print('picrust2')
" 2>/dev/null || echo "picrust2")

  if conda env list | grep -q "^${PICRUST_ENV} "; then
    if conda run -n "$PICRUST_ENV" picrust2_pipeline.py --version &>/dev/null 2>&1; then
      ok "PICRUSt2: env '$PICRUST_ENV' already exists"
    else
      warn "PICRUSt2 env exists but picrust2 not working — reinstalling…"
      mamba install -n "$PICRUST_ENV" -c conda-forge -c bioconda picrust2 -y
      ok "PICRUSt2: reinstalled"
    fi
  else
    warn "Creating PICRUSt2 conda env '$PICRUST_ENV'…"
    mamba create -n "$PICRUST_ENV" -c conda-forge -c bioconda picrust2 python=3.8 -y
    ok "PICRUSt2: env '$PICRUST_ENV' created"
  fi
fi

# ─────────────────────────────────────────────────────────────────────────────
# 6. CONFIG VALIDATION
# ─────────────────────────────────────────────────────────────────────────────
section "6 / 7  Config validation"

ALL_OK=true
for f in config/config.yaml config/samples.tsv config/metadata.tsv; do
  if [[ -f "$f" ]]; then
    ok "Found: $f"
  else
    warn "Missing: $f"
    ALL_OK=false
  fi
done

CLF="refs/silva-138-99-nb-classifier.qza"
if [[ -f "$CLF" ]]; then
  ok "Classifier: $CLF ($(du -sh "$CLF" | cut -f1))"
else
  warn "Classifier not found: $CLF"
  warn "Download command:"
  warn "  mkdir -p refs"
  warn "  wget https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza \\"
  warn "       -O refs/silva-138-99-nb-classifier.qza"
  ALL_OK=false
fi

# Check that sample FASTQ files exist
if [[ -f "config/samples.tsv" ]]; then
  python3 -c "
import sys
missing = []
with open('config/samples.tsv') as f:
    for i, line in enumerate(f):
        line = line.strip()
        if i == 0 or not line or line.startswith('#'): continue
        parts = line.split('\t')
        if len(parts) >= 2:
            import os
            r1 = parts[1]
            if not os.path.exists(r1):
                missing.append(f'  {parts[0]}: {r1}')
if missing:
    print('WARNING: Missing FASTQ files:')
    for m in missing: print(m)
else:
    print('All FASTQ files found.')
" 2>&1 | while read -r line; do
    [[ "$line" == WARNING* ]] && warn "$line" || ok "$line"
  done
fi

# ─────────────────────────────────────────────────────────────────────────────
# 7. DRY RUN
# ─────────────────────────────────────────────────────────────────────────────
section "7 / 7  Pipeline dry-run"

if $ALL_OK; then
  echo "Running base pipeline dry-run…"
  if conda run -n "$QIIME2_ENV" snakemake -n --snakefile workflow/Snakefile 2>&1 | \
     grep -E "(Job stats|total|Error)" | tail -5; then
    ok "Base pipeline dry-run PASSED"
  else
    warn "Base pipeline dry-run reported issues (see output above)"
  fi

  echo ""
  echo "Running extended analysis dry-run…"
  if conda run -n "$QIIME2_ENV" snakemake all_extended -n --snakefile workflow/Snakefile 2>&1 | \
     grep -E "(Job stats|total|Error)" | tail -5; then
    ok "Extended pipeline dry-run PASSED"
  else
    warn "Extended dry-run reported issues (see output above)"
  fi
else
  warn "Skipping dry-run — fix missing files above first"
fi

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo -e "${GREEN}  Setup complete!${NC}"
echo ""
echo "  Activate env:   mamba activate $QIIME2_ENV"
echo ""
echo "  Base pipeline:  snakemake --cores all --snakefile workflow/Snakefile"
echo "  Full analysis:  snakemake all_extended --cores all --snakefile workflow/Snakefile"
echo "  Dry-run:        snakemake -n --snakefile workflow/Snakefile"
echo ""
echo "  Flags:"
echo "    --skip-r         skip R packages"
echo "    --skip-krona     skip KRONA install"
echo "    --skip-picrust   skip PICRUSt2 env"
echo "============================================================"
