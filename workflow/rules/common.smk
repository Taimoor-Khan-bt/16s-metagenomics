# workflow/rules/common.smk

import os
import subprocess

# Get IDs
UID = subprocess.check_output(["id", "-u"]).decode().strip()
GID = subprocess.check_output(["id", "-g"]).decode().strip()

WDIR = os.path.abspath(".")
HOME = os.path.expanduser("~")
IMG  = config["docker"]["image"]

# THE NEW STRATEGY: 
# 1. We mount the project folder (WDIR) and HOME.
# 2. We do NOT mount the host's /tmp to the container's /tmp.
# 3. We let Docker create a fresh internal /tmp for every run.
DOCKER = (
    f'docker run --rm '
    f'--user {UID}:{GID} ' 
    f'--volume "{HOME}:{HOME}" '
    f'--workdir "{WDIR}" '
    f'"{IMG}"'
)

OUT     = config["output_dir"]
OUT_VIZ = config.get("viz_dir", f"{OUT}/visualizations")  # all QZVs, PDFs, HTML plots

# Resolve Rscript from the same conda env as this snakemake process.
# Falls back to bare 'Rscript' if not found (e.g. Rscript already on PATH).
import shutil, sys as _sys
_snakemake_bin = os.path.dirname(os.path.abspath(_sys.executable))  # e.g. .../envs/qiime2/bin
_r_candidate   = os.path.join(_snakemake_bin, "Rscript")
RSCRIPT = _r_candidate if os.path.isfile(_r_candidate) else (shutil.which("Rscript") or "Rscript")