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