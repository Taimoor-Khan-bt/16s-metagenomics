# =============================================================================
# rules/common.smk — Shared helpers and Docker configuration
# =============================================================================

import os

# ── Resolved paths (evaluated once at parse time) ─────────────────────────────
WDIR = os.path.abspath(".")          # absolute project root (quoted in docker)
HOME = os.path.expanduser("~")       # user home, mounted at same path
IMG  = config["docker"]["image"]     # e.g. quay.io/qiime2/amplicon:2024.10

# ── Docker run prefix ──────────────────────────────────────────────────────────
# Mounted paths must be absolute so they work identically inside the container.
# --user preserves file ownership; --volume mounts home + tmp.
# Quoting WDIR handles paths with spaces (e.g. "16s metagenomics").
DOCKER = (
    f'docker run --rm '
    f'--user "$(id -u):$(id -g)" '
    f'--volume "{HOME}:{HOME}" '
    f'--volume "/tmp:/tmp" '
    f'--workdir "{WDIR}" '
    f'"{IMG}"'
)

# ── Output root (shortcut used in all rule files) ──────────────────────────────
OUT = config["output_dir"]
