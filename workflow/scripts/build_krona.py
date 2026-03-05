#!/usr/bin/env python3
"""
build_krona.py — Publication-ready multi-level Krona HTML generator
=====================================================================
Produces:
  krona_all_samples.html              — all samples merged (+ per-sample tabs)
  krona_level{L}_{name}.html          — one HTML per --collapse entry

Usage
-----
python3 build_krona.py \
    --taxonomy  exported/taxonomy/taxonomy.tsv \
    --table     exported/feature_table_filtered/feature-table.tsv \
    --metadata  config/metadata.tsv \
    --alpha-dir exported/alpha_diversity \
    --lefse     differential/lefse_results.tsv \
    --out-dir   visualizations/composition/krona \
    --krona-env "" \
    --per-sample \
    --collapse 2:phylum:composition/2_relfreq.tsv \
    --collapse 3:class:composition/3_relfreq.tsv \
    --collapse 6:genus:composition/6_relfreq.tsv

Requirements
------------
  pip install pandas
  mamba install -n qiime2 -c bioconda krona
  ktUpdateTaxonomy.sh   # run once after installing krona
"""

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

# ── pandas is the only external dependency ───────────────────────────────────
try:
    import pandas as pd
except ImportError:
    sys.exit("ERROR: pandas is required — pip install pandas")


# =============================================================================
# Helper: run ktImportText with graceful fallback
# =============================================================================

def run_ktImportText(input_files: list[tuple[str, str]],
                     output_html: str,
                     krona_env: str = "") -> bool:
    """
    Call ktImportText to generate a Krona HTML.

    Parameters
    ----------
    input_files : list of (tab_file, label) tuples
    output_html : destination HTML path
    krona_env   : conda env that contains ktImportText; "" = try PATH first

    Returns True on success, False if ktImportText is not available
    (a placeholder HTML is written in that case so Snakemake targets are met).
    """
    # Build argument list:  tab_file,label  (Krona tab-separated input format)
    kt_args = []
    for tab_file, label in input_files:
        kt_args.append(f"{tab_file},{label}")

    cmd_base = ["ktImportText", "-o", output_html] + kt_args

    def _run(cmd):
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"[Krona] WARN: {result.stderr[:500]}", file=sys.stderr)
        return result.returncode == 0

    # 1) Try PATH
    try:
        if _run(cmd_base):
            return True
    except FileNotFoundError:
        pass

    # 2) Try mamba run -n krona_env
    if krona_env:
        mamba = _find_mamba()
        if mamba:
            try:
                if _run([mamba, "run", "-n", krona_env] + cmd_base):
                    return True
            except FileNotFoundError:
                pass

    # 3) Placeholder HTML so downstream Snakemake targets are satisfied
    print(
        "[Krona] WARNING: ktImportText not found. Writing placeholder HTML.\n"
        "Install: mamba install -n qiime2 -c bioconda krona && ktUpdateTaxonomy.sh",
        file=sys.stderr,
    )
    _write_placeholder_html(output_html, kt_args)
    return False


def _find_mamba():
    for candidate in ["mamba", "micromamba", "conda"]:
        r = subprocess.run(["which", candidate], capture_output=True, text=True)
        if r.returncode == 0:
            return r.stdout.strip()
    return None


def _write_placeholder_html(path: str, labels: list[str]):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as fh:
        fh.write(
            "<!DOCTYPE html><html><body>"
            "<h2>Krona not available</h2>"
            "<p>Install: <code>mamba install -n qiime2 -c bioconda krona</code> "
            "then re-run the pipeline.</p>"
            f"<p>Expected datasets: {', '.join(labels)}</p>"
            "</body></html>"
        )


# =============================================================================
# Loaders
# =============================================================================

def load_taxonomy(tsv_path: str) -> dict[str, list[str]]:
    """
    Returns {feature_id: [rank_name, ...]} from a QIIME 2 taxonomy TSV.
    Ranks are split on ';' and stripped of SILVA prefixes (d__, p__, etc.).
    """
    df = pd.read_csv(tsv_path, sep="\t", index_col=0, comment="#")
    taxonomy = {}
    for fid, row in df.iterrows():
        raw = str(row.get("Taxon", ""))
        ranks = [r.strip() for r in raw.split(";")]
        # Strip SILVA-style prefixes d__, p__, c__, o__, f__, g__, s__
        clean = []
        for r in ranks:
            if len(r) >= 3 and r[1:3] == "__":
                r = r[3:]
            if r and r not in ("", "unidentified", "uncultured"):
                clean.append(r)
        taxonomy[str(fid)] = clean if clean else ["Unclassified"]
    return taxonomy


def load_feature_table(tsv_path: str) -> tuple[list[str], dict[str, dict[str, float]]]:
    """
    Returns (sample_ids, {feature_id: {sample_id: count}}) from a BIOM-derived TSV.
    The first row in QIIME 2 BIOM exports is '# Constructed from biom file'.
    """
    df = pd.read_csv(tsv_path, sep="\t", skiprows=1, index_col=0)
    samples = list(df.columns)
    table: dict[str, dict[str, float]] = {}
    for fid, row in df.iterrows():
        table[str(fid)] = {sid: float(row[sid]) for sid in samples if row[sid] > 0}
    return samples, table


def load_alpha_diversity(alpha_dir: str,
                         metrics: list[str] | None = None) -> dict[str, dict[str, float]]:
    """
    Returns {sample_id: {metric: value}} by scanning subdirectories of alpha_dir.
    If metrics is None, scans all sub-directories named <metric>/alpha-diversity.tsv.
    """
    base = Path(alpha_dir)
    if not base.exists():
        return {}
    data: dict[str, dict[str, float]] = {}
    for sub in base.iterdir():
        if not sub.is_dir():
            continue
        f = sub / "alpha-diversity.tsv"
        if not f.exists():
            continue
        if metrics and sub.name not in metrics:
            continue
        try:
            df = pd.read_csv(f, sep="\t", index_col=0)
            col = df.columns[0]
            for sid, val in df[col].items():
                data.setdefault(str(sid), {})[sub.name] = float(val)
        except Exception as exc:
            print(f"[Krona] WARN loading {f}: {exc}", file=sys.stderr)
    return data


def load_lefse(lefse_tsv: str | None) -> dict[str, tuple[float, str]]:
    """
    Returns {taxon_name: (lda_score, direction)} from a LEfSe results TSV.
    Expected columns: feature | lda_score | direction | p_value   (or similar)
    Gracefully returns {} if file does not exist or cannot be parsed.
    """
    if not lefse_tsv or not os.path.exists(lefse_tsv):
        return {}
    try:
        df = pd.read_csv(lefse_tsv, sep="\t")
        lda_map: dict[str, tuple[float, str]] = {}
        # Flexible column detection
        feat_col  = next((c for c in df.columns if c.lower() in ("feature", "taxon", "taxa")), df.columns[0])
        lda_col   = next((c for c in df.columns if "lda" in c.lower()), None)
        dir_col   = next((c for c in df.columns if "dir" in c.lower() or "group" in c.lower()), None)
        if lda_col is None:
            return {}
        for _, row in df.iterrows():
            name  = str(row[feat_col])
            score = float(row[lda_col]) if pd.notna(row[lda_col]) else 0.0
            direc = str(row[dir_col]) if (dir_col and pd.notna(row[dir_col])) else "?"
            lda_map[name] = (score, direc)
        return lda_map
    except Exception as exc:
        print(f"[Krona] WARN loading LEfSe {lefse_tsv}: {exc}", file=sys.stderr)
        return {}


def load_collapsed_tsv(tsv_path: str) -> tuple[list[str], dict[str, dict[str, float]]]:
    """
    Load a relative-frequency TSV produced by QIIME 2 biom convert.
    Returns (sample_ids, {taxon_name: {sample_id: rel_abundance}}).
    """
    df = pd.read_csv(tsv_path, sep="\t", skiprows=1, index_col=0)
    samples = list(df.columns)
    data: dict[str, dict[str, float]] = {}
    for taxon, row in df.iterrows():
        data[str(taxon)] = {sid: float(row[sid]) for sid in samples if row[sid] > 0}
    return samples, data


def load_metadata_groups(metadata_path: str,
                          group_col: str | None = None) -> dict[str, str]:
    """Returns {sample_id: group_value} for the first non-#SampleID column (or group_col)."""
    try:
        df = pd.read_csv(metadata_path, sep="\t", index_col=0, comment="#")
        if not group_col:
            group_col = df.columns[0] if len(df.columns) > 0 else None
        if group_col and group_col in df.columns:
            return {str(sid): str(df.loc[sid, group_col]) for sid in df.index
                    if sid in df.index}
    except Exception:
        pass
    return {}


# =============================================================================
# Writers
# =============================================================================

def write_krona_text(fpath: str,
                     feature_table: dict[str, dict[str, float]],
                     taxonomy: dict[str, list[str]],
                     sample_filter: list[str] | None = None):
    """
    Write a ktImportText-compatible tab file from raw ASV feature table.
    Format: count [TAB] rank1 [TAB] rank2 ...  (one row per feature × sample)
    sample_filter=None → aggregate all samples.
    """
    os.makedirs(os.path.dirname(fpath) or ".", exist_ok=True)
    with open(fpath, "w") as fh:
        for fid, sample_counts in feature_table.items():
            ranks = taxonomy.get(fid, ["Unclassified"])
            if sample_filter is not None:
                total = sum(sample_counts.get(sid, 0) for sid in sample_filter)
            else:
                total = sum(sample_counts.values())
            if total <= 0:
                continue
            fh.write("\t".join([str(int(total))] + ranks) + "\n")


def write_krona_collapsed(fpath: str,
                           collapsed: dict[str, dict[str, float]],
                           name: str,
                           sample_filter: list[str] | None = None,
                           lda_map: dict[str, tuple[float, str]] | None = None):
    """
    Write a ktImportText tab file from a collapsed relative-frequency table.
    Ranks are obtained by splitting the taxon string on ';'.
    If lda_map is provided, the last rank leaf is annotated with the LDA score.
    """
    os.makedirs(os.path.dirname(fpath) or ".", exist_ok=True)
    lda_map = lda_map or {}
    with open(fpath, "w") as fh:
        for taxon, sample_vals in collapsed.items():
            if sample_filter is not None:
                total = sum(sample_vals.get(sid, 0.0) for sid in sample_filter)
            else:
                total = sum(sample_vals.values())
            if total <= 0:
                continue

            # Split taxon into ranks and clean SILVA prefixes
            parts = [p.strip() for p in str(taxon).split(";")]
            clean = []
            for p in parts:
                if len(p) >= 3 and p[1:3] == "__":
                    p = p[3:]
                if p:
                    clean.append(p)
            if not clean:
                clean = ["Unclassified"]

            # Annotate last leaf with LDA score if available
            last = clean[-1]
            if last in lda_map:
                score, direction = lda_map[last]
                arrow = "↑" if direction not in ("", "?") else "~"
                clean[-1] = f"{last} [{arrow} LDA={score:.2f}]"

            # Use relative abundance × 10000 to avoid float rounding in Krona
            fh.write("\t".join([str(int(total * 10000))] + clean) + "\n")


def make_sample_title(sid: str,
                      alpha_data: dict[str, dict[str, float]],
                      group: str | None = None) -> str:
    """Compose a descriptive tab label: SampleID (Group=X, H'=4.21, ASVs=312)."""
    parts = []
    if group:
        parts.append(f"Group={group}")
    a = alpha_data.get(sid, {})
    if "shannon" in a:
        parts.append(f"H'={a['shannon']:.2f}")
    if "observed_features" in a:
        parts.append(f"ASVs={int(a['observed_features'])}")
    if parts:
        return f"{sid} ({', '.join(parts)})"
    return sid


# =============================================================================
# Main
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(
        description="Build multi-level annotated Krona HTML visualizations."
    )
    p.add_argument("--taxonomy",   required=True,  help="Taxonomy TSV (QIIME 2 export)")
    p.add_argument("--table",      required=True,  help="Feature table TSV (BIOM export)")
    p.add_argument("--metadata",   required=True,  help="Sample metadata TSV")
    p.add_argument("--alpha-dir",  default="",     help="Directory with alpha diversity subdirs")
    p.add_argument("--lefse",      default="",     help="LEfSe results TSV (optional)")
    p.add_argument("--out-dir",    required=True,  help="Output directory for HTML files")
    p.add_argument("--krona-env",  default="",     help="Conda env with ktImportText")
    p.add_argument("--per-sample", action="store_true",
                   help="Generate one Krona tab per sample in merged HTML")
    p.add_argument("--collapse",   action="append", default=[],
                   metavar="LEVEL:NAME:TSV",
                   help="Collapsed table: level:human_name:path.tsv (repeatable)")
    return p.parse_args()


def main():
    args = parse_args()
    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)

    print("[Krona] Loading taxonomy ...", flush=True)
    taxonomy = load_taxonomy(args.taxonomy)

    print("[Krona] Loading feature table ...", flush=True)
    samples, feature_table = load_feature_table(args.table)

    print("[Krona] Loading alpha diversity ...", flush=True)
    alpha_data = load_alpha_diversity(args.alpha_dir) if args.alpha_dir else {}

    print("[Krona] Loading LEfSe results ...", flush=True)
    lda_map = load_lefse(args.lefse)
    if lda_map:
        print(f"[Krona]   → {len(lda_map)} LEfSe taxa loaded")

    print("[Krona] Loading metadata groups ...", flush=True)
    groups = load_metadata_groups(args.metadata)

    # ── 1. All-samples Krona (merged + optional per-sample tabs) ──────────────
    print("[Krona] Building all-samples Krona ...", flush=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        tab_inputs: list[tuple[str, str]] = []

        # Merged tab (sum of all samples)
        merged_tab = os.path.join(tmpdir, "all_merged.txt")
        write_krona_text(merged_tab, feature_table, taxonomy)
        tab_inputs.append((merged_tab, "All Samples"))

        # Per-sample tabs
        if args.per_sample:
            for sid in samples:
                tab_file = os.path.join(tmpdir, f"{sid}.txt")
                write_krona_text(tab_file, feature_table, taxonomy,
                                 sample_filter=[sid])
                label = make_sample_title(sid, alpha_data, groups.get(sid))
                tab_inputs.append((tab_file, label))

        all_html = os.path.join(out_dir, "krona_all_samples.html")
        run_ktImportText(tab_inputs, all_html, args.krona_env)
        print(f"[Krona]   → {all_html}")

    # ── 2. Per-level collapsed Krona ──────────────────────────────────────────
    for entry in args.collapse:
        parts = entry.split(":", 2)
        if len(parts) != 3:
            print(f"[Krona] WARN: invalid --collapse value '{entry}' (expected LEVEL:NAME:TSV)",
                  file=sys.stderr)
            continue
        level, name, tsv_path = parts
        if not os.path.exists(tsv_path):
            print(f"[Krona] WARN: collapsed TSV not found: {tsv_path}", file=sys.stderr)
            # Write placeholder so Snakemake output is satisfied
            placeholder = os.path.join(out_dir, f"krona_level{level}_{name}.html")
            _write_placeholder_html(placeholder, [name])
            continue

        print(f"[Krona] Building level {level} ({name}) Krona ...", flush=True)
        coll_samples, collapsed = load_collapsed_tsv(tsv_path)

        with tempfile.TemporaryDirectory() as tmpdir:
            tab_inputs: list[tuple[str, str]] = []

            merged_tab = os.path.join(tmpdir, "merged.txt")
            write_krona_collapsed(merged_tab, collapsed, name,
                                   lda_map=lda_map)
            tab_inputs.append((merged_tab, f"All — {name}"))

            if args.per_sample:
                for sid in coll_samples:
                    tab_file = os.path.join(tmpdir, f"{sid}.txt")
                    write_krona_collapsed(tab_file, collapsed, name,
                                           sample_filter=[sid],
                                           lda_map=lda_map)
                    label = make_sample_title(sid, alpha_data, groups.get(sid))
                    tab_inputs.append((tab_file, label))

            level_html = os.path.join(out_dir, f"krona_level{level}_{name}.html")
            run_ktImportText(tab_inputs, level_html, args.krona_env)
            print(f"[Krona]   → {level_html}")

    print("[Krona] Done.", flush=True)


if __name__ == "__main__":
    main()
