#!/usr/bin/env python3
"""
merge_ancombc_csvs.py — Merge ANCOM-BC CSV exports into a single TSV
Usage: python3 merge_ancombc_csvs.py <export_dir> <output_tsv>
"""
import os, csv, sys

export_dir = sys.argv[1]
output_tsv = sys.argv[2]

rows = []
header = None

csv_files = sorted(f for f in os.listdir(export_dir) if f.endswith(".csv"))
if not csv_files:
    print(f"WARNING: No CSV files found in {export_dir}")
    open(output_tsv, "w").close()
    sys.exit(0)

for fn in csv_files:
    with open(os.path.join(export_dir, fn)) as f:
        reader = csv.DictReader(f)
        if header is None:
            header = list(reader.fieldnames) + ["source_file"]
            rows.append(header)
        for row in reader:
            row["source_file"] = fn
            rows.append([row.get(h, "") for h in header])

with open(output_tsv, "w") as f:
    for r in rows:
        f.write("\t".join(r) + "\n")

n_rows = len(rows) - 1
print(f"Exported {n_rows} rows from {len(csv_files)} CSV files → {output_tsv}")
