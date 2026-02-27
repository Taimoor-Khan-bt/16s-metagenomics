#!/usr/bin/env python3
"""
make_krona_input.py — Convert taxonomy TSV + feature table TSV to KRONA input
Usage: python3 make_krona_input.py <taxonomy_tsv> <feature_table_tsv> <out_txt>

The KRONA input format is tab-separated: count  [taxon1]  [taxon2]  ...
"""
import sys, csv

taxonomy_file = sys.argv[1]
table_file    = sys.argv[2]
out_file      = sys.argv[3]

# Load taxonomy: Feature ID → list of taxon levels
taxonomy = {}
with open(taxonomy_file) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        fid  = row["Feature ID"]
        tax  = row["Taxon"].replace("_", " ")
        # Split on "; " and strip d__/p__/c__ prefixes
        levels = []
        for t in tax.split("; "):
            t = t.strip()
            if "__" in t:
                t = t.split("__", 1)[1].strip()
            if t and t != "":
                levels.append(t)
        taxonomy[fid] = levels

# Load feature table (skip first comment line "#Constructed from ...")
with open(table_file) as f:
    lines = f.readlines()

# Skip lines starting with #
data_lines = [l for l in lines if not l.startswith("#")]
reader = csv.DictReader(data_lines, delimiter="\t")
samples = [c for c in reader.fieldnames if c != "OTU ID" and c != "#OTU ID"]

# Aggregate counts per sample and write one KRONA input file (summed over samples)
with open(out_file, "w") as out:
    for row in reader:
        fid   = row.get("OTU ID") or row.get("#OTU ID") or list(row.values())[0]
        total = sum(float(row.get(s, 0)) for s in samples)
        count = int(total)
        if count == 0:
            continue
        tax_levels = taxonomy.get(fid, [fid])
        out.write(str(count) + "\t" + "\t".join(tax_levels) + "\n")

print(f"KRONA input written to: {out_file}")
