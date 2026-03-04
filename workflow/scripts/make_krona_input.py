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

# 1. Load Taxonomy
taxonomy = {}
with open(taxonomy_file) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        fid  = row["Feature ID"]
        tax  = row["Taxon"].replace("_", " ")
        levels = [t.split("__")[-1].strip() for t in tax.split(";") if t.strip()]
        taxonomy[fid] = levels

# 2. Process Feature Table
with open(table_file) as f:
    # Filter out only the metadata comments, but KEEP the line starting with #OTU ID
    # We do this by dropping any line that starts with # but NOT with #OTU ID
    clean_lines = []
    for line in f:
        if line.startswith("#") and not line.startswith("#OTU ID"):
            continue
        clean_lines.append(line)

# Use the cleaned lines for DictReader
reader = csv.DictReader(clean_lines, delimiter="\t")

# Identify Sample Names (everything except the ID column)
id_col = "#OTU ID" if "#OTU ID" in reader.fieldnames else "OTU ID"
samples = [c for c in reader.fieldnames if c != id_col]

# 3. Write Krona Input
with open(out_file, "w") as out:
    for row in reader:
        fid = row[id_col]
        try:
            # Only sum columns we know are samples
            total = sum(float(row[s]) for s in samples if row[s])
            count = int(total)
        except ValueError:
            # This skips lines that might still be misformatted
            continue

        if count <= 0:
            continue
            
        tax_levels = taxonomy.get(fid, [fid])
        out.write(f"{count}\t" + "\t".join(tax_levels) + "\n")

print(f"Success! KRONA input written to: {out_file}")
