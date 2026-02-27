#!/usr/bin/env python3
"""
build_sample_sheet.py — Automatically generate samples.tsv by matching metadata IDs to FASTQ files.
Usage: python3 build_sample_sheet.py <metadata.tsv> <fastq_dir> <output_samples.tsv>
"""
import os
import sys
import csv
import re

def main():
    if len(sys.argv) != 4:
        print(f"Usage: python3 {sys.argv[0]} <metadata.tsv> <fastq_dir> <output_samples.tsv>")
        sys.exit(1)

    metadata_path = sys.argv[1]
    fastq_dir = sys.argv[2]
    out_sheet = sys.argv[3]

    # 1. Read true SampleIDs from metadata
    sample_ids = []
    with open(metadata_path, "r") as f:
        # Sniff delimiter
        first_line = f.readline()
        f.seek(0)
        delim = "\t" if "\t" in first_line else ","
        
        reader = csv.DictReader(f, delimiter=delim)
        header = reader.fieldnames
        
        # Guess the ID column
        id_cols = ["#SampleID", "SampleID", "sample_id", "id", "Code"]
        col_name = None
        for c in id_cols:
            if c in header:
                col_name = c
                break
        if not col_name:
            col_name = header[0] # Fallback to first column
            print(f"Warning: Did not find standard ID column, using first column '{col_name}'")

        for row in reader:
            sid = row.get(col_name, "").strip()
            if sid and not sid.startswith("#"):
                sample_ids.append(sid)
    
    print(f"Loaded {len(sample_ids)} Sample IDs from {metadata_path}")

    # 2. Get all fastq files in the directory
    try:
        files = os.listdir(fastq_dir)
    except FileNotFoundError:
        print(f"Error: Directory {fastq_dir} not found.")
        sys.exit(1)

    r1_files = sorted([f for f in files if f.endswith((".fq.gz", ".fastq.gz", ".fq", ".fastq")) and ("_1." in f or "_R1" in f)])
    
    # 3. Match SampleIDs to R1 FASTQ files
    # A simple greedy match: if the FASTQ starts with the SampleID or SampleID is in the FASTQ name
    matches = {}
    missing = []
    
    for sid in sample_ids:
        # Sort by length ascending so we pick the cleanest match, or just find the first that starts with it
        matched_file = None
        for r1 in r1_files:
            # e.g., KMUN-001 matches KMUN-001-LFM4957_L2_1.fq.gz
            # we check if r1 starts with sid + a non-alphanumeric character (like _, -, .) to avoid partial matches
            if r1.startswith(sid + "_") or r1.startswith(sid + "-") or r1.startswith(sid + "."):
                matched_file = r1
                break
        
        if not matched_file:
            # fallback: just is it in the filename anywhere?
            for r1 in r1_files:
                if sid in r1:
                    matched_file = r1
                    break

        if matched_file:
            r1_path = os.path.join(fastq_dir, matched_file)
            # Try to infer R2 based on R1 name
            r2_name = matched_file.replace("_1.", "_2.").replace("_R1", "_R2")
            r2_path = os.path.join(fastq_dir, r2_name)
            
            if not os.path.exists(r2_path):
                r2_path = ""
            
            matches[sid] = {"r1": r1_path, "r2": r2_path}
        else:
            missing.append(sid)

    # 4. Write samples.tsv
    with open(out_sheet, "w") as f:
        f.write("# =============================================================================\n")
        f.write(f"# Auto-generated Paired-End sample sheet from {fastq_dir}\n")
        f.write("# =============================================================================\n")
        f.write("sample_id\tr1\tr2\n")
        for sid in sample_ids:
            if sid in matches:
                r1 = matches[sid]["r1"]
                r2 = matches[sid]["r2"]
                f.write(f"{sid}\t{r1}\t{r2}\n")

    print(f"\nSuccessfully wrote {len(matches)} matches to {out_sheet}")
    if missing:
        print(f"WARNING: Could not find FASTQ files for {len(missing)} samples:")
        for m in missing:
            print(f"  - {m}")

if __name__ == "__main__":
    main()
