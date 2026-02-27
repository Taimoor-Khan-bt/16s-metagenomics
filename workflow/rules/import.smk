# =============================================================================
# rules/import.smk — Step 1–2: Build manifest + import reads into QIIME 2
# =============================================================================

rule build_manifest:
    """
    Build a QIIME 2 PairedEndFastqManifestPhred33V2 TSV from samples.tsv.
    Each row maps a sample_id to its absolute R1 and R2 FASTQ paths.
    """
    input:
        samples = config["samples"],
    output:
        manifest = f"{OUT}/manifest.tsv",
    run:
        import pandas as pd, os
        samples = pd.read_table(input.samples, comment="#").set_index("sample_id")
        os.makedirs(os.path.dirname(output.manifest), exist_ok=True)
        with open(output.manifest, "w") as fh:
            fh.write("sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")
            for sid, row in samples.iterrows():
                r1 = os.path.abspath(row["r1"])
                r2 = os.path.abspath(row.get("r2", ""))
                if not os.path.exists(r1):
                    raise FileNotFoundError(f"R1 not found for {sid}: {r1}")
                if not r2 or not os.path.exists(r2):
                    raise FileNotFoundError(f"R2 not found for {sid}: {r2}")
                fh.write(f"{sid}\t{r1}\t{r2}\n")
        print(f"Manifest written: {len(samples)} samples → {output.manifest}")


rule import_reads:
    """
    Import FASTQ reads into a QIIME 2 artifact (SampleData[PairedEndSequencesWithQuality]).
    Uses PairedEnd manifest format (Phred33 quality scores).
    """
    input:
        manifest = f"{OUT}/manifest.tsv",
    output:
        seqs = f"{OUT}/sequences.qza",
    params:
        docker = DOCKER,
    log:
        f"{OUT}/logs/import_reads.log",
    shell:
        """
        mkdir -p $(dirname {log})
        {params.docker} qiime tools import \
            --type 'SampleData[PairedEndSequencesWithQuality]' \
            --input-path '{input.manifest}' \
            --input-format PairedEndFastqManifestPhred33V2 \
            --output-path '{output.seqs}' \
            2>&1 | tee {log}
        """
