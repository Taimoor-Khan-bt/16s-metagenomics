# =============================================================================
# rules/taxonomy.smk — Step 4: Taxonomic classification (SILVA NB classifier)
# =============================================================================

rule classify_taxonomy:
    """
    Assign taxonomy to ASVs using a pre-trained naive-Bayes classifier.
    The classifier file is configured in config/config.yaml (classifier key).

    Recommended classifiers (download once):
      Full-length:  refs/silva-138-99-nb-classifier.qza          (208 MB) ← default
      V4 specific:  refs/silva-138-99-515-806-nb-classifier.qza  (125 MB)

    Download:
      wget https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza \
           -O refs/silva-138-99-nb-classifier.qza
    """
    input:
        rep_seqs   = f"{OUT}/rep_seqs.qza",
        classifier = config["classifier"],
    output:
        taxonomy = f"{OUT}/taxonomy.qza",
    params:
        docker  = DOCKER,
        threads = config["dada2"]["threads"],   # reuse thread count
    log:
        f"{OUT}/logs/classify_taxonomy.log",
    shell:
        """
        mkdir -p $(dirname {log})
        {params.docker} qiime feature-classifier classify-sklearn \
            --i-classifier      '{input.classifier}' \
            --i-reads           '{input.rep_seqs}' \
            --p-n-jobs          {params.threads} \
            --o-classification  '{output.taxonomy}' \
            --verbose \
            2>&1 | tee {log}
        """


# =============================================================================
# filter_taxa — Remove mitochondria, chloroplast and unassigned features.
# Creates table_filtered.qza which is used by ALL downstream analyses. The
# original table.qza is preserved for raw QC comparisons only.
# Exclude list is driven by config taxa_processing.exclude_list (case-insensitive).
# =============================================================================

rule filter_taxa:
    """
    Filter out unwanted taxa (mitochondria, chloroplast, unassigned) from the
    raw feature table using QIIME 2 taxa filter-table.

    Output: table_filtered.qza — used by diversity, composition, and
    differential abundance rules. table.qza remains untouched for QC.
    """
    input:
        table    = f"{OUT}/table.qza",
        taxonomy = f"{OUT}/taxonomy.qza",
    output:
        filtered = f"{OUT}/table_filtered.qza",
    params:
        docker      = DOCKER,
        exclude_str = ",".join(config.get("taxa_processing", {}).get(
                          "exclude_list", ["mitochondria", "chloroplast", "unassigned"]
                      )),
    log:
        f"{OUT}/logs/filter_taxa.log",
    shell:
        """
        mkdir -p $(dirname {log})
        {params.docker} qiime taxa filter-table \
            --i-table         '{input.table}' \
            --i-taxonomy      '{input.taxonomy}' \
            --p-exclude       '{params.exclude_str}' \
            --p-mode          contains \
            --o-filtered-table '{output.filtered}' \
            --verbose \
            2>&1 | tee {log}
        """
