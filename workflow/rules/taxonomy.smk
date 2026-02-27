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
