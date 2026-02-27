# QIIME 2 Manual Pipeline Run — Example Data

**Samples:** KMUN-001, KMUN-004 | **Image:** `quay.io/qiime2/amplicon:2024.10` (q2cli 2024.10.1)  
**Input:** `data/example/` — 250 bp PE reads, R1 used as single-end  
**Classifier:** `refs/silva-138-99-nb-classifier.qza` (SILVA 138, full-length NB)  
**Output:** `output/example_qiime2/qiime2/`

---

## Shell Variables (set once per session)

```bash
export IMG="quay.io/qiime2/amplicon:2024.10"
export WDIR="/home/taimoor/taimoor-data/genomics/16s metagenomics"
export OUT="$WDIR/output/example_qiime2/qiime2"
export META="$WDIR/output/example_qiime2/metadata_qiime2.tsv"
export CLF="$WDIR/refs/silva-138-99-nb-classifier.qza"
export DIV="$OUT/core_diversity"
export THREADS=$(nproc)

# Docker helper (mounts $HOME at same path so all absolute paths work as-is)
dqiime() {
  docker run --rm \
    --user "$(id -u):$(id -g)" \
    --volume "${HOME}:${HOME}" \
    --volume "/tmp:/tmp" \
    --workdir "$WDIR" \
    "$IMG" qiime "$@"
}
```

---

## Step 0 — Verify Docker + QIIME 2

```bash
docker run --rm \
  --user "$(id -u):$(id -g)" \
  --volume "${HOME}:${HOME}" \
  --workdir "$WDIR" \
  "$IMG" qiime --version 2>&1
```
✅ `q2cli version 2024.10.1`

---

## Step 1 — Build Sample Manifest

```bash
mkdir -p "$OUT"
MANIFEST="$OUT/manifest.tsv"
printf "sample-id\tabsolute-filepath\n" > "$MANIFEST"
printf "KMUN-001\t%s/data/example/KMUN-001_1.fq.gz\n" "$WDIR" >> "$MANIFEST"
printf "KMUN-004\t%s/data/example/KMUN-004_1.fq.gz\n" "$WDIR" >> "$MANIFEST"
```
✅ `manifest.tsv` — 2 samples

---

## Step 2 — Import Reads

```bash
dqiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-path "$OUT/manifest.tsv" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --output-path "$OUT/sequences.qza"
```
✅ `sequences.qza` — 86 MB

---

## Step 3 — DADA2 Denoising

Truncates reads to 220 bp (from 250) to remove low-quality 3′ ends.

```bash
dqiime dada2 denoise-single \
  --i-demultiplexed-seqs "$OUT/sequences.qza" \
  --p-trunc-len 220 \
  --p-trim-left 0 \
  --p-n-threads "$THREADS" \
  --o-table                    "$OUT/table.qza" \
  --o-representative-sequences "$OUT/rep_seqs.qza" \
  --o-denoising-stats          "$OUT/dada2_stats.qza" \
  --verbose
```
✅ `table.qza` 35 KB | `rep_seqs.qza` 47 KB | `dada2_stats.qza` 13 KB | **637 ASVs**

---

## Step 4 — Taxonomic Classification

```bash
dqiime feature-classifier classify-sklearn \
  --i-classifier   "$CLF" \
  --i-reads        "$OUT/rep_seqs.qza" \
  --p-n-jobs       "$THREADS" \
  --o-classification "$OUT/taxonomy.qza" \
  --verbose
```
✅ `taxonomy.qza` — 118 KB

---

## Step 5 — Phylogenetic Tree

```bash
dqiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences        "$OUT/rep_seqs.qza" \
  --p-n-threads        "$THREADS" \
  --o-alignment        "$OUT/aligned_rep_seqs.qza" \
  --o-masked-alignment "$OUT/masked_aligned_rep_seqs.qza" \
  --o-tree             "$OUT/unrooted_tree.qza" \
  --o-rooted-tree      "$OUT/rooted_tree.qza" \
  --verbose
```
✅ `rooted_tree.qza` — 62 KB

---

## Step 6 — Summary Visualizations

```bash
# DADA2 statistics
dqiime metadata tabulate \
  --m-input-file "$OUT/dada2_stats.qza" \
  --o-visualization "$OUT/dada2_stats.qzv"

# Feature table summary
dqiime feature-table summarize \
  --i-table "$OUT/table.qza" \
  --m-sample-metadata-file "$META" \
  --o-visualization "$OUT/table_summary.qzv"

# Representative sequences
dqiime feature-table tabulate-seqs \
  --i-data "$OUT/rep_seqs.qza" \
  --o-visualization "$OUT/rep_seqs.qzv"
```
✅ `dada2_stats.qzv` 1.2 MB | `table_summary.qzv` 421 KB | `rep_seqs.qzv` 280 KB

---

## Step 7 — Check Feature Table Counts

```bash
mkdir -p "$OUT/exported/feature_table"
dqiime tools export \
  --input-path "$OUT/table.qza" \
  --output-path "$OUT/exported/feature_table"

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --volume "${HOME}:${HOME}" \
  --workdir "$WDIR" \
  "$IMG" biom summarize-table \
    -i "$OUT/exported/feature_table/feature-table.biom"
```
✅ Results:
```
Num samples: 2   Num observations: 637   Total count: 509,678
KMUN-001: 276,490 reads
KMUN-004: 233,188 reads
Min: 233,188  →  Rarefaction depth: 200,000
```

---

## Step 8 — Core Diversity Metrics

> ⚠️ Remove `core_diversity/` if it exists — QIIME 2 refuses to overwrite existing output directories.

```bash
rm -rf "$OUT/core_diversity"

dqiime diversity core-metrics-phylogenetic \
  --i-phylogeny         "$OUT/rooted_tree.qza" \
  --i-table             "$OUT/table.qza" \
  --p-sampling-depth    200000 \
  --m-metadata-file     "$META" \
  --p-n-jobs-or-threads "$THREADS" \
  --output-dir          "$OUT/core_diversity"
```
✅ 17 files: alpha vectors (faith_pd, shannon, evenness, observed_features), beta distance matrices, rarefied table, and 4 Emperor PCoA QZVs

> ⚠️ **Emperor QZVs appear blank with 2 samples.** PCoA produces at most `n-1` principal coordinates. With `n=2`, there is only 1 axis — the 3D plot has no meaningful content. This is expected and correct behaviour. Requires ≥3 samples for a useful Emperor visualization.

---

## Step 9 — Taxonomy Barplot

```bash
dqiime taxa barplot \
  --i-table         "$OUT/table.qza" \
  --i-taxonomy      "$OUT/taxonomy.qza" \
  --m-metadata-file "$META" \
  --o-visualization "$OUT/taxa_barplot.qzv"
```
✅ `taxa_barplot.qzv` — interactive composition view at all 7 taxonomic levels ⭐

---

## Step 10 — Alpha Diversity Correlation

`alpha-group-significance` requires **≥2 samples per group** and cannot run with 2 total samples.  
`alpha-correlation` (Spearman) works with any sample size and correlates diversity against continuous metadata variables (e.g. AgeMonths).

```bash
for METRIC in faith_pd shannon evenness observed_features; do
  dqiime diversity alpha-correlation \
    --i-alpha-diversity "$DIV/${METRIC}_vector.qza" \
    --m-metadata-file   "$META" \
    --p-method          spearman \
    --o-visualization   "$DIV/${METRIC}_correlation.qzv"
done
```
✅ 4 QZVs: `faith_pd_correlation.qzv`, `shannon_correlation.qzv`, `evenness_correlation.qzv`, `observed_features_correlation.qzv` (~334 KB each)

---

## Step 11 — Export Alpha Diversity Values

```bash
mkdir -p "$OUT/exported/alpha_diversity"
for METRIC in faith_pd shannon evenness observed_features; do
  dqiime tools export \
    --input-path  "$DIV/${METRIC}_vector.qza" \
    --output-path "$OUT/exported/alpha_diversity/${METRIC}"
done
```
✅ TSV per metric in `exported/alpha_diversity/`

**Alpha diversity values:**

| Sample | Faith's PD | Shannon | Pielou Evenness |
|--------|-----------|---------|-----------------|
| KMUN-001 (Sujad M 42 mo) | 27.98 | 3.63 | 0.465 |
| KMUN-004 (Jalwa F 48 mo) | 24.12 | 6.00 | 0.680 |

---

## Step 12 — Export All Artifacts to Flat Files

```bash
# Taxonomy TSV
dqiime tools export \
  --input-path "$OUT/taxonomy.qza" \
  --output-path "$OUT/exported/taxonomy"

# Phylogenetic tree (Newick)
dqiime tools export \
  --input-path "$OUT/rooted_tree.qza" \
  --output-path "$OUT/exported/tree"

# Representative sequences (FASTA)
dqiime tools export \
  --input-path "$OUT/rep_seqs.qza" \
  --output-path "$OUT/exported/rep_seqs"

# DADA2 stats TSV
dqiime tools export \
  --input-path "$OUT/dada2_stats.qza" \
  --output-path "$OUT/exported/dada2_stats"

# Convert BIOM to human-readable TSV
docker run --rm \
  --user "$(id -u):$(id -g)" \
  --volume "${HOME}:${HOME}" \
  --workdir "$WDIR" \
  "$IMG" biom convert \
    -i "$OUT/exported/feature_table/feature-table.biom" \
    -o "$OUT/exported/feature_table/feature-table.tsv" \
    --to-tsv
```
✅ Flat files:
```
exported/feature_table/feature-table.tsv   ← ASV counts per sample
exported/feature_table/feature-table.biom
exported/taxonomy/taxonomy.tsv             ← Feature ID | Taxon | Confidence
exported/tree/tree.nwk                     ← Newick format
exported/rep_seqs/dna-sequences.fasta
exported/dada2_stats/stats.tsv
exported/alpha_diversity/*/alpha-diversity.tsv
```

---

## Viewing Visualizations

Upload any `.qzv` to **[https://view.qiime2.org](https://view.qiime2.org)**

| QZV | Description | Works with 2 samples? |
|-----|-------------|----------------------|
| `taxa_barplot.qzv` | Interactive taxonomy by level | ✅ Yes |
| `dada2_stats.qzv` | Denoising filter summary | ✅ Yes |
| `table_summary.qzv` | Read depth distribution | ✅ Yes |
| `rep_seqs.qzv` | BLAST-searchable ASV sequences | ✅ Yes |
| `core_diversity/*_correlation.qzv` | Alpha diversity vs metadata | ✅ Yes |
| `core_diversity/*_emperor.qzv` | 3D beta diversity PCoA | ❌ Blank (need ≥3 samples) |
| `core_diversity/*_significance.qzv` | Group diversity tests | ❌ Need ≥2 samples/group |

---

## Output Summary

```
output/example_qiime2/qiime2/
├── sequences.qza          86 MB
├── table.qza              35 KB   637 ASVs
├── rep_seqs.qza           47 KB
├── dada2_stats.qza        13 KB
├── taxonomy.qza          118 KB
├── rooted_tree.qza        62 KB
├── dada2_stats.qzv         1.2 MB  ⭐
├── table_summary.qzv     421 KB   ⭐
├── rep_seqs.qzv          280 KB   ⭐
├── taxa_barplot.qzv              ⭐ main visualization
├── core_diversity/
│   ├── *_emperor.qzv ×4         (blank with 2 samples)
│   ├── *_correlation.qzv ×4     ⭐ alpha diversity plots
│   ├── *_vector.qza ×4
│   └── *_distance_matrix.qza ×4
└── exported/
    ├── feature_table/feature-table.tsv
    ├── taxonomy/taxonomy.tsv
    ├── tree/tree.nwk
    ├── rep_seqs/dna-sequences.fasta
    ├── dada2_stats/stats.tsv
    └── alpha_diversity/*/alpha-diversity.tsv
```

---

*Next: wrap into `scripts/run_qiime2_example.sh` automation*
