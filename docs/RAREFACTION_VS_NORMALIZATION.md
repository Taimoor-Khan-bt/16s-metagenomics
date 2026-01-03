# Rarefaction vs Alternative Normalization Methods

## Overview

This document discusses rarefaction and alternative normalization approaches for 16S rRNA gene sequencing data, helping you make informed decisions about which method is appropriate for your analysis.

## What is Rarefaction?

**Rarefaction** is a normalization method that randomly subsamples all samples to the same sequencing depth (number of reads). This addresses the issue of unequal library sizes across samples.

### How It Works

1. Determine rarefaction depth (e.g., 5,000 reads)
2. For each sample with ≥ 5,000 reads:
   - Randomly sample exactly 5,000 reads without replacement
   - Discard samples below the threshold
3. All retained samples now have equal depth

### Advantages

✓ **Simple and intuitive** - Easy to understand and explain  
✓ **Removes library size bias** - All samples have equal weight  
✓ **Works with any diversity metric** - Compatible with phylogenetic diversity  
✓ **No distributional assumptions** - Non-parametric  
✓ **Conservative** - Reduces false positives from depth differences  

### Disadvantages

✗ **Discards data** - Throws away reads from deeper-sequenced samples  
✗ **Loss of power** - Reduces statistical power by reducing information  
✗ **Arbitrary threshold** - Choice of depth can be subjective  
✗ **Not ideal for differential abundance** - Better methods exist for this specific task  

## Alternative Normalization Methods

### 1. DESeq2 (Size Factor Normalization)

**Best for:** Differential abundance testing of individual taxa

**Method:** Uses median-of-ratios to calculate sample-specific size factors

**Advantages:**
- Retains all data
- Accounts for compositionality
- Robust to outliers
- Built-in statistical testing

**Disadvantages:**
- Assumes most taxa are not differentially abundant
- Requires count data (not relative abundance)
- More complex to explain

**When to use:**
- Identifying differentially abundant taxa
- When you have good sequencing depth across samples
- For count-based statistical modeling

### 2. CSS (Cumulative Sum Scaling)

**Best for:** Samples with highly variable library sizes

**Method:** Normalizes to a quantile of the count distribution

**Advantages:**
- Robust to differential abundance
- Works with sparse data
- Preserves information

**Disadvantages:**
- Less intuitive than rarefaction
- Implementation primarily in metagenomeSeq package

**When to use:**
- Highly uneven sequencing depths
- Sparse microbiome data
- When many taxa are differentially abundant

### 3. Relative Abundance (Proportions)

**Best for:** Compositional data analysis

**Method:** Convert counts to proportions (divide by total reads per sample)

**Advantages:**
- Simple
- Accounts for different library sizes
- Natural for compositional data

**Disadvantages:**
- Compositionality issues (parts of a whole)
- Not appropriate for diversity metrics
- Can amplify noise in low-abundance taxa

**When to use:**
- Visualizing taxonomic composition
- Compositional data analysis frameworks (ALDEx2, ANCOM)
- Not recommended for diversity metrics

### 4. TMM/RLE (edgeR methods)

**Best for:** RNA-seq-inspired normalizations

**Method:** Similar philosophy to DESeq2 but different calculations

**Advantages:**
- Well-tested in RNA-seq context
- Handles compositionality
- Statistical framework included

**Disadvantages:**
- Designed for RNA-seq, not microbiome
- May not handle sparsity as well

## Decision Guide

### Use Rarefaction When:

✓ Calculating **alpha diversity** metrics (Shannon, Simpson, Observed)  
✓ Calculating **beta diversity** / community dissimilarity  
✓ You want a **simple, defensible** approach  
✓ Library sizes are **reasonably similar** (within 2-3x)  
✓ You have **adequate sequencing depth** (>5,000 reads/sample after QC)  
✓ **Interpretability** is important for your audience  

### Use DESeq2/CSS When:

✓ Performing **differential abundance testing**  
✓ Library sizes are **highly variable** (>5-10x difference)  
✓ You want to **retain all data**  
✓ You're testing **specific hypotheses** about individual taxa  
✓ You have **statistical modeling expertise**  

### Don't Use Either When:

✗ Samples have **very low read counts** (<1,000 reads)  
✗ Data is **highly sparse** (>90% zeros)  
✗ You haven't **removed contaminants** yet  

## Rarefaction in This Pipeline

### Current Implementation

The pipeline uses **rarefaction** as the default normalization for:
- Alpha diversity calculations
- Beta diversity / ordination
- General community comparisons

**Rarefaction depth** is determined by:
1. Examining read depth distribution after filtering
2. Choosing depth that retains most samples (typically 90-95%)
3. Balancing depth vs sample retention trade-off

### Configuration

```yaml
diversity:
  rarefaction:
    enabled: true
    depth: 5000  # Adjust based on your data
```

### Decision Workflow

The pipeline logs read depths after filtering:

```
Sample Read Depths (after filtering):
  Minimum: 2,431
  Q1: 8,234
  Median: 12,456
  Q3: 18,234
  Maximum: 45,123
```

**Recommendations:**
- If all samples > 5,000: Use 5,000 or higher
- If 90% > 10,000: Use 10,000
- If highly variable: Consider CSS or subsetting samples

## Justifying Your Choice

### For Rarefaction

**Justification template:**

> "We normalized samples by rarefying to X reads per sample to account for 
> differences in sequencing depth while maintaining simplicity and 
> interpretability. This depth was chosen to retain Y% of samples while 
> ensuring adequate coverage for diversity estimates. Rarefaction is 
> appropriate for our diversity analyses as it provides equal weight to all 
> samples and is compatible with phylogenetic diversity metrics."

**Citations:**
- McMurdie PJ & Holmes S (2014) Waste not, want not: why rarefying microbiome data is inadmissible. PLoS Comput Biol 10(4):e1003531
- Willis AD (2019) Rarefaction, alpha diversity, and statistics. Front Microbiol 10:2407
- Weiss S et al. (2017) Normalization and microbial differential abundance strategies depend upon data characteristics. Microbiome 5:27

### For Alternative Methods

**DESeq2 justification:**

> "For differential abundance testing, we used DESeq2 size factor normalization 
> rather than rarefaction to retain all sequencing information and properly model 
> count data with appropriate distributional assumptions (negative binomial)."

**Citation:**
- Love MI, Huber W, Anders S (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol 15:550

## Practical Recommendations

### 1. Use Multiple Approaches

Consider running key analyses with different normalizations:
- Rarefaction for diversity
- DESeq2 for differential abundance
- Relative abundance for visualization

### 2. Report Your Decision

Always report:
- Which normalization was used
- Rarefaction depth chosen
- How many samples were discarded
- Justification for the choice

### 3. Sensitivity Analysis

If controversial, test sensitivity:
- Try multiple rarefaction depths
- Compare rarefaction vs DESeq2 results
- Report concordance

### 4. Don't Mix Methods

⚠️ **Warning:** Don't mix normalizations in the same analysis
- Use rarefied data for diversity metrics
- Use non-rarefied counts for DESeq2
- Don't compare diversity from rarefied to differential abundance from DESeq2

## Pipeline Recommendations

### Default: Rarefaction for Diversity

The pipeline uses rarefaction as default because:

1. **Most users need diversity metrics** - Rarefaction is appropriate
2. **Simplicity** - Easier to explain to collaborators
3. **Conservative** - Reduces false positives from depth bias
4. **Widely accepted** - Standard in microbiome field
5. **Compatible** - Works with all diversity metrics

### When to Override

Consider alternative approaches if:
- Performing focused differential abundance analysis
- Library sizes vary >10-fold
- You have statistical modeling expertise
- Reviewers request specific methods

## Further Reading

**Pro-rarefaction:**
- Weiss S et al. (2017) Microbiome 5:27 - "Rarefying is valid for richness estimation"
- Willis AD (2019) Front Microbiol - "Rarefaction is appropriate for certain questions"

**Anti-rarefaction:**
- McMurdie & Holmes (2014) PLoS Comput Biol - "Rarefying is inadmissible"
- Fung TC et al. (2023) - Alternative methods can have more power

**Nuanced perspective:**
- Cameron ES et al. (2021) mSystems - "Context-dependent; different methods for different questions"

## Summary

| Approach | Best For | Advantages | Disadvantages |
|----------|----------|------------|---------------|
| **Rarefaction** | Diversity metrics | Simple, conservative | Discards data |
| **DESeq2** | Differential abundance | Retains data, statistical framework | Complex, assumes |
| **CSS** | Variable library sizes | Robust | Less common |
| **Proportions** | Visualization | Intuitive | Not for diversity |

**Bottom line:** Rarefaction is appropriate and defensible for alpha/beta diversity in most cases. Use DESeq2 or similar for differential abundance testing of specific taxa.
