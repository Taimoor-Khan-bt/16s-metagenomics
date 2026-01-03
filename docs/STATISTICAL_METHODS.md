# Statistical Methods Documentation

## Overview

This document provides comprehensive documentation of all statistical methods used in the 16S metagenomics pipeline, including their assumptions, appropriate use cases, and interpretation guidelines.

---

## Table of Contents

1. [Alpha Diversity Statistics](#alpha-diversity-statistics)
2. [Beta Diversity Statistics](#beta-diversity-statistics)
3. [Multiple Testing Correction](#multiple-testing-correction)
4. [Assumptions and Diagnostics](#assumptions-and-diagnostics)
5. [Effect Sizes and Power](#effect-sizes-and-power)
6. [Reporting Guidelines](#reporting-guidelines)

---

## Alpha Diversity Statistics

### Wilcoxon Rank-Sum Test (Mann-Whitney U)

**Purpose:** Compare alpha diversity between two groups

**When to use:**
- Exactly 2 groups to compare
- Non-normally distributed data (common for diversity metrics)
- Independent samples

**Assumptions:**
✓ Independent observations  
✓ Ordinal or continuous data  
✓ Similar distribution shapes (for location interpretation)  

**Null hypothesis:** The two groups have the same distribution

**Implementation:**
```r
wilcox.test(Shannon ~ Group, data = alpha_div)
```

**Interpretation:**
- p < 0.05: Evidence that groups differ in alpha diversity
- Report median and IQR for each group
- Consider effect size (e.g., rank-biserial correlation)

**Citations:**
- Mann HB, Whitney DR (1947) Ann Math Stat 18:50-60
- Hollander M, Wolfe DA (1999) Nonparametric Statistical Methods, 2nd ed

### Kruskal-Wallis Test

**Purpose:** Compare alpha diversity among ≥3 groups

**When to use:**
- Three or more groups to compare
- Non-normally distributed data
- Independent samples

**Assumptions:**
✓ Independent observations  
✓ Ordinal or continuous data  
✓ Similar distribution shapes across groups  

**Null hypothesis:** All groups have the same distribution

**Implementation:**
```r
kruskal.test(Shannon ~ Group, data = alpha_div)
```

**Post-hoc testing:**
If significant, perform pairwise Wilcoxon tests with multiple testing correction:
```r
pairwise.wilcox.test(alpha_div$Shannon, alpha_div$Group, 
                    p.adjust.method = "fdr")
```

**Interpretation:**
- p < 0.05: At least one group differs from others
- Follow up with pairwise comparisons
- Report median ± IQR for each group

**Citations:**
- Kruskal WH, Wallis WA (1952) J Am Stat Assoc 47:583-621

---

## Beta Diversity Statistics

### PERMANOVA (Permutational Multivariate Analysis of Variance)

**Purpose:** Test if community composition differs among groups

**When to use:**
- Comparing microbial community structure
- Any distance metric (Bray-Curtis, UniFrac, etc.)
- ≥2 groups

**Assumptions:**
✓ Independent observations  
✓ **Homogeneity of multivariate dispersions** (critical!)  
✓ Exchangeability under null hypothesis  
✓ Adequate sample size (n ≥ 3-5 per group)  

**Null hypothesis:** Centroids of groups are equivalent for all groups

**Implementation:**
```r
library(vegan)
dist_mat <- distance(ps, method = "bray")
adonis2(dist_mat ~ Group, data = metadata, permutations = 999)
```

**Critical: Check Assumptions!**

PERMANOVA is sensitive to dispersion differences. **Always run betadisper first:**

```r
# Check homogeneity of dispersions
bd <- betadisper(dist_mat, metadata$Group)
permutest(bd, permutations = 999)
```

**Interpretation of betadisper:**
- p > 0.05: Assumption met → PERMANOVA results are reliable
- p ≤ 0.05: **Assumption violated** → PERMANOVA may detect dispersion, not location

**What to do if assumption is violated:**
1. Report both location (PERMANOVA) and dispersion (betadisper) effects
2. Consider alternative tests (ANOSIM, MRPP)
3. Try data transformation
4. Interpret PERMANOVA cautiously

**PERMANOVA Interpretation:**
- **R²**: Proportion of variance explained (effect size)
  - < 0.01: Very small effect
  - 0.01-0.06: Small effect
  - 0.06-0.14: Medium effect
  - > 0.14: Large effect
- **p-value**: Significance of group differences
- Report: F-statistic, R², p-value, and betadisper result

**Citations:**
- Anderson MJ (2001) Austral Ecol 26:32-46
- Anderson MJ, Walsh DCI (2013) Ecol Monogr 83:557-574

### Pairwise PERMANOVA

**Purpose:** Identify which specific group pairs differ

**When to use:**
- After significant overall PERMANOVA
- ≥3 groups

**Multiple testing correction required!**

**Implementation:**
```r
# This pipeline's implementation includes FDR correction
pairwise_results <- pairwise_permanova(ps, "Group", 
                                      distance = "bray",
                                      correction = "fdr")
```

**Interpretation:**
- Use adjusted p-values for conclusions
- Report both raw and adjusted p-values
- Consider effect sizes (R²) alongside p-values

### Alternative Tests

#### ANOSIM (Analysis of Similarities)

**When to use:**
- Alternative to PERMANOVA
- When dispersions are unequal
- Focuses on rank dissimilarities

**Advantages:**
- Less sensitive to dispersion differences
- Easier interpretation (R statistic -1 to +1)

**Disadvantages:**
- Lower power than PERMANOVA
- Doesn't handle complex designs well

**Implementation:**
```r
anosim(dist_mat, metadata$Group, permutations = 999)
```

#### MRPP (Multi-Response Permutation Procedure)

**When to use:**
- Alternative to PERMANOVA
- Tests if groups are more similar within than between

**Interpretation:**
- δ (delta): Within-group average distance
- A: Effect size (chance-corrected within-group agreement)

---

## Multiple Testing Correction

### Why Correction is Necessary

When performing multiple statistical tests, the probability of false positives (Type I errors) increases:

- 1 test at α=0.05: 5% false positive rate
- 20 tests: Expected 1 false positive
- 100 tests: Expected 5 false positives

**Solution:** Adjust p-values to control family-wise error rate (FWER) or false discovery rate (FDR).

### Correction Methods

#### 1. False Discovery Rate (FDR) - Benjamini-Hochberg

**Recommended for most microbiome analyses**

**Controls:** Proportion of false positives among rejected hypotheses

**Advantages:**
- More powerful than Bonferroni
- Appropriate for exploratory research
- Balances discovery and rigor

**When to use:**
- Multiple alpha diversity tests (e.g., Shannon, Simpson, Observed)
- Pairwise PERMANOVA comparisons
- Differential abundance testing
- Most microbiome applications

**Implementation:**
```r
p.adjust(p_values, method = "fdr")
```

**Interpretation:**
- Adjusted p < 0.05: < 5% of discoveries expected to be false positives

**Citation:**
- Benjamini Y, Hochberg Y (1995) J R Stat Soc B 57:289-300

#### 2. Bonferroni Correction

**Most conservative correction**

**Controls:** Family-wise error rate (probability of ≥1 false positive)

**Advantages:**
- Strict control of Type I error
- Simple to understand
- Appropriate for confirmatory research

**Disadvantages:**
- Very conservative (low power)
- May miss real effects

**When to use:**
- Small number of pre-planned comparisons
- Confirmatory studies
- When false positives are very costly

**Implementation:**
```r
p.adjust(p_values, method = "bonferroni")
# Or manually: p_values * number_of_tests
```

**Interpretation:**
- Adjusted p < 0.05: Very strong evidence

#### 3. Holm's Method

**Stepwise correction (compromise)**

**Controls:** Family-wise error rate (like Bonferroni)

**Advantages:**
- More powerful than Bonferroni
- Still controls FWER
- No assumptions about test independence

**When to use:**
- Intermediate between FDR and Bonferroni
- Confirmatory research with multiple tests

**Implementation:**
```r
p.adjust(p_values, method = "holm")
```

#### 4. No Correction

**When appropriate:**
- Single pre-planned comparison
- Pilot/exploratory studies (report as such)
- Already using methods with built-in correction

**⚠️ Warning:** Always report if no correction was applied and justify why.

### Pipeline Implementation

The pipeline provides automatic correction:

```r
# Alpha diversity tests with FDR correction
alpha_results <- test_alpha_diversity(ps, "Group", 
                                     measures = c("Shannon", "Observed", "Simpson"),
                                     correction = "fdr")

# Pairwise PERMANOVA with FDR correction
pairwise_results <- pairwise_permanova(ps, "Group",
                                      correction = "fdr")
```

### Reporting Corrected Results

**Always report:**
1. Number of tests performed
2. Correction method used
3. Both raw and adjusted p-values (in tables/supplements)
4. Significance threshold for adjusted values

**Example:**
> "Alpha diversity was compared among three groups using Kruskal-Wallis tests  
> for Shannon, Simpson, and Observed richness (3 tests). P-values were adjusted  
> for multiple testing using the Benjamini-Hochberg FDR procedure (adjusted  
> p < 0.05 considered significant)."

---

## Assumptions and Diagnostics

### Testing Assumptions

#### 1. Independence

**Assumption:** Samples are independent

**Check:**
- Experimental design review
- No repeated measures without accounting for it
- No pseudoreplication

**Violations:**
- Multiple samples from same individual
- Time series without appropriate modeling
- Spatial autocorrelation

**Solutions:**
- Use mixed-effects models
- Block by individual/time/location
- Account for correlation structure

#### 2. Normality (for parametric tests)

**Check:**
```r
# Visual: Q-Q plot
qqnorm(alpha_div$Shannon)
qqline(alpha_div$Shannon)

# Statistical: Shapiro-Wilk test
shapiro.test(alpha_div$Shannon)
```

**Note:** Most diversity metrics are NOT normally distributed → use non-parametric tests

#### 3. Homogeneity of Variances

**For parametric tests:**
```r
# Levene's test
car::leveneTest(Shannon ~ Group, data = alpha_div)
```

**For PERMANOVA (multivariate dispersions):**
```r
# Betadisper (CRITICAL!)
bd <- betadisper(dist_mat, Group)
permutest(bd)
```

### What to Do When Assumptions Fail

| Assumption Violated | Solution |
|---------------------|----------|
| Normality | Use non-parametric tests (Wilcoxon, Kruskal-Wallis) |
| Equal variances | Use Welch's t-test, or report with caution |
| Homogeneity of dispersions | Report both location and dispersion, use ANOSIM |
| Independence | Use mixed models, block designs |

---

## Effect Sizes and Power

### Why Effect Sizes Matter

P-values tell you **if** an effect exists.  
Effect sizes tell you **how large** it is.

**Always report effect sizes alongside p-values!**

### Alpha Diversity Effect Sizes

#### Rank-Biserial Correlation (2 groups)

**Range:** -1 to +1

**Interpretation:**
- 0.1: Small effect
- 0.3: Medium effect
- 0.5: Large effect

**Implementation:**
```r
library(effectsize)
rank_biserial(Shannon ~ Group, data = alpha_div)
```

#### Epsilon-squared (≥3 groups, Kruskal-Wallis)

**Range:** 0 to 1

**Interpretation:**
- 0.01: Small effect
- 0.06: Medium effect
- 0.14: Large effect

### Beta Diversity Effect Sizes

#### R² from PERMANOVA

**Interpretation:**
- Proportion of variance explained by grouping
- 0.01-0.06: Small
- 0.06-0.14: Medium
- >0.14: Large

**Example:**
> "Group explained 12% of the variance in community composition  
> (PERMANOVA: F=3.45, R²=0.12, p=0.001), indicating a medium effect size."

### Sample Size and Power

**General guidelines:**
- Minimum 3-5 samples per group for basic tests
- 10-15 samples per group for adequate power
- 20-30 samples per group for detecting small effects

**Post-hoc power analysis:**
```r
library(pwr)
# For t-tests/Wilcoxon equivalent
pwr.t.test(n = 15, d = 0.5, sig.level = 0.05, type = "two.sample")
```

**For PERMANOVA:**
- Use simulation or pilot data
- Consider RVAideMemoire::pwr.Adonis2()

---

## Reporting Guidelines

### Minimum Reporting Standards

For each statistical test, report:

1. **Test name** (e.g., "Kruskal-Wallis test")
2. **What was tested** (e.g., "Shannon diversity among three diet groups")
3. **Test statistic** (e.g., "H = 8.23")
4. **P-value** (e.g., "p = 0.016")
5. **Sample sizes** (e.g., "n = 15, 18, 12 per group")
6. **Effect size** (e.g., "ε² = 0.08")
7. **Multiple testing correction** (if applicable)

### Example Reporting

#### Alpha Diversity

> "Alpha diversity (Shannon index) differed significantly among the three treatment  
> groups (Kruskal-Wallis: H=12.45, p=0.002, ε²=0.14, indicating a medium-to-large  
> effect size). Post-hoc pairwise Wilcoxon tests with FDR correction revealed that  
> the high-fat diet group had significantly lower diversity than both the control  
> (W=45, p=0.003, adjusted p=0.009) and low-fat groups (W=38, p=0.001, adjusted  
> p=0.003), while control and low-fat groups did not differ significantly  
> (W=112, p=0.234, adjusted p=0.234)."

#### Beta Diversity

> "Microbial community composition differed significantly among treatment groups  
> (PERMANOVA on Bray-Curtis distances: F=4.32, R²=0.15, p=0.001, 999 permutations).  
> Prior to PERMANOVA, we verified homogeneity of multivariate dispersions  
> (betadisper: F=1.23, p=0.31), confirming that observed differences reflected  
> location (community composition) rather than dispersion. Treatment explained  
> 15% of the variance in community composition, representing a large effect size  
> (Anderson 2001)."

### What NOT to Do

❌ Report p-values without test statistics  
❌ Report "p < 0.05" without exact value  
❌ Ignore failed assumptions  
❌ Perform multiple tests without correction  
❌ Report significance without effect sizes  
❌ Use parametric tests on non-normal data without justification  

### Checklist

- [ ] Test assumptions checked and reported
- [ ] Appropriate test chosen for data type
- [ ] Sample sizes reported
- [ ] Test statistics reported (not just p-values)
- [ ] Effect sizes reported
- [ ] Multiple testing correction applied (if needed)
- [ ] Raw and adjusted p-values provided
- [ ] Methods section includes software/package versions

---

## Software and Versions

Always report software versions:

```r
# Get package versions
packageVersion("phyloseq")
packageVersion("vegan")
packageVersion("ggplot2")

# Session info
sessionInfo()
```

**Example:**
> "Statistical analyses were performed in R version 4.3.0 (R Core Team 2023)  
> using the vegan package (v2.6-4; Oksanen et al. 2022) for PERMANOVA,  
> phyloseq (v1.44.0; McMurdie & Holmes 2013) for diversity calculations,  
> and base R for Kruskal-Wallis and Wilcoxon tests. Multiple testing  
> correction used the Benjamini-Hochberg FDR procedure."

---

## Key References

### General Microbiome Statistics
- Knight R et al. (2018) Best practices for analysing microbiomes. Nat Rev Microbiol 16:410-422

### PERMANOVA
- Anderson MJ (2001) A new method for non-parametric multivariate analysis of variance. Austral Ecol 26:32-46
- Anderson MJ, Walsh DCI (2013) PERMANOVA, ANOSIM, and the Mantel test in the face of heterogeneous dispersions. Ecol Monogr 83:557-574

### Multiple Testing
- Benjamini Y, Hochberg Y (1995) Controlling the false discovery rate. J R Stat Soc B 57:289-300
- Noble WS (2009) How does multiple testing correction work? Nat Biotechnol 27:1135-1137

### Effect Sizes
- Tomczak M, Tomczak E (2014) The need to report effect size estimates revisited. Trends Sport Sci 21:19-25

### Normalization
- McMurdie PJ, Holmes S (2014) Waste not, want not: why rarefying microbiome data is inadmissible. PLoS Comput Biol 10:e1003531
- Weiss S et al. (2017) Normalization and microbial differential abundance strategies depend upon data characteristics. Microbiome 5:27

---

## Pipeline Implementation

This pipeline implements these methods in:
- `scripts/modules/statistical_tests.R` - Multiple testing corrections
- `scripts/modules/beta_diversity_plots.R` - PERMANOVA with assumption checking
- `scripts/modules/alpha_diversity_plots.R` - Alpha diversity tests

Access these functions:
```r
# Alpha tests with correction
source('scripts/modules/statistical_tests.R')
results <- test_alpha_diversity(ps, "Group", correction = "fdr")
print_statistical_results(results, "alpha")

# Pairwise PERMANOVA with correction
pairwise <- pairwise_permanova(ps, "Group", correction = "fdr")
print_statistical_results(pairwise, "pairwise_permanova")

# PERMANOVA with assumption checking
source('scripts/modules/beta_diversity_plots.R')
result <- run_permanova(ps, cfg, "~ Group", check_assumptions = TRUE)
```

---

## Questions?

For method-specific questions:
1. Check assumptions first
2. Consult primary literature (citations above)
3. Consider consulting a statistician for complex designs
4. Report what you did transparently

**Remember:** No method is perfect. The key is choosing appropriate methods, checking assumptions, and reporting transparently.
