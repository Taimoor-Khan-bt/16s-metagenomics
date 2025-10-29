# Release v2.1.0 - Enhanced Visualizations & Statistical Analysis

**Release Date**: October 30, 2025

## 🎯 Overview

This release focuses on **publication-ready visualizations** with comprehensive bug fixes, statistical enhancements, and improved plot quality. All visualizations now meet journal standards with 600 DPI resolution, proper statistical annotations, and professional styling.

## ✨ Major Features

### Statistical Annotations
- **Alpha Diversity Plots**: Added automated statistical testing with ggpubr
  - Wilcoxon test for 2-group comparisons
  - Kruskal-Wallis test for 3+ groups
  - P-value significance symbols (ns/*/\*\*/\*\*\*/\*\*\*\*)
  - Adaptive comparison display (pairwise for ≤4 groups, overall for >4)
  - Dynamic captions showing test method and significance levels

### Enhanced Plot Quality
- **Resolution**: Increased from 300 to 600 DPI for publication quality
- **Typography**: Base font size set to 14pt for readability
- **Composition Plots**: 
  - Reduced to top 10 taxa (from 15) for clarity
  - Increased legend key size to 0.7cm
  - Added sample counts to subtitles
  - Fixed deprecated `aes_string()` syntax
- **Heatmaps**: Reduced from 50 to 25 taxa with 90° x-axis labels
- **Phylogenetic Trees**: 
  - Simplified to 100 tips (from 500) for clarity
  - Removed cluttered tip labels
  - Fixed deprecated `size` parameter (now uses `linewidth`)
  - Improved subtitles with ASV and phyla counts

## 🐛 Critical Bug Fixes

### Metadata Handling
- **Fixed metadata column mismatch** causing pipeline crashes
  - Changed `metadata.id_column` from "Sample" to "SampleID"
  - Updated validation logic in `validate_and_align_metadata()`
  - Added automatic Sample→SampleID column renaming
  - Improved type coercion and error messages

### Alpha Diversity
- **Fixed variable handling** preventing plot generation
  - Corrected `Observed_log1p` calculation (now properly creates column)
  - Replaced string keys with actual transformed columns
  - Added robust metadata merging with validation
  - Fixed NA handling in diversity calculations

### Deprecated Syntax
- Replaced deprecated `aes_string()` with modern `aes()` + `.data[[]]` syntax
- Updated ggtree `size=` to `linewidth=` for line aesthetics
- Eliminated all controllable deprecation warnings

## 📊 Visualization Improvements

### Before & After
| Aspect | Before | After |
|--------|--------|-------|
| DPI | 300 | 600 |
| Base Font | 12pt | 14pt |
| Top Taxa (Composition) | 15 | 10 |
| Heatmap Taxa | 50 | 25 |
| Tree Tips | 500 | 100 |
| Tree Labels | Cluttered | Clean (removed) |
| Statistical Tests | None | Automated p-values |
| Legend Size | 0.4cm | 0.7cm |

### New Metadata Schema
- Created `docs/metadata_schema.csv` defining required/optional columns
- Documents validation rules and data types
- Supports extensible metadata with 8 optional columns

## 📁 File Changes

### Modified Files
- `config/config.yaml` - Updated plot settings, metadata column, beta metrics
- `scripts/visualization.R` - Major enhancements to all plot functions
- `scripts/analysis_16s.R` - Improved metadata validation and alignment
- `README.md` - Streamlined documentation

### New Files
- `docs/metadata_schema.csv` - Formal metadata specification

## 🔧 Technical Details

### Dependencies
- Added `ggpubr` for statistical annotations (optional, graceful fallback)
- All other dependencies unchanged

### Configuration Changes
```yaml
# Updated settings
metadata:
  id_column: "SampleID"  # Changed from "Sample"

plots:
  dpi: 600              # Increased from 300
  base_size: 14         # New: explicit font size
  top_taxa: 10          # Reduced from 15
  heatmap_top: 25       # New: explicit heatmap limit

analysis:
  beta_metrics:         # Removed unsupported 'aitchison'
    - bray
    - jaccard
    - unifrac
    - wunifrac

amplicon:
  phylogeny:
    max_tips: 100       # Reduced from 500
```

## 📈 Impact

### Resolved Issues
- ✅ Alpha diversity plots not generating
- ✅ Metadata mismatch causing crashes  
- ✅ Cluttered tree visualizations
- ✅ Oversized heatmaps (unreadable labels)
- ✅ Composition plots using deprecated syntax
- ✅ Missing statistical comparisons
- ✅ Inconsistent plot styling

### Quality Metrics
- **All plots generating successfully**: ✓
- **Publication-ready quality**: ✓ (600 DPI)
- **Statistical rigor**: ✓ (p-values on alpha plots)
- **Professional appearance**: ✓ (consistent themes, proper sizing)
- **Pipeline stability**: ✓ (no crashes, robust validation)

## 🚀 Upgrade Guide

### From v2.0.0

1. **Pull latest changes**:
   ```bash
   git pull origin master
   ```

2. **Install ggpubr** (optional but recommended):
   ```r
   install.packages("ggpubr")
   ```

3. **Update your metadata**:
   - Ensure sample ID column is named "SampleID" (not "Sample")
   - Or keep "Sample" - the pipeline now auto-converts

4. **Review config** (optional):
   - Your existing config will work, but consider adopting new defaults
   - See `config/config.yaml` for updated settings

5. **Re-run pipeline**:
   ```bash
   Rscript scripts/runner.R --config config/config.yaml
   ```

## 🔗 Commits

This release includes 2 commits since v2.0.0:

1. `56ae261` - Fix visualization integration and improve plot quality
2. `06e2a89` - Fix deprecated size parameter in ggtree - use linewidth instead

## 🙏 Notes

### External Package Warnings
Some deprecation warnings remain from external packages (ggtree, phyloseq) that use outdated ggplot2 syntax internally. These don't affect functionality and will be resolved when those packages update.

### Breaking Changes
**None** - This release is backward compatible. The pipeline auto-converts "Sample" to "SampleID" if needed.

---

**Full Changelog**: https://github.com/Taimoor-Khan-bt/16s-metagenomics/compare/v2.0.0...v2.1.0
