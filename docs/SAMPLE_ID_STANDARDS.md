# Sample ID Naming Standards

## Overview
This document defines the strict naming conventions for sample identifiers in the 16S metagenomics pipeline. Consistent naming prevents regex matching errors, ensures proper file pairing, and maintains data integrity throughout the analysis.

## Standard Format

### Basic Pattern
```
{PROJECT}[-]{SAMPLE_NUM}[-]{REPLICATE}
```

**Components:**
- `PROJECT`: Project/cohort identifier (alphanumeric, no spaces)
- `SAMPLE_NUM`: Sample number (zero-padded to 3 digits recommended)
- `REPLICATE`: Optional replicate identifier (e.g., A, B, R1)

### Valid Examples
```
KMUN-001
KMUS-023
STUDY-042-A
EXP1-001
CTRL-001-R1
```

### Invalid Examples
```
Sample 1          # Contains spaces
KMUN_001          # Uses underscore instead of hyphen
KM-1              # Not zero-padded
sample-001        # Lowercase (not recommended)
001               # Missing project identifier
```

## File Naming Convention

### Paired-End FASTQ Files

**Format:**
```
{SAMPLE_ID}_{READ}.{extension}
```

**Components:**
- `SAMPLE_ID`: Sample identifier following the standard format above
- `READ`: Read direction indicator (see patterns below)
- `extension`: File extension (see patterns below)

**Supported READ Patterns:**
- `R1` / `R2` (Illumina standard)
- `1` / `2` (simplified)
- `_1` / `_2` (alternative)

**Supported Extensions:**
- `.fastq.gz` (preferred)
- `.fq.gz` (alternative)
- `.fastq` (uncompressed, not recommended)

**Valid Filename Examples:**
```
KMUN-001_R1.fastq.gz
KMUN-001_R2.fastq.gz
KMUS-023_1.fq.gz
KMUS-023_2.fq.gz
STUDY-042-A_R1.fastq.gz
STUDY-042-A_R2.fastq.gz
```

**Invalid Filename Examples:**
```
KMUN-001.R1.fastq.gz     # Uses dot instead of underscore
KMUN-001-R1.fastq.gz     # Uses hyphen instead of underscore
KMUN_001_R1.fastq.gz     # Inconsistent delimiter (underscore in ID)
Sample 1_R1.fastq.gz     # Contains spaces in ID
KMUN-001_forward.fastq.gz # Non-standard read indicator
```

## Metadata Requirements

### Sample ID Column
The metadata CSV file must contain a column matching the sample IDs exactly.

**Recommended column name:** `sample_id`

**Requirements:**
- Must match the sample IDs in FASTQ filenames (excluding read suffix and extension)
- No duplicate IDs
- No missing values
- Case-sensitive matching

**Example metadata.csv:**
```csv
sample_id,group,condition,time_point
KMUN-001,control,healthy,0
KMUN-002,control,healthy,0
KMUS-001,treatment,disease,0
KMUS-002,treatment,disease,14
```

## Validation Rules

### Automated Checks
The pipeline performs the following validation checks:

1. **Format Validation:**
   - Sample IDs match the pattern: `^[A-Z0-9]+-[0-9]{3,}(-[A-Z0-9]+)?$`
   - No spaces, special characters (except hyphens)
   - Consistent delimiters

2. **Pairing Validation:**
   - Each R1 file has a corresponding R2 file
   - Identical sample IDs (excluding read suffix)
   - Matching file sizes (within reasonable tolerance)

3. **Metadata Alignment:**
   - All FASTQ sample IDs present in metadata
   - All metadata sample IDs have corresponding FASTQ files
   - No orphaned samples

4. **Uniqueness:**
   - No duplicate sample IDs
   - Each sample has unique filename

## Migration Guide

### For Existing Datasets
If your existing dataset doesn't follow these conventions:

1. **Rename files systematically:**
   ```bash
   # Example: Convert underscores to hyphens in sample IDs
   for f in *_R1.fastq.gz; do
     new=$(echo "$f" | sed 's/_/-/g')
     mv "$f" "$new"
   done
   ```

2. **Update metadata:**
   - Ensure sample_id column matches new filenames
   - Update any documentation referencing old names

3. **Validate:**
   ```bash
   ./scripts/preflight_check.sh --input-dir <your_data> --output-dir output
   ```

## Troubleshooting

### Common Issues

**Issue:** "Unequal R1 and R2 files"
- **Cause:** Missing paired files or naming mismatch
- **Solution:** Check all R1 files have R2 counterparts with identical base names

**Issue:** "Sample not found in metadata"
- **Cause:** Sample ID in filename doesn't match metadata
- **Solution:** Ensure exact match between FASTQ filenames and metadata sample_id column

**Issue:** "Invalid sample ID format"
- **Cause:** Sample ID contains invalid characters or wrong pattern
- **Solution:** Rename files following the standard format above

## Best Practices

1. **Zero-pad numbers:** Use `001` instead of `1` for better sorting
2. **Use uppercase:** `KMUN-001` instead of `kmun-001`
3. **Be consistent:** Pick one pattern and stick with it
4. **Document exceptions:** If you must deviate, document why
5. **Validate early:** Run preflight checks before starting long pipelines
6. **Version control:** Keep metadata in version control alongside configs

## References

- NCBI SRA naming guidelines: https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/
- Illumina naming conventions: https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm
