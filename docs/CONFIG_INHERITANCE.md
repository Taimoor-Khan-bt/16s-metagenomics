# Configuration Inheritance Guide

## Overview
Configuration inheritance allows you to create project-specific configs that inherit default settings from a base configuration. This reduces duplication, ensures consistency, and makes configs easier to maintain.

## Basic Concept

Instead of repeating all settings in every project config, you can:
1. Define common defaults in `base_config.yaml`
2. Create project configs that only specify what's different
3. The pipeline automatically merges them

## Usage

### Simple Example

**Base config** (`config/base_config.yaml`):
```yaml
project:
  threads: 4
  random_seed: 1234

quality:
  maxEE_fwd: 2
  maxEE_rev: 2
```

**Project config** (`config/my_project.yaml`):
```yaml
inherit_from: "base_config.yaml"

project:
  name: "MyProject"
  threads: 8  # Override threads, keep random_seed

quality:
  maxEE_fwd: 3  # Override forward only
```

**Result after merging**:
```yaml
project:
  name: "MyProject"
  threads: 8          # From my_project.yaml
  random_seed: 1234   # From base_config.yaml

quality:
  maxEE_fwd: 3        # From my_project.yaml
  maxEE_rev: 2        # From base_config.yaml
```

## How It Works

### Inheritance Directive
Add this to your project config:
```yaml
inherit_from: "base_config.yaml"
```

- **Relative paths** are resolved from the config directory
- **Absolute paths** work too: `inherit_from: "/path/to/base.yaml"`

### Merging Rules

1. **Deep merge**: Nested structures are merged recursively
2. **Override wins**: Project values override base values
3. **Additive**: Missing fields are filled from base

**Example:**
```yaml
# Base
analysis:
  mode: "standard"
  compare_by: "group"

# Override
analysis:
  mode: "comprehensive"

# Result
analysis:
  mode: "comprehensive"  # Overridden
  compare_by: "group"    # Inherited
```

## Common Patterns

### Pattern 1: Minimal Project Config

Only specify what's essential:

```yaml
inherit_from: "base_config.yaml"

project:
  name: "QuickAnalysis"

io:
  input_dir: "data/samples"
  metadata_csv: "metadata.csv"
```

Everything else comes from base_config.yaml.

### Pattern 2: Override Analysis Parameters

```yaml
inherit_from: "base_config.yaml"

project:
  name: "HighQualityRun"

quality:
  truncLen_fwd: 250
  truncLen_rev: 230
  maxEE_fwd: 1  # Stricter quality

diversity:
  rarefaction:
    depth: 10000  # Specific depth
```

### Pattern 3: Multi-Environment Setup

**Production config:**
```yaml
inherit_from: "base_config.yaml"

project:
  name: "Production-Analysis"
  threads: 16
  output_dir: "/data/production/output"

io:
  input_dir: "/data/production/fastq"
```

**Test config:**
```yaml
inherit_from: "base_config.yaml"

project:
  name: "Test-Analysis"
  threads: 2
  output_dir: "test_output"

io:
  input_dir: "test_data"

filtering:
  min_reads_per_sample: 100  # Relaxed for testing
```

## Validation

The pipeline validates configs **after** inheritance:

```bash
./run_pipeline.sh config/my_project.yaml
```

Output:
```
[validate] Loading config with inheritance...
[config] Inheriting from: config/base_config.yaml
[validate] Validating configuration...
✓ All required fields present
✓ Input directory exists
✓ Metadata file found
```

## Troubleshooting

### Issue: "Base config not found"

**Cause:** Relative path doesn't resolve correctly

**Solution:** Use absolute path or verify file exists:
```yaml
# Relative from config directory
inherit_from: "base_config.yaml"

# Or absolute
inherit_from: "/full/path/to/base_config.yaml"
```

### Issue: Values not overriding

**Cause:** Incorrect YAML structure (indentation)

**Check:**
```yaml
# WRONG - creates new top-level key
project:
quality:
  maxEE_fwd: 3

# CORRECT - overrides nested value
quality:
  maxEE_fwd: 3
```

### Issue: Unexpected merge behavior

**Debug:** Print merged config:
```r
# In R console
source("scripts/config_inheritance.R")
cfg <- load_config_with_inheritance("config/my_project.yaml")
str(cfg)
```

## Best Practices

### 1. Keep Base Generic
Base config should contain:
- Sensible defaults for most projects
- Common tool parameters
- Standard workflows

### 2. Override Sparingly
Only override what's truly different:
```yaml
# BAD - repeating unchanged values
quality:
  truncLen_fwd: 240
  truncLen_rev: 200
  maxN: 0
  maxEE_fwd: 2
  maxEE_rev: 2

# GOOD - only what changed
quality:
  truncLen_fwd: 250
```

### 3. Document Overrides
Add comments explaining why you override:
```yaml
quality:
  maxEE_fwd: 3  # Higher error threshold due to poor library quality
```

### 4. Version Control
- Keep `base_config.yaml` in version control
- Track changes with clear commit messages
- Tag stable versions: `v1.0`, `v2.0`

### 5. Testing
Test new base configs with multiple projects:
```bash
# Test with different project configs
./run_pipeline.sh config/project_a.yaml
./run_pipeline.sh config/project_b.yaml
```

## Migration Guide

### Converting Existing Configs

**Before** - Standalone config:
```yaml
project:
  name: "OldProject"
  threads: 4
  random_seed: 1234

io:
  input_dir: "data"
  metadata_csv: "metadata.csv"

amplicon:
  region: "V3-V4"
  fwd_primer: "CCTACGGGNGGCWGCAG"
  # ... many more lines
```

**After** - Inherited config:
```yaml
inherit_from: "base_config.yaml"

project:
  name: "OldProject"

io:
  input_dir: "data"
  metadata_csv: "metadata.csv"

# That's it! Everything else from base
```

### Migration Steps

1. **Create base_config.yaml** with common settings
2. **Identify project-specific values** in existing config
3. **Create minimal override** with just those values
4. **Test** to ensure identical behavior
5. **Document** any intentional changes

## Advanced Usage

### Chained Inheritance

Base → Intermediate → Project:

**base_config.yaml:**
```yaml
project:
  threads: 4
```

**illumina_defaults.yaml:**
```yaml
inherit_from: "base_config.yaml"

amplicon:
  region: "V3-V4"
```

**my_project.yaml:**
```yaml
inherit_from: "illumina_defaults.yaml"

project:
  name: "MyIlluminaProject"
```

### Conditional Inheritance

Different bases for different platforms:

```yaml
# For Illumina data
inherit_from: "base_illumina.yaml"

# For PacBio data
# inherit_from: "base_pacbio.yaml"
```

## Examples

### Example 1: Multi-Site Study
```yaml
inherit_from: "base_config.yaml"

project:
  name: "MultiSite-Gut-2025"

io:
  input_dir: "data/all_sites"
  metadata_csv: "metadata_multisite.csv"

analysis:
  compare_by: "site"
  adjust_by: ["age", "sex", "bmi"]
```

### Example 2: Quick QC Run
```yaml
inherit_from: "base_config.yaml"

project:
  name: "QC-Check"

diversity:
  rarefaction:
    enabled: false  # Skip for QC

visualization:
  figure_format: "pdf"  # Quick review
```

## References

- Base configuration: `config/base_config.yaml`
- Inheritance implementation: `scripts/config_inheritance.R`
- Validation: `scripts/validate_config.R`
