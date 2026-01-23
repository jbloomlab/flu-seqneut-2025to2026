# CLAUDE.md - Context for flu-seqneut-2025to2026

## Repository Overview

This is `flu-seqneut-2025to2026`, a project repository for measuring human neutralizing antibody responses to influenza viruses circulating in late 2025 to early 2026. The goal is to provide near real-time data to inform vaccine strain selection for the **February 2026 Northern Hemisphere vaccine formulation**.

This study uses **sequencing-based neutralization assays** to measure titers to influenza viruses (H1N1 and H3N2) against human sera collected in mid to late 2025.

**Important**: This repository uses the [seqneut-pipeline](https://github.com/jbloomlab/seqneut-pipeline) as a git submodule for the main analysis workflow. See `seqneut-pipeline/CLAUDE.md` for detailed documentation of the pipeline itself.

## Project Status

### Currently Implemented
- ✅ **Viral Library Design**: Influenza strains (H1N1 and H3N2) from late 2025 circulation with multiple barcodes per strain for redundancy (see `data/viral_libraries/`)
- ✅ **Viral Library Validation**: Automated validation of library CSVs against README specifications (see `scripts/validate_viral_library.py`)
- ✅ **Library Pooling QC**: Equal volume pool tested on 2026-01-08 to optimize pooling proportions
- ✅ **Pipeline Configuration**: `config.yml` configured with viral library, neutralization standards, and QC thresholds
- ✅ **Cluster Execution**: Slurm submission script ready (`run_Hutch_cluster.bash`)
- ✅ **Nextstrain Trees**: Interactive protein trees via [nextstrain-prot-titers-tree](https://github.com/jbloomlab/nextstrain-prot-titers-tree) submodule (titers will be added when available)

### Pending
- ⏳ Human serum sample collection and processing
- ⏳ Neutralization assay plates with serum dilutions
- ⏳ Main pipeline execution for titer calculations
- ⏳ Aggregated results and visualization

## Critical Scientific Coding Principles

**This is scientific research code.** Data integrity and reproducibility are paramount. Follow these principles:

### 1. Fail Fast - No Silent Errors
- **NEVER allow silent failures** or default to placeholder values
- All data processing should raise explicit exceptions when issues are encountered
- Validate inputs at entry points (file loading, configuration parsing)
- Use assertions for critical assumptions
- Log warnings for unexpected but non-fatal conditions

**Example - Good**:
```python
if barcode_counts.sum() < min_counts:
    raise ValueError(f"Barcode counts {barcode_counts.sum()} below minimum {min_counts}")
```

**Example - Bad**:
```python
if barcode_counts.sum() < min_counts:
    print("Warning: low counts")  # Silent failure - DO NOT DO THIS
    barcode_counts = None  # Might cause issues downstream
```

### 2. Single Source of Truth (DRY Principle)
- **Parameters should be specified in exactly ONE place** (typically `config.yml`)
- Never duplicate parameter values in code, documentation, or multiple config sections
- If a parameter exists in config, reference it - don't redefine it
- This prevents inconsistencies and improves maintainability

**Example - Good**:
```python
min_counts = config["qc_thresholds"]["min_counts"]
```

**Example - Bad**:
```python
min_counts = 500  # Duplicates value from config - DO NOT DO THIS
```

**Documentation Principle**:
- Both code documentation and README should **reference where** values are set (e.g., "see `config.yml`", "configured in `data/viral_libraries/`")
- Do NOT repeat current configuration values in Markdown text
- Describe WHAT parameters control and HOW to set them, not their current values
- This keeps documentation maintainable as configuration changes

**Example - Good**:
```markdown
QC thresholds are configured in `config.yml` under `default_serum_qc_thresholds`.
Key parameters include `min_replicates` and `max_fold_change_from_median`.
```

**Example - Bad**:
```markdown
The minimum replicates threshold is set to 1, and outliers are flagged at 3-fold change.
```

### 3. Explicit Over Implicit
- Be explicit about data transformations and filtering
- Document QC drops in YAML files (already implemented in pipeline)
- Avoid "magic numbers" - use named configuration parameters
- Type hints and docstrings for complex functions

### 4. Reproducibility
- All analysis controlled by `config.yml`
- Random seeds set where stochastic methods used
- Track exact versions (conda environment, submodule versions)
- Results committed to git for key QC files

## Directory Structure and What to Focus On

### Main Pipeline Components (PRIMARY FOCUS)
```
flu-seqneut-2025to2026/
├── Snakefile                    # Top-level workflow orchestration
├── config.yml                   # Main configuration (CRITICAL FILE)
├── run_Hutch_cluster.bash       # Slurm cluster submission
├── scripts/
│   ├── validate_viral_library.py  # Viral library validation script
│   └── nextstrain_prot_titers_tree_alignment_and_metadata.py  # Tree input prep
├── data/
│   ├── viral_libraries/         # Barcode-to-strain mappings
│   ├── neut_standard_sets/      # Neutralization control barcodes
│   ├── miscellaneous_plates/    # QC and pooling plate data
│   └── nextstrain-prot-titers-tree_data/  # Outgroup and site numbering for trees
├── results/                     # Pipeline outputs (partially tracked in git)
├── auspice/                     # Nextstrain auspice JSON files for tree visualization
├── seqneut-pipeline/            # Analysis pipeline (git submodule)
│   └── CLAUDE.md                # Detailed pipeline documentation
└── nextstrain-prot-titers-tree/ # Tree building pipeline (git submodule)
```

### Non-Pipeline Analyses (GENERALLY IGNORE)

The `non-pipeline_analyses/` directory contains one-off analyses for library design and pooling optimization. These are **NOT part of the main neutralization assay pipeline** and are documented separately in that directory. Ignore unless the user specifically asks about them.

## Code Style and Quality Requirements

All code must pass these checks before committing:

### Python
```bash
ruff check .        # Linting (fast, comprehensive)
black .             # Code formatting (auto-fix)
```

### Snakemake
```bash
snakefmt .          # Snakemake formatting
snakemake --lint    # Snakemake validation
```

### Testing
```bash
# Run test example (when modifying pipeline)
cd seqneut-pipeline/test_example
snakemake -j 4 --software-deployment-method conda
```

**Important**: The pipeline already has these configured in `seqneut-pipeline/pyproject.toml`. Run checks from the repository root.

## Key Data Files

### Viral Library
**Location**: `data/viral_libraries/`

The viral library contains influenza strains (H1N1 pdm09 and H3N2) from late 2025 circulation with multiple barcodes per strain for redundancy. Files may include:
- `*-designed.csv`: Designed library from strain selection
- `*-actual.csv`: Actual library after synthesis and QC (if applicable)

Each row represents one barcode-to-strain mapping with HA ectodomain sequences and metadata.

**Key Columns**:
- `strain`: Strain designation (e.g., "A/Wisconsin/67/2022_H1N1")
- `barcode`: Unique barcode sequence (typically 16-nt)
- `subtype`: H1N1 or H3N2
- `subclade`: Phylogenetic subclade
- `num_date`: Decimal collection date
- `nt_sequence_HA_ectodomain`: Circulating HA ectodomain nucleotide sequence
- `protein_sequence_HA_ectodomain`: Translated protein sequence
- Plus: accession, name, strain_type, vaccine_type

### Neutralization Standards
**File**: `data/neut_standard_sets/loes2023_neut_standards.csv`
- Control barcodes used to normalize plate-to-plate variation
- From Loes et al. (2024) publication
- Used in QC thresholds for minimum fraction per well

### Configuration
**File**: `config.yml` (2965 bytes)
- Single source of truth for all pipeline parameters
- See "Configuration Guide" section below for details

## Configuration Guide

The `config.yml` file controls all analysis parameters. Key sections:

### Viral Library Configuration
```yaml
viral_libraries:
  designed:
    barcode_to_strain_csv: data/viral_libraries/flu-seqneut-2025to2026-barcode-to-strain-designed.csv
```

Which library is used may differ across plates; we may eventually have an actual as well as designed library if some strains are dropped during QC.

### QC Thresholds (Plate-Level)
Configured per plate in `plates.{plate_name}.qc_thresholds`. Critical parameters:
- `avg_barcode_counts_per_well`: Minimum average counts (ensure adequate sequencing depth)
- `min_neut_standard_frac_per_well`: Minimum fraction from neut standards (QC for proper normalization)
- `min_neut_standard_count_per_well`: Absolute minimum neut standard counts

### QC Thresholds (Serum-Level)
Configured in `default_serum_qc_thresholds`. Key parameters:
- `min_replicates`: Minimum replicates required to report a titer
- `max_fold_change_from_median`: Threshold for flagging outlier replicates
- `viruses_ignore_qc`: List of viruses exempt from QC (all sera)

Can be overridden per serum or group in `sera_override_defaults`.

### Plate Configuration Structure
- **Regular plates** (`plates:`): Serum neutralization assay plates with full curve fitting
  - Each plate requires: group, date, viral_library, neut_standard_set, samples_csv
  - Plus: qc_thresholds, curvefit_params, manual_drops (optional)
- **Miscellaneous plates** (`miscellaneous_plates:`): QC/pooling plates with barcode counting only
  - Simpler configuration: viral_library, neut_standard_set, samples_csv

See `config.yml` for current plate configuration.

### Nextstrain Tree Configuration
The `nextstrain-prot-titers-tree_config` section in `config.yml` controls tree building for each subtype (H3N2, H1N1). Key parameters:
- `alignment`, `metadata`: Built by pipeline from viral library
- `outgroup`: Reference sequence for tree rooting (in `data/nextstrain-prot-titers-tree_data/`)
- `site_numbering_map`: Maps alignment positions to HA1/HA2 numbering
- `color_by_metadata`: Metadata columns to color tree by (strain_type, subclade, titers when available)
- `titers`: Set to `null` until titer data available; later points to titers TSV
- `auspice_json`: Output path for Nextstrain visualization

The `recent_vaccine_strains` config maps strain names to display labels for vaccine strains to include alongside circulating strains.

## Common Tasks for Claude

### When User Asks About Pipeline Details
**First**: Direct them to `seqneut-pipeline/CLAUDE.md` for comprehensive pipeline documentation. It covers:
- Complete workflow (barcode counting → plate processing → curve fitting → aggregation)
- QC system details
- Output file descriptions
- Troubleshooting guides

**This file** focuses on project-specific context (viral library, current status, strain selection).

### When Adding Neutralization Assay Plates
1. Obtain plate metadata:
   - Plate name (no underscores or pipes in group names!)
   - Date (YYYY-MM-DD format)
   - Sample CSV with columns: `well`, `serum`, `dilution_factor`, `replicate`, `fastq`
2. Add to `config.yml` under `plates:`
3. Configure QC thresholds (start with defaults from seqneut-pipeline docs)
4. Run pipeline: `bash run_Hutch_cluster.bash` or `snakemake -j <cores>`
5. Examine QC drops: `results/qc_drops/plate_qc_drops.yml`
6. Review plate notebook: `docs/plates/{plate_name}/process_plate.html`

### When Troubleshooting Pipeline Issues
1. Check Snakemake error messages (often very informative)
2. Look at `results/qc_drops/*.yml` for data quality issues
3. Examine plate notebooks in `docs/` for visual diagnostics
4. Verify `config.yml` structure (run `snakemake --lint`)
5. Check that FASTQ files exist and are readable

### When Asked About Specific Strains
1. The viral library is in `data/viral_libraries/` (check config.yml for active library file)
2. Each barcode maps to one strain (but strains can have multiple barcodes)
3. Metadata includes: subtype (H1N1/H3N2), subclade, collection date, vaccine type

### When Modifying Non-Pipeline Analyses
**Warning**: Ask user to confirm - these analyses are typically one-off and complete. They are documented separately in the `non-pipeline_analyses/` directory.

## Relationship to Similar Projects

This repository follows the same structure as [flu-seqneut-2025](https://github.com/jbloomlab/flu-seqneut-2025), which:
- Informed September 2025 vaccine selection
- Used seqneut-pipeline v5.0.0+ (this project uses v6.1.0)
- Processed human sera from 5 cohorts
- Published as Kikawa et al. (2025) in Virus Evolution

**Key Differences**:
- Different viral library (2025-2026 circulation vs 2024-2025)
- Different timeframe (late 2025 - early 2026 vs late 2024 - spring 2025)
- Different target (February 2026 vs September 2025 vaccine decision)

## Key Publications and References

1. **Pipeline Methods**: [Loes et al. (2024)](https://doi.org/10.1128/jvi.00689-24) - Journal of Virology
   - Describes sequencing-based neutralization assay methodology
2. **Previous Application**: [Kikawa et al. (2025)](https://academic.oup.com/ve/article/11/1/veaf086/8313343) - Virus Evolution
   - flu-seqneut-2025 project for September 2025 vaccine selection
3. **Pipeline Repository**: [seqneut-pipeline on GitHub](https://github.com/jbloomlab/seqneut-pipeline)
   - Modular Snakemake pipeline used as submodule

## Getting Help

- **Pipeline issues**: See `seqneut-pipeline/CLAUDE.md` and pipeline documentation
- **Configuration questions**: Check `config.yml` and seqneut-pipeline docs
- **Non-pipeline analyses**: See `non-pipeline_analyses/` (documented separately in that directory)
- **Git submodule issues**: `git submodule update --init --recursive`

## Quick Reference

**Run pipeline locally**:
```bash
conda activate seqneut-pipeline
snakemake -j 4 --software-deployment-method conda
```

**Run on Slurm cluster**:
```bash
bash run_Hutch_cluster.bash
```

**Check code quality**:
```bash
ruff check .
black .
snakefmt .
snakemake --lint
```

**Update submodule**:
```bash
git submodule update --remote seqneut-pipeline
```
