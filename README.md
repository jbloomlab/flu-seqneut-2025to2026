# Near real-time data on the human neutralizing antibody landscape to influenza virus in late 2025 to early 2026 to inform vaccine-strain selection in February 2026

This study uses sequencing-based neutralization assays to measure titers to influenza viruses with HAs from seasonal H3N2 and H1N1 viruses representative of those circulating in late 2025 and early 2026 against human sera collected in mid to late 2025.

The goal is to provide near real-time neutralization data to inform the **February 2026 Northern Hemisphere vaccine-strain selection** decision.

The data here are described in *Add link to preprint*.

## Quick links to key results
Here are quick links with key data/results:

- [results/final_titer_data](results/final_titer_data) has the final titer data and information on the viruses and sera for which these final data were obtained. Specifically:
  + For human sera:
    * [results/final_titer_data/human_titers.csv](results/final_titer_data/human_titers.csv): final set of titers for each virus/serum pair (keeping only viruses for which titers measured against most sera, and sera for which titers measured against most viruses).
    * [results/final_titer_data/human_sera.csv](results/final_titer_data/human_sera.csv): detailed information about the sera for which these titers were measured.
    * [results/final_titer_data/human_sera_multicohort.csv](results/final_titer_data/human_sera_multicohort.csv): sera assigned to multiple cohorts (original, "All", and days-post-vax cohorts); differs from [results/final_titer_data/human_sera.csv](results/final_titer_data/human_sera.csv) in that each serum may appear in multiple rows.
    * [results/final_titer_data/human_viruses.csv](results/final_titer_data/human_viruses.csv): detailed information about the viruses for which these titers were measured.
    * [results/final_titer_data/human_titers_summarized_by_virus.csv](results/final_titer_data/human_titers_summarized_by_virus.csv): summary statistics about the titers against each virus.

- Interactive Nextstrain protein trees (can be colored by subclade, median titer, or fraction below titer cutoff for each cohort):
  - [H3N2](https://nextstrain.org/community/jbloomlab/flu-seqneut-2025to2026@main/H3N2)
  - [H1N1](https://nextstrain.org/community/jbloomlab/flu-seqneut-2025to2026@main/H1N1)

## Study Overview

### Assay Method
Sequencing-based neutralization assays as described in [Loes et al. (2024), *Journal of Virology*](https://doi.org/10.1128/jvi.00689-24) and [Kikawa et al. (2025), *Virus Evolution*](https://academic.oup.com/ve/article/11/1/veaf086/8313343).

### Final QC-ed titer data
If you don't care to understand the overall repo structure and are just looking for final QC-ed titer data and information on the sera / viruses for which those titers were generated, look in [results/final_titer_data/](results/final_titer_data/).

### Viral Library
- Influenza strains (H3N2 and H1N1 pdm09) from late 2025 to early 2026 circulation
- Multiple barcodes per strain for redundancy
- Each viral construct contains the circulating HA ectodomain with signal peptide, endodomain, and cytoplasmid tail constant across strains for a given subtype
- Library defined in in CSVs in [data/viral_libraries/](data/viral_libraries) (see [config.yml](config.yml) for active library file for each experiment)
- Each viral library CSV must have the following required columns (additional columns are allowed):
  - *strain*: strain name; there can be multiple rows for each strain (as we may have several barcodes for a strain), but each strain name must be uniquely paired with a single value for all other required columns except *barcode* (where a strain can be paired with multiple *barcode* values). Note that non-required columns (like *Twist_name*) may vary per barcode. Strain name must also end in "_H3N2" or "_H1N1"
  - *subtype*: H3N2 or H1N1
  - *strain_type*: should be "vaccine" or the *circulating_strain_type* defined in `config.yml` (eg "circulating_2025to2026")
  - *vaccine_type*: if *strain_type* is "vaccine", this should be "egg" or "cell"; otherwise it should be null
  - *barcode*: should be 16-nt string, all barcode entries should be unique for each row
  - *accession*: Genbank accession if available
  - *subclade*: can only be null for *strain_type* of "vaccine"
  - *derived_haplotype*: description of derived haplotype from clade (eg, clade plus additional HA1 mutations), can only be null for *strain_type* of "vaccine". Note that HA2 mutations are not included in the derived_haplotype naming, so multiple strains may share the same derived_haplotype if they differ only in HA2.
  - *collection_date*: collection date as float (eg, 2025.5); in many cases this may refer to the date that HA1 haplotype was last identified rather than the actual collection date of the particular named strain.
  - *nt_sequence_HA_ectodomain*: nucleotide sequence of HA ectodomain, must be all A, T, C, or G (either case) and length multiple of three. Note that this is not the same as the full Twist synthesized insert, rather it is the part of the HA that is taken from that strain, and does not include flanking constant regions not derived from the strain (recall signal peptide and transmembrane / cytoplasmic tail are constant in our barcoded constructs), so it is not the full-length HA.
  - *protein_sequence_HA_ectodomain*: protein sequence, must meet these criteria:
      + must not contain any stop codons ("*" characters) and must be the result of translating *nt_sequence_HA_ectodomain*.
      + must start with "CIGY" if a H1N1 *subtype*, and with "Q[KNR][IL]P" if a H3N2 *subtype* (note this is start of HA ectodomain for H3N2, for H1N1 a "DTL" must be added to get HA ectodomain since the plasmid construct used in experiments provides the first 3 ectodomain amino acids from WSN strain).
      + must end with "NNRFQ" if a H3N2 *subtype*, and "[EK]IDG[VI]" if a *H1N1* subtype.
      + the length should be 500 for H1N1 *subtype* and 501 for *H3N2* subtype (note this might change of H1N1 is modified to include pre-2009 seasonal H1N1 as well as pdmH1N1 lineage).
      + Each unique protein sequence should only be associated with one *strain*.

There are two CSVs in [data/viral_libraries](data/viral_libraries), one giving the originally designed library and the other giving the library that had strains that passed various QC and are used for the actual titer measurements.
The strain used for the actual measurements are in [data/viral_libraries/flu-seqneut-2025to2026-barcode-to-strain-actual.csv](data/viral_libraries/flu-seqneut-2025to2026-barcode-to-strain-actual.csv).

Details about the viruses in these files that are actually used to generate final titer data are in [results/final_titer_data/](results/final_titer_data/).

### Sera
- Human sera collected in mid to late 2025 from multiple cohorts
- Per-cohort metadata CSVs stored in [data/sera_metadata/](data/sera_metadata); note that sera may be listed here for which no titers were measured
- Each sera metadata CSV must have the following required columns (additional columns are allowed):
  - *bloom_lab_id*: unique identifier for each serum sample; must be non-null and unique across all files
  - *cohort*: cohort identifier (e.g., "HKU", "NIID", "PENN"); must be non-null
  - *species*: species of serum donor (e.g., "human"); must be non-null
  - *age*: age of donor; can be numeric (e.g., "45"), a range with or without 'y' suffix (e.g., "10-19y", "18-29"), or open-ended (e.g., "75+")
  - *sex*: sex of donor; accepts various formats ("M"/"F", "male"/"female", "Male"/"Female") which are normalized during aggregation
  - *collection_date*: date of serum collection; accepts "Mon-YYYY" (e.g., "Aug-2025") or "Mon-YY" (e.g., "Nov-25") formats

The per-cohort metadata files are aggregated and validated by `scripts/aggregate_sera_metadata.py` (Snakemake rule `aggregate_sera_metadata`) into [results/sera_metadata/all_sera_metadata.csv](results/sera_metadata/all_sera_metadata.csv).
The aggregation performs these standardizations:
  - Renames *bloom_lab_id* to *serum*
  - Normalizes *sex* to "Male", "Female", or "Unknown"
  - Parses *age* to create *age_numeric* column (midpoint in years for ranges, lower bound for open-ended)
  - Standardizes *collection_date* to "YYYY-MM" format
  - Validates uniqueness of serum IDs across all cohorts

Details about the sera in these files for which we actually obtained high-quality titer data are in [results/final_titer_data/](results/final_titer_data/).

### Aggregated titers
The titers aggregated across plates are in [results/aggregated_titers/](results/aggregated_titers/).

The titers after classifying by species from which the sera came and subsetting for just titers for viruses and sera which are well measured are in [results/final_titer_data/](results/final_titer_data/).

## Repository Structure

### Main Pipeline Components

```
flu-seqneut-2025to2026/
├── Snakefile                    # Top-level Snakemake workflow
├── config.yml                   # Main pipeline configuration
├── run_Hutch_cluster.bash       # SLURM cluster submission script
│
├── data/                        # Input data
│   ├── viral_libraries/         # Barcode-to-strain mappings
│   ├── neut_standard_sets/      # Neutralization control barcodes (Loes et al. 2024)
│   ├── plates/                  # Per-plate sample metadata CSVs
│   ├── miscellaneous_plates/    # QC and pooling validation data
│   ├── sera_metadata/           # Per-cohort serum metadata CSVs
│   └── nextstrain-prot-titers-tree_data/  # Outgroup and site numbering for trees
│
├── results/                     # Pipeline outputs (key files tracked in git)
│   ├── barcode_counts/          # Per-sample barcode counts
│   ├── barcode_fates/           # Read processing statistics
│   ├── plates/                  # Per-plate analyses
│   ├── sera/                    # Per-serum titer aggregations
│   ├── aggregated_titers/       # Final titer matrices
│   ├── qc_drops/                # QC filtering documentation
│   └── miscellaneous_plates/    # Library pooling QC results
│
├── auspice/                     # Nextstrain auspice JSON files for tree visualization
│
├── seqneut-pipeline/            # Analysis pipeline git submodule
└── nextstrain-prot-titers-tree/ # Tree building git submodule
```

Note all input data for the pipeline are in [./data/](data).
All output generated by pipeline goes in [./results/](results), although only some results are tracked in the GitHub repo (see [.gitignore](.gitignore)).

### Snakemake pipeline
The pipeline is specified in [Snakefile](Snakefile), which is a snakemake file that includes analyses run by the submodules as well as additional custom rules.

### Submodules
This repository uses the following git submodules:

 - [seqneut-pipeline](https://github.com/jbloomlab/seqneut-pipeline) for the core analysis workflow.
 - [nextstrain-prot-titers-tree](https://github.com/jbloomlab/nextstrain-prot-titers-tree) for building interactive Nextstrain protein trees with titer data.

### Non-Pipeline Analyses

The [`non-pipeline_analyses/`](non-pipeline_analyses/) subdirectory contains auxiliary workflows (library design and pooling optimization) that are separate from the main neutralization assay pipeline. See the [README in non-pipeline_analyses](non-pipeline_analyses/README.md) for details.

## Configuration

All pipeline parameters are specified in [`config.yml`](config.yml). Key configuration sections include:

- **Viral Libraries**: Barcode-to-strain CSV file paths
- **Neutralization Standards**: Control barcode definitions (Loes et al. 2024)
- **Illumina Barcode Parser**: Parameters for extracting barcodes from FASTQ files
- **Plates**: Per-plate configuration including sample metadata, QC thresholds, and curve fitting parameters
- **Miscellaneous Plates**: Plates for QC/pooling (barcode counting only, no curve fitting)
- **Serum QC**: Thresholds for replicate consistency and minimum replicates
- **Titer Calculation**: Method (midpoint or NT50) and display options

See [`seqneut-pipeline/README.md`](seqneut-pipeline/README.md) for detailed configuration documentation.

## Running the Pipeline

### Environment Setup

```bash
# Create conda environment
conda env create -f seqneut-pipeline/environment.yml

# Activate environment
conda activate seqneut-pipeline
```

### Local Execution

```bash
snakemake -j <num_cores> --software-deployment-method conda
```

### Cluster Execution (Fred Hutch)

```bash
bash run_Hutch_cluster.bash
```

## Documentation on GitHub Pages
As described in [seqneut-pipeline/README.md](seqneut-pipeline/README.md), the pipeline generates HTML documentation that can be pushed to a *gh-pages* branch for publishing on GitHub Pages with:
```bash
./seqneut-pipeline/publish_docs_gh-pages.sh
```

This is then displayed on GitHub Pages via the *gh-pages* branch at [https://jbloomlab.github.io/flu-seqneut-2025to2026](https://jbloomlab.github.io/flu-seqneut-2025to2026). 

The command to publish to the *gh-pages* branch must be run manually, and GitHub must be manually configured to display the Pages from this branch.

## Detailed Pipeline Steps

The complete analysis workflow proceeds through these stages:

### 1. Viral Library Design and Production

The viral library (barcode-to-strain mapping CSV in `data/viral_libraries/`) is designed and produced prior to running the main pipeline. See `non-pipeline_analyses/` for library design and pooling optimization workflows (documented separately in that directory).

The pipeline validates viral library CSVs against the specifications above using `scripts/validate_viral_library.py`. Validation output is written to `results/validate_viral_library/`.

### 2. Neutralization Assay Setup

**Goal**: Incubate viral library with serially diluted human sera to measure neutralization.

**Process** (when serum samples available):
- Sera are heat inactivated and serially diluted
- Diluted sera incubated with viral library
- Non-neutralized viruses infect target cells and express barcodes
- Each plate includes serum dilution series, replicates, and no-serum controls

**Configuration**: Plates configured in `config.yml` with sample metadata CSV specifying wells, sera, dilutions, replicates, and FASTQ file paths.

### 3. High-Throughput Sequencing

**Goal**: Generate FASTQ files with barcode counts for each well.

**Process**:
- RNA extraction from infected cells
- PCR amplification of barcode region with well-specific indexes
- Illumina sequencing
- Demultiplexing by well index

**Output**: FASTQ files listed in plate sample CSVs.

### 4. Barcode Counting

**Pipeline Rule**: `count_barcodes`

**Process**:
- Parses FASTQ files using `dms_variants.illuminabarcodeparser.IlluminaBarcodeParser`
- Extracts barcodes and validates against viral library
- Classifies read fates (valid, invalid, unparseable, low quality)
- Aggregates counts per barcode per well

**Key Outputs**:
- `results/barcode_counts/{sample}.csv` - Barcode counts
- `results/barcode_fates/{sample}.csv` - Read fate statistics
- `results/barcode_invalid/{sample}.csv` - Counts of invalid barcodes; not tracked in GitHub repo but can be useful for debugging

### 5. Plate Processing and QC

**Pipeline Rule**: `process_plate`

**Process**:
- Applies multi-stage QC filters (thresholds specified in `config.yml`)
- Computes fraction infectivity by normalizing to neutralization standards and no-serum wells
- Applies manual drops specified in configuration
- Generates interactive marimo notebook with QC visualizations

**Key Outputs**:
- `results/plates/{plate}/frac_infectivity.csv` - Normalized infectivity matrix
- `results/plates/{plate}/qc_drops.yml` - QC filtering documentation
- `results/plates/{plate}/process_plate.html` - Interactive QC report

### 6. Neutralization Curve Fitting

**Pipeline Rule**: Integrated into `process_plate`

**Process**:
- Groups fraction infectivity by barcode, serum replicate, and dilution
- Fits Hill curves using `neutcurve.CurveFits` (parameters specified in `config.yml`)
- Computes neutralization titers as midpoint or NT50 (configurable)
- Applies goodness-of-fit QC (thresholds in `config.yml`)

**Key Outputs**:
- `results/plates/{plate}/curvefits.csv` - Fitted parameters and titers
- `results/plates/{plate}/curvefits.pickle` - Pickled `neutcurve.CurveFits` objects
- Curve plots in `process_plate.html`

### 7. Serum-Level Titer Aggregation

**Pipeline Rule**: `group_serum_titers`

**Process**:
- Groups replicates of each serum across plates
- Applies serum-level QC (thresholds in `config.yml`)
- Computes median titers across passing replicates
- Generates per-serum neutralization curve visualizations

**Key Outputs**:
- `results/sera/{group}_{serum}/titers.csv` - Median titers per virus
- `results/sera/{group}_{serum}/titers_per_replicate.csv` - Individual replicate titers
- `results/sera/{group}_{serum}/curves.pdf` - Neutralization curves
- `results/sera/{group}_{serum}/qc_drops.yml` - Serum-level QC drops

### 8. Final Titer Aggregation

**Pipeline Rule**: `aggregate_titers`

**Process**:
- Combines all sera titers into matrices (sera × viruses)
- Creates summary statistics and visualizations
- Generates interactive marimo notebook for data exploration

**Key Outputs**:
- `results/aggregated_titers/titers_{group}.csv` - Final titer matrix (tracked in git)
- `results/aggregated_titers/curvefits_{group}.pickle` - All curve fits for custom analysis

### 9. QC Summary and Documentation

**Pipeline Rules**: `aggregate_qc_drops`, `build_docs`

**Process**:
- Aggregates QC filtering across all plates and sera
- Builds HTML documentation with navigation index
- Collects all marimo notebook outputs into `docs/` for GitHub Pages

**Key Outputs**:
- `results/qc_drops/*.yml` - Comprehensive QC summaries (tracked in git)
- `docs/` - Complete HTML documentation

### 10. Nextstrain Protein Trees

**Pipeline Rules**: `nextstrain_prot_titers_tree_alignment_and_metadata`, plus rules from [nextstrain-prot-titers-tree](https://github.com/jbloomlab/nextstrain-prot-titers-tree) submodule

**Goal**: Build interactive phylogenetic trees of library sequences that can be colored by neutralization titers.

**Process**:
- Extracts protein sequences from viral library for circulating strains and recent vaccine strains
- Prepends prefix amino acids if needed (H1N1 sequences need "DTL" added to match full ectodomain)
- Generates alignment FASTA and metadata TSV for each subtype
- When titers configured: processes summarized titers to add median titer and fraction below cutoff columns to metadata, and creates per-serum titers TSV for measurements panel
- Builds phylogenetic tree using IQ-TREE
- Refines tree and extracts amino acid mutations using TreeTime
- Exports to Auspice JSON format for Nextstrain visualization
- Overlays per-serum titer measurements in interactive measurements panel

**Key Inputs** (in `data/nextstrain-prot-titers-tree_data/`):
- `{subtype}_outgroup.fa` - Reference sequence for rooting the tree
- `{subtype}_site_numbering_map.tsv` - Maps alignment positions to HA1/HA2 numbering

**Key Outputs**:
- `results/nextstrain-prot-titers-tree/{subtype}/alignment.fa` - Protein alignment
- `results/nextstrain-prot-titers-tree/{subtype}/metadata.tsv` - Strain metadata with titer summary columns
- `results/nextstrain-prot-titers-tree/{subtype}/titers.tsv` - Per-serum titers for tree overlay
- `auspice/flu-seqneut-2025to2026_{subtype}.json` - Auspice tree JSON (tracked in git)
- `auspice/flu-seqneut-2025to2026_{subtype}_measurements.json` - Auspice titer measurements JSON (tracked in git), created only if titers specified for tree

**Visualization**: Trees are viewable as [Nextstrain Community Builds](https://docs.nextstrain.org/en/latest/guides/share/community-builds.html):
- H3N2: https://nextstrain.org/community/jbloomlab/flu-seqneut-2025to2026@main/H3N2
- H1N1: https://nextstrain.org/community/jbloomlab/flu-seqneut-2025to2026@main/H1N1

**Configuration**: Tree parameters are in `config.yml` under `nextstrain-prot-titers-tree_config`. Key settings:
- `nextstrain-prot-titers-tree_viral_library`: Which viral library key to use for building tree alignment and metadata (e.g., `designed` or `actual`).
- `nextstrain-prot-titers-tree_titers_from`: Species for titer data (e.g., `human`). When set, enables titer overlay on trees with input files derived from `results/final_titer_data/{species}_*.csv`.
- `serum_cohorts_for_tree`: List of cohorts to include in tree coloring. Generates `median_titer_{cohort}_sera` and `frac_w_titer_below_{cutoff}_{cohort}_sera` columns.
- `color_by_metadata`: Columns available for tree coloring (strain_type, subclade, and titer columns).
- `titers`: Configuration for the per-serum titers measurement panel.
- `display_defaults`: Default tree display options (color_by, distance_measure, branch_label).

### 11. Final Titer Data Processing

**Pipeline Rule**: `process_final_titer_data`

**Goal**: Produce QC'd, publication-ready titer datasets with associated metadata.

**Process**:
- Validates that all sera and viruses in titers have corresponding metadata
- Drops explicitly specified sera and viruses (configured in `process_final_titer_data`)
- Filters sera/viruses that don't meet minimum completeness thresholds
- Generates separate CSV files for titers, sera metadata, and virus metadata
- Computes summary statistics by virus and cohort

**Key Outputs** (in `results/final_titer_data/`):
- `{group}_titers.csv` - One row per serum-virus titer measurement
- `{group}_sera.csv` - Metadata for sera included in final dataset (one row per serum)
- `{group}_sera_multicohort.csv` - Sera with multi-cohort assignments; each serum appears in multiple rows (once for "All", once for its original cohort, and once for each days-post-vax cohort if applicable)
- `{group}_viruses.csv` - Metadata for viruses included in final dataset
- `{group}_titers_summarized_by_virus.csv` - Per-virus statistics across sera, computed for each cohort
- `{group}_summary.txt` - Processing log with counts at each step

**Configuration**: Parameters in `config.yml` under `process_final_titer_data`.

## Quality Control Philosophy

The pipeline implements **multi-stage QC**:

1. **Technical QC (Plate-Level)**: Removes wells and barcodes with insufficient sequencing depth or technical issues
2. **Curve Fit QC**: Removes poorly-fitted neutralization curves
3. **Biological QC (Serum-Level)**: Removes replicate outliers and requires minimum replicates

All QC thresholds are configurable in `config.yml`. QC decisions are logged in YAML files and visualized in interactive notebooks.

## References

1. **Assay Methodology**: Loes, A. N. et al. (2024). High-throughput sequencing-based neutralization assay reveals how repeated vaccinations impact titers to recent human H1N1 influenza strains. *Journal of Virology* 98(5):e00689-24. [https://doi.org/10.1128/jvi.00689-24](https://doi.org/10.1128/jvi.00689-24)

2. **Previous Application**:
 - Kikawa, C. et al. (2025). Near real-time data on the human neutralizing antibody landscape to influenza virus to inform vaccine-strain selection in September 2025. *Virus Evolution*. [https://academic.oup.com/ve/article/11/1/veaf086/8313343]([https://academic.oup.com/ve/article/11/1/veaf086/8313343])
 - Kikawa, C. et al. (2025). High-throughput neutralization measurements correlate strongly with evolutionary success of human influenza strains. *eLife*. [https://elifesciences.org/reviewed-preprints/106811](https://elifesciences.org/reviewed-preprints/106811)

3. **Pipeline Repository**: [https://github.com/jbloomlab/seqneut-pipeline](https://github.com/jbloomlab/seqneut-pipeline)

4. **Curve Fitting Package**: [https://jbloomlab.github.io/neutcurve/](https://jbloomlab.github.io/neutcurve/)

## Related Projects

- **flu-seqneut-2025**: Data for September 2025 vaccine selection - [https://github.com/jbloomlab/flu-seqneut-2025](https://github.com/jbloomlab/flu-seqneut-2025)

## Contact

For questions about this study:
- **Caroline Kikawa** (ckikawa@fredhutch.org)
- **Jesse Bloom** (jbloom@fredhutch.org)

For pipeline issues: Open an issue at [seqneut-pipeline GitHub](https://github.com/jbloomlab/seqneut-pipeline/issues)

## License

[MIT License](LICENSE)
