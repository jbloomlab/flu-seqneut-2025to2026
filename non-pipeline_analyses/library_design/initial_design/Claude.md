# Claude.md - Pipeline Context

## Overview
This is a Snakemake pipeline for selecting influenza HA (hemagglutinin) strains for neutralization assays, specifically for the Feb-2026 vaccine strain selection. The pipeline helps curate ~100 diverse H3N2 and H1N1 strains based on genetic distances, epitope mutations, and metadata.

## Key Concepts

### Haplotypes
- A haplotype is a variant of the HA protein sequence
- Recent haplotypes are those observed since July-1-2025
- Each haplotype has:
  - `derived_haplotype`: A name (usually subclade + mutations)
  - `representative_strain`: An exemplar strain name
  - `representative_strain_ha_sequence`: The full-length HA protein sequence
  - `count`: Number of sequences with this haplotype
  - Metadata: days since latest sequence, median genetic advancement (ga), etc.

### Selection Process
- Manual curation is done by adding values to selection columns (e.g., `John_Huddleston_selection`, `Sam_Turner_selection`)
- The pipeline validates sequences, checks for ambiguous/missing amino acids
- Invalid haplotypes can be manually fixed and added back via `additional_haplotypes`
- Haplotypes can be excluded via `override_select_recent_haplotypes` in config

### Distance Metrics
- **Hamming distance**: Number of amino acid differences across entire HA sequence
- **Epitope distance**: Number of differences at epitope sites only
  - H3N2: Wolf and Koel epitope sites
  - H1N1: Caton epitope sites
- Used to ensure library has diverse coverage while avoiding redundant strains

## Workflow Steps

1. **validate_recent_haplotypes**: Validates HA sequences, checks for ambiguous/missing amino acids, applies selection criteria
2. **create_site_annotations**: Downloads GFF files and epitope definitions, creates site numbering map
3. **curate_library_sequences**: Combines validated + additional haplotypes, computes Hamming and epitope distances
4. **match_genbank**: Matches selected HA proteins to closest GenBank nucleotide sequences
5. **prepare_nextstrain_tree**: Prepares alignment and metadata for visualization
6. **nextstrain-prot-titers-tree module**: Builds protein-based phylogenetic trees (via git submodule)
7. **aggregate_library_strains**: Creates final TSV with all selected strains and GenBank matches

## Key Files

### Inputs
- `config.yaml`: All configuration including input files, selection columns, epitope definitions
- `data/2025-NH-VCM-neutralization-library-strain-selection-{H3N2,H1N1}.tsv`: Recent haplotypes with metadata
- `data/manually_fixed_recent_{H3N2,H1N1}_haplotypes.tsv`: Manually corrected sequences
- `data/{H3N2,H1N1}_outgroup.fa`: Outgroup sequences for tree rooting

### Outputs
- `results/recent_haplotype_validation/`: Validation reports and plots
- `results/curated_library/`: Selected vs non-selected haplotypes with distance metrics
- `results/tree/`: Nextstrain JSON files for interactive visualization (upload to auspice.us)
- `results/aggregated_library_strains/library_strains.tsv`: **Final output** with GenBank matches

### Scripts
- `scripts/validate_recent_haplotypes.py`: Sequence validation and selection logic
- `scripts/curate_library_sequences.py`: Distance computation and library curation
- `scripts/create_site_annotations.py`: Parses GFF and epitope JSONs
- `scripts/prepare_nextstrain_tree.py`: Formats data for tree building
- `scripts/aggregate_library_strains.py`: Combines results across subtypes

## Running the Pipeline

```bash
# Build conda environment
conda env create -f environment.yml
conda activate library_design

# Run locally
snakemake --use-conda -j <ncpus>

# Or on Hutch cluster
sbatch run_Hutch_cluster.bash
```

## Iterative Workflow

1. Run pipeline to see initial selections
2. Review validation reports in `results/recent_haplotype_validation/`
3. Examine curated library TSVs in `results/curated_library/`
4. Upload trees from `results/tree/` to https://auspice.us/ for visualization
5. Identify:
   - Selected strains too similar (low Hamming/epitope distance to other library strains)
   - Non-selected strains too distant (high distance, indicating coverage gaps)
6. Adjust selections in config:
   - Modify selection columns in input TSVs
   - Add to `override_select_recent_haplotypes` for exclusions
   - Add to `additional_haplotypes` for inclusions
7. Re-run pipeline until satisfied with library composition

## Important Notes

- All selected haplotypes MUST have valid amino acids (no ambiguous or missing residues)
- Sequences must all be the same length for alignment
- The pipeline validates for duplicate strain names
- GenBank matching can take a while to run
- Tree JSON files are interactive and show selection status, distances, and metadata
