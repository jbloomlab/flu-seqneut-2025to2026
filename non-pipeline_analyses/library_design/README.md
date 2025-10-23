# Design of the library

The [./initial_design/](initial_design) contains the initial library design of strains designed to cover existing HA haplotypes as of late October 2025.
The key final file with the selected haplotypes is [initial_design/results/aggregated_library_strains/library_strains.tsv](initial_design/results/aggregated_library_strains/library_strains.tsv).

**To-do**:
- @ckikawa to add new subdirectory with designing the actual constructs to order from [initial_design/results/aggregated_library_strains/library_strains.tsv](initial_design/results/aggregated_library_strains/library_strains.tsv)
- In mid to late November, look at any new HA haplotypes with respect to those in [./initial_design/](initial_design) to see if more should be added; if so they would be added via a new subdirectory called something like `design_updates`. Note that for visualization, similar code to in [./initial_design/](initial_design) could be used specifying the already designed strains as *additional_haplotypes* in the `config.yaml`
