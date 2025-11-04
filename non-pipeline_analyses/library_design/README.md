# Design of the library

The [./initial_design/](initial_design) contains the initial library design of strains designed to cover existing HA haplotypes as of late October 2025.
The key final file with the selected haplotypes is [initial_design/results/aggregated_library_strains/library_strains.tsv](initial_design/results/aggregated_library_strains/library_strains.tsv).

The [./construct_order/](construct_order) contains the code generating the actual insert sequences corresponding to each designed HA gene in [initial_design/results/aggregated_library_strains/library_strains.tsv](initial_design/results/aggregated_library_strains/library_strains.tsv) and a unique 16-nucleotide barcode. The key output files are [construct_order/results/H3_order1.csv](construct_order/results/H3_order1.csv) and [construct_order/results/H1_order1.csv](construct_order/results/H1_order1.csv) which are the ordersheets submitted to Twist Biosciences for synthesis. 

**To-do**:
- In mid to late November, look at any new HA haplotypes with respect to those in [./initial_design/](initial_design) to see if more should be added; if so they would be added via a new subdirectory called something like `design_updates`. Note that for visualization, similar code to in [./initial_design/](initial_design) could be used specifying the already designed strains as *additional_haplotypes* in the `config.yaml`
