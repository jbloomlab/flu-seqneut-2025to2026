# Design of the library

The [./initial_design/](initial_design) subdirectory contains the initial library design of strains designed to cover existing HA haplotypes as of late October 2025.
The key final file with the selected haplotypes is [initial_design/results/aggregated_library_strains/library_strains.tsv](initial_design/results/aggregated_library_strains/library_strains.tsv).

The [./construct_order/](construct_order) subdirectory contains the code generating the actual insert sequences corresponding to each designed HA gene in [initial_design/results/aggregated_library_strains/library_strains.tsv](initial_design/results/aggregated_library_strains/library_strains.tsv) and a unique 16-nucleotide barcode. The key output files are [construct_order/results/H3_order1.csv](construct_order/results/H3_order1.csv) and [construct_order/results/H1_order1.csv](construct_order/results/H1_order1.csv) which are the ordersheets submitted to Twist Biosciences for synthesis. 

Possible design updates:
  - [./design_updates_jdb](design_updates_jdb) contains possible new haplotypes to add based on recent frequent HAs not in the library as suggested by Jesse Bloom on Nov-26-2025. The key file is [./design_updates_jdb/to_add/strains_to_add.tsv](design_updates_jdb/to_add/strains_to_add.tsv) which contains possible HAs to consider adding.
  - [./design_updates_jlhudd](design_updates_jlhudd) contains possible new haplotypes to add based on recent frequent HAs not in the library as suggested by John Huddleston on Nov-26-2025. The key file is [design_updates_jlhudd/strains_to_add.tsv](design_updates_jlhudd/strains_to_add.tsv) which contains H3N2 strains to add.
