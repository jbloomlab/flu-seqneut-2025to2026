# Update to library design suggested by John Huddleston, Nov-26-2025
This subdirectory contains a Snakemake pipeline for checking if we need updates to the initial library design described in [../initial_design](../initial_design).
It does this by using the latest sequence data to see if there are any new strains in current sequencing data not well covered in that initial design.

The suggested strains to possibly add are in [to_add/strains_to_add.tsv](to_add/strains_to_add.tsv).

## How to use and look at results
This pipeline works like that in [../initial_design](../initial_design) with the following changes (see [config.yaml](config.yaml)):

 - the *recent_haplotypes* are a more recent set generated on Nov-26-2025:
 - the *additional_haplotypes* are everything selected in [../initial_design](../initial_design)

This is basically a copy of [../initial_design](../initial_design) with new input data and the designed library we already ordered specified as additional haplotypes in [config.yaml](config.yaml).

So run the regular pipeline in [Snakefile](Snakefile) by building the `conda` environment in [environment.yml](environment.yml) and then running:

    snakemake --use-conda -j <ncpus> -s Snakefile

This will create [./results/](results).

Then post-process the results with [get_new_strains.py](get_new_strains.py) to create [strains_to_add.tsv](strains_to_add.tsv) which has the new H3N2 strains to add.

### Configuration
All pipeline configuration is in [config.yaml](config.yaml), which specifies input data and various configurations and options, and should be largely self explanatory.

### Input data
All input data are in [./data/](data):

 - Recent set of HA haplotypes with various scores provided by John Huddleston on Nov-26-2025, and selections by Jesse Bloom or John Huddleston indicated, and haplotypes already in library removed:
   + [data/H1N1_haplotypes.tsv](data/H1N1_haplotypes.tsv)
   + [data/H3N2_haplotypes.tsv](data/H3N2_haplotypes.tsv)
 - Outgroup sequences for phylogenetic tree visualization:
   + [H1N1_outgroup.fa](data/H1N1_outgroup.fa)
   + [H3N2_outgroup.fa](data/H3N2_outgroup.fa)

### Submodules
Two submodules are included; note these are included via `git submodule` at the top-level:
 - [match_prot_to_genbank_nt](match_prot_to_genbank_nt): Matches protein sequences to closest GenBank nucleotide sequences
 - [nextstrain-prot-titers-tree](nextstrain-prot-titers-tree): Builds nextstrain phylogenetic trees from protein alignments
