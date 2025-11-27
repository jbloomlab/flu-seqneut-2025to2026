# Update to library design suggested by Jesse Bloom, Nov-26-2025
This subdirectory contains a Snakemake pipeline for checking if we need updates to the initial library design described in [../initial_design](../initial_design).
It does this by using the latest sequence data to see if there are any new strains in current sequencing data not well covered in that initial design.

The suggested strains to possibly add are in [to_add/strains_to_add.tsv](to_add/strains_to_add.tsv).

## How to use and look at results
This pipeline works like that in [../initial_design](../initial_design) with the following changes (see [config.yaml](config.yaml)):

 - the *recent_haplotypes* are a more recent set generated on Nov-25-2025:
   + complete HA from A H3N2 human influenza with submission date after Oct-15-2025 and collection date after Aug-15-2025.
   + complete HA from A H1N1 pdm09 human influenza with submission date after Oct-15-2025 and collection date after Aug-15-2025.
 - the *additional_haplotypes* are everything selected in [../initial_design](../initial_design)

This is basically a copy of [../initial_design](../initial_design) with new input data and the designed library we already ordered specified as additional haplotypes in [config.yaml](config.yaml).

So run the regular pipeline in [Snakefile](Snakefile) by building the `conda` environment in [environment.yml](environment.yml) and then running:

    snakemake --use-conda -j <ncpus> -s Snakefile

This will create [./results/](results).

The post-process those results to select strains to add.
This post-processing is done primarily by the marimo script [select_additional_strains.py](select_additional_strains.py), which has the details of how this is done hardcoded into it in a self-explanatory way.
Then to process the results of that to get nucleotide sequences, run [Snakefile_addtl_strains](Snakefile_addtl_strains) with:

    snakemake --use-conda -j <ncpus> -s Snakefile_addtl_strains

The result of all of this is the file [to_add/strains_to_add.tsv](to_add/strains_to_add.tsv) which contains the strains to maybe add.

Admittedly, this post-processing part is a bit hacky, but right now it works as we just need to do some simple post-processing of the regular pipeline results.

### Configuration
All pipeline configuration is in [config.yaml](config.yaml), which specifies input data and various configurations and options, and should be largely self explanatory.

### Input data
All input data are in [./data/](data):

 - Recent HA haplotypes, generated on Nov-25-2025 by downloading from GISAID all complete HA proteins from H3N2 or H1N1 pdm09 from humans with a submission date after Oct-15-2025 and a collection date after Aug-15-2025, and then processing to create the haplotype CSVs:
   + [data/recent_H3N2_haplotypes.tsv](data/recent_H3N2_haplotypes.tsv)
   + [data/recent_H1N1_haplotypes.tsv](data/recent_H1N1_haplotypes.tsv)
 - All strains matching each HA haplotype in the above TSVs:
   + [data/haplotype_strains_H3N2.tsv](data/haplotype_strains_H3N2.tsv)
   + [data/haplotype_strains_H1N1.tsv](data/haplotype_strains_H1N1.tsv)
 - Outgroup sequences for phylogenetic tree visualization:
   + [H1N1_outgroup.fa](data/H1N1_outgroup.fa)
   + [H3N2_outgroup.fa](data/H3N2_outgroup.fa)

### Submodules
Two submodules are included; note these are included via `git submodule` at the top-level:
 - [match_prot_to_genbank_nt](match_prot_to_genbank_nt): Matches protein sequences to closest GenBank nucleotide sequences
 - [nextstrain-prot-titers-tree](nextstrain-prot-titers-tree): Builds nextstrain phylogenetic trees from protein alignments
