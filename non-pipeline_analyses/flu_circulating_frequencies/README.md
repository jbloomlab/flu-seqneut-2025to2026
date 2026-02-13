# Analyze how well circulating sequences are represented in library

## Overview 
The goal of this analysis is to determine how well all known recent HA sequences are represented in the library.
To do this, we obtain all recent available sequeences.
Then we compare protein sequences beteween the set of sequences and those in our library, and identify exact and close matches. 

## Input data 
The configuration for the analysis is in [config.yaml](config.yaml), which as more explanation about what is in each input file.

Note that the input data in `./data/` are **not** tracked in this repo due to data sharing rules.

## Workflow
First, build and activate the conda environment with:

    conda env create -f environment.yml
    conda activate flu_circulating_frequencies

Configure the analysis in [config.yml](config.yml).
Then run the pipeline with:
        
    snakemake -j 16 --software-deployment-method conda

### Output
All the results are placed in [results](results) and are organized by the `group` flag designated in the analysis configuration. 
