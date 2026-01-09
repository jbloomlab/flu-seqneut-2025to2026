"""Top-level ``snakemake`` file that runs analysis."""

import pandas as pd
from os.path import join

configfile: "config.yml"

include: "seqneut-pipeline/seqneut-pipeline.smk"

rule all:
    input:
        # output from seqneut-pipeline
        seqneut_pipeline_outputs


