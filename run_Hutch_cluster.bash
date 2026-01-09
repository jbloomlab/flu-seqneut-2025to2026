#!/bin/bash
#
#SBATCH -c 16

snakemake \
    --use-conda \
    -j 16 \
    --rerun-incomplete \
    --software-deployment-method conda \
    --keep-going
