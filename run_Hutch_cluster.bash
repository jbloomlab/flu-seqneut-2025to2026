#!/bin/bash
#
#SBATCH -c 32

snakemake \
    -j 32 \
    --rerun-incomplete \
    --software-deployment-method conda \
    --keep-going
