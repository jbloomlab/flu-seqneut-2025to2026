"""Top-level ``snakemake`` file that runs analysis."""

import pandas as pd
from os.path import join

configfile: "config.yml"

include: "seqneut-pipeline/seqneut-pipeline.smk"

rule all:
    input:
        # output from seqneut-pipeline
        seqneut_pipeline_outputs,
        # validation of viral libraries
        expand(
            "results/validate_viral_library/{viral_library}_validation.txt",
            viral_library=config["viral_libraries"],
        ),


# =======================================================================================
# Additional rules outside make seqneut-pipeline
# =======================================================================================

# Validate and process viral libraries --------------------------------------------------

rule validate_viral_library:
    """Validate a viral library."""
    input:
        csv=lambda wc: config["viral_libraries"][wc.viral_library],
    output:
        validation="results/validate_viral_library/{viral_library}_validation.txt",
    log:
        "results/logs/validate_viral_library_{viral_library}.txt",
    params:
        circulating_strain_type=config["circulating_strain_type"],
    conda:
        "seqneut-pipeline/environment.yml"
    script:
        "scripts/validate_viral_library.py"
