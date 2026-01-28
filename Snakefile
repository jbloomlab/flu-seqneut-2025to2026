"""Top-level ``snakemake`` file that runs analysis."""

import os

import pandas as pd
from os.path import join


configfile: "config.yml"


include: "seqneut-pipeline/seqneut-pipeline.smk"


# auspice JSONs as target output
auspice_jsons = []
if "nextstrain-prot-titers-tree_config" in config:
    for d in config["nextstrain-prot-titers-tree_config"].values():
        auspice_jsons.append(d["auspice_json"])
        if d.get("titers"):  # only add measurements.json if titers configured
            auspice_jsons.append(
                os.path.splitext(d["auspice_json"])[0] + "_measurements.json"
            )


rule all:
    input:
        # output from seqneut-pipeline
        seqneut_pipeline_outputs,
        # validation of viral libraries
        expand(
            "results/validate_viral_library/{viral_library}_validation.txt",
            viral_library=config["viral_libraries"],
        ),
        # auspice JSONs from nextstrain-prot-titers-tree
        auspice_jsons,
        # aggregated sera metadata
        "results/sera_metadata/all_sera_metadata.csv",
        # final processed titer data
        expand(
            "results/final_titer_data/{group}_{output_type}.csv",
            group=groups,
            output_type=["titers", "sera", "viruses", "titers_summarized_by_virus"],
        ),
        expand("results/final_titer_data/{group}_summary.txt", group=groups),


# =======================================================================================
# Additional rules outside seqneut-pipeline
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


# Aggregate and validate sera metadata --------------------------------------------------


rule aggregate_sera_metadata:
    """Aggregate and validate sera metadata from multiple cohorts."""
    input:
        csvs=config["sera_metadata"],
    output:
        csv="results/sera_metadata/all_sera_metadata.csv",
    log:
        "results/logs/aggregate_sera_metadata.txt",
    conda:
        "seqneut-pipeline/environment.yml"
    script:
        "scripts/aggregate_sera_metadata.py"


# Process and QC final titer data -----------------------------------------------


rule process_final_titer_data:
    """Process and QC final titer data for each group."""
    input:
        sera_metadata="results/sera_metadata/all_sera_metadata.csv",
        viral_library=lambda wc: config["viral_libraries"][
            config["process_final_titer_data"]["viral_library"]
        ],
        titers="results/aggregated_titers/titers_{group}.csv",
    output:
        titers="results/final_titer_data/{group}_titers.csv",
        sera="results/final_titer_data/{group}_sera.csv",
        viruses="results/final_titer_data/{group}_viruses.csv",
        titers_summarized="results/final_titer_data/{group}_titers_summarized_by_virus.csv",
        summary="results/final_titer_data/{group}_summary.txt",
    params:
        config=config["process_final_titer_data"],
    log:
        "results/logs/process_final_titer_data_{group}.txt",
    conda:
        "seqneut-pipeline/environment.yml"
    script:
        "scripts/process_final_titer_data.py"


# Build nextstrain-prot-titers-tree inputs ----------------------------------------------


rule nextstrain_prot_titers_tree_alignment_and_metadata:
    """Build alignment and metadata used by `nextstrain-prot-titers-tree`."""
    input:
        viral_libraries_csv=config["viral_libraries"]["designed"],
    output:
        **{
            f"alignment_{subtype}": f"results/nextstrain-prot-titers-tree/{subtype}/alignment.fa"
            for subtype in config["subtypes"]
        },
        **{
            f"metadata_{subtype}": f"results/nextstrain-prot-titers-tree/{subtype}/metadata.tsv"
            for subtype in config["subtypes"]
        },
    params:
        subtypes=config["subtypes"],
        circulating_strain_type=config["circulating_strain_type"],
        recent_vaccine_strains=config["recent_vaccine_strains"],
        prefix_alignment=config["nextstrain-prot-titers-tree_prefix_alignment"],
    conda:
        "seqneut-pipeline/environment.yml"
    log:
        "results/logs/nextstrain_prot_titers_tree_alignment_and_metadata.txt",
    script:
        "scripts/nextstrain_prot_titers_tree_alignment_and_metadata.py"


# run the nextstrain-prot-titers-tree submodule on each lineage
for subtype in config["subtypes"]:
    module_name = f"nextstrain-prot-titers-tree_{subtype}"

    module:
        name: module_name
        snakefile:
            "nextstrain-prot-titers-tree/Snakefile"
        config:
            config["nextstrain-prot-titers-tree_config"][subtype]

    use rule * from module_name as module_name*
