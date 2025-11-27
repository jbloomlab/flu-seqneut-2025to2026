"""``snakemake`` file for ``nextstrain-prot-titers-tree``."""

import os.path
import shlex


# get some variables from config
results_subdir = config["results_subdir"]
log_subdir = os.path.join(results_subdir, "logs")
titers = config["titers"]
if titers == "null":
    titers = None

# if we are getting titers, they are stored in file similar to the
# auspice_json but with _measurements suffix as a sidecar JSON
auspice_json = config["auspice_json"]
auspice_titers_json = os.path.splitext(auspice_json)[0] + "_measurements.json"


rule all:
    input:
        [auspice_json, auspice_titers_json] if titers else [auspice_json],


rule prep_and_sanitize_data:
    """Prep alignment and remove special characters from strain names."""
    input:
        **({"titers": titers["titers_tsv"]} if titers else {}),
        alignment=config["alignment"],
        outgroup=config["outgroup"],
        metadata=config["metadata"],
    output:
        **({"titers": os.path.join(results_subdir, "titers.tsv")} if titers else {}),
        alignment=os.path.join(results_subdir, "alignment_w_outgroup.fa"),
        metadata=os.path.join(results_subdir, "metadata.tsv"),
    params:
        **(
            {
                "titer_cols": [
                    "strain",
                    "serum",
                    titers["titer_col"],
                    *[
                        c
                        for c in titers["grouping_columns"]
                        if c not in {"strain", "serum", titers["titer_col"]}
                    ],
                ]
            }
            if titers
            else {}
        ),
        color_by_metadata=list(config["color_by_metadata"]),
        have_titers=bool(titers),
    log:
        os.path.join(log_subdir, "prep_and_sanitize_data.txt"),
    conda:
        "environment.yml"
    script:
        "scripts/prep_and_sanitize_data.py"


rule tree:
    """Build tree using Poisson model that weights all amino-acid mutations equally."""
    input:
        alignment=rules.prep_and_sanitize_data.output.alignment,
    output:
        outdir=directory(os.path.join(results_subdir, "iqtree")),
        raw_tree=os.path.join(results_subdir, "raw_tree.nwk"),
    params:
        prefix=lambda _, output: os.path.join(output.outdir, "tree"),
    threads: 8
    log:
        os.path.join(log_subdir, "tree.txt"),
    conda:
        "environment.yml"
    shell:
        """
        mkdir -p {output.outdir} &> {log}
        iqtree \
            -s {input.alignment} \
            --seqtype AA \
            -m Poisson \
            -fast \
            --seed 1 \
            -T {threads} \
            --prefix {params.prefix} \
            &>> {log}
        cp {params.prefix}.treefile {output.raw_tree} &>> {log}
        """


rule treetime:
    """Use `treetime` to refine tree and get ancestral states / inferred mutations.
    
    We use `treetime` directly rather than `augur refine` and `augur ancestral`
    in order to be able to process amino-acid sequences.
    """
    input:
        raw_tree=rules.tree.output.raw_tree,
        metadata=rules.prep_and_sanitize_data.output.metadata,
        alignment=rules.prep_and_sanitize_data.output.alignment,
    output:
        timetree_w_outgroup=os.path.join(results_subdir, "treetime/timetree.nexus"),
        divtree_w_outgroup=os.path.join(
            results_subdir, "treetime/divergence_tree.nexus"
        ),
        ancestral_sequences=os.path.join(
            results_subdir, "treetime/ancestral_sequences.fasta"
        ),
        dates=os.path.join(results_subdir, "treetime/dates.tsv"),
    params:
        outdir=lambda _, output: os.path.dirname(output.timetree_w_outgroup),
    log:
        os.path.join(log_subdir, "treetime.txt"),
    conda:
        "environment.yml"
    shell:
        """
        treetime \
            --tree {input.raw_tree} \
            --aln {input.alignment} \
            --dates {input.metadata} \
            --name-column strain \
            --date-column date \
            --reroot outgroup \
            --branch-length-mode input \
            --aa \
            --gtr jtt92 \
            --rng-seed 1 \
            --outdir {params.outdir} \
            &> {log}
        """


rule process_treetime_output:
    """Process the output of `treetime` to get data suitable for augur."""
    input:
        timetree_w_outgroup=rules.treetime.output.timetree_w_outgroup,
        divtree_w_outgroup=rules.treetime.output.divtree_w_outgroup,
        ancestral_sequences=rules.treetime.output.ancestral_sequences,
        dates=rules.treetime.output.dates,
        site_numbering_map=config["site_numbering_map"],
    output:
        timetree=os.path.join(results_subdir, "timetree.nwk"),
        divtree=os.path.join(results_subdir, "divergence_tree.nwk"),
        brlens=os.path.join(results_subdir, "brlens.json"),
        aa_muts=os.path.join(results_subdir, "aa_muts.json"),
    log:
        os.path.join(log_subdir, "process_treetime_output.txt"),
    conda:
        "environment.yml"
    script:
        "scripts/process_treetime_output.py"


rule auspice_config:
    """Build auspice configuration file."""
    input:
        metadata=rules.prep_and_sanitize_data.output.metadata,
    output:
        auspice_config=os.path.join(results_subdir, "auspice_config.json"),
    params:
        display_defaults=config["display_defaults"],
        color_by_metadata=config["color_by_metadata"],
        has_titers=bool(titers),
    log:
        os.path.join(log_subdir, "auspice_config.txt"),
    conda:
        "environment.yml"
    script:
        "scripts/auspice_config.py"


rule export:
    """Export the JSON for visualization with `auspice`."""
    input:
        tree=rules.process_treetime_output.output.divtree,
        brlens=rules.process_treetime_output.output.brlens,
        aa_muts=rules.process_treetime_output.output.aa_muts,
        metadata=rules.prep_and_sanitize_data.output.metadata,
        auspice_config=rules.auspice_config.output.auspice_config,
    output:
        auspice_json=auspice_json,
    params:
        color_by_metadata_args=(
            "--color-by-metadata "
            + " ".join(shlex.quote(col) for col in config["color_by_metadata"])
            if config["color_by_metadata"]
            else ""
        ),
        addtl_export_args=" ".join(
            f"--{k} {shlex.quote(v)}" for k, v in config["addtl_export_args"].items()
        ),
    log:
        os.path.join(log_subdir, "export.txt"),
    conda:
        "environment.yml"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --node-data {input.brlens} {input.aa_muts} \
            --include-root-sequence-inline \
            --auspice-config {input.auspice_config} \
            --metadata {input.metadata} \
            {params.color_by_metadata_args} \
            {params.addtl_export_args} \
            --output {output.auspice_json} \
            &> {log}
        """


if titers:

    rule export_measurements:
        """Export the JSON of titer measurements."""
        input:
            titers=rules.prep_and_sanitize_data.output.titers,
        output:
            auspice_titers_json=auspice_titers_json,
        params:
            titer_col=titers["titer_col"],
            grouping_columns=titers["grouping_columns"],
            addtl_args=" ".join(
                f"--{k} {shlex.quote(v)}" for k, v in titers["addtl_args"].items()
            ),
        log:
            os.path.join(log_subdir, "export_measurements.txt"),
        conda:
            "environment.yml"
        shell:
            """
            augur measurements export \
                --collection {input.titers} \
                --output-json {output.auspice_titers_json} \
                --value-column {params.titer_col} \
                --grouping-column {params.grouping_columns:q} \
                {params.addtl_args} \
                &> {log}
            """
