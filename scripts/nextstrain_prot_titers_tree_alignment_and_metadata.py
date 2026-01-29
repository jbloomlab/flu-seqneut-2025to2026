"""Build alignment, metadata, and titers for nextstrain-prot-titers-tree.

This script prepares input files for the nextstrain-prot-titers-tree pipeline:
- Alignment FASTA with protein sequences (with optional prefix for H1N1)
- Metadata TSV with strain information and titer summary columns
- Titers TSV with per-serum titer data for tree overlay (when titers configured)

For titers, uses sera_multicohort.csv so each serum appears once per cohort it
belongs to (original cohort, "All", and any days-post-vax cohorts), enabling
grouping/filtering by all cohorts in the tree measurements panel.
"""

import datetime
import sys

import numpy as np
import pandas as pd


sys.stdout = sys.stderr = open(snakemake.log[0], "w")

subtypes = snakemake.params.subtypes
circulating_strain_type = snakemake.params.circulating_strain_type
recent_vaccine_strains = snakemake.params.recent_vaccine_strains
prefix_alignment = snakemake.params.prefix_alignment
frac_below_cols = snakemake.params.frac_below_cols
serum_cohorts_for_tree = snakemake.params.serum_cohorts_for_tree
has_titers = snakemake.params.has_titers

viruses = pd.read_csv(snakemake.input.viral_libraries_csv)[
    [
        "strain",
        "subtype",
        "derived_haplotype",
        "strain_type",
        "protein_sequence_HA_ectodomain",
        "subclade",
        "collection_date",
    ]
].drop_duplicates()

assert len(viruses) == viruses["strain"].nunique(), "Duplicate strain entries found"

# Validate that recent_vaccine_strains are in the viral library
if recent_vaccine_strains:
    assert set(recent_vaccine_strains).issubset(
        viruses["strain"]
    ), f"recent_vaccine_strains not found in viral library: {set(recent_vaccine_strains) - set(viruses['strain'])}"

# Filter to circulating strains and recent vaccine strains
df = viruses[
    (viruses["strain_type"] == circulating_strain_type)
    | viruses["strain"].isin(recent_vaccine_strains)
].copy()

# Relabel strain_type for vaccine strains using the label from recent_vaccine_strains dict
if recent_vaccine_strains:
    df["strain_type"] = df.apply(
        lambda x: recent_vaccine_strains.get(x["strain"], x["strain_type"]),
        axis=1,
    )

print(
    f"{len(df)=} of {len(viruses)} are {circulating_strain_type=} or in {recent_vaccine_strains=}"
)

# Ensure collection_date is in valid format (numerical year)
year = datetime.datetime.now().year
if all((df["collection_date"] > year - 100) & (df["collection_date"] < year + 1)):
    df = df.rename(columns={"collection_date": "date"})
else:
    raise ValueError(f"Not valid numerical dates in {df['collection_date'].tolist()}")

# Process titer data if available
if has_titers:
    print("\nProcessing titer data...")

    # Read summarized titers (for metadata columns)
    summarized_titers = pd.read_csv(snakemake.input.summarized_titers_csv).rename(
        columns={"virus": "strain"}
    )
    print(f"Read {len(summarized_titers)=} summarized titer rows")

    # Validate frac_below_cols exist
    assert set(frac_below_cols).issubset(
        summarized_titers.columns
    ), f"{frac_below_cols=} not all in {summarized_titers.columns.tolist()=}"

    # Filter summarized titers to only include strains in tree (may exclude historical vaccines)
    strains_not_in_tree = set(summarized_titers["strain"]) - set(df["strain"])
    if strains_not_in_tree:
        print(
            f"Filtering out {len(strains_not_in_tree)} strains not in tree: {strains_not_in_tree}"
        )
        summarized_titers = summarized_titers[
            summarized_titers["strain"].isin(df["strain"])
        ]
        print(f"After filtering: {len(summarized_titers)=} rows")

    # Pivot median_titer and frac_below columns by cohort and merge into metadata
    for col in ["median_titer"] + frac_below_cols:
        assert col not in df.columns, f"{col} already in df columns"
        pivoted = (
            summarized_titers.assign(
                cohort_col=lambda x: f"{col}_" + x["cohort"] + "_sera"
            )
            .pivot_table(index="strain", values=col, columns="cohort_col")
            .reset_index()
        )
        df = df.merge(pivoted, on="strain", how="left", validate="one_to_one")

    # Read individual titers (for per-serum titers TSV)
    titers = (
        pd.read_csv(snakemake.input.titers_csv)[["serum", "virus", "titer"]]
        .rename(columns={"virus": "strain"})
        .assign(log2_titer=lambda x: np.log2(x["titer"]))
    )
    print(f"Read {len(titers)=} individual titer rows")

    # Filter titers to only include strains in tree
    titers = titers[titers["strain"].isin(df["strain"])]
    print(f"After filtering to tree strains: {len(titers)=} rows")

    # Read sera metadata (multicohort format: each serum appears once per cohort)
    sera_multicohort = pd.read_csv(snakemake.input.sera_metadata_csv)[
        ["serum", "cohort", "serum_collection_date", "age_numeric", "sex"]
    ]
    # Validate each (serum, cohort) pair is unique
    assert len(sera_multicohort) == len(
        sera_multicohort[["serum", "cohort"]].drop_duplicates()
    ), "Duplicate (serum, cohort) pairs"
    assert set(titers["serum"]).issubset(
        sera_multicohort["serum"]
    ), f"Titers have sera not in metadata: {set(titers['serum']) - set(sera_multicohort['serum'])}"

    # Validate unique serum-strain pairs (before join expands rows)
    assert len(titers) == len(
        titers[["serum", "strain"]].drop_duplicates()
    ), "Duplicate serum-strain pairs in titers"

    # Join titers with sera multicohort (one-to-many: each titer row expands to
    # multiple rows, one per cohort the serum belongs to)
    n_titers_before = len(titers)
    titers = titers.merge(sera_multicohort, on="serum", validate="many_to_many")
    print(
        f"Joined titers with sera multicohort: {n_titers_before} -> {len(titers)} rows "
        f"({titers['serum'].nunique()} sera, {titers['cohort'].nunique()} cohorts)"
    )

    # Create unique serum identifier for tree (serum + cohort) since
    # nextstrain-prot-titers-tree requires unique serum identifiers
    titers["serum_for_tree"] = titers["serum"] + "_" + titers["cohort"]
    assert len(titers) == len(
        titers[["serum_for_tree", "strain"]].drop_duplicates()
    ), "Duplicate serum_for_tree-strain pairs after cohort expansion"

# Process each subtype
for subtype in subtypes:
    print(f"\nProcessing {subtype=}")
    subtype_df = df[df["subtype"] == subtype].drop(columns="subtype")
    print(f"{len(subtype_df)=} of {len(df)=} are {subtype=}")

    if len(subtype_df) == 0:
        raise ValueError(f"No strains found for {subtype=}")

    # Remove subtype suffix from strain names if present (e.g., "_H3N2")
    # Keep original names for filtering titers before renaming
    original_strain_names = set(subtype_df["strain"])
    strain_rename = {
        s: (s[: -len(subtype) - 1] if s.endswith(f"_{subtype}") else s)
        for s in subtype_df["strain"]
    }
    subtype_df["strain"] = subtype_df["strain"].map(strain_rename)
    assert len(subtype_df) == subtype_df["strain"].nunique()

    alignment_file = snakemake.output[f"alignment_{subtype}"]
    metadata_file = snakemake.output[f"metadata_{subtype}"]

    print(f"Writing alignment to {alignment_file=}")
    with open(alignment_file, "w") as f:
        for tup in subtype_df.itertuples():
            seq = prefix_alignment[subtype] + tup.protein_sequence_HA_ectodomain
            f.write(f">{tup.strain}\n{seq}\n")

    print(f"Writing metadata to {metadata_file=}")
    (
        subtype_df.drop(columns=["protein_sequence_HA_ectodomain"]).to_csv(
            metadata_file, index=False, sep="\t", float_format="%.6g"
        )
    )

    # Write titers TSV if titer data available
    if has_titers:
        titers_file = snakemake.output[f"titers_{subtype}"]
        # Filter using original strain names (before suffix removal)
        subtype_titers = titers[titers["strain"].isin(original_strain_names)]
        print(
            f"{len(subtype_titers)} of {len(titers)} titer rows are for {subtype=} strains "
            f"({subtype_titers[['serum', 'strain']].drop_duplicates().shape[0]} unique serum-strain pairs)"
        )

        print(f"Writing titers to {titers_file=}")
        (
            subtype_titers.assign(
                strain=lambda x: x["strain"].map(strain_rename),
                serum=lambda x: x["serum_for_tree"],  # Use cohort-specific serum ID
            ).to_csv(titers_file, sep="\t", index=False, float_format="%.6g")
        )

print("\nDone!")
