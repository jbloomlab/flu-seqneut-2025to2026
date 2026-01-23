"""Build alignment and metadata for nextstrain-prot-titers-tree.

This script prepares input files for the nextstrain-prot-titers-tree pipeline:
- Alignment FASTA with protein sequences (with optional prefix for H1N1)
- Metadata TSV with strain information

Currently builds trees without titer data. When titer data becomes available,
this script can be extended to include titer columns in the metadata.
"""

import datetime
import sys

import pandas as pd


sys.stdout = sys.stderr = open(snakemake.log[0], "w")

subtypes = snakemake.params.subtypes
circulating_strain_type = snakemake.params.circulating_strain_type
recent_vaccine_strains = snakemake.params.recent_vaccine_strains
prefix_alignment = snakemake.params.prefix_alignment

viruses = pd.read_csv(snakemake.input.viral_libraries_csv)[
    [
        "strain",
        "subtype",
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

# Process each subtype
for subtype in subtypes:
    print(f"\nProcessing {subtype=}")
    subtype_df = df[df["subtype"] == subtype].drop(columns="subtype")
    print(f"{len(subtype_df)=} of {len(df)=} are {subtype=}")

    if len(subtype_df) == 0:
        raise ValueError(f"No strains found for {subtype=}")

    # Remove subtype suffix from strain names if present (e.g., "_H3N2")
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

print("\nDone!")
