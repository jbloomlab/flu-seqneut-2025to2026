"""Check input metadata and alignment and remove special charactes from strain names."""

import collections
import re
import sys

import Bio.SeqIO

import pandas as pd


sys.stdout = sys.stderr = open(snakemake.log[0], "w")

alignment = list(Bio.SeqIO.parse(snakemake.input.alignment, format="fasta"))
assert len(alignment), "empty alignment in {snakemake.input.alignment=}"

req_metadata_cols = {"strain", "date"}
metadata = pd.read_csv(snakemake.input.metadata, sep="\t")

if not req_metadata_cols.issubset(metadata.columns):
    raise ValueError(f"{metadata.columns=} lacks {req_metadata_cols=}")

color_by_metadata = snakemake.params.color_by_metadata
assert len(color_by_metadata) == len(set(color_by_metadata))
assert not req_metadata_cols.intersection(color_by_metadata)
if not set(color_by_metadata).issubset(metadata.columns):
    raise ValueError(f"{metadata.columns=} lacks {color_by_metadata=}")

# make sure dates all numeric
dates = pd.to_numeric(metadata["date"], errors="raise")
assert dates.notnull().all(), "date not all non-null"

assert len(metadata) == len(alignment)
assert len(metadata) == metadata["strain"].nunique()

strain_renames = {
    orig: re.sub(r"[^A-Za-z0-9._\-/]", "_", orig.replace("'", ""))
    for orig in metadata["strain"].unique()
}
if len(strain_renames) != len(set(strain_renames.values())):
    raise ValueError(f"re-named strains not unique:\n{strain_renames=}")

metadata = metadata.assign(strain=lambda x: x["strain"].map(strain_renames))
metadata.to_csv(snakemake.output.metadata, sep="\t")

outgroup = Bio.SeqIO.read(snakemake.input.outgroup, "fasta")
assert "outgroup" not in set(strain_renames.values())

seqlengths = collections.defaultdict(int)
with open(snakemake.output.alignment, "w") as f:
    for seq in alignment:
        strain = seq.id
        assert strain in strain_renames, strain
        strain = strain_renames[strain]
        seq = str(seq.seq)
        seqlengths[len(seq)] += 1
        f.write(f">{strain}\n{seq}\n")

    if len(seqlengths) != 1:
        raise ValueError(f"Not all sequences same length in alignment:\n{seqlengths=}")
    assert len(outgroup) == list(seqlengths)[0], f"{seqlengths=}, {len(outgroup)=}"
    f.write(f">outgroup\n{str(outgroup.seq)}\n")

if snakemake.params.have_titers:
    titers = pd.read_csv(snakemake.input.titers, sep="\t")
    titer_cols = snakemake.params.titer_cols
    if not set(titer_cols).issubset(titers.columns):
        raise ValueError(f"{titer_cols=} not all in {titers.columns=}")
    extra_strains = set(titers["strain"]) - set(strain_renames)
    if extra_strains:
        raise ValueError(
            f"strains w titers but not specified as part of tree: {extra_strains}"
        )
    titers = titers[titer_cols].assign(strain=lambda x: x["strain"].map(strain_renames))
    assert set(titers["strain"]).issubset(metadata["strain"])
    titers_per_strain_serum = (
        titers.groupby(["strain", "serum"])
        .aggregate(
            n_titers=pd.NamedAgg(titer_cols[2], "size"),
            titers=pd.NamedAgg(titer_cols[2], "unique"),
        )
        .sort_values("n_titers", ascending=False)
        .query("n_titers > 1")
    )
    if len(titers_per_strain_serum):
        raise ValueError(
            f"multiple titers for some strains/sera:\n{titers_per_strain_serum}"
        )
    titers.to_csv(snakemake.output.titers, sep="\t", float_format="%.5g")
