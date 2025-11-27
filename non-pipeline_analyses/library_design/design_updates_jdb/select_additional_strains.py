import marimo

__generated_with = "0.18.0"
app = marimo.App(width="full")


@app.cell
def _(mo):
    mo.md(r"""
    # Select additional strains with high counts and substantial difference to existing library strains
    This notebook assumes the results of running the rest of the pipeline, and then is designed to be run manually.

    The selection of which strains to add is just done manually in the code below, based on:
     - counts of haplotype
     - Hamming distance from other strains
     - epitope distance from other strains
     - presence of HA1 mutation(s)

    The results are written to a TSV file and FASTA file.
    """)
    return


@app.cell
def _():
    import os

    import marimo as mo

    import pandas as pd

    subtypes = ["H3N2", "H1N1"]

    all_to_add = []

    for subtype in subtypes:
        mo.output.append(mo.md(f"## Analyzing {subtype=}"))

        # read unselected haplotypes to consider
        all_unselected = (
            pd.read_csv(f"results/curated_library/{subtype}_nonselected_haplotypes.tsv", sep="\t")
            .sort_values("count", ascending=False)
            .query("hamming_distance_nearest_library_strain > 0")
            .query("not representative_strain_ha_sequence.str.contains('X')")
        )

        # annotate haplotypes with strain names
        all_strains = pd.read_csv(f"data/haplotype_strains_{subtype}.tsv", sep="\t").set_index("haplotype")["strains"].to_dict()
        all_unselected = all_unselected.assign(
            strains=lambda x: x["representative_strain_ha_sequence"].map(all_strains)
        )

        # find differences from closest library strain
        all_selected = (
            pd.read_csv(f"results/curated_library/{subtype}_selected_haplotypes.tsv", sep="\t")
            .set_index("representative_strain_ha_sequence")["representative_strain"].to_dict()
        )
        site_annotations = (
            pd.read_csv(f"results/site_annotations/{subtype}_site_annotations.tsv", sep="\t")
            .set_index("sequential_site")
            [["protein", "protein_site"]]
            .to_dict(orient="index")
        )
        def get_closest(seq):
            closest = None
            for (s, n) in all_selected.items():
                diffs = " ".join(
                    [
                        f"{site_annotations[r]['protein']}: {wt}{site_annotations[r]['protein_site']}{mut}"
                        for r, (wt, mut) in enumerate(zip(s, seq, strict=True), start=1)
                        if wt != mut
                    ]
                )
                if closest is None:
                    closest = [diffs, n]
                elif len(diffs) < len(closest[0]):
                    closest = [diffs, n]
            return closest

        all_unselected = all_unselected.assign(
            closest_library_strain=lambda x: x["representative_strain_ha_sequence"].map(lambda s: get_closest(s)[1]),
            diffs_from_closest_library_strain=lambda x: x["representative_strain_ha_sequence"].map(lambda s: get_closest(s)[0])
        )

        # filter which unselected haplotypes to add
        epitope_cols = [c for c in all_unselected.columns if "epitope_distance" in c]
        to_add = (
            all_unselected
            [
                # high overall sequence count
                (all_unselected["count"] > 20)
                # modest sequence count and large Hamming distance
                | ((all_unselected["count"] > 5) & (all_unselected["hamming_distance_nearest_library_strain"] > 2))
                # modest sequence count and any epitope distance
                | ((all_unselected["count"] > 5) & all_unselected[epitope_cols].any(axis=1))
            ]
            # require at least one HA1 mutation
            .query("diffs_from_closest_library_strain.str.contains('HA1')")
            .assign(
                subtype=subtype,
                representative_strain=lambda x: x["representative_strain"] + "_" + x["subtype"],
            )
            [
                [
                    "subtype",
                    "representative_strain",
                    "count",
                    "latest_sequence",
                    "hamming_distance_nearest_library_strain",
                    *epitope_cols,
                    "closest_library_strain",
                    "diffs_from_closest_library_strain",
                    "strains",
                    "representative_strain_ha_sequence",
                ]
            ]
            .reset_index(drop=True)
        )
        all_to_add.append(to_add.drop(columns=epitope_cols))

        mo.output.append(to_add)

    all_to_add = pd.concat(all_to_add)
    to_add_tsv = "to_add/haplotypes_to_add.tsv"
    to_add_fasta = "to_add/haplotypes_to_add.fa"
    os.makedirs(os.path.dirname(to_add_tsv), exist_ok=True)
    print(f"Writing all sequences to add to {to_add_tsv} and {to_add_fasta}")
    all_to_add.to_csv(to_add_tsv, index=False, sep="\t")
    with open(to_add_fasta, "w") as f:
        for tup in all_to_add.itertuples(index=False):
            f.write(f">{tup.representative_strain}\n{tup.representative_strain_ha_sequence}\n")
    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
 
    """)
    return


if __name__ == "__main__":
    app.run()
