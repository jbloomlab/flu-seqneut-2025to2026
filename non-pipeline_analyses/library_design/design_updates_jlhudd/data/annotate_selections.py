import marimo

__generated_with = "0.18.0"
app = marimo.App(width="full")


@app.cell
def _():
    # Annotate selections

    import pandas as pd

    import marimo as mo

    additions = {
        "H3N2": {
            "JohnHuddleston_selections": {
                "derived_haplotype": [
                    # previousy selected "K:V88I",
                    # previously selected "K:V88I,I214T",
                    "K:S145N",
                    "K:K278E",
                    # previously selected "K:S145N,A272T,R299K",
                ],
            },
            "JesseBloom_selections": {
                "derived_haplotype": ["K", "K:V102M"],
            }
        },
        "H1N1": {
            "JohnHuddleston_selections": {
                "derived_haplotype": ["D.3.1:R113K,A139D,S157L,E283K,K302E"],
            },
            "JesseBloom_selections": {
                "derived_haplotype": ["D.3.1:R113K,A139D,R205K,E283K,K302E", "C.1.9.3:T82K,P137S,A141V,I166V"],
            },
        },
    }

    already_selected = set(
        pd.read_csv("../../initial_design/results/aggregated_library_strains/library_strains.tsv", sep="\t")
        ["prot_sequence"]
    )

    for subtype, add in additions.items():

        haplotypes = pd.read_csv(f"{subtype}_haplotypes_original.tsv", sep="\t")
        mo.output.append(mo.md(f"## {subtype=}, {len(haplotypes)=}"))
        selection_names = []
        for selection_name, selections in add.items():
            selection_names.append(selection_name)
            haplotypes[selection_name] = False
            for col, vals in selections.items():
                assert col in haplotypes.columns
                for val in vals:
                    assert sum(haplotypes[col] == val) == 1, f"No {val=} in {col=}"
                    mo.output.append(f"annotating {col=} {val=} as {selection_name=}")
                    haplotypes[selection_name] = haplotypes[selection_name].where(haplotypes[col] != val, True)
            mo.output.append(f"Overall, for {subtype=} {selection_name=}, {haplotypes[selection_name].sum()=}")

        jdb_counts = pd.read_csv(f"../../design_updates_jdb/data/recent_{subtype}_haplotypes.tsv", sep="\t")[[
            "representative_strain_ha_sequence", "count"
        ]].rename(columns={"count": "count_jdb"})

        mo.output.append(
            haplotypes
            [haplotypes[selection_names].any(axis=1)]
            .merge(jdb_counts, on="representative_strain_ha_sequence", how="left", validate="one_to_one")
        )

        previously_selected = (
            haplotypes
            [haplotypes[selection_names].any(axis=1)]
            .query("representative_strain_ha_sequence.isin(@already_selected)")
        )
        if len(previously_selected):
            raise ValueError(previously_selected)

        (
            haplotypes
            .query("not representative_strain_ha_sequence.isin(@already_selected)")
            .to_csv(f"{subtype}_haplotypes.tsv", sep="\t", index=False)
        )


    return


if __name__ == "__main__":
    app.run()
