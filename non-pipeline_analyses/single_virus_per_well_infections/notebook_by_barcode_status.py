import marimo

__generated_with = "0.17.6"
app = marimo.App(width="medium")


@app.cell
def _(mo):
    mo.md(r"""
    # Analysis notebook for single virus per well infections
    Author: Caroline Kikawa.

    See README and inline annotation for details.
    """)
    return


@app.cell
def _():
    import marimo as mo
    import os
    from pathlib import Path

    import sys
    import pandas as pd

    import altair as alt
    _ = alt.data_transformers.disable_max_rows()

    output_max_bytes = 13_000_000
    return Path, alt, mo, os, pd


@app.cell
def _(Path, mo, os):
    # Marimo path to notebook
    notebook_directory: Path = mo.notebook_dir()

    # Id input and output
    countsdir = './results/miscellaneous_plates/20260115_single_well/'
    samplesfile = './data/miscellaneous_plates/2026-01-15-single-well.csv'
    viral_library_file = './data/viral_libraries/flu-seqneut-2025to2026-barcode-to-strain-designed.csv'
    neut_standard_file = 'data/neut_standard_sets/loes2023_neut_standards.csv'

    # Identify all counts and fates CSVs
    count_csvs = []
    fate_csvs = []
    file_list = os.listdir(countsdir)
    for f in file_list:
        location = countsdir + f
        if "_counts" in f:
            count_csvs.append(location)
        elif "_fates" in f:
            fate_csvs.append(location)
    return (
        count_csvs,
        countsdir,
        fate_csvs,
        neut_standard_file,
        samplesfile,
        viral_library_file,
    )


@app.cell
def _(pd, samplesfile):
    # Define a samples dataframe using the samples file
    samples_df = pd.read_csv(samplesfile)
    samples_df.drop(columns=['fastq'], inplace=True)
    samples_df['sample'] = samples_df.apply(
        lambda x: '-'.join(x.astype(str)), axis=1
    )

    samples = samples_df["sample"].unique().tolist()
    print(f"There are {len(samples)} barcode runs.")
    return samples, samples_df


@app.cell
def _(alt, countsdir, fate_csvs, pd, samples, samples_df):
    # Aggregate and plot barcode fates across wells
    fates = (
        pd.concat([pd.read_csv(f).assign(well=f.replace(countsdir,'').replace('_fates.csv','')) for f, s in zip(fate_csvs, samples)])
        .merge(samples_df, validate="many_to_one", on="well")
        .assign(
            fate_counts=lambda x: x.groupby("fate")["count"].transform("sum"),
            # sample_well=lambda x: x["sample"] + " (" + x["well"] + ")",
        )
        .query("fate_counts > 0")[  # only keep fates with at least one count
            ["fate", "count", "well", "dilution_factor"]
        ]
    )

    assert len(fates) == len(fates.drop_duplicates())


    fates_chart = (
        alt.Chart(fates)
        .encode(
            alt.X("count", scale=alt.Scale(nice=False, padding=3)),
            alt.Y(
                "well",
                title=None,
            ),
            alt.Color("fate", sort=sorted(fates["fate"].unique(), reverse=True)),
            alt.Order("fate", sort="descending"),
            tooltip=fates.columns.tolist(),
        )
        .mark_bar(height={"band": 0.85})
        .properties(
            height=alt.Step(10),
            width=200,
            title=f"Barcode parsing for initial titering plate",
        )
        .configure_axis(grid=False)
    )

    fates_chart
    return


@app.cell
def _(
    count_csvs,
    countsdir,
    neut_standard_file,
    pd,
    samples,
    samples_df,
    viral_library_file,
):
    # Aggregate and plot barcode counts across wells, classifying each barcode as expected
    # in the well, from a virus in an adjacent well, a distant well, or not in plate
    counts = (
        pd.concat([pd.read_csv(c).assign(well=c.replace(countsdir,'').replace('_counts.csv','')) for c, s in zip(count_csvs, samples)])
        .merge(samples_df, validate="many_to_one", on="well")
        .drop(columns=["replicate"])
    )

    # classify barcodes as viral or neut standard
    viral_library = pd.read_csv(viral_library_file)
    barcode_class = pd.concat(
        [
            viral_library[["barcode", "strain"]],
            pd.read_csv(neut_standard_file)[["barcode"]].assign(strain="neut_standard"),
        ],
        ignore_index=True,
    )

    # Merge counts and classification of barcodes
    assert set(counts["barcode"]) == set(barcode_class["barcode"])
    counts = (
        counts
        .merge(barcode_class, on="barcode", validate="m:1")
        [["barcode", "count", "well", "strain"]]
        .rename(columns={"strain": "barcode_strain"})
    )

    well_to_strain = (
        viral_library
        .query("name.notnull()")
        .assign(well=lambda x: x["name"].str[: -4])
        [["strain", "well"]]
        .drop_duplicates()
    )
    assert set(counts["well"]).issubset(well_to_strain["well"])
    assert len(well_to_strain["well"]) == well_to_strain["well"].nunique()

    well_strain_number = (
        samples_df[["well"]].reset_index(names="well_number")
        .assign(well_number=lambda x: x["well_number"] + 1)
        .merge(well_to_strain, on="well", validate="1:1")
    )
    assert set(well_strain_number["well"]) == set(counts["well"])
    assert len(well_strain_number) == counts["well"].nunique()

    counts = (
        counts
        .merge(well_to_strain, on="well", validate="m:1")
        .rename(columns={"strain": "well_strain"})
        .merge(well_strain_number[["well", "well_number"]], on="well", validate="m:1")
        .merge(
            well_strain_number[["strain", "well_number"]].rename(
                columns={"strain": "barcode_strain", "well_number": "barcode_well_number"}
            ),
            on="barcode_strain",
            validate="many_to_one",
            how="left",
        )
    )

    def well_distance(
        w1: int,
        w2: int,
        n_cols: int = 12,
        n_rows: int = 8,
        row_first: bool = True,
    ) -> str:
        """
        Classify the relationship between two 1-indexed well numbers in a plate layout.

        If row_first=True (default): number left-to-right across columns, then top-to-bottom across rows.
          i.e., row-major: (row, col) -> idx = row*n_cols + col + 1

        If row_first=False: number top-to-bottom down rows, then left-to-right across columns.
          i.e., column-major: (row, col) -> idx = col*n_rows + row + 1

        Returns: "same", "adjacent" (including diagonals), or "distant".
        """
        total = n_cols * n_rows
        if not (1 <= w1 <= total) or not (1 <= w2 <= total):
            raise ValueError(f"Well numbers must be in [1, {total}]: {w1=}, {w2=}")

        if w1 == w2:
            return "same"

        def to_rc(w: int) -> tuple[int, int]:
            w0 = w - 1  # 0-index
            if row_first:
                r, c = divmod(w0, n_cols)   # row-major
            else:
                c, r = divmod(w0, n_rows)   # column-major
            return r, c

        r1, c1 = to_rc(w1)
        r2, c2 = to_rc(w2)

        return "adjacent" if (abs(r1 - r2) <= 1 and abs(c1 - c2) <= 1) else "distant"

    counts["barcode_status"] = counts.apply(
        lambda row: (
            "not in plate"
            if pd.isnull(row["barcode_well_number"])
            else well_distance(row["well_number"], row["barcode_well_number"])
        ),
        axis=1,
    )

    counts
    return (counts,)


@app.cell
def _(alt, counts, pd):
    # Plot how many counts are of each type (same, adjacent, etc)

    count_by_barcode_status = (
        counts
        .sort_values("count", ascending=False)
        .groupby(["well", "well_strain", "well_number", "barcode_status"], as_index=False)
        .aggregate(
            count=pd.NamedAgg("count", "sum"),
            top_strain_w_status=pd.NamedAgg("barcode_strain", "first"),
            top_strain_w_status_counts=pd.NamedAgg("count", "first"),
        )
        .assign(
            well_counts=lambda x: x.groupby("well")["count"].transform("sum"),
            frac=lambda x: x["count"] / x["well_counts"],
        )
        .melt(
            id_vars=["well", "well_number", "well_strain", "barcode_status", "top_strain_w_status", "top_strain_w_status_counts", "well_counts"],
            value_vars=["count", "frac"],
            var_name="stat_type",
            value_name="stat",
        )
    )

    count_by_barcode_status_chart = (
        alt.Chart(count_by_barcode_status)
        .encode(
            alt.X("stat"),
            alt.Y("well", sort=alt.SortField("well_number")),
            alt.Color("barcode_status"),
            alt.Column("stat_type", title=None),
            tooltip=["well", "well_number",
                     "well_strain", "barcode_status", "top_strain_w_status", "top_strain_w_status_counts", "well_counts"],
        )
        .mark_bar()
        .resolve_scale(x="independent")
        .properties(
            height=alt.Step(10),
            width=200,
            title="Barcode counts from virus in same well, adjacent well, distant well, or not in plate",
        )
    )

    count_by_barcode_status_chart.save("count_by_barcode_status.html")

    count_by_barcode_status_chart
    return


@app.cell
def _(counts):
    # Look at counts for A/Bangkok/P2277/2025
    counts.query("barcode_strain.str.contains('A/Bangkok/P2277/2025')").sort_values("count", ascending=False)
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
