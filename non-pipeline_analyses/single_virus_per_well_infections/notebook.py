import marimo

__generated_with = "0.17.2"
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
        notebook_directory,
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
    # Aggregate and plot barcode counts across wells
    counts = (
        pd.concat([pd.read_csv(c).assign(well=c.replace(countsdir,'').replace('_counts.csv','')) for c, s in zip(count_csvs, samples)])
        .merge(samples_df, validate="many_to_one", on="well")
        .drop(columns=["replicate"])
        # .assign(sample_well=lambda x: x["sample"] + " (" + x["well"] + ")")
    )

    # classify barcodes as viral or neut standard
    barcode_class = pd.concat(
        [
            pd.read_csv(viral_library_file)[["barcode", "strain", 'name']].assign(
                neut_standard=False
            ),
            pd.read_csv(neut_standard_file)[["barcode"]].assign(
                neut_standard=True,
                strain=pd.NA,
                name=pd.NA,
            ),
        ],
        ignore_index=True,
    )

    # Merge counts and classification of barcodes
    assert set(counts["barcode"]) == set(barcode_class["barcode"])
    counts = (counts
        .merge(barcode_class, on="barcode")
              .rename(columns={'name': 'barcode_strain_link'})
             )

    # Map barcodes to see if they are expected, unexpected (linked to other strains), or unlinked to any strain at all
    linked_barcode_class_list = []

    for index, row in counts.iterrows():
        expected_strain_in_well = (row.well)

        if pd.notna(row.barcode_strain_link):
            strain_that_barcode_links = (row.barcode_strain_link.split('_bc')[0])
        else:
            strain_that_barcode_links = 'no linked barcode'

        if expected_strain_in_well == strain_that_barcode_links:
            linked_barcode_class_list.append([row.barcode, row.well, 'expected linked barcode'])

        elif expected_strain_in_well != strain_that_barcode_links:
            if strain_that_barcode_links == 'no linked barcode':
                linked_barcode_class_list.append([row.barcode, row.well, 'no linked barcode'])
            else:
                linked_barcode_class_list.append([row.barcode, row.well, f'unexpected linked barcode, {strain_that_barcode_links}'])

    linked_barcode_class=pd.DataFrame(linked_barcode_class_list, columns=['barcode', 'well', 'linked_barcode_class'])

    # Merge counts and classification of barcodes
    counts = (counts
        .merge(linked_barcode_class, on=["barcode", 'well'])
             )
    return (counts,)


@app.cell
def _(alt, counts, notebook_directory: "Path"):
    # Custom color scheme
    custom_colors = [
        '#B8B8B8',  # Muted gray for first category
        '#D4C4B0',  # Muted beige for second category
    ]

    # Subset counts data to just those barcodes with at least X counts
    counts_for_barchart = counts.query("count > 10") # At least 10 counts

    # Plot as colored bar chart
    bar_chart = alt.Chart(counts_for_barchart).mark_bar(height={"band": 0.85}).encode(
        x = 'fraction_all_valid_and_invalid_counts',
        y = 'well',
        color=alt.Color('linked_barcode_class', 
                        scale=alt.Scale(range=['#B8B8B8', '#D4C4B0'] +  # 2 muted colors
                                        ['#E63946', '#457B9D', '#2A9D8F', '#F77F00', '#9B59B6',  # vibrant
                                         '#06D6A0', '#EF476F', '#118AB2', '#FFD60A', '#7209B7'] * 15  # repeat as needed
                                       ),
                        legend=alt.Legend(labelLimit=500, columns=3, symbolLimit=0)),
        tooltip = ['strain', 'fraction_all_valid_and_invalid_counts', 'barcode','linked_barcode_class'],
    ).configure_axis(grid=False).properties(
            height = alt.Step(10),
            width = 200,
            title = "Barcode fraction with expected vs. unexpected categorization for each strain in flu-seqneut-2025to2026 library")

    # Save chart to HTML file instead of displaying inline
    bar_chart.save(notebook_directory / './barcode_chart.html')
    return


@app.cell
def _(counts):
    # Print some stats
    wells_with_low_level_contamination = []
    wells_with_measurable_contamination = []
    for well in counts.well.unique():
        low_level_df = counts.query(f'well == "{well}"').query("count > 10")
        measurable_df = counts.query(f'well == "{well}"').query("fraction_all_valid_and_invalid_counts > 0.001")
        if len(measurable_df.linked_barcode_class.unique()) > 2:
            wells_with_measurable_contamination.append(well)
        if len(low_level_df.linked_barcode_class.unique()) > 2:
            wells_with_low_level_contamination.append(well)

    print('Number of strains with barcode contamination 10 counts or greater:', len(wells_with_low_level_contamination))
    print('Number of strains with barcode contamination 0.001 fraction of well reads or greater:', len(wells_with_measurable_contamination))

    print('Displaying the wells with measurable contamination...')
    counts[counts['well'].isin(wells_with_measurable_contamination)].query('linked_barcode_class != "no linked barcode"').query("fraction_all_valid_and_invalid_counts > 0.001")
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
