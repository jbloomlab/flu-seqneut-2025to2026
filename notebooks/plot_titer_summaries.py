import marimo

__generated_with = "0.17.6"
app = marimo.App(width="full")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    # Titer summary plots

    Interactive Altair plots showing median titers, individual serum titers,
    interquartile ranges, and fraction of sera below titer cutoffs.
    """
    )
    return


@app.cell
def _():
    # Load context from pickled file.
    #
    # This cell supports multiple ways to provide context:
    # 1. Via command-line: marimo export html notebook.py -- --context-pickle path/to/context.pickle
    # 2. Via saved pickle: Manually save a context pickle to a dev location
    # 3. Stub context: If no pickle available, creates minimal empty context for exploration
    #
    # For interactive development with `marimo edit`, you can:
    # - Run the pipeline once to generate a real context pickle, then copy it to a dev location
    # - Or work with the stub context (downstream cells will show warnings/empty data)

    import argparse
    import os
    import pathlib
    import pickle
    import sys

    import marimo as mo

    # Check if context-pickle argument is provided (run by driver script)
    from_cmdline = "--context-pickle" in sys.argv

    if from_cmdline:
        # Running via driver script - parse args
        print("Loading context from command-line argument")
        p = argparse.ArgumentParser()
        p.add_argument("--context-pickle", required=True)
        args = p.parse_args()
        context_pickle_path = pathlib.Path(args.context_pickle)
    else:
        # Running in marimo edit - try to use development pickle
        print("Running in marimo edit mode")
        # set `context_pickle_path` to valid pickle if running via marimo edit
        context_pickle_path = None
        context_pickle_path = pathlib.Path(
            "results/titer_plots/human_plot_titer_summaries_vertical_context.pickle"
        )

    # Load context if pickle path exists and is valid
    if context_pickle_path and context_pickle_path.exists():
        print(f"Reading context from {context_pickle_path}")
        with open(context_pickle_path, "rb") as f_context:
            context = pickle.load(f_context)

        # Handle working directory
        context_workdir = context["workdir"]
        current_workdir = os.getcwd()

        if from_cmdline:
            # Running via snakemake - verify workdir matches
            if context_workdir != current_workdir:
                raise RuntimeError(
                    f"Context workdir mismatch!\n"
                    f"  Context was created in: {context_workdir}\n"
                    f"  Currently running in:   {current_workdir}\n"
                    f"This should not happen when running via Snakemake."
                )
            print(f"Verified working directory: {current_workdir}")
        else:
            # Running in marimo edit - change to context workdir
            if context_workdir and context_workdir != current_workdir:
                print(f"Changing directory from {current_workdir} to {context_workdir}")
                os.chdir(context_workdir)
            elif context_workdir:
                print(f"Already in correct working directory: {context_workdir}")
    else:
        # Create a minimal stub context for interactive development
        print("Creating minimal stub context that you need to complete")
        context = {
            "input": {},
            "output": {},
            "params": {},
            "wildcards": {},
            "threads": 1,
            "resources": {},
        }
    return context, mo


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Read data
    """
    )
    return


@app.cell
def _():
    import itertools
    import json

    import altair as alt
    import pandas as pd

    _ = alt.data_transformers.disable_max_rows()
    return alt, itertools, json, pd


@app.cell
def _(context):
    # Extract variables from context - raises KeyError if required keys missing
    titers_csv = context["input"]["titers_csv"]
    sera_csv = context["input"]["sera_csv"]
    sera_multicohort_csv = context["input"]["sera_multicohort_csv"]
    viruses_csv = context["input"]["viruses_csv"]
    recent_vaccine_strains = context["params"]["recent_vaccine_strains"]
    circulating_strain_type = context["params"]["circulating_strain_type"]
    plot_titer_summaries_params = context["params"]["plot_titer_summaries_params"]
    subtypes = context["params"]["subtypes"]
    facet_orientation = context["params"]["facet_orientation"]
    chart_htmls = context["output"]["chart_htmls"]
    return (
        chart_htmls,
        circulating_strain_type,
        facet_orientation,
        plot_titer_summaries_params,
        recent_vaccine_strains,
        sera_csv,
        sera_multicohort_csv,
        subtypes,
        titers_csv,
        viruses_csv,
    )


@app.cell
def _(
    circulating_strain_type,
    mo,
    pd,
    recent_vaccine_strains,
    sera_csv,
    sera_multicohort_csv,
    subtypes,
    titers_csv,
    viruses_csv,
):
    # Read titers
    titers = pd.read_csv(titers_csv)
    mo.output.append(mo.md(f"Read {len(titers)=} titers from {titers_csv=}"))

    # Validate required columns in titers
    required_titer_cols = {"serum", "virus", "titer"}
    missing_titer_cols = required_titer_cols - set(titers.columns)
    if missing_titer_cols:
        raise ValueError(f"titers_csv missing required columns: {missing_titer_cols}")

    # Read sera metadata (one row per serum)
    metadata = pd.read_csv(sera_csv)
    mo.output.append(mo.md(f"\nRead {len(metadata)=} sera from {sera_csv=}"))

    # Validate required columns in sera metadata
    required_sera_cols = {
        "serum",
        "cohort",
        "age_numeric",
        "serum_collection_date",
        "age",
        "sex",
    }
    missing_sera_cols = required_sera_cols - set(metadata.columns)
    if missing_sera_cols:
        raise ValueError(f"sera_csv missing required columns: {missing_sera_cols}")

    # Validate sera match between titers and metadata
    titers_sera = set(titers["serum"])
    metadata_sera = set(metadata["serum"])
    if titers_sera != metadata_sera:
        raise ValueError(
            f"Serum mismatch between titers and sera metadata.\n"
            f"  In titers but not metadata: {titers_sera - metadata_sera}\n"
            f"  In metadata but not titers: {metadata_sera - titers_sera}"
        )

    # Read sera multicohort (multiple rows per serum, one per cohort assignment)
    sera_multicohort = pd.read_csv(sera_multicohort_csv)
    mo.output.append(
        mo.md(f"\nRead {len(sera_multicohort)=} rows from {sera_multicohort_csv=}")
    )

    # Validate sera match between titers and multicohort
    multicohort_sera = set(sera_multicohort["serum"])
    if titers_sera != multicohort_sera:
        raise ValueError(
            f"Serum mismatch between titers and sera_multicohort.\n"
            f"  In titers but not multicohort: {titers_sera - multicohort_sera}\n"
            f"  In multicohort but not titers: {multicohort_sera - titers_sera}"
        )

    # Validate "All" cohort exists in multicohort
    if "All" not in sera_multicohort["cohort"].values:
        raise ValueError(
            "Expected 'All' cohort in sera_multicohort but not found. "
            f"Available cohorts: {sera_multicohort['cohort'].unique().tolist()}"
        )

    # Get list of all cohorts (for legend), with "All" first
    all_cohorts = sera_multicohort["cohort"].unique().tolist()
    all_cohorts = ["All"] + sorted([c for c in all_cohorts if c != "All"])
    mo.output.append(mo.md(f"\nCohorts: {all_cohorts}"))

    # Aggregate cohorts into a list per serum and add to metadata
    cohorts_per_serum = sera_multicohort.groupby("serum")["cohort"].apply(list)
    metadata = metadata.merge(
        cohorts_per_serum.rename("cohorts").reset_index(),
        on="serum",
        how="left",
    )
    sera_missing_cohorts = metadata.loc[metadata["cohorts"].isna(), "serum"].tolist()
    if sera_missing_cohorts:
        raise ValueError(
            f"Some sera missing cohort assignments: {sera_missing_cohorts[:10]}"
            + (
                f" (and {len(sera_missing_cohorts) - 10} more)"
                if len(sera_missing_cohorts) > 10
                else ""
            )
        )

    # Read viruses
    viruses = pd.read_csv(viruses_csv)
    mo.output.append(mo.md(f"\nRead {len(viruses)=} viruses from {viruses_csv=}"))

    # Validate required columns in viruses
    required_virus_cols = {
        "virus",
        "subtype",
        "strain_type",
        "subclade",
        "derived_haplotype",
        "vaccine_type",
    }
    missing_virus_cols = required_virus_cols - set(viruses.columns)
    if missing_virus_cols:
        raise ValueError(f"viruses_csv missing required columns: {missing_virus_cols}")

    # Validate recent_vaccine_strains are in viruses
    missing_vaccine_strains = set(recent_vaccine_strains) - set(viruses["virus"])
    if missing_vaccine_strains:
        raise ValueError(
            f"recent_vaccine_strains not in viruses: {missing_vaccine_strains}"
        )

    # Validate strain_type values
    valid_strain_types = {circulating_strain_type, "vaccine"}
    invalid_strain_types = set(viruses["strain_type"]) - valid_strain_types
    if invalid_strain_types:
        raise ValueError(
            f"Invalid strain_type values: {invalid_strain_types}. "
            f"Expected: {valid_strain_types}"
        )

    # Mark recent vaccine strains as "recent_vaccine" strain_type
    viruses["strain_type"] = viruses["strain_type"].where(
        ~viruses["virus"].isin(recent_vaccine_strains), "recent_vaccine"
    )

    # Validate subtypes param matches virus data
    data_subtypes = set(viruses["subtype"].unique())
    param_subtypes = set(subtypes)
    if not param_subtypes.issubset(data_subtypes):
        raise ValueError(
            f"subtypes param {param_subtypes} not all in viruses data. "
            f"Available subtypes: {data_subtypes}"
        )

    # Compute strain plot order: sort by subclade then alphabetically
    # Strains without subclade (vaccine strains) sort to end
    viruses_sorted = viruses.sort_values(
        by=["subclade", "virus"],
        key=lambda col: col.fillna("zzz") if col.name == "subclade" else col,
    )
    viral_strain_plot_order = viruses_sorted["virus"].tolist()
    assert set(viruses["virus"]) == set(viral_strain_plot_order)
    mo.output.append(
        mo.md(f"\nComputed strain plot order: {len(viral_strain_plot_order)} viruses")
    )
    return all_cohorts, metadata, titers, viral_strain_plot_order, viruses


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Plot all the titers

    ### Assign label colors by subclade / vaccine type
    Define color mapping from subclade (circulating strains) or vaccine type (vaccine strains) to colors for label coloring, then create an expression that can be passed to `altair` *labelColor*:
    """
    )
    return


@app.cell
def _(
    circulating_strain_type,
    json,
    mo,
    pd,
    plot_titer_summaries_params,
    viral_strain_plot_order,
    viruses,
):
    strain_color_prop = viruses.assign(
        strain=lambda x: pd.Categorical(
            x["virus"], viral_strain_plot_order, ordered=True
        ),
        color_prop=lambda x: x["subclade"].where(
            x["strain_type"] == circulating_strain_type, x["vaccine_type"] + " vaccine"
        ),
    ).sort_values("strain")

    assert strain_color_prop["color_prop"].notnull().all()
    assert set(viruses["virus"]) == set(strain_color_prop["strain"])

    viruses["color_prop"] = viruses["virus"].map(
        strain_color_prop.set_index("strain")["color_prop"].to_dict()
    )

    prop_colors = dict(plot_titer_summaries_params["prop_colors"])

    other_prop_colors = plot_titer_summaries_params["other_prop_colors"]

    for _subtype in strain_color_prop["subtype"].unique():
        subtype_color_props = (
            strain_color_prop[strain_color_prop["subtype"] == _subtype]["color_prop"]
            .unique()
            .tolist()
        )
        props_not_yet_colored = [p for p in subtype_color_props if p not in prop_colors]
        if len(props_not_yet_colored) > len(other_prop_colors):
            raise ValueError(
                f"props_not_yet_colored={props_not_yet_colored!r} longer than "
                f"other_prop_colors={other_prop_colors!r}"
            )
        prop_colors.update(dict(zip(props_not_yet_colored, other_prop_colors)))

    mo.output.append(
        pd.Series(prop_colors).rename("color").rename_axis("property").to_frame()
    )
    assert set(strain_color_prop["color_prop"]).issubset(prop_colors)

    strain_color_prop = strain_color_prop.assign(
        color=lambda x: x["color_prop"].map(prop_colors)
    )

    color_mapping = strain_color_prop.set_index("virus")["color"].to_dict()

    # make a different color map for each subtype as they are plotted separately
    labelColor_expr = f"({json.dumps(color_mapping)})[datum.label] || 'black'"

    # Create label mapping: derived_haplotype if exists, else strain (without subtype suffix)
    label_mapping = {
        row["virus"]: (
            row["derived_haplotype"]
            if pd.notna(row["derived_haplotype"])
            else row["virus"].rsplit("_", 1)[0]  # remove _H1N1 or _H3N2 suffix
        )
        for _, row in viruses.iterrows()
    }
    labelText_expr = f"({json.dumps(label_mapping)})[datum.label] || datum.label"
    return labelColor_expr, labelText_expr, strain_color_prop


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ### Now make nicely formatted charts
    """
    )
    return


@app.cell
def _(
    all_cohorts,
    alt,
    facet_orientation,
    labelColor_expr,
    labelText_expr,
    metadata,
    mo,
    pd,
    plot_titer_summaries_params,
    titers,
    viral_strain_plot_order,
    viruses,
):
    # First set up the base chart and selections

    facet_size = plot_titer_summaries_params["facet_size"]
    if facet_orientation not in {"vertical", "horizontal"}:
        raise ValueError(
            f"facet_orientation must be 'vertical' or 'horizontal', got {facet_orientation!r}"
        )
    titer_encoding = "x" if facet_orientation == "vertical" else "y"

    # Validate no duplicate serum-virus pairs in titers
    duplicate_pairs = titers.groupby(["serum", "virus"]).size()
    duplicate_pairs = duplicate_pairs[duplicate_pairs > 1]
    if len(duplicate_pairs) > 0:
        raise ValueError(
            f"Found {len(duplicate_pairs)} duplicate serum-virus pairs in titers:\n"
            f"{duplicate_pairs.head(10).to_string()}"
        )

    # Validate each virus has unique subtype and strain_type
    virus_groups = viruses.groupby("virus")[["subtype", "strain_type"]].nunique()
    inconsistent = virus_groups[(virus_groups > 1).any(axis=1)]
    if len(inconsistent) > 0:
        raise ValueError(
            f"Found viruses with inconsistent subtype/strain_type:\n{inconsistent}"
        )

    mo.output.append(mo.md(f"Plotting cohorts={all_cohorts} with {facet_orientation=}"))

    virus_selection = alt.selection_point(
        fields=["virus"], on="mouseover", empty=False, clear="mouseout", nearest=False
    )

    serum_selection = alt.selection_point(
        fields=["serum"], on="mouseover", empty=False, clear="mouseout", nearest=False
    )

    cohort_selection = alt.selection_point(
        fields=["cohort"], bind="legend", empty="all", toggle="true", clear=False
    )

    # select by color used to color strain labels
    color_prop_selection = alt.selection_point(
        fields=["color_prop"], bind="legend", empty="all", toggle="true", clear=False
    )

    max_age = 5 * int(metadata["age_numeric"].max() // 5) + 5
    assert all(metadata["age_numeric"] <= max_age)
    min_age_slider = alt.param(
        value=0,
        bind=alt.binding_range(
            min=0, max=max_age, step=5, name="minimum subject age (years)"
        ),
    )
    max_age_slider = alt.param(
        value=max_age,
        bind=alt.binding_range(
            min=0, max=max_age, step=5, name="maximum subject age (years)"
        ),
    )

    # make the chart base, using transform_lookup to make it as small as possible
    # by looking up serum-specific and virus-specific annotations
    titers_base_nolookup = (
        alt.Chart(titers[["serum", "virus", "titer"]])
        .add_params(
            virus_selection,
            serum_selection,
            cohort_selection,
            color_prop_selection,
            min_age_slider,
            max_age_slider,
        )
        .encode(
            **{
                ("y" if facet_orientation == "vertical" else "x"): alt.Y(
                    "virus",
                    sort=list(reversed(viral_strain_plot_order)),
                    axis=alt.Axis(
                        labelLimit=500,
                        labelColor={"expr": labelColor_expr},
                        labelFontWeight=600,  # make a bit bolder so colors show
                        labelExpr=labelText_expr,
                    ),
                ),
            },
        )
        .properties(
            **(
                {"height": alt.Step(11), "width": facet_size}
                if facet_orientation == "vertical"
                else {"width": alt.Step(11), "height": facet_size}
            )
        )
    )

    # dummy chart to bind the selectable legend for serum cohort
    dummy_cohort_chart = (
        alt.Chart(pd.DataFrame({"cohort": all_cohorts}))
        .add_params(cohort_selection)
        .mark_point(opacity=0)
        .encode(
            fill=alt.Fill(
                "cohort",
                title="serum cohort (click to select)",
                scale=alt.Scale(domain=all_cohorts, range=["gray"]),
                legend=alt.Legend(
                    symbolStrokeColor="black", symbolOpacity=1, columns=6
                ),
            )
        )
        .properties(width=1, height=1)  # tiny plot; legend renders outside
    )

    # because of scoping issues when layering and faceting charts with
    # transform_lookups (faceting must be done before lookups), we add
    # this function to do the faceting and lookups
    def facet_and_add_lookups(chart):
        return (
            chart
            # facet
            .facet(
                {
                    (
                        "column" if facet_orientation == "vertical" else "row"
                    ): alt.Column("cohort_n:N", title=None)
                }
            )
            # lookup additional data
            .transform_lookup(
                lookup="serum",
                from_=alt.LookupData(
                    data=metadata,
                    key="serum",
                    fields=[
                        "cohorts",
                        "serum_collection_date",
                        "age",
                        "age_numeric",
                        "sex",
                    ],
                ),
            )
            .transform_lookup(
                lookup="virus",
                from_=alt.LookupData(
                    data=viruses,
                    key="virus",
                    fields=[
                        "subtype",
                        "strain_type",
                        "subclade",
                        "color_prop",
                        "derived_haplotype",
                    ],
                ),
            )
            # flatten cohorts list (from sera_multicohort) to one row per cohort
            .transform_flatten(["cohorts"], as_=["cohort"])
            # filter by property used to color strain labels
            .transform_filter(color_prop_selection)
            # filter by cohort and age
            .transform_filter(cohort_selection)
            .transform_filter(alt.datum["age_numeric"] >= min_age_slider)
            .transform_filter(alt.datum["age_numeric"] <= max_age_slider)
            # make facet labels w n per cohort
            .transform_joinaggregate(n_per_cohort="distinct(serum)", groupby=["cohort"])
            .transform_calculate(
                cohort_n="datum.cohort + ' (n=' + datum.n_per_cohort + ')'"
            )
        )

    return (
        color_prop_selection,
        dummy_cohort_chart,
        facet_and_add_lookups,
        serum_selection,
        titer_encoding,
        titers_base_nolookup,
        virus_selection,
    )


@app.cell
def _(alt, mo, plot_titer_summaries_params):
    # set titer scale
    titer_lower_limit = plot_titer_summaries_params["titer_lower_limit"]
    mo.output.append(mo.md(f"Using {titer_lower_limit=}"))
    titer_scale = alt.Scale(
        type="log", nice=False, domainMin=titer_lower_limit, padding=4
    )
    return titer_lower_limit, titer_scale


@app.cell
def _(alt, titer_encoding, titer_scale, titers_base_nolookup, virus_selection):
    # make median titer point chart
    median_points = (
        titers_base_nolookup.transform_aggregate(
            median_titer="median(titer)",
            groupby=[
                "virus",
                "subtype",
                "derived_haplotype",
                "strain_type",
                "subclade",
                "cohort",
            ],
        )
        .encode(
            **{
                titer_encoding: alt.X(
                    "median_titer:Q", title="titer", scale=titer_scale
                )
            },
            tooltip=[
                "virus",
                "derived_haplotype:N",
                alt.Tooltip("median_titer:Q", format=".1f"),
                "strain_type:N",
                "subclade:N",
            ],
            color=alt.condition(virus_selection, alt.value("red"), alt.value("black")),
            size=alt.condition(virus_selection, alt.value(80), alt.value(40)),
        )
        .mark_circle(opacity=1)
    )

    # facet_and_add_lookups(median_points)
    return (median_points,)


@app.cell
def _(alt, serum_selection, titer_encoding, titer_scale, titers_base_nolookup):
    # make per-serum lines
    serum_lines = titers_base_nolookup.encode(
        **{titer_encoding: alt.X("titer", scale=titer_scale)},
        detail=alt.Detail("serum"),
        tooltip=[
            "virus",
            "derived_haplotype:N",
            "serum",
            alt.Tooltip("titer", format=".1f"),
            alt.Tooltip("serum_collection_date:N", title="serum date"),
            alt.Tooltip("age:N", title="age"),
            "sex:N",
        ],
        size=alt.condition(serum_selection, alt.value(3), alt.value(1.5)),
        opacity=alt.condition(serum_selection, alt.value(1), alt.value(0.2)),
    ).mark_line()

    # facet_and_add_lookups(serum_lines + median_points)
    return (serum_lines,)


@app.cell
def _(alt, titer_encoding, titer_scale, titers_base_nolookup):
    # make interquartile range chart

    interquartile_range = (
        titers_base_nolookup.transform_joinaggregate(
            median_titer="median(titer)",
            titer_q1="q1(titer)",
            titer_q3="q3(titer)",
            groupby=["virus"],
        )
        .encode(
            **{titer_encoding: alt.X("titer", scale=titer_scale)},
            tooltip=[
                "virus",
                "derived_haplotype:N",
                alt.Tooltip("median_titer:Q", format=".1f"),
                alt.Tooltip("titer_q1:Q", format=".1f"),
                alt.Tooltip("titer_q3:Q", format=".1f"),
                "strain_type:N",
                "subclade:N",
            ],
        )
        .mark_errorband(extent="iqr", opacity=0.5, interpolate="linear")
    )

    # facet_and_add_lookups(interquartile_range + median_points)
    return (interquartile_range,)


@app.cell
def _(
    alt,
    mo,
    plot_titer_summaries_params,
    titer_encoding,
    titer_lower_limit,
    titers_base_nolookup,
    virus_selection,
):
    # make fraction below titer cutoff chart

    titer_cutoff = plot_titer_summaries_params["titer_cutoff"]
    mo.output.append(mo.md(f"Setting initial {titer_cutoff=}"))

    titer_cutoff_slider = alt.param(
        value=titer_cutoff,
        bind=alt.binding_range(
            min=titer_lower_limit,
            max=1000,
            step=5,
            name="fraction sera below this cutoff",
        ),
    )

    # make titer cutoff chart
    frac_below_cutoff = (
        titers_base_nolookup.add_params(titer_cutoff_slider)
        .transform_calculate(below_cutoff=alt.datum["titer"] < titer_cutoff_slider)
        .transform_aggregate(
            n_below_cutoff="sum(below_cutoff)",
            n_total="distinct(serum)",
            groupby=[
                "virus",
                "subtype",
                "derived_haplotype",
                "strain_type",
                "subclade",
                "cohort",
            ],
        )
        .transform_calculate(
            frac_below_cutoff=alt.datum["n_below_cutoff"] / alt.datum["n_total"]
        )
        .encode(
            **{
                titer_encoding: alt.X(
                    "frac_below_cutoff:Q", title="fraction below cutoff"
                )
            },
            tooltip=[
                "virus",
                "derived_haplotype:N",
                alt.Tooltip("frac_below_cutoff:Q", format=".2f"),
                "strain_type:N",
                "subclade:N",
            ],
            color=alt.condition(virus_selection, alt.value("red"), alt.value("black")),
        )
        .mark_bar(opacity=0.8)
    )

    # facet_and_add_lookups(frac_below_cutoff)
    return (frac_below_cutoff,)


@app.cell
def _(
    alt,
    chart_htmls,
    circulating_strain_type,
    color_prop_selection,
    dummy_cohort_chart,
    facet_and_add_lookups,
    facet_orientation,
    frac_below_cutoff,
    interquartile_range,
    itertools,
    median_points,
    mo,
    serum_lines,
    strain_color_prop,
    subtypes,
):
    made_chart = {c: False for c in chart_htmls}

    for _subtype, strain_type, (chart_obj, chart_desc, title) in itertools.product(
        subtypes,
        ["recent", "vaccine"],
        [
            (
                serum_lines + median_points,
                "individual_sera",
                "median (points) and per-serum (lines) titers",
            ),
            (
                interquartile_range + median_points,
                "interquartile_range",
                "median (points) and interquartile range titers",
            ),
            (
                frac_below_cutoff,
                "frac_below_cutoff",
                "fraction sera below titer cutoff",
            ),
        ],
    ):
        filepattern = f"{_subtype}_{strain_type}_{chart_desc}"
        filename = [c for c in chart_htmls if filepattern in c]
        assert (
            len(filename) == 1
        ), f"did not find one filepattern={filepattern!r} in chart_htmls={chart_htmls!r}"
        filename = filename[0]

        # strain types to plot
        strain_types = {
            "recent": [circulating_strain_type, "recent_vaccine"],
            "vaccine": ["vaccine", "recent_vaccine"],
        }[strain_type]

        # ---- Make the legend for the colored strain labels ------------------------
        # get the virus colors plotted for the labels
        plotted_colors = strain_color_prop[
            (strain_color_prop["subtype"] == _subtype)
            & (strain_color_prop["strain_type"].isin(strain_types))
        ][["color_prop", "color"]].drop_duplicates()

        label_color_legend = (
            alt.Chart(plotted_colors)
            .add_params(color_prop_selection)
            .mark_point(opacity=0)  # invisible mark; we just want the legend
            .encode(
                fill=alt.Fill(
                    "color_prop",
                    title="virus type (click to select)",
                    scale=alt.Scale(
                        domain=list(reversed(plotted_colors["color_prop"].tolist())),
                        range=list(reversed(plotted_colors["color"].tolist())),
                    ),
                    legend=alt.Legend(symbolType="square"),
                )
            )
            .properties(width=1, height=1)  # tiny plot; legend renders outside
        )
        # ---- Finished making the legend for the colored strain labels -------------

        chart = (
            alt.vconcat(
                (
                    facet_and_add_lookups(chart_obj)
                    .transform_filter(alt.datum["subtype"] == _subtype)
                    .transform_filter(
                        alt.FieldOneOfPredicate("strain_type", strain_types)
                    )
                ),
                label_color_legend,
                dummy_cohort_chart,
                spacing=1,
            )
            .resolve_scale(fill="independent")
            .configure_axis(
                grid=False,
                titleFontWeight="normal",
                titleFontSize=13,
                labelOverlap=True,
            )
            .configure_header(
                title=None,
                labelOrient="top" if facet_orientation == "vertical" else "right",
                labelFontSize=13,
                labelPadding=2,
            )
            .configure_view(stroke="black")
            .configure_facet(spacing=8)
            .configure_legend(
                labelFontSize=12,
                titleFontSize=13,
                symbolStrokeWidth=1,
                symbolOpacity=1,
                symbolStrokeColor="black",
                columns=12,
                orient="bottom",
            )
            .properties(
                title=alt.TitleParams(
                    f"{title} for {_subtype} {strain_type} strains",
                    anchor="middle",
                    fontSize=13,
                )
            )
        )
        if not any(made_chart.values()):
            mo.output.append(
                mo.md("Displaying just the first chart here (since they are large).")
            )
            mo.output.append(chart)

        mo.output.append(mo.md(f"Saving to filename={filename!r}\n"))
        chart.save(filename)

        made_chart[filename] = True

    assert all(made_chart.values()), f"made_chart={made_chart!r}"
    return


if __name__ == "__main__":
    app.run()
