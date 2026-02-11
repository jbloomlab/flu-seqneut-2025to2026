import marimo

__generated_with = "0.17.6"
app = marimo.App(width="full")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    # Fold-change titer plots

    Interactive Altair plots comparing titers between paired sera cohorts
    (e.g., pre- vs post-vaccination), with fold-change panels.
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

    import argparse
    import os
    import pathlib
    import pickle
    import sys

    import marimo as mo

    from_cmdline = "--context-pickle" in sys.argv

    if from_cmdline:
        print("Loading context from command-line argument")
        p = argparse.ArgumentParser()
        p.add_argument("--context-pickle", required=True)
        args = p.parse_args()
        context_pickle_path = pathlib.Path(args.context_pickle)
    else:
        print("Running in marimo edit mode")
        context_pickle_path = None
        context_pickle_path = pathlib.Path(
            "results/titer_plots/PENN_pre_post_vax_plot_fold_changes_vertical_context.pickle"
        )

    if context_pickle_path and context_pickle_path.exists():
        print(f"Reading context from {context_pickle_path}")
        with open(context_pickle_path, "rb") as f_context:
            context = pickle.load(f_context)

        context_workdir = context["workdir"]
        current_workdir = os.getcwd()

        if from_cmdline:
            if context_workdir != current_workdir:
                raise RuntimeError(
                    f"Context workdir mismatch!\n"
                    f"  Context was created in: {context_workdir}\n"
                    f"  Currently running in:   {current_workdir}\n"
                    f"This should not happen when running via Snakemake."
                )
            print(f"Verified working directory: {current_workdir}")
        else:
            if context_workdir and context_workdir != current_workdir:
                print(f"Changing directory from {current_workdir} to {context_workdir}")
                os.chdir(context_workdir)
            elif context_workdir:
                print(f"Already in correct working directory: {context_workdir}")
    else:
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
    ## Read data and pair sera across cohorts
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
    # Extract variables from context
    titers_csv = context["input"]["titers_csv"]
    sera_multicohort_csv = context["input"]["sera_multicohort_csv"]
    viruses_csv = context["input"]["viruses_csv"]
    recent_vaccine_strains = context["params"]["recent_vaccine_strains"]
    circulating_strain_type = context["params"]["circulating_strain_type"]
    plot_titer_summaries_params = context["params"]["plot_titer_summaries_params"]
    subtypes = context["params"]["subtypes"]
    fold_change_config = context["params"]["fold_change_config"]
    facet_orientation = context["params"]["facet_orientation"]
    chart_htmls = context["output"]["chart_htmls"]
    return (
        chart_htmls,
        circulating_strain_type,
        facet_orientation,
        fold_change_config,
        plot_titer_summaries_params,
        recent_vaccine_strains,
        sera_multicohort_csv,
        subtypes,
        titers_csv,
        viruses_csv,
    )


@app.cell
def _(fold_change_config, mo, pd, sera_multicohort_csv, titers_csv):
    # Read titers
    titers = pd.read_csv(titers_csv)
    mo.output.append(mo.md(f"Read {len(titers)=} titers from {titers_csv=}"))

    required_titer_cols = {"serum", "virus", "titer"}
    missing_titer_cols = required_titer_cols - set(titers.columns)
    if missing_titer_cols:
        raise ValueError(f"titers_csv missing required columns: {missing_titer_cols}")

    # Read sera multicohort
    sera_mc = pd.read_csv(sera_multicohort_csv)
    mo.output.append(mo.md(f"\nRead {len(sera_mc)=} rows from {sera_multicohort_csv=}"))

    # Extract fold-change config
    cohorts_config = fold_change_config["cohorts"]
    condition_colors = fold_change_config["condition_colors"]
    fold_change_title = fold_change_config["title"]

    if len(condition_colors) != len(cohorts_config):
        raise ValueError(
            f"condition_colors has {len(condition_colors)} entries but "
            f"cohorts has {len(cohorts_config)} entries; must match"
        )

    mo.output.append(mo.md(f"\nCohorts config: {cohorts_config}"))
    mo.output.append(mo.md(f"Condition colors: {condition_colors}"))

    # Validate all specified cohorts exist
    available_cohorts = set(sera_mc["cohort"].unique())
    for cohort_name in cohorts_config:
        if cohort_name not in available_cohorts:
            raise ValueError(
                f"Cohort {cohort_name!r} not found in sera_multicohort. "
                f"Available cohorts: {sorted(available_cohorts)}"
            )

    # Filter to specified cohorts with non-null subject_id
    sera_filtered = sera_mc[
        sera_mc["cohort"].isin(cohorts_config)
        & sera_mc["subject_id"].notna()
        & (sera_mc["subject_id"] != "")
    ].copy()
    mo.output.append(
        mo.md(
            f"\nFiltered to {len(sera_filtered)} rows with specified cohorts "
            f"and non-null subject_id"
        )
    )

    # Find subject_ids present in ALL specified cohorts
    subjects_per_cohort = sera_filtered.groupby("cohort")["subject_id"].apply(set)
    paired_subjects = set.intersection(*subjects_per_cohort.values)
    if not paired_subjects:
        raise ValueError(
            "No subject_ids found in all specified cohorts. "
            "Subjects per cohort: "
            + ", ".join(
                f"{c}: {len(s)} subjects" for c, s in subjects_per_cohort.items()
            )
        )
    mo.output.append(mo.md(f"\nFound {len(paired_subjects)} paired subjects"))

    if "all subjects" in paired_subjects:
        raise ValueError(
            "Found subject_id 'all subjects' in data; this conflicts with "
            "the dropdown label used in individual_sera plots"
        )

    # Filter to paired subjects and add condition labels
    sera_paired = sera_filtered[
        sera_filtered["subject_id"].isin(paired_subjects)
    ].copy()
    sera_paired["condition"] = sera_paired["cohort"].map(cohorts_config)

    # Validate each subject has exactly one serum per cohort
    serum_counts = sera_paired.groupby(["subject_id", "cohort"])["serum"].nunique()
    multi_sera = serum_counts[serum_counts > 1]
    if len(multi_sera) > 0:
        raise ValueError(f"Some subjects have multiple sera in a cohort:\n{multi_sera}")

    # Join with titers
    paired_titers = sera_paired[
        [
            "serum",
            "subject_id",
            "condition",
            "cohort",
            "age",
            "sex",
            "age_numeric",
            "serum_collection_date",
        ]
    ].merge(titers[["serum", "virus", "titer"]], on="serum")
    mo.output.append(mo.md(f"\nPaired titers: {len(paired_titers)} rows"))

    # Compute fold change relative to baseline (first cohort)
    baseline_cohort = list(cohorts_config.keys())[0]
    baseline_label = cohorts_config[baseline_cohort]
    baseline = paired_titers[paired_titers["cohort"] == baseline_cohort][
        ["subject_id", "virus", "titer"]
    ].rename(columns={"titer": "baseline_titer"})

    paired_titers = paired_titers.merge(baseline, on=["subject_id", "virus"])
    paired_titers["fold_change"] = (
        paired_titers["titer"] / paired_titers["baseline_titer"]
    )
    mo.output.append(
        mo.md(
            f"\nComputed fold changes relative to baseline: "
            f"{baseline_label!r} ({baseline_cohort})"
        )
    )

    return (
        baseline_label,
        cohorts_config,
        condition_colors,
        fold_change_title,
        paired_titers,
        titers,
    )


@app.cell
def _(
    circulating_strain_type,
    mo,
    paired_titers,
    pd,
    recent_vaccine_strains,
    subtypes,
    viruses_csv,
):
    # Read viruses and process same as plot_titer_summaries
    viruses = pd.read_csv(viruses_csv)
    mo.output.append(mo.md(f"Read {len(viruses)=} viruses from {viruses_csv=}"))

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

    # Mark recent vaccine strains
    viruses["strain_type"] = viruses["strain_type"].where(
        ~viruses["virus"].isin(recent_vaccine_strains), "recent_vaccine"
    )

    # Validate subtypes
    data_subtypes = set(viruses["subtype"].unique())
    param_subtypes = set(subtypes)
    if not param_subtypes.issubset(data_subtypes):
        raise ValueError(
            f"subtypes param {param_subtypes} not all in viruses data. "
            f"Available subtypes: {data_subtypes}"
        )

    # Compute strain plot order
    viruses_sorted = viruses.sort_values(
        by=["subclade", "virus"],
        key=lambda col: col.fillna("zzz") if col.name == "subclade" else col,
    )
    viral_strain_plot_order = viruses_sorted["virus"].tolist()
    assert set(viruses["virus"]) == set(viral_strain_plot_order)

    # Merge virus metadata into paired_titers
    virus_cols = [
        "virus",
        "subtype",
        "strain_type",
        "subclade",
        "derived_haplotype",
        "vaccine_type",
    ]
    paired_titers_full = paired_titers.merge(
        viruses[virus_cols], on="virus", how="left"
    )

    # Validate all viruses have metadata
    missing_virus_meta = paired_titers_full["subtype"].isna().sum()
    if missing_virus_meta > 0:
        missing = set(paired_titers["virus"]) - set(viruses["virus"])
        raise ValueError(f"Viruses in titers but not in viruses_csv: {missing}")

    mo.output.append(
        mo.md(f"\nComputed strain plot order: {len(viral_strain_plot_order)} viruses")
    )
    return paired_titers_full, viral_strain_plot_order, viruses


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Build charts

    ### Assign label colors by subclade / vaccine type
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
    labelColor_expr = f"({json.dumps(color_mapping)})[datum.label] || 'black'"

    label_mapping = {
        row["virus"]: (
            row["derived_haplotype"]
            if pd.notna(row["derived_haplotype"])
            else row["virus"].rsplit("_", 1)[0]
        )
        for _, row in viruses.iterrows()
    }
    labelText_expr = f"({json.dumps(label_mapping)})[datum.label] || datum.label"
    return labelColor_expr, labelText_expr, strain_color_prop


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ### Make charts with titer + fold-change panels
    """
    )
    return


@app.cell
def _(
    alt,
    baseline_label,
    chart_htmls,
    circulating_strain_type,
    cohorts_config,
    condition_colors,
    facet_orientation,
    fold_change_title,
    itertools,
    labelColor_expr,
    labelText_expr,
    mo,
    paired_titers_full,
    pd,
    plot_titer_summaries_params,
    strain_color_prop,
    subtypes,
    viral_strain_plot_order,
    viruses,
):
    facet_size = plot_titer_summaries_params["facet_size"]
    titer_lower_limit = plot_titer_summaries_params["titer_lower_limit"]

    if facet_orientation not in {"vertical", "horizontal"}:
        raise ValueError(
            f"facet_orientation must be 'vertical' or 'horizontal', "
            f"got {facet_orientation!r}"
        )

    # Condition color scale
    condition_domain = list(cohorts_config.values())
    condition_range = list(condition_colors)

    # Titer scale (log)
    titer_scale = alt.Scale(
        type="log", nice=False, domainMin=titer_lower_limit, padding=4
    )

    # Fold-change scale (log)
    fold_change_scale = alt.Scale(type="log", nice=False, padding=4)

    # Shared virus axis encoding
    virus_sort = list(reversed(viral_strain_plot_order))
    virus_axis = alt.Axis(
        labelLimit=500,
        labelColor={"expr": labelColor_expr},
        labelFontWeight=600,
        labelExpr=labelText_expr,
    )

    # Mouseover selections for cross-panel highlighting (added to concatenated chart)
    virus_selection = alt.selection_point(
        fields=["virus"], on="mouseover", empty=False, clear="mouseout", nearest=False
    )
    subject_selection = alt.selection_point(
        fields=["subject_id"],
        on="mouseover",
        empty=False,
        clear="mouseout",
        nearest=False,
    )
    color_prop_selection = alt.selection_point(
        fields=["color_prop"], bind="legend", empty="all", toggle="true", clear=False
    )

    # Age sliders
    max_age = 5 * int(paired_titers_full["age_numeric"].max() // 5) + 5
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

    # Condition legend (clickable)
    condition_selection = alt.selection_point(
        fields=["condition"], bind="legend", empty="all", toggle="true", clear=False
    )

    # --- Build and save charts ---
    made_chart = {c: False for c in chart_htmls}

    for _subtype, strain_type, (chart_desc, chart_title) in itertools.product(
        subtypes,
        ["recent", "vaccine"],
        [
            ("individual_sera", "per-subject (lines) and median (points)"),
            ("interquartile_range", "median (points) and interquartile range"),
        ],
    ):
        filepattern = f"{_subtype}_{strain_type}_{chart_desc}"
        filename = [c for c in chart_htmls if filepattern in c]
        assert (
            len(filename) == 1
        ), f"did not find one filepattern={filepattern!r} in chart_htmls={chart_htmls!r}"
        filename = filename[0]

        strain_types = {
            "recent": [circulating_strain_type, "recent_vaccine"],
            "vaccine": ["vaccine", "recent_vaccine"],
        }[strain_type]

        # Filter data for this subtype and strain_type
        chart_data = paired_titers_full[
            (paired_titers_full["subtype"] == _subtype)
            & (paired_titers_full["strain_type"].isin(strain_types))
        ].copy()

        # Minimal chart data for Altair (reduces serialized spec size via transform_lookup)
        chart_data_minimal = chart_data[
            ["serum", "virus", "titer", "baseline_titer", "fold_change"]
        ]
        chart_viruses = viruses[viruses["virus"].isin(chart_data["virus"].unique())][
            [
                "virus",
                "strain_type",
                "subclade",
                "derived_haplotype",
                "vaccine_type",
                "color_prop",
            ]
        ]
        chart_sera = chart_data[
            [
                "serum",
                "subject_id",
                "condition",
                "age",
                "age_numeric",
                "sex",
                "serum_collection_date",
            ]
        ].drop_duplicates(subset=["serum"])

        if len(chart_data) == 0:
            mo.output.append(
                mo.md(f"**Skipping** {_subtype} {strain_type}: no data after filtering")
            )
            chart_data_empty = (
                alt.Chart(pd.DataFrame({"x": [0], "y": [0], "text": ["No data"]}))
                .mark_text()
                .encode(x="x:Q", y="y:Q", text="text:N")
            )
            chart_data_empty.save(filename)
            made_chart[filename] = True
            continue

        # Virus axis encoding kwargs (shared between panels)
        if facet_orientation == "vertical":
            virus_enc_key = "y"
            titer_enc_key = "x"
        else:
            virus_enc_key = "x"
            titer_enc_key = "y"

        virus_enc_labeled = {
            virus_enc_key: alt.Y("virus:N", sort=virus_sort, axis=virus_axis)
        }
        virus_enc_no_labels = {
            virus_enc_key: alt.Y(
                "virus:N",
                sort=virus_sort,
                axis=alt.Axis(labels=False, ticks=False, title=None),
            )
        }
        step_props = (
            {"height": alt.Step(11), "width": facet_size}
            if facet_orientation == "vertical"
            else {"width": alt.Step(11), "height": facet_size}
        )

        # vertical (hconcat): titer left with labels, fc right without
        # horizontal (vconcat): titer top without labels, fc bottom with
        titer_virus_enc = (
            virus_enc_labeled
            if facet_orientation == "vertical"
            else virus_enc_no_labels
        )
        fc_virus_enc = (
            virus_enc_no_labels
            if facet_orientation == "vertical"
            else virus_enc_labeled
        )

        # Color encoding for conditions
        condition_color_enc = alt.Color(
            "condition:N",
            title="condition",
            scale=alt.Scale(domain=condition_domain, range=condition_range),
            legend=alt.Legend(orient="bottom", columns=6),
        )

        # Base chart for titer panel
        titer_base = (
            alt.Chart(chart_data_minimal)
            .transform_lookup(
                lookup="virus",
                from_=alt.LookupData(
                    data=chart_viruses,
                    key="virus",
                    fields=[
                        "strain_type",
                        "subclade",
                        "derived_haplotype",
                        "vaccine_type",
                        "color_prop",
                    ],
                ),
            )
            .transform_lookup(
                lookup="serum",
                from_=alt.LookupData(
                    data=chart_sera,
                    key="serum",
                    fields=[
                        "subject_id",
                        "condition",
                        "age",
                        "age_numeric",
                        "sex",
                        "serum_collection_date",
                    ],
                ),
            )
            .add_params(
                color_prop_selection,
                condition_selection,
                min_age_slider,
                max_age_slider,
            )
            .transform_filter(color_prop_selection)
            .transform_filter(condition_selection)
            .transform_filter(alt.datum["age_numeric"] >= min_age_slider)
            .transform_filter(alt.datum["age_numeric"] <= max_age_slider)
            .encode(**titer_virus_enc)
            .properties(**step_props)
        )

        # Base chart for fold-change panel
        fc_base = (
            alt.Chart(chart_data_minimal)
            .transform_lookup(
                lookup="virus",
                from_=alt.LookupData(
                    data=chart_viruses,
                    key="virus",
                    fields=[
                        "strain_type",
                        "subclade",
                        "derived_haplotype",
                        "vaccine_type",
                        "color_prop",
                    ],
                ),
            )
            .transform_lookup(
                lookup="serum",
                from_=alt.LookupData(
                    data=chart_sera,
                    key="serum",
                    fields=[
                        "subject_id",
                        "condition",
                        "age",
                        "age_numeric",
                        "sex",
                        "serum_collection_date",
                    ],
                ),
            )
            .transform_filter(f"datum.condition !== '{baseline_label}'")
            .transform_filter(color_prop_selection)
            .transform_filter(condition_selection)
            .transform_filter(alt.datum["age_numeric"] >= min_age_slider)
            .transform_filter(alt.datum["age_numeric"] <= max_age_slider)
            .encode(**fc_virus_enc)
            .properties(**step_props)
        )

        # Subject dropdown filter (only active for individual_sera charts)
        if chart_desc == "individual_sera":
            chart_subject_ids = sorted(chart_data["subject_id"].unique())
            subject_dropdown = alt.param(
                name="subject_dropdown",
                value="all subjects",
                bind=alt.binding_select(
                    options=["all subjects"] + chart_subject_ids,
                    name="subject ",
                ),
            )
            subject_filter_expr = (
                "subject_dropdown == 'all subjects' "
                "|| datum.subject_id == subject_dropdown"
            )
            titer_base_filtered = titer_base.transform_filter(subject_filter_expr)
            fc_base_filtered = fc_base.transform_filter(subject_filter_expr)
        else:
            subject_dropdown = None
            titer_base_filtered = titer_base
            fc_base_filtered = fc_base

        # --- Titer panel median points ---
        titer_median = (
            titer_base_filtered.transform_aggregate(
                median_titer="median(titer)",
                groupby=[
                    "virus",
                    "condition",
                    "derived_haplotype",
                    "subclade",
                    "strain_type",
                ],
            )
            .encode(
                **{
                    titer_enc_key: alt.X(
                        "median_titer:Q", title="titer", scale=titer_scale
                    )
                },
                color=condition_color_enc,
                tooltip=[
                    "virus:N",
                    alt.Tooltip("derived_haplotype:N"),
                    "subclade:N",
                    "strain_type:N",
                    "condition:N",
                    alt.Tooltip("median_titer:Q", format=".1f", title="median titer"),
                ],
                size=alt.condition(virus_selection, alt.value(80), alt.value(40)),
            )
            .mark_circle(opacity=1)
        )

        # --- Fold-change panel median points ---
        fc_median = (
            fc_base_filtered.transform_aggregate(
                median_fold_change="median(fold_change)",
                groupby=[
                    "virus",
                    "condition",
                    "derived_haplotype",
                    "subclade",
                    "strain_type",
                ],
            )
            .encode(
                **{
                    titer_enc_key: alt.X(
                        "median_fold_change:Q",
                        title="fold change in titer",
                        scale=fold_change_scale,
                    )
                },
                color=condition_color_enc,
                tooltip=[
                    "virus:N",
                    alt.Tooltip("derived_haplotype:N"),
                    "subclade:N",
                    "strain_type:N",
                    "condition:N",
                    alt.Tooltip(
                        "median_fold_change:Q", format=".2f", title="median fold change"
                    ),
                ],
                size=alt.condition(virus_selection, alt.value(80), alt.value(40)),
            )
            .mark_circle(opacity=1)
        )

        # --- Fold-change reference line at 1.0 ---
        if facet_orientation == "vertical":
            fc_ref_line = (
                alt.Chart(pd.DataFrame({"val": [1.0]}))
                .mark_rule(strokeDash=[4, 4], color="gray")
                .encode(x="val:Q")
            )
        else:
            fc_ref_line = (
                alt.Chart(pd.DataFrame({"val": [1.0]}))
                .mark_rule(strokeDash=[4, 4], color="gray")
                .encode(y="val:Q")
            )

        if chart_desc == "individual_sera":
            # --- Titer panel: individual lines ---
            titer_lines = titer_base_filtered.encode(
                **{titer_enc_key: alt.X("titer:Q", scale=titer_scale)},
                detail="subject_id:N",
                color=condition_color_enc,
                tooltip=[
                    "virus:N",
                    alt.Tooltip("derived_haplotype:N"),
                    "subclade:N",
                    "strain_type:N",
                    "subject_id:N",
                    "condition:N",
                    "serum:N",
                    alt.Tooltip("titer:Q", format=".1f"),
                    alt.Tooltip("serum_collection_date:N", title="serum date"),
                    alt.Tooltip("age:N"),
                    "sex:N",
                ],
                size=alt.condition(subject_selection, alt.value(3), alt.value(1.5)),
                opacity=alt.condition(subject_selection, alt.value(1), alt.value(0.2)),
            ).mark_line()

            titer_panel = titer_lines + titer_median

            # --- Fold-change panel: individual lines ---
            fc_lines = fc_base_filtered.encode(
                **{
                    titer_enc_key: alt.X(
                        "fold_change:Q",
                        title="fold change in titer",
                        scale=fold_change_scale,
                    )
                },
                detail="subject_id:N",
                color=condition_color_enc,
                tooltip=[
                    "virus:N",
                    alt.Tooltip("derived_haplotype:N"),
                    "subclade:N",
                    "strain_type:N",
                    "subject_id:N",
                    "condition:N",
                    alt.Tooltip("fold_change:Q", format=".2f"),
                    alt.Tooltip("titer:Q", format=".1f"),
                    alt.Tooltip("baseline_titer:Q", format=".1f"),
                ],
                size=alt.condition(subject_selection, alt.value(3), alt.value(1.5)),
                opacity=alt.condition(subject_selection, alt.value(1), alt.value(0.2)),
            ).mark_line()

            fc_panel = fc_ref_line + fc_lines + fc_median

        else:
            # interquartile_range
            # --- Titer panel: IQR bands ---
            titer_iqr = titer_base.encode(
                **{titer_enc_key: alt.X("titer:Q", scale=titer_scale)},
                color=condition_color_enc,
                tooltip=[
                    "virus:N",
                    alt.Tooltip("derived_haplotype:N"),
                    "subclade:N",
                    "strain_type:N",
                    "condition:N",
                ],
            ).mark_errorband(extent="iqr", opacity=0.3, interpolate="linear")

            titer_panel = titer_iqr + titer_median

            # --- Fold-change panel: IQR bands ---
            fc_iqr = fc_base.encode(
                **{
                    titer_enc_key: alt.X(
                        "fold_change:Q",
                        title="fold change in titer",
                        scale=fold_change_scale,
                    )
                },
                color=condition_color_enc,
                tooltip=[
                    "virus:N",
                    alt.Tooltip("derived_haplotype:N"),
                    "subclade:N",
                    "strain_type:N",
                    "condition:N",
                ],
            ).mark_errorband(extent="iqr", opacity=0.3, interpolate="linear")

            fc_panel = fc_ref_line + fc_iqr + fc_median

        # --- Virus type color legend ---
        plotted_colors = strain_color_prop[
            (strain_color_prop["subtype"] == _subtype)
            & (strain_color_prop["strain_type"].isin(strain_types))
        ][["color_prop", "color"]].drop_duplicates()

        label_color_legend = (
            alt.Chart(plotted_colors)
            .add_params(color_prop_selection)
            .mark_point(opacity=0)
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
            .properties(width=1, height=1)
        )

        # --- Assemble full chart ---
        # Concatenate titer and fold-change panels sharing the virus axis.
        # vertical orientation: viruses on Y → hconcat (side-by-side) shares Y axis
        # horizontal orientation: viruses on X → vconcat (stacked) shares X axis
        # Mouseover selections added here so they work bidirectionally across panels.
        concat_params = [virus_selection, subject_selection]
        if subject_dropdown is not None:
            concat_params.append(subject_dropdown)

        if facet_orientation == "vertical":
            panels_concat = alt.hconcat(titer_panel, fc_panel, spacing=5).add_params(
                *concat_params
            )
        else:
            panels_concat = alt.vconcat(titer_panel, fc_panel, spacing=5).add_params(
                *concat_params
            )

        chart = (
            alt.vconcat(
                panels_concat,
                label_color_legend,
                spacing=5,
            )
            .resolve_scale(
                color="shared",
                fill="independent",
            )
            .configure_axis(
                grid=False,
                titleFontWeight="normal",
                titleFontSize=13,
                labelOverlap=True,
            )
            .configure_view(stroke="black")
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
                    f"{fold_change_title}: {_subtype} {strain_type} strains",
                    subtitle=chart_title,
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
