"""Process and QC final titer data for each serum group.

Reads aggregated titers, sera metadata, and viral library, performs validation
and filtering, and outputs cleaned datasets for downstream analysis.
"""

import sys

import numpy as np
import pandas as pd


# =============================================================================
# Setup logging - redirect stdout/stderr to log file
# =============================================================================
sys.stderr = sys.stdout = log = open(snakemake.log[0], "w")

# Summary messages are written to both log and summary output file
summary_messages = []


def log_message(msg):
    """Log message to stdout and accumulate for summary file."""
    print(msg)
    summary_messages.append(msg)


# =============================================================================
# Read snakemake context
# =============================================================================
sera_metadata_csv = snakemake.input.sera_metadata
viral_library_csv = snakemake.input.viral_library
titers_csv = snakemake.input.titers

output_titers_csv = snakemake.output.titers
output_sera_csv = snakemake.output.sera
output_sera_multicohort_csv = snakemake.output.sera_multicohort
output_viruses_csv = snakemake.output.viruses
output_titers_summarized_csv = snakemake.output.titers_summarized
output_summary_txt = snakemake.output.summary

group = snakemake.wildcards.group
config = snakemake.params.config

# Extract config parameters
min_frac_viruses = config["min_frac_viruses"]
min_frac_sera = config["min_frac_sera"]
min_frac_action = config["min_frac_action"]
sera_to_drop = config["sera_to_drop"]
viruses_to_drop = config["viruses_to_drop"]
titer_cutoffs = config["titer_cutoffs"]

# Validate min_frac_action
if min_frac_action not in ("raise", "drop"):
    raise ValueError(
        f"min_frac_action must be 'raise' or 'drop', got: {min_frac_action!r}"
    )

log_message("=" * 70)
log_message(f"Processing final titer data for group: {group}")
log_message("=" * 70)
log_message("")

# =============================================================================
# Load input data
# =============================================================================
log_message("Loading input data...")

titers_df = pd.read_csv(titers_csv)
log_message(f"  Loaded {len(titers_df)} titer measurements from {titers_csv}")

sera_metadata_df = pd.read_csv(sera_metadata_csv)
log_message(f"  Loaded {len(sera_metadata_df)} sera from {sera_metadata_csv}")

viral_library_df = pd.read_csv(viral_library_csv)
log_message(f"  Loaded {len(viral_library_df)} barcodes from {viral_library_csv}")
log_message("")


# =============================================================================
# Validation functions
# =============================================================================
def validate_group_column(titers_df, expected_group):
    """Validate that all rows have the expected group value."""
    invalid_groups = titers_df[titers_df["group"] != expected_group]
    if len(invalid_groups) > 0:
        invalid_values = invalid_groups["group"].unique().tolist()
        raise ValueError(
            f"Titer data contains rows with group != '{expected_group}': "
            f"{invalid_values}"
        )


def validate_unique_serum_virus_pairs(titers_df):
    """Validate that serum-virus pairs are unique."""
    duplicates = titers_df.duplicated(subset=["serum", "virus"], keep=False)
    if duplicates.any():
        dup_pairs = titers_df.loc[duplicates, ["serum", "virus"]].drop_duplicates()
        raise ValueError(
            f"Duplicate serum-virus pairs in titer data:\n{dup_pairs.head(10)}"
        )


def validate_sera_in_metadata(titers_df, sera_metadata_df, group):
    """Validate that all sera in titers have matching metadata with species=group."""
    # Filter metadata to only include rows where species matches group
    group_sera = sera_metadata_df[sera_metadata_df["species"] == group][
        "serum"
    ].unique()
    titers_sera = titers_df["serum"].unique()

    missing_sera = set(titers_sera) - set(group_sera)
    if missing_sera:
        raise ValueError(
            f"Sera in titers not found in metadata with species='{group}': "
            f"{sorted(missing_sera)[:10]}{'...' if len(missing_sera) > 10 else ''}"
        )


def validate_viruses_in_library(titers_df, viral_library_df):
    """Validate that all viruses in titers are in the viral library."""
    library_strains = viral_library_df["strain"].unique()
    titers_viruses = titers_df["virus"].unique()

    missing_viruses = set(titers_viruses) - set(library_strains)
    if missing_viruses:
        ellipsis = "..." if len(missing_viruses) > 10 else ""
        raise ValueError(
            f"Viruses in titers not found in viral library: "
            f"{sorted(missing_viruses)[:10]}{ellipsis}"
        )


# Run validations
log_message("Validating input data...")
validate_group_column(titers_df, group)
log_message("  All rows have correct group")

validate_unique_serum_virus_pairs(titers_df)
log_message("  All serum-virus pairs are unique")

validate_sera_in_metadata(titers_df, sera_metadata_df, group)
log_message("  All sera found in metadata")

validate_viruses_in_library(titers_df, viral_library_df)
log_message("  All viruses found in viral library")
log_message("")

# Extract viral library strains for filtering
viral_library_strains = viral_library_df["strain"].unique()


# =============================================================================
# Processing functions
# =============================================================================
def drop_explicit_items(titers_df, sera_to_drop, viruses_to_drop):
    """Drop explicitly specified sera and viruses."""
    n_before = len(titers_df)

    # Drop viruses
    if viruses_to_drop:
        titers_df = titers_df[~titers_df["virus"].isin(viruses_to_drop)]
        n_after_viruses = len(titers_df)
        log_message(
            f"  Dropped {n_before - n_after_viruses} rows for "
            f"{len(viruses_to_drop)} viruses in viruses_to_drop"
        )
    else:
        n_after_viruses = n_before
        log_message("  No viruses specified in viruses_to_drop")

    # Drop sera
    if sera_to_drop:
        titers_df = titers_df[~titers_df["serum"].isin(sera_to_drop)]
        n_after_sera = len(titers_df)
        log_message(
            f"  Dropped {n_after_viruses - n_after_sera} rows for "
            f"{len(sera_to_drop)} sera in sera_to_drop"
        )
    else:
        log_message("  No sera specified in sera_to_drop")

    return titers_df


def apply_min_frac_filters(
    titers_df, min_frac_viruses, min_frac_sera, action, viral_library_strains
):
    """Apply minimum fraction filters iteratively until stable.

    Args:
        titers_df: DataFrame with titer data
        min_frac_viruses: Minimum fraction of viruses a serum must have titers for
        min_frac_sera: Minimum fraction of sera a virus must have titers from
        action: "drop" to remove failing items, "raise" to error
        viral_library_strains: Set of all virus names from viral library to check

    Returns:
        Tuple of (filtered DataFrame, set of excluded viruses)
    """
    iteration = 0
    max_iterations = 100  # Safety limit
    all_dropped_sera = []  # Track all dropped sera with their fractions
    all_dropped_viruses = []  # Track all dropped viruses with their fractions
    # Track all viruses failing min_frac_sera across iterations
    excluded_viruses = set()

    while iteration < max_iterations:
        iteration += 1
        n_sera_before = titers_df["serum"].nunique()
        n_viruses_before = titers_df["virus"].nunique()

        # Check min_frac_viruses - each serum must have titers for X% of viruses
        total_viruses = titers_df["virus"].nunique()
        sera_virus_counts = titers_df.groupby("serum")["virus"].nunique()
        sera_frac = sera_virus_counts / total_viruses
        failing_sera_frac = sera_frac[sera_frac < min_frac_viruses]
        failing_sera = failing_sera_frac.index.tolist()

        if failing_sera:
            if action == "raise":
                raise ValueError(
                    f"Sera below min_frac_viruses={min_frac_viruses}: {failing_sera}"
                )
            # Record dropped sera with their fractions
            for serum in failing_sera:
                all_dropped_sera.append((serum, failing_sera_frac[serum]))
            titers_df = titers_df[~titers_df["serum"].isin(failing_sera)]

        # Check min_frac_sera - each virus must have titers from X% of sera
        # This includes ALL viruses from viral library, including those with 0 titers
        total_sera = titers_df["serum"].nunique()
        if total_sera > 0:
            # Only check viruses not yet excluded
            viruses_to_check_iter = viral_library_strains - excluded_viruses

            # Get counts from titers_df (missing viruses have count=0)
            virus_sera_counts_from_titers = titers_df.groupby("virus")[
                "serum"
            ].nunique()

            # Create complete series for all viruses to check (missing = 0)
            virus_sera_counts = pd.Series(
                {
                    v: virus_sera_counts_from_titers.get(v, 0)
                    for v in viruses_to_check_iter
                }
            )

            virus_frac = virus_sera_counts / total_sera
            failing_virus_frac = virus_frac[virus_frac < min_frac_sera]
            failing_viruses = failing_virus_frac.index.tolist()

            if failing_viruses:
                if action == "raise":
                    raise ValueError(
                        f"Viruses below min_frac_sera={min_frac_sera}: "
                        f"{failing_viruses}"
                    )
                # Record newly dropped viruses (avoid duplicates)
                for virus in failing_viruses:
                    all_dropped_viruses.append((virus, failing_virus_frac[virus]))
                    excluded_viruses.add(virus)

                # Remove failing viruses from titers_df
                # (no-op for viruses with 0 titers)
                titers_df = titers_df[~titers_df["virus"].isin(failing_viruses)]

        # Check for convergence
        n_sera_after = titers_df["serum"].nunique()
        n_viruses_after = titers_df["virus"].nunique()

        if n_sera_after == n_sera_before and n_viruses_after == n_viruses_before:
            log_message(f"  min_frac filters converged after {iteration} iterations")
            break

        log_message(
            f"  Iteration {iteration}: dropped {n_sera_before - n_sera_after} sera, "
            f"{n_viruses_before - n_viruses_after} viruses"
        )

    if iteration >= max_iterations:
        raise RuntimeError(
            f"min_frac filters did not converge after {max_iterations} iterations"
        )

    # Log details of dropped sera
    if all_dropped_sera:
        log_message("")
        log_message(
            f"  Sera dropped for min_frac_viruses < {min_frac_viruses} "
            f"({len(all_dropped_sera)} total):"
        )
        for serum, frac in sorted(all_dropped_sera, key=lambda x: x[1]):
            log_message(f"    {serum}: frac_viruses={frac:.4f}")

    # Log details of dropped viruses
    if all_dropped_viruses:
        log_message("")
        log_message(
            f"  Viruses dropped for min_frac_sera < {min_frac_sera} "
            f"({len(all_dropped_viruses)} total):"
        )
        for virus, frac in sorted(all_dropped_viruses, key=lambda x: (x[1], x[0])):
            log_message(f"    {virus}: frac_sera={frac:.4f}")

    return titers_df, excluded_viruses


# =============================================================================
# Apply processing steps
# =============================================================================
log_message("Processing titer data...")
initial_sera = titers_df["serum"].nunique()
initial_viruses = titers_df["virus"].nunique()
initial_rows = len(titers_df)
log_message(
    f"  Initial: {initial_rows} rows, {initial_sera} sera, {initial_viruses} viruses"
)
log_message("")

# Step 1: Drop explicit items
log_message("Step 1: Dropping explicitly specified sera and viruses...")
titers_df = drop_explicit_items(titers_df, sera_to_drop, viruses_to_drop)
log_message("")

# Determine viruses to check against min_frac_sera threshold
# (all library viruses except those explicitly dropped)
viruses_to_check = set(viral_library_strains) - set(viruses_to_drop)

# Step 2: Apply min_frac filters
log_message(
    f"Step 2: Applying min_frac_viruses={min_frac_viruses}, "
    f"min_frac_sera={min_frac_sera} (action={min_frac_action})..."
)
titers_df, excluded_viruses_set = apply_min_frac_filters(
    titers_df, min_frac_viruses, min_frac_sera, min_frac_action, viruses_to_check
)
log_message("")

# Final counts
final_sera = titers_df["serum"].nunique()
final_viruses = titers_df["virus"].nunique()
final_rows = len(titers_df)
log_message("Processing complete:")
log_message(f"  Final: {final_rows} rows, {final_sera} sera, {final_viruses} viruses")
log_message(
    f"  Dropped: {initial_rows - final_rows} rows, "
    f"{initial_sera - final_sera} sera, {initial_viruses - final_viruses} viruses"
)
log_message("")


# =============================================================================
# Generate output files
# =============================================================================
log_message("Generating output files...")

# Get final list of sera and viruses
final_sera_list = titers_df["serum"].unique()
final_viruses_list = titers_df["virus"].unique()

# --- Output 1: Titers CSV ---
titers_output = titers_df[
    ["serum", "virus", "titer", "titer_bound", "titer_sem", "n_replicates", "titer_as"]
].copy()
titers_output.to_csv(output_titers_csv, index=False, float_format="%.4g")
log_message(f"  Written: {output_titers_csv} ({len(titers_output)} rows)")

# --- Output 2: Sera CSV ---
# Filter to sera in final titers, drop species, rename collection_date
sera_output = sera_metadata_df[sera_metadata_df["serum"].isin(final_sera_list)].copy()

# Verify no existing serum_collection_date column
if "serum_collection_date" in sera_output.columns:
    raise ValueError(
        "Sera metadata already has 'serum_collection_date' column; "
        "cannot rename 'collection_date'"
    )

sera_output = sera_output.drop(columns=["species"])
sera_output = sera_output.rename(columns={"collection_date": "serum_collection_date"})

# Reorder columns
sera_cols = ["serum", "cohort", "age", "sex", "serum_collection_date", "age_numeric"]
other_sera_cols = [c for c in sera_output.columns if c not in sera_cols]
sera_output = sera_output[sera_cols + other_sera_cols]

sera_output.to_csv(output_sera_csv, index=False)
log_message(f"  Written: {output_sera_csv} ({len(sera_output)} rows)")

# --- Output 3: Sera multicohort CSV ---
# First validate no conflicting cohort names (fail fast)
original_cohorts = sera_output["cohort"].unique()

# Check for "All" conflict (case-insensitive)
if any(c.lower() == "all" for c in original_cohorts if pd.notna(c)):
    raise ValueError(
        "Cannot have a cohort named 'all' (case-insensitive) as it conflicts with "
        "the derived 'All' cohort"
    )

# Build list of derived days-post-vax cohort names that will be created
derived_dpv_cohorts = set()
if "days_post_vax" in sera_output.columns:
    for _, row in sera_output.iterrows():
        if pd.notna(row.get("days_post_vax")):
            days = int(row["days_post_vax"])
            dpv_cohort = f"{row['cohort']}_{days}d-post-vax"
            derived_dpv_cohorts.add(dpv_cohort)

# Check for conflicts with derived days-post-vax cohort names
conflicting_cohorts = set(original_cohorts) & derived_dpv_cohorts
if conflicting_cohorts:
    raise ValueError(
        f"Original cohort names conflict with derived days-post-vax cohort names: "
        f"{sorted(conflicting_cohorts)}"
    )

# Build cohort assignments: each serum gets "All", their original cohort,
# and optionally a days-post-vax cohort
cohort_assignments = []
for _, row in sera_output.iterrows():
    serum = row["serum"]
    original_cohort = row["cohort"]

    # All sera belong to "All" cohort
    cohort_assignments.append({"serum": serum, "cohort": "All"})

    # Original cohort
    cohort_assignments.append({"serum": serum, "cohort": original_cohort})

    # Days-post-vax cohort if applicable
    if "days_post_vax" in sera_output.columns and pd.notna(row.get("days_post_vax")):
        days = int(row["days_post_vax"])
        dpv_cohort = f"{original_cohort}_{days}d-post-vax"
        cohort_assignments.append({"serum": serum, "cohort": dpv_cohort})

cohort_assignments_df = pd.DataFrame(cohort_assignments)

# Merge with sera metadata (drop original cohort column first to avoid conflict)
sera_multicohort = cohort_assignments_df.merge(
    sera_output.drop(columns=["cohort"]),
    on="serum",
    how="left",
    validate="many_to_one",
)

# Reorder columns
sera_mc_cols = ["serum", "cohort", "age", "sex", "serum_collection_date", "age_numeric"]
other_mc_cols = [c for c in sera_multicohort.columns if c not in sera_mc_cols]
sera_multicohort = sera_multicohort[sera_mc_cols + other_mc_cols]

sera_multicohort.to_csv(output_sera_multicohort_csv, index=False)
log_message(f"  Written: {output_sera_multicohort_csv} ({len(sera_multicohort)} rows)")

# --- Output 4: Viruses CSV ---
# Filter to viruses that passed all filtering
# (all library viruses except explicitly dropped and failing min_frac_sera)
viruses_that_passed = (
    set(viral_library_strains) - set(viruses_to_drop) - excluded_viruses_set
)
viruses_output = viral_library_df[
    viral_library_df["strain"].isin(viruses_that_passed)
].copy()

# Verify no existing virus or virus_collection_date column
if "virus" in viruses_output.columns:
    raise ValueError("Viral library already has 'virus' column; cannot rename 'strain'")
if "virus_collection_date" in viruses_output.columns:
    raise ValueError(
        "Viral library already has 'virus_collection_date' column; "
        "cannot rename 'collection_date'"
    )

# Rename columns
viruses_output = viruses_output.rename(
    columns={"strain": "virus", "collection_date": "virus_collection_date"}
)

# Select only the columns we want (excludes barcode, Twist_name, etc.)
virus_cols = [
    "virus",
    "subtype",
    "strain_type",
    "vaccine_type",
    "subclade",
    "derived_haplotype",
    "virus_collection_date",
    "protein_sequence_HA_ectodomain",
    "nt_sequence_HA_ectodomain",
]
# Only include columns that exist
virus_cols = [c for c in virus_cols if c in viruses_output.columns]
viruses_output = viruses_output[virus_cols]

# Drop duplicates (multiple barcodes per strain have identical metadata)
viruses_output = viruses_output.drop_duplicates()

# Verify single row per virus
virus_counts = viruses_output["virus"].value_counts()
if (virus_counts > 1).any():
    multi_row_viruses = virus_counts[virus_counts > 1].index.tolist()
    raise ValueError(
        f"Multiple rows per virus after dropping duplicates: " f"{multi_row_viruses}"
    )

viruses_output.to_csv(output_viruses_csv, index=False)
log_message(f"  Written: {output_viruses_csv} ({len(viruses_output)} rows)")

# --- Output 5: Titers summarized by virus ---
# Merge titers with multi-cohort sera assignments (one-to-many: each titer row
# is duplicated for each cohort the serum belongs to)
titers_multicohort = titers_output.merge(
    sera_multicohort[["serum", "cohort"]],
    on="serum",
    how="left",
    validate="many_to_many",  # many titers per serum, many cohorts per serum
)


def compute_virus_summary(df):
    """Compute summary statistics for each virus-cohort combination."""
    summary = (
        df.groupby(["virus", "cohort"])
        .agg(
            n_sera=pd.NamedAgg("serum", "nunique"),
            median_titer=pd.NamedAgg("titer", "median"),
            geomean_titer=pd.NamedAgg("titer", lambda x: np.exp(np.mean(np.log(x)))),
            titer_q1=pd.NamedAgg("titer", lambda x: x.quantile(0.25)),
            titer_q3=pd.NamedAgg("titer", lambda x: x.quantile(0.75)),
        )
        .reset_index()
    )

    # Add fraction below each cutoff
    for cutoff in titer_cutoffs:
        col_name = f"frac_w_titer_below_{cutoff}"
        frac_below = (
            df.groupby(["virus", "cohort"])["titer"]
            .apply(lambda x: (x < cutoff).mean())
            .reset_index(name=col_name)
        )
        summary = summary.merge(frac_below, on=["virus", "cohort"])

    return summary


titers_summarized = compute_virus_summary(titers_multicohort)


# Sort by virus, then cohort (with "All" first, original cohorts second, dpv last)
def cohort_sort_key(cohort):
    if cohort == "All":
        return (0, "")
    elif "d-post-vax" in cohort:
        return (2, cohort)
    else:
        return (1, cohort)


titers_summarized["_sort_key"] = titers_summarized["cohort"].apply(cohort_sort_key)
titers_summarized = titers_summarized.sort_values(["virus", "_sort_key"]).drop(
    columns=["_sort_key"]
)

# Get cohort counts for summary logging
cohort_info = (
    sera_multicohort.groupby("cohort")["serum"]
    .nunique()
    .reset_index()
    .rename(columns={"serum": "n_sera"})
)
# Add titer counts
titer_counts = titers_multicohort.groupby("cohort").size().reset_index(name="n_titers")
cohort_info = cohort_info.merge(titer_counts, on="cohort")
# Sort cohort info same way as summary
cohort_info["_sort_key"] = cohort_info["cohort"].apply(cohort_sort_key)
cohort_info = cohort_info.sort_values("_sort_key").drop(columns=["_sort_key"])

titers_summarized.to_csv(output_titers_summarized_csv, index=False, float_format="%.4g")
log_message(
    f"  Written: {output_titers_summarized_csv} ({len(titers_summarized)} rows)"
)
log_message("")

# --- Output 6: Summary text file ---
log_message("=" * 70)
log_message("FINAL SUMMARY")
log_message("=" * 70)
log_message(f"Group: {group}")
log_message(f"Sera: {final_sera}")
log_message(f"Viruses: {final_viruses}")
log_message(f"Total titers: {final_rows}")
log_message("")
log_message("Cohorts in titers_summarized_by_virus.csv:")
for _, row in cohort_info.iterrows():
    log_message(f"  {row['cohort']}: {row['n_sera']} sera, {row['n_titers']} titers")
log_message("")

log_message("Viruses per subtype:")
subtype_counts = viruses_output["subtype"].value_counts().sort_index()
for subtype, count in subtype_counts.items():
    log_message(f"  {subtype}: {count}")
log_message("")

# Write summary file
with open(output_summary_txt, "w") as f:
    f.write("\n".join(summary_messages))

log_message(f"Written: {output_summary_txt}")
log_message("")
log_message("Done!")

# Flush log
log.flush()
