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
            f"Titer data contains rows with group != '{expected_group}': {invalid_values}"
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
        raise ValueError(
            f"Viruses in titers not found in viral library: "
            f"{sorted(missing_viruses)[:10]}{'...' if len(missing_viruses) > 10 else ''}"
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


def apply_min_frac_filters(titers_df, min_frac_viruses, min_frac_sera, action):
    """Apply minimum fraction filters iteratively until stable.

    Args:
        titers_df: DataFrame with titer data
        min_frac_viruses: Minimum fraction of viruses a serum must have titers for
        min_frac_sera: Minimum fraction of sera a virus must have titers from
        action: "drop" to remove failing items, "raise" to error

    Returns:
        Filtered DataFrame
    """
    iteration = 0
    max_iterations = 100  # Safety limit
    all_dropped_sera = []  # Track all dropped sera with their fractions
    all_dropped_viruses = []  # Track all dropped viruses with their fractions

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
        total_sera = titers_df["serum"].nunique()
        if total_sera > 0:
            virus_sera_counts = titers_df.groupby("virus")["serum"].nunique()
            virus_frac = virus_sera_counts / total_sera
            failing_virus_frac = virus_frac[virus_frac < min_frac_sera]
            failing_viruses = failing_virus_frac.index.tolist()

            if failing_viruses:
                if action == "raise":
                    raise ValueError(
                        f"Viruses below min_frac_sera={min_frac_sera}: {failing_viruses}"
                    )
                # Record dropped viruses with their fractions
                for virus in failing_viruses:
                    all_dropped_viruses.append((virus, failing_virus_frac[virus]))
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
        for virus, frac in sorted(all_dropped_viruses, key=lambda x: x[1]):
            log_message(f"    {virus}: frac_sera={frac:.4f}")

    return titers_df


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

# Step 2: Drop any sera/viruses with zero remaining titers (shouldn't happen, but safety)
sera_with_titers = titers_df["serum"].unique()
viruses_with_titers = titers_df["virus"].unique()

# Step 3: Apply min_frac filters
log_message(
    f"Step 2: Applying min_frac_viruses={min_frac_viruses}, "
    f"min_frac_sera={min_frac_sera} (action={min_frac_action})..."
)
titers_df = apply_min_frac_filters(
    titers_df, min_frac_viruses, min_frac_sera, min_frac_action
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

# --- Output 3: Viruses CSV ---
# Filter to viruses in final titers
viruses_output = viral_library_df[
    viral_library_df["strain"].isin(final_viruses_list)
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

# Select only the columns we want (this excludes barcode, Twist_name, etc.)
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

# Drop duplicates (multiple barcodes per strain have identical metadata after selecting columns)
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

# --- Output 4: Titers summarized by virus ---
# Merge titers with sera metadata to get cohort and days_post_vax
titers_with_cohort = titers_df.merge(
    sera_output[
        ["serum", "cohort"]
        + (["days_post_vax"] if "days_post_vax" in sera_output.columns else [])
    ],
    on="serum",
    how="left",
    validate="many_to_one",
)

# Check for cohort named "all" (case-insensitive)
cohorts = titers_with_cohort["cohort"].unique()
if any(c.lower() == "all" for c in cohorts if pd.notna(c)):
    raise ValueError(
        "Cannot have a cohort named 'all' (case-insensitive) as it conflicts with "
        "the 'All' summary cohort"
    )


def compute_virus_summary(df, cohort_name):
    """Compute summary statistics for each virus."""
    summary = (
        df.groupby("virus")
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
            df.groupby("virus")["titer"]
            .apply(lambda x: (x < cutoff).mean())
            .reset_index(name=col_name)
        )
        summary = summary.merge(frac_below, on="virus")

    summary.insert(1, "cohort", cohort_name)
    return summary


# Track all cohorts info for logging (cohort_name, n_sera, n_titers)
all_cohorts_info = []

# Compute summary for "All" cohorts
summaries = [compute_virus_summary(titers_with_cohort, "All")]
all_cohorts_info.append(
    ("All", titers_with_cohort["serum"].nunique(), len(titers_with_cohort))
)

# Compute summary for each cohort
for cohort in sorted(cohorts):
    cohort_df = titers_with_cohort[titers_with_cohort["cohort"] == cohort]
    summaries.append(compute_virus_summary(cohort_df, cohort))
    all_cohorts_info.append((cohort, cohort_df["serum"].nunique(), len(cohort_df)))

# Compute summary for cohort + days_post_vax combinations if days_post_vax exists
if "days_post_vax" in titers_with_cohort.columns:
    # Get unique cohort-days combinations with non-null days_post_vax
    cohort_days = titers_with_cohort[titers_with_cohort["days_post_vax"].notna()][
        ["cohort", "days_post_vax"]
    ].drop_duplicates()

    for _, row in cohort_days.iterrows():
        cohort = row["cohort"]
        days = int(row["days_post_vax"])
        cohort_name = f"{cohort}_{days}-days-post-vax"
        cohort_df = titers_with_cohort[
            (titers_with_cohort["cohort"] == cohort)
            & (titers_with_cohort["days_post_vax"] == days)
        ]
        if len(cohort_df) > 0:
            summaries.append(compute_virus_summary(cohort_df, cohort_name))
            all_cohorts_info.append(
                (cohort_name, cohort_df["serum"].nunique(), len(cohort_df))
            )

titers_summarized = pd.concat(summaries, ignore_index=True)
titers_summarized.to_csv(output_titers_summarized_csv, index=False, float_format="%.4g")
log_message(
    f"  Written: {output_titers_summarized_csv} ({len(titers_summarized)} rows)"
)
log_message("")

# --- Output 5: Summary text file ---
log_message("=" * 70)
log_message("FINAL SUMMARY")
log_message("=" * 70)
log_message(f"Group: {group}")
log_message(f"Sera: {final_sera}")
log_message(f"Viruses: {final_viruses}")
log_message(f"Total titers: {final_rows}")
log_message("")
log_message("Cohorts in titers_summarized_by_virus.csv:")
for cohort_name, n_sera, n_titers in all_cohorts_info:
    log_message(f"  {cohort_name}: {n_sera} sera, {n_titers} titers")
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
