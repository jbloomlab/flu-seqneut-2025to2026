"""Validate a viral library CSV file.

Performs comprehensive validation of viral library files according to the
specifications in README.md, outputs an informative summary, and raises
ValueError if any validations fail.
"""

import re
import sys

from Bio.Seq import Seq

import pandas as pd


sys.stderr = sys.stdout = log = open(snakemake.log[0], "w")

# Get inputs and params from snakemake
csv_path = snakemake.input.csv
output_path = snakemake.output.validation
circulating_strain_type = snakemake.params.circulating_strain_type

# Read the viral library CSV
print(f"Reading viral library from: {csv_path}")
df = pd.read_csv(csv_path)
print(f"Loaded {len(df)} rows")

# Track all validation results
validation_results = []


def log_validation(name, passed, details=""):
    """Log a validation result and track it."""
    status = "PASS" if passed else "FAIL"
    msg = f"[{status}] {name}"
    if details:
        msg += f": {details}"
    print(msg)
    validation_results.append((name, passed, details))
    if not passed:
        raise ValueError(f"Validation failed: {name} - {details}")


# =============================================================================
# Validation 1: Required columns exist
# =============================================================================
required_columns = [
    "strain",
    "subtype",
    "strain_type",
    "vaccine_type",
    "barcode",
    "accession",
    "subclade",
    "derived_haplotype",
    "collection_date",
    "nt_sequence_HA_ectodomain",
    "protein_sequence_HA_ectodomain",
]

missing_columns = [col for col in required_columns if col not in df.columns]
log_validation(
    "Required columns exist",
    len(missing_columns) == 0,
    f"Missing: {missing_columns}" if missing_columns else "All present",
)

# =============================================================================
# Validation 2: strain - must end in "_H3N2" or "_H1N1"
# =============================================================================
invalid_strain_suffix = df[
    ~df["strain"].str.endswith("_H3N2") & ~df["strain"].str.endswith("_H1N1")
]
log_validation(
    "Strain names end with _H3N2 or _H1N1",
    len(invalid_strain_suffix) == 0,
    (
        f"Invalid: {invalid_strain_suffix['strain'].tolist()}"
        if len(invalid_strain_suffix) > 0
        else "All valid"
    ),
)

# =============================================================================
# Validation 3: strain consistency - each strain has same values for all
# required columns except barcode
# =============================================================================
# Only check consistency for required columns (excluding barcode)
# Other columns like Twist_name may intentionally vary per barcode
consistency_check_cols = [col for col in required_columns if col != "barcode"]
strain_consistency_issues = []
for strain, group in df.groupby("strain"):
    for col in consistency_check_cols:
        unique_vals = group[col].dropna().unique()
        # Handle case where column has null and non-null values
        has_null = group[col].isna().any()
        has_non_null = group[col].notna().any()
        if len(unique_vals) > 1 or (has_null and has_non_null):
            strain_consistency_issues.append(
                f"{strain}: column '{col}' has inconsistent values"
            )

log_validation(
    "Strain values consistent across barcodes",
    len(strain_consistency_issues) == 0,
    (
        f"Issues: {strain_consistency_issues[:5]}"
        if strain_consistency_issues
        else "All consistent"
    ),
)

# =============================================================================
# Validation 4: subtype must be "H3N2" or "H1N1"
# =============================================================================
invalid_subtype = df[~df["subtype"].isin(["H3N2", "H1N1"])]
log_validation(
    "Subtype is H3N2 or H1N1",
    len(invalid_subtype) == 0,
    (
        f"Invalid: {invalid_subtype['subtype'].unique().tolist()}"
        if len(invalid_subtype) > 0
        else "All valid"
    ),
)

# =============================================================================
# Validation 5: strain_type must be "vaccine" or circulating_strain_type
# =============================================================================
valid_strain_types = ["vaccine", circulating_strain_type]
invalid_strain_type = df[~df["strain_type"].isin(valid_strain_types)]
log_validation(
    f"Strain type is 'vaccine' or '{circulating_strain_type}'",
    len(invalid_strain_type) == 0,
    (
        f"Invalid: {invalid_strain_type['strain_type'].unique().tolist()}"
        if len(invalid_strain_type) > 0
        else "All valid"
    ),
)

# =============================================================================
# Validation 6: vaccine_type - "egg" or "cell" for vaccines, null otherwise
# =============================================================================
vaccine_rows = df[df["strain_type"] == "vaccine"]
non_vaccine_rows = df[df["strain_type"] != "vaccine"]

# Vaccines must have vaccine_type of "egg" or "cell"
invalid_vaccine_type = vaccine_rows[~vaccine_rows["vaccine_type"].isin(["egg", "cell"])]

# Non-vaccines must have null vaccine_type
non_vaccine_with_type = non_vaccine_rows[non_vaccine_rows["vaccine_type"].notna()]

vaccine_type_issues = []
if len(invalid_vaccine_type) > 0:
    vaccine_type_issues.append(
        f"Vaccines with invalid vaccine_type: "
        f"{invalid_vaccine_type['strain'].unique().tolist()}"
    )
if len(non_vaccine_with_type) > 0:
    vaccine_type_issues.append(
        f"Non-vaccines with vaccine_type set: "
        f"{non_vaccine_with_type['strain'].unique().tolist()}"
    )

log_validation(
    "Vaccine type valid",
    len(vaccine_type_issues) == 0,
    "; ".join(vaccine_type_issues) if vaccine_type_issues else "All valid",
)

# =============================================================================
# Validation 7: barcode - 16-nt string, all unique
# =============================================================================
# Check barcode format (16-nt, all ATCG)
barcode_pattern = re.compile(r"^[ATCG]{16}$", re.IGNORECASE)
invalid_barcode_format = df[~df["barcode"].str.match(barcode_pattern, na=False)]
log_validation(
    "Barcode format (16-nt ATCG)",
    len(invalid_barcode_format) == 0,
    (
        f"Invalid: {invalid_barcode_format['barcode'].tolist()[:5]}"
        if len(invalid_barcode_format) > 0
        else "All valid"
    ),
)

# Check barcode uniqueness
duplicate_barcodes = df[df["barcode"].duplicated(keep=False)]
if len(duplicate_barcodes) > 0:
    dup_list = duplicate_barcodes["barcode"].unique().tolist()[:5]
else:
    dup_list = []
log_validation(
    "Barcodes unique",
    len(duplicate_barcodes) == 0,
    f"Duplicates: {dup_list}" if dup_list else "All unique",
)

# =============================================================================
# Validation 8: accession - no validation (can be null)
# =============================================================================
print("[INFO] Accession: no validation required (can be null)")

# =============================================================================
# Validation 9: subclade - can only be null for strain_type == "vaccine"
# =============================================================================
non_vaccine_null_subclade = non_vaccine_rows[non_vaccine_rows["subclade"].isna()]
log_validation(
    "Subclade present for non-vaccines",
    len(non_vaccine_null_subclade) == 0,
    (
        f"Missing subclade: {non_vaccine_null_subclade['strain'].unique().tolist()[:5]}"
        if len(non_vaccine_null_subclade) > 0
        else "All non-vaccines have subclade"
    ),
)

# =============================================================================
# Validation 10: derived_haplotype - null only for vaccines, unique per strain
# =============================================================================
# Check non-vaccines have derived_haplotype
non_vaccine_null_haplotype = non_vaccine_rows[
    non_vaccine_rows["derived_haplotype"].isna()
]
log_validation(
    "Derived haplotype present for non-vaccines",
    len(non_vaccine_null_haplotype) == 0,
    (
        f"Missing: {non_vaccine_null_haplotype['strain'].unique().tolist()[:5]}"
        if len(non_vaccine_null_haplotype) > 0
        else "All non-vaccines have derived_haplotype"
    ),
)

# Note: derived_haplotype does not need to be unique per strain, as HA2 mutations
# are not included in the derived_haplotype naming convention

# =============================================================================
# Validation 11: collection_date - must be a valid float
# =============================================================================
# Check if collection_date can be converted to float
try:
    dates = pd.to_numeric(df["collection_date"], errors="coerce")
    invalid_dates = df[dates.isna()]
    log_validation(
        "Collection date is valid float",
        len(invalid_dates) == 0,
        (
            f"Invalid: {invalid_dates[['strain', 'collection_date']].values.tolist()[:5]}"
            if len(invalid_dates) > 0
            else "All valid"
        ),
    )
except Exception as e:
    log_validation("Collection date is valid float", False, str(e))

# =============================================================================
# Validation 12: nt_sequence_HA_ectodomain - all ATCG, length multiple of 3
# =============================================================================
nt_pattern = re.compile(r"^[ATCG]+$", re.IGNORECASE)
invalid_nt_chars = df[~df["nt_sequence_HA_ectodomain"].str.match(nt_pattern, na=False)]
log_validation(
    "NT sequence contains only ATCG",
    len(invalid_nt_chars) == 0,
    (
        f"Invalid strains: {invalid_nt_chars['strain'].tolist()[:5]}"
        if len(invalid_nt_chars) > 0
        else "All valid"
    ),
)

invalid_nt_len = df[df["nt_sequence_HA_ectodomain"].str.len() % 3 != 0]
log_validation(
    "NT sequence length multiple of 3",
    len(invalid_nt_len) == 0,
    (
        f"Invalid: {invalid_nt_len[['strain']].values.tolist()[:5]}"
        if len(invalid_nt_len) > 0
        else "All valid"
    ),
)

# =============================================================================
# Validation 13: protein_sequence_HA_ectodomain - comprehensive checks
# =============================================================================

# Check for stop codons
has_stop = df[df["protein_sequence_HA_ectodomain"].str.contains(r"\*", na=False)]
log_validation(
    "Protein sequence has no stop codons",
    len(has_stop) == 0,
    (
        f"Has stop: {has_stop['strain'].tolist()[:5]}"
        if len(has_stop) > 0
        else "No stop codons"
    ),
)

# Check translation matches
translation_mismatches = []
for idx, row in df.iterrows():
    nt_seq = row["nt_sequence_HA_ectodomain"].upper()
    protein_seq = row["protein_sequence_HA_ectodomain"]
    translated = str(Seq(nt_seq).translate())
    # Remove trailing stop codon if present in translation
    if translated.endswith("*"):
        translated = translated[:-1]
    if translated != protein_seq:
        translation_mismatches.append(row["strain"])

log_validation(
    "Protein equals translation of NT sequence",
    len(translation_mismatches) == 0,
    (
        f"Mismatches: {translation_mismatches[:5]}"
        if translation_mismatches
        else "All match"
    ),
)

# Check H1N1 protein sequences
h1n1_rows = df[df["subtype"] == "H1N1"]
h1n1_start_invalid = h1n1_rows[
    ~h1n1_rows["protein_sequence_HA_ectodomain"].str.startswith("CIGY")
]
log_validation(
    "H1N1 protein starts with CIGY",
    len(h1n1_start_invalid) == 0,
    (
        f"Invalid: {h1n1_start_invalid['strain'].tolist()[:5]}"
        if len(h1n1_start_invalid) > 0
        else "All valid"
    ),
)

h1n1_end_pattern = re.compile(r"[EK]IDG[VI]$")
h1n1_end_invalid = h1n1_rows[
    ~h1n1_rows["protein_sequence_HA_ectodomain"].str.match(r".*[EK]IDG[VI]$", na=False)
]
log_validation(
    "H1N1 protein ends with [EK]IDG[VI]",
    len(h1n1_end_invalid) == 0,
    (
        f"Invalid: {h1n1_end_invalid['strain'].tolist()[:5]}"
        if len(h1n1_end_invalid) > 0
        else "All valid"
    ),
)

h1n1_len_invalid = h1n1_rows[
    h1n1_rows["protein_sequence_HA_ectodomain"].str.len() != 500
]
log_validation(
    "H1N1 protein length is 500",
    len(h1n1_len_invalid) == 0,
    (
        f"Invalid lengths: "
        f"{h1n1_len_invalid[['strain']].assign(len=h1n1_len_invalid['protein_sequence_HA_ectodomain'].str.len()).values.tolist()[:5]}"
        if len(h1n1_len_invalid) > 0
        else "All valid"
    ),
)

# Check H3N2 protein sequences
h3n2_rows = df[df["subtype"] == "H3N2"]
h3n2_start_pattern = re.compile(r"^Q[KNR][IL]P")
h3n2_start_invalid = h3n2_rows[
    ~h3n2_rows["protein_sequence_HA_ectodomain"].str.match(r"^Q[KNR][IL]P", na=False)
]
log_validation(
    "H3N2 protein starts with Q[KNR][IL]P",
    len(h3n2_start_invalid) == 0,
    (
        f"Invalid: {h3n2_start_invalid['strain'].tolist()[:5]}"
        if len(h3n2_start_invalid) > 0
        else "All valid"
    ),
)

h3n2_end_invalid = h3n2_rows[
    ~h3n2_rows["protein_sequence_HA_ectodomain"].str.endswith("NNRFQ")
]
log_validation(
    "H3N2 protein ends with NNRFQ",
    len(h3n2_end_invalid) == 0,
    (
        f"Invalid: {h3n2_end_invalid['strain'].tolist()[:5]}"
        if len(h3n2_end_invalid) > 0
        else "All valid"
    ),
)

h3n2_len_invalid = h3n2_rows[
    h3n2_rows["protein_sequence_HA_ectodomain"].str.len() != 501
]
log_validation(
    "H3N2 protein length is 501",
    len(h3n2_len_invalid) == 0,
    (
        f"Invalid lengths: "
        f"{h3n2_len_invalid[['strain']].assign(len=h3n2_len_invalid['protein_sequence_HA_ectodomain'].str.len()).values.tolist()[:5]}"
        if len(h3n2_len_invalid) > 0
        else "All valid"
    ),
)

# Check unique protein sequence per strain
protein_to_strains = (
    df.drop_duplicates(["strain", "protein_sequence_HA_ectodomain"])
    .groupby("protein_sequence_HA_ectodomain")["strain"]
    .apply(list)
)
dup_protein_strains = protein_to_strains[protein_to_strains.apply(len) > 1]
log_validation(
    "Each unique protein sequence maps to one strain",
    len(dup_protein_strains) == 0,
    (
        f"Shared proteins: {dup_protein_strains.tolist()[:3]}"
        if len(dup_protein_strains) > 0
        else "All unique"
    ),
)

# =============================================================================
# Summary
# =============================================================================
print("\n" + "=" * 70)
print("VALIDATION SUMMARY")
print("=" * 70)

# Library summary
n_strains = df["strain"].nunique()
n_barcodes = len(df)
n_h1n1 = len(df[df["subtype"] == "H1N1"])
n_h3n2 = len(df[df["subtype"] == "H3N2"])
n_vaccine = len(df[df["strain_type"] == "vaccine"])
n_circulating = len(df[df["strain_type"] == circulating_strain_type])
n_vaccine_strains = df[df["strain_type"] == "vaccine"]["strain"].nunique()
n_circulating_strains = df[df["strain_type"] == circulating_strain_type][
    "strain"
].nunique()

print(f"\nLibrary file: {csv_path}")
print(f"Total rows (barcodes): {n_barcodes}")
print(f"Unique strains: {n_strains}")
print(f"  - Vaccine strains: {n_vaccine_strains} ({n_vaccine} barcodes)")
print(f"  - Circulating strains: {n_circulating_strains} ({n_circulating} barcodes)")
print("\nBy subtype:")
print(f"  - H1N1: {n_h1n1} barcodes")
print(f"  - H3N2: {n_h3n2} barcodes")

print(f"\nAll {len(validation_results)} validations PASSED")

# Flush log before reading it back
log.flush()

# Re-read log file and write to output
with open(snakemake.log[0]) as log_file:
    log_content = log_file.read()

with open(output_path, "w") as f:
    f.write(log_content)

print(f"\nValidation output written to: {output_path}")
