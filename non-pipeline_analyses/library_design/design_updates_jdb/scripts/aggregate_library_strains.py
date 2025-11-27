"""
Aggregate library strains from multiple subtypes into a single TSV.

This script combines the match_genbank results from multiple subtypes (H1N1, H3N2)
into a single output file with validation checks.
"""

import pandas as pd
import sys


def aggregate_library_strains(input_csvs, output_tsv, subtypes):
    """
    Aggregate library strains from multiple subtype-specific CSVs.

    Parameters
    ----------
    input_csvs : list of str
        List of paths to input CSV files from match_genbank rule
    output_tsv : str
        Path to output TSV file
    subtypes : list of str
        List of subtypes corresponding to each input CSV
    """
    # Read and process each subtype
    dfs = []
    for csv_path, subtype in zip(input_csvs, subtypes):
        df = pd.read_csv(csv_path)

        # Add subtype column
        df['subtype'] = subtype

        # Add subtype suffix to strain names
        df['strain'] = df['strain'] + f'_{subtype}'

        # Select only required columns
        df = df[['subtype', 'strain', 'accession_w_aa_muts_added', 'prot_sequence', 'nt_sequence']]

        dfs.append(df)

    # Combine all dataframes
    combined_df = pd.concat(dfs, ignore_index=True)

    # Validation 1: Check for null values
    null_counts = combined_df.isnull().sum()
    if null_counts.any():
        null_info = null_counts[null_counts > 0]
        error_msg = "ERROR: Found null values in the following columns:\n"
        for col, count in null_info.items():
            error_msg += f"  - {col}: {count} null values\n"

        # Show rows with null values for debugging
        error_msg += "\nRows with null values:\n"
        null_rows = combined_df[combined_df.isnull().any(axis=1)]
        error_msg += null_rows.to_string()

        raise ValueError(error_msg)

    # Validation 2: Check for duplicate strain names
    duplicate_strains = combined_df['strain'].duplicated()
    if duplicate_strains.any():
        duplicate_names = combined_df.loc[duplicate_strains, 'strain'].unique()
        error_msg = f"ERROR: Found {len(duplicate_names)} duplicate strain names:\n"
        for name in duplicate_names:
            count = (combined_df['strain'] == name).sum()
            error_msg += f"  - {name} (appears {count} times)\n"

        raise ValueError(error_msg)

    # Write output
    combined_df.to_csv(output_tsv, sep='\t', index=False)

    # Print summary statistics
    print(f"Successfully aggregated library strains:")
    print(f"  - Total strains: {len(combined_df)}")
    for subtype in subtypes:
        count = (combined_df['subtype'] == subtype).sum()
        print(f"  - {subtype}: {count} strains")
    print(f"  - Output written to: {output_tsv}")


if __name__ == "__main__":
    # Redirect stdout and stderr to log file
    sys.stderr = sys.stdout = open(snakemake.log[0], "w")

    # Snakemake workflow will provide these variables
    aggregate_library_strains(
        input_csvs=snakemake.input.csvs,
        output_tsv=snakemake.output.tsv,
        subtypes=snakemake.params.subtypes
    )
