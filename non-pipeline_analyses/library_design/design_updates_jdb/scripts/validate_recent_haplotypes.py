"""
Validate HA protein sequences in haplotype data.

This script:
1. Ensures all HA sequences are of equal length
2. Checks which sequences contain only valid amino acids
3. Generates a histogram of haplotype counts colored by sequence validity
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend

# Valid amino acid characters (20 standard amino acids)
VALID_AA = set('ACDEFGHIKLMNPQRSTVWY')


def validate_sequences(tsv_path, report_path, scatter_path, valid_path, invalid_path, select_cols, override_haplotypes):
    """
    Validate HA sequences and generate reports.

    Parameters
    ----------
    tsv_path : str
        Path to input TSV file with haplotype data
    report_path : str
        Path to output validation report
    scatter_path : str
        Path to output scatter plot PDF
    valid_path : str
        Path to output TSV of valid haplotypes
    invalid_path : str
        Path to output TSV of invalid haplotypes
    select_cols : list
        List of column names to use for selecting haplotypes
    override_haplotypes : list
        List of haplotype names (from derived_haplotype column) to exclude from selection
    """
    # Read data
    df = pd.read_csv(tsv_path, sep='\t')

    # Check for duplicate values in key columns
    duplicate_errors = []

    for col_name in ['representative_strain', 'derived_haplotype']:
        if col_name not in df.columns:
            raise ValueError(f"Required column '{col_name}' not found in input data")

        duplicates = df[col_name].duplicated(keep=False)
        if duplicates.any():
            dup_values = df[duplicates][col_name].value_counts()
            examples = []
            for val, count in dup_values.head(5).items():
                # Show which rows contain this duplicate
                dup_rows = df[df[col_name] == val].index.tolist()
                if col_name == 'representative_strain_ha_sequence':
                    # For sequences, truncate to show first/last part
                    display_val = f"{val[:30]}...{val[-30:]}" if len(val) > 65 else val
                else:
                    display_val = val
                examples.append(f"    '{display_val}' appears {count} times (rows: {dup_rows[:5]})")

            duplicate_errors.append(
                f"  Column '{col_name}' contains {len(dup_values)} duplicate value(s):\n" +
                "\n".join(examples) +
                (f"\n    ... and {len(dup_values) - 5} more" if len(dup_values) > 5 else "")
            )

    if duplicate_errors:
        raise ValueError(
            "Input data contains duplicate values in columns that must be unique:\n" +
            "\n\n".join(duplicate_errors)
        )

    # Check if valid_sequence column already exists
    if 'valid_sequence' in df.columns:
        raise ValueError(
            "Input data already contains a 'valid_sequence' column. "
            "This column will be created by this script and must not already exist."
        )

    # Check if selected_haplotype column already exists
    if 'selected_haplotype' in df.columns:
        raise ValueError(
            "Input data already contains a 'selected_haplotype' column. "
            "This column will be created by this script and must not already exist."
        )

    # Validate that all selection columns exist in the DataFrame
    missing_cols = [col for col in select_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(
            f"The following selection columns are missing from the input TSV: {missing_cols}\n"
            f"Available columns: {list(df.columns)}"
        )

    # Convert count column to numeric, handling comma separators
    # Some count values may be formatted with commas (e.g., "1,333")
    if df['count'].dtype == 'object':
        df['count'] = df['count'].str.replace(',', '').astype(int)
    else:
        df['count'] = df['count'].astype(int)

    # Verify all counts are valid
    if df['count'].isna().any():
        raise ValueError(
            f"Some 'count' values could not be converted to numeric.\n"
            f"Rows with invalid counts: {df[df['count'].isna()][['derived_haplotype']].head()}"
        )

    # Extract sequences
    sequences = df['representative_strain_ha_sequence']

    # Check sequence lengths
    seq_lengths = sequences.str.len()
    unique_lengths = seq_lengths.unique()

    if len(unique_lengths) != 1:
        raise ValueError(
            f"Sequences have different lengths: {sorted(unique_lengths)}\n"
            f"All sequences must be the same length for proper alignment."
        )

    seq_length = unique_lengths[0]

    # Check amino acid validity
    def has_valid_amino_acids(seq):
        """Check if sequence contains only valid amino acids."""
        return all(aa in VALID_AA for aa in seq)

    df['valid_sequence'] = sequences.apply(has_valid_amino_acids)

    # Create selected_haplotype column
    # A haplotype is selected if ANY of the selection columns evaluate to True
    # NA values should be treated as False
    df['selected_haplotype'] = df[select_cols].fillna(False).any(axis=1)

    # Validate that all override haplotypes exist in the data
    if override_haplotypes:
        existing_haplotypes = set(df['derived_haplotype'])
        missing_overrides = [h for h in override_haplotypes if h not in existing_haplotypes]
        if missing_overrides:
            raise ValueError(
                f"The following override haplotypes were not found in the input data:\n" +
                "\n".join(f"  - {h}" for h in missing_overrides) +
                f"\n\nAvailable haplotypes: {sorted(existing_haplotypes)[:10]}..."
            )

    # Apply override exclusions
    # Log which haplotypes would have been selected but are being excluded
    if override_haplotypes:
        initially_selected = df['selected_haplotype'].copy()
        excluded_mask = df['derived_haplotype'].isin(override_haplotypes)
        excluded_haplotypes = df[initially_selected & excluded_mask]

        if len(excluded_haplotypes) > 0:
            print(f"\nExcluding {len(excluded_haplotypes)} haplotype(s) via override_select_recent_haplotypes:")
            for _, row in excluded_haplotypes.iterrows():
                print(f"  - {row['derived_haplotype']} (strain: {row['representative_strain']}, "
                      f"count: {row['count']}, valid_sequence: {row['valid_sequence']})")

        # Apply the override by setting selected_haplotype to False for override entries
        df.loc[excluded_mask, 'selected_haplotype'] = False

    # Validate that all selected haplotypes have valid sequences
    selected_invalid = df[df['selected_haplotype'] & ~df['valid_sequence']]
    if len(selected_invalid) > 0:
        examples = []
        for _, row in selected_invalid.head().iterrows():
            examples.append(
                f"  - {row['derived_haplotype']} (strain: {row['representative_strain']})\n"
                f"    Sequence: {row['representative_strain_ha_sequence']}"
            )
        raise ValueError(
            f"Found {len(selected_invalid)} selected haplotype(s) with invalid sequences.\n"
            f"All selected haplotypes must have valid amino acid sequences.\n"
            f"Examples:\n" +
            "\n".join(examples)
        )

    # Count valid vs invalid
    n_valid = df['valid_sequence'].sum()
    n_invalid = (~df['valid_sequence']).sum()
    n_selected = df['selected_haplotype'].sum()

    # Calculate total sequence counts
    total_sequences_valid = df[df['valid_sequence']]['count'].sum()
    total_sequences_invalid = df[~df['valid_sequence']]['count'].sum()

    # Identify invalid characters if any
    invalid_chars = set()
    if n_invalid > 0:
        for seq in sequences[~df['valid_sequence']]:
            invalid_chars.update(set(seq) - VALID_AA)

    # Save valid and invalid haplotypes to separate files
    valid_df = df[df['valid_sequence']]
    invalid_df = df[~df['valid_sequence']]

    valid_df.to_csv(valid_path, sep='\t', index=False)
    invalid_df.to_csv(invalid_path, sep='\t', index=False)

    # Write validation report
    with open(report_path, 'w') as f:
        f.write("HA Sequence Validation Report\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Input file: {tsv_path}\n\n")

        f.write("Sequence Length Check:\n")
        f.write(f"  All sequences have length: {seq_length}\n")
        f.write(f"  ✓ All sequences are the same length (properly aligned)\n\n")

        f.write("Amino Acid Validity Check:\n")
        f.write(f"  Total haplotypes: {len(df)}\n")
        f.write(f"  Haplotypes with valid amino acids only: {n_valid} ({100*n_valid/len(df):.1f}%)\n")
        f.write(f"  Haplotypes with invalid characters: {n_invalid} ({100*n_invalid/len(df):.1f}%)\n")
        f.write(f"\n  Total sequences (valid haplotypes): {total_sequences_valid:,}\n")
        f.write(f"  Total sequences (invalid haplotypes): {total_sequences_invalid:,}\n")
        f.write(f"  Total sequences (all): {total_sequences_valid + total_sequences_invalid:,}\n")

        f.write("\nHaplotype Selection:\n")
        f.write(f"  Selection columns: {select_cols}\n")
        if override_haplotypes:
            f.write(f"  Override exclusions: {len(override_haplotypes)} haplotype(s) specified\n")
            excluded_mask = df['derived_haplotype'].isin(override_haplotypes)
            n_excluded = excluded_mask.sum()
            if n_excluded > 0:
                f.write(f"  Excluded via override: {n_excluded} haplotype(s)\n")
                for haplotype in override_haplotypes:
                    if haplotype in df['derived_haplotype'].values:
                        row = df[df['derived_haplotype'] == haplotype].iloc[0]
                        f.write(f"    - {haplotype} (strain: {row['representative_strain']}, count: {row['count']})\n")
        f.write(f"  Selected haplotypes: {n_selected} ({100*n_selected/len(df):.1f}%)\n")
        total_sequences_selected = df[df['selected_haplotype']]['count'].sum()
        f.write(f"  Total sequences (selected haplotypes): {total_sequences_selected:,}\n")
        f.write(f"  ✓ All selected haplotypes have valid sequences\n")

        if n_invalid > 0:
            f.write(f"\n  Invalid characters found: {sorted(invalid_chars)}\n")
            f.write("\n  Examples of haplotypes with invalid characters:\n")
            invalid_examples = df[~df['valid_sequence']].head(5)
            for idx, row in invalid_examples.iterrows():
                invalid_in_seq = set(row['representative_strain_ha_sequence']) - VALID_AA
                f.write(f"    {row['derived_haplotype']}: contains {sorted(invalid_in_seq)}\n")
        else:
            f.write(f"  ✓ All sequences contain only valid amino acids\n")

        f.write("\nSequence Count Distribution:\n")
        f.write(f"  Total sequences represented: {df['count'].sum()}\n")
        f.write(f"  Mean sequences per haplotype: {df['count'].mean():.1f}\n")
        f.write(f"  Median sequences per haplotype: {df['count'].median():.0f}\n")
        f.write(f"  Max sequences for a haplotype: {df['count'].max()}\n")

    # Create scatter plot
    fig, ax = plt.subplots(figsize=(12, 8))

    # Group by count to get number of haplotypes at each count level
    # Separate valid and invalid
    valid_counts = df[df['valid_sequence']].groupby('count').size()
    invalid_counts = df[~df['valid_sequence']].groupby('count').size()

    # Plot scatter points
    if len(valid_counts) > 0:
        ax.scatter(valid_counts.index, valid_counts.values,
                  alpha=0.6, s=60, marker='o',
                  color='steelblue', edgecolors='darkblue', linewidth=0.8,
                  label=f'Valid AA only (n={n_valid}, seqs={total_sequences_valid:,})',
                  zorder=2)

    if len(invalid_counts) > 0:
        ax.scatter(invalid_counts.index, invalid_counts.values,
                  alpha=0.7, s=100, marker='^',
                  color='coral', edgecolors='darkred', linewidth=0.8,
                  label=f'Contains invalid AA (n={n_invalid}, seqs={total_sequences_invalid:,})',
                  zorder=3)

    ax.set_xlabel('Number of sequences with haplotype (count)', fontsize=12)
    ax.set_ylabel('Number of haplotypes', fontsize=12)
    ax.set_title(f'Distribution of Haplotype Counts\n{tsv_path.split("/")[-1]}', fontsize=14)
    ax.legend(loc='upper right')
    ax.grid(axis='both', alpha=0.3, linestyle='--', zorder=1)

    # Use log scale for x-axis if counts span large range
    max_count = df['count'].max()
    if max_count > 100:
        ax.set_xscale('log')
        ax.set_xlabel('Number of sequences with haplotype (log scale)', fontsize=12)

    plt.tight_layout()
    plt.savefig(scatter_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Validation complete:")
    print(f"  Sequence length: {seq_length}")
    print(f"  Valid haplotypes: {n_valid}/{len(df)} ({100*n_valid/len(df):.1f}%)")
    print(f"  Invalid haplotypes: {n_invalid}/{len(df)} ({100*n_invalid/len(df):.1f}%)")
    print(f"  Selected haplotypes: {n_selected}/{len(df)} ({100*n_selected/len(df):.1f}%)")
    print(f"  Total sequences (valid): {total_sequences_valid:,}")
    print(f"  Total sequences (selected): {total_sequences_selected:,}")
    print(f"  Report: {report_path}")
    print(f"  Scatter plot: {scatter_path}")
    print(f"  Valid haplotypes TSV: {valid_path}")
    print(f"  Invalid haplotypes TSV: {invalid_path}")


if __name__ == '__main__':
    # Snakemake provides these variables
    sys.stdout = sys.stderr = open(snakemake.log[0], "w")
    validate_sequences(
        tsv_path=snakemake.input.tsv,
        report_path=snakemake.output.report,
        scatter_path=snakemake.output.scatter,
        valid_path=snakemake.output.valid_haplotypes,
        invalid_path=snakemake.output.invalid_haplotypes,
        select_cols=snakemake.params.select_cols,
        override_haplotypes=snakemake.params.override_haplotypes
    )
