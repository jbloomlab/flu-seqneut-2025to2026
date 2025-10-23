"""
Curate library sequences by combining validated recent haplotypes and additional haplotypes.

This script:
1. Reads validated recent haplotypes (with selected_haplotype already defined)
2. Reads any additional haplotypes (and sets selected_haplotype=True for them)
3. Combines into a single dataframe
4. Validates that representative_strain is unique
5. Validates that all HA sequences are the same length and contain only valid amino acids
6. Computes Hamming distance to nearest library strain (excluding self)
7. Computes epitope-specific distances to nearest library strain for each epitope set
8. Outputs two TSVs: one for selected haplotypes (sorted by distance ascending),
   one for non-selected haplotypes (sorted by distance descending)
"""

import sys
import pandas as pd

# Valid amino acid characters (20 standard amino acids)
VALID_AA = set('ACDEFGHIKLMNPQRSTVWY')


def hamming_distance(seq1, seq2):
    """
    Calculate Hamming distance between two sequences of equal length.

    Parameters
    ----------
    seq1 : str
        First sequence
    seq2 : str
        Second sequence

    Returns
    -------
    int
        Number of positions where sequences differ
    """
    if len(seq1) != len(seq2):
        raise ValueError(f"Sequences must be same length: {len(seq1)} vs {len(seq2)}")
    return sum(a != b for a, b in zip(seq1, seq2))


def curate_library_sequences(valid_haplotypes_path, additional_haplotypes_paths,
                              site_annotations_path, selected_output_path, nonselected_output_path,
                              library_fasta_path):
    """
    Curate library sequences by combining validated and additional haplotypes.

    Parameters
    ----------
    valid_haplotypes_path : str
        Path to TSV file with validated recent haplotypes (from validate_recent_haplotypes)
    additional_haplotypes_paths : list of str
        List of paths to TSV files with additional haplotypes to include
    site_annotations_path : str
        Path to TSV file with site annotations including epitope definitions
    selected_output_path : str
        Path to output TSV of selected haplotypes
    nonselected_output_path : str
        Path to output TSV of non-selected haplotypes
    library_fasta_path : str
        Path to output FASTA file with library (selected) haplotype sequences
    """
    # Read validated recent haplotypes
    df_validated = pd.read_csv(valid_haplotypes_path, sep='\t')

    # Check that selected_haplotype column exists
    if 'selected_haplotype' not in df_validated.columns:
        raise ValueError(
            f"Input file '{valid_haplotypes_path}' does not contain 'selected_haplotype' column.\n"
            f"Available columns: {list(df_validated.columns)}"
        )

    # Check that required columns exist
    required_cols = ['representative_strain', 'representative_strain_ha_sequence']
    for col in required_cols:
        if col not in df_validated.columns:
            raise ValueError(
                f"Input file '{valid_haplotypes_path}' does not contain required column '{col}'.\n"
                f"Available columns: {list(df_validated.columns)}"
            )

    print(f"Read {len(df_validated)} validated haplotypes from {valid_haplotypes_path}")
    print(f"  Selected: {df_validated['selected_haplotype'].sum()}")
    print(f"  Non-selected: {(~df_validated['selected_haplotype']).sum()}")

    # Read additional haplotypes if any
    dfs_to_combine = [df_validated]

    if additional_haplotypes_paths:
        for path in additional_haplotypes_paths:
            df_additional = pd.read_csv(path, sep='\t')

            # Check required columns
            for col in required_cols:
                if col not in df_additional.columns:
                    raise ValueError(
                        f"Additional haplotypes file '{path}' does not contain required column '{col}'.\n"
                        f"Available columns: {list(df_additional.columns)}"
                    )

            # Set selected_haplotype to True for all additional haplotypes
            df_additional['selected_haplotype'] = True

            print(f"\nRead {len(df_additional)} additional haplotypes from {path}")
            print(f"  All additional haplotypes are marked as selected")

            dfs_to_combine.append(df_additional)

    # Combine all dataframes
    df_combined = pd.concat(dfs_to_combine, ignore_index=True)

    print(f"\nCombined dataframe contains {len(df_combined)} total haplotypes")
    print(f"  Selected: {df_combined['selected_haplotype'].sum()}")
    print(f"  Non-selected: {(~df_combined['selected_haplotype']).sum()}")

    # Check for duplicate representative_strain values
    duplicates = df_combined['representative_strain'].duplicated(keep=False)
    if duplicates.any():
        dup_values = df_combined[duplicates]['representative_strain'].value_counts()
        examples = []
        for val, count in dup_values.head(10).items():
            # Show which rows contain this duplicate
            dup_rows = df_combined[df_combined['representative_strain'] == val].index.tolist()
            examples.append(f"    '{val}' appears {count} times (rows: {dup_rows})")

        raise ValueError(
            f"Combined data contains {len(dup_values)} duplicate representative_strain value(s).\n"
            f"Each strain name must be unique.\n"
            f"Examples:\n" +
            "\n".join(examples) +
            (f"\n    ... and {len(dup_values) - 10} more" if len(dup_values) > 10 else "")
        )

    print("  ✓ All representative_strain values are unique")

    # Check sequence lengths
    sequences = df_combined['representative_strain_ha_sequence']
    seq_lengths = sequences.str.len()
    unique_lengths = seq_lengths.unique()

    if len(unique_lengths) != 1:
        # Show examples of different lengths
        length_examples = []
        for length in sorted(unique_lengths)[:5]:
            strains = df_combined[seq_lengths == length]['representative_strain'].head(3).tolist()
            length_examples.append(f"    Length {length}: {strains}")

        raise ValueError(
            f"HA sequences have different lengths: {sorted(unique_lengths)}\n"
            f"All sequences must be the same length for proper alignment.\n"
            f"Examples:\n" +
            "\n".join(length_examples) +
            (f"\n    ... and {len(unique_lengths) - 5} more lengths" if len(unique_lengths) > 5 else "")
        )

    seq_length = unique_lengths[0]
    print(f"  ✓ All sequences have length {seq_length}")

    # Check amino acid validity
    def has_valid_amino_acids(seq):
        """Check if sequence contains only valid amino acids."""
        return all(aa in VALID_AA for aa in seq)

    invalid_seqs = ~sequences.apply(has_valid_amino_acids)

    if invalid_seqs.any():
        # Identify invalid characters and show examples
        invalid_examples = []
        for idx in df_combined[invalid_seqs].index[:5]:
            row = df_combined.iloc[idx]
            seq = row['representative_strain_ha_sequence']
            invalid_chars = sorted(set(seq) - VALID_AA)
            invalid_examples.append(
                f"    {row['representative_strain']}: contains invalid character(s) {invalid_chars}"
            )

        n_invalid = invalid_seqs.sum()
        raise ValueError(
            f"Found {n_invalid} sequence(s) with invalid amino acids.\n"
            f"All sequences must contain only the 20 standard amino acids: {sorted(VALID_AA)}\n"
            f"Examples:\n" +
            "\n".join(invalid_examples) +
            (f"\n    ... and {n_invalid - 5} more" if n_invalid > 5 else "")
        )

    print(f"  ✓ All sequences contain only valid amino acids")

    # Split into selected and non-selected
    df_selected = df_combined[df_combined['selected_haplotype']].copy()
    df_nonselected = df_combined[~df_combined['selected_haplotype']].copy()

    print(f"\nComputing Hamming distances to nearest library strain...")

    # Get sequences and strain names for selected haplotypes (library)
    library_sequences = df_selected['representative_strain_ha_sequence'].values
    library_strains = df_selected['representative_strain'].values

    # For each haplotype (both selected and non-selected), find the nearest library strain
    def find_nearest_library_strain(sequence, current_strain_name):
        """
        Find the nearest selected (library) strain to a given sequence.

        Parameters
        ----------
        sequence : str
            The sequence to compare
        current_strain_name : str
            The name of the current strain (to exclude from comparison)

        Returns
        -------
        tuple of (int, str)
            Hamming distance and strain name of nearest library strain
        """
        distances = []
        valid_indices = []

        for idx, (lib_seq, lib_strain) in enumerate(zip(library_sequences, library_strains)):
            # Skip if this is the same strain (comparing to itself)
            if lib_strain == current_strain_name:
                continue
            distances.append(hamming_distance(sequence, lib_seq))
            valid_indices.append(idx)

        if not distances:
            # This should only happen if there's only one library strain
            raise ValueError(
                f"Cannot find nearest library strain for {current_strain_name}: "
                f"no other library strains available for comparison"
            )

        min_distance = min(distances)
        # Get the first (earliest in dataframe order) strain with minimum distance
        min_idx_in_valid = distances.index(min_distance)
        actual_idx = valid_indices[min_idx_in_valid]
        return min_distance, library_strains[actual_idx]

    # Compute for all haplotypes in the combined dataframe
    hamming_distances = []
    nearest_strains = []

    for _, row in df_combined.iterrows():
        sequence = row['representative_strain_ha_sequence']
        strain_name = row['representative_strain']
        dist, nearest_strain = find_nearest_library_strain(sequence, strain_name)
        hamming_distances.append(dist)
        nearest_strains.append(nearest_strain)

    # Add the new columns to the combined dataframe
    df_combined['hamming_distance_nearest_library_strain'] = hamming_distances
    df_combined['nearest_library_strain'] = nearest_strains

    # Update the selected and non-selected dataframes
    df_selected = df_combined[df_combined['selected_haplotype']].copy()
    df_nonselected = df_combined[~df_combined['selected_haplotype']].copy()

    print(f"  Computed Hamming distances for {len(df_combined)} haplotypes")
    print(f"  Selected haplotypes: distance range {df_selected['hamming_distance_nearest_library_strain'].min()}-{df_selected['hamming_distance_nearest_library_strain'].max()}")
    if len(df_nonselected) > 0:
        print(f"  Non-selected haplotypes: distance range {df_nonselected['hamming_distance_nearest_library_strain'].min()}-{df_nonselected['hamming_distance_nearest_library_strain'].max()}")

    # Read site annotations to identify epitope sites
    print(f"\nReading site annotations from {site_annotations_path}")
    df_sites = pd.read_csv(site_annotations_path, sep='\t')

    # Identify epitope columns (columns ending with "_epitope")
    epitope_cols = [col for col in df_sites.columns if col.endswith('_epitope')]

    if not epitope_cols:
        raise ValueError(
            f"No epitope columns found in site annotations file {site_annotations_path}.\n"
            f"Available columns: {df_sites.columns.tolist()}"
        )

    print(f"  Found {len(epitope_cols)} epitope sets: {epitope_cols}")

    # For each epitope set, extract the positions (0-indexed) where epitope = 1
    epitope_positions = {}
    for epitope_col in epitope_cols:
        # Get sequential_site values where this epitope is 1 (convert to 0-indexed for sequence access)
        positions = (df_sites[df_sites[epitope_col] == 1]['sequential_site'] - 1).tolist()
        epitope_positions[epitope_col] = positions
        print(f"    {epitope_col}: {len(positions)} sites")

    # Compute epitope-specific distances
    print(f"\nComputing epitope-specific distances to nearest library strain...")

    # For each epitope set, compute distances
    epitope_distance_cols = {}
    for epitope_col in epitope_cols:
        # Create distance column name (e.g., "Wolf_epitope" -> "Wolf_epitope_distance")
        distance_col_name = epitope_col.replace('_epitope', '_epitope_distance')
        epitope_distance_cols[epitope_col] = distance_col_name

        positions = epitope_positions[epitope_col]

        if not positions:
            # No epitope sites for this epitope set
            print(f"  Warning: {epitope_col} has no epitope sites, setting all distances to 0")
            df_combined[distance_col_name] = 0
            continue

        # Function to compute distance at epitope sites only
        def epitope_distance_at_positions(seq1, seq2, positions):
            """Compute Hamming distance at specific positions."""
            return sum(seq1[pos] != seq2[pos] for pos in positions)

        # Compute for each haplotype
        epitope_distances = []
        for _, row in df_combined.iterrows():
            sequence = row['representative_strain_ha_sequence']
            strain_name = row['representative_strain']

            # Find nearest library strain (excluding self)
            min_distance = float('inf')
            for lib_seq, lib_strain in zip(library_sequences, library_strains):
                if lib_strain == strain_name:
                    continue
                dist = epitope_distance_at_positions(sequence, lib_seq, positions)
                if dist < min_distance:
                    min_distance = dist

            if min_distance == float('inf'):
                # Only one library strain, can't compute distance
                raise ValueError(
                    f"Cannot compute epitope distance for {strain_name}: "
                    f"no other library strains available"
                )

            epitope_distances.append(min_distance)

        # Add column to dataframe
        df_combined[distance_col_name] = epitope_distances

        print(f"  {distance_col_name}: computed for {len(df_combined)} haplotypes")

    # Update the selected and non-selected dataframes
    df_selected = df_combined[df_combined['selected_haplotype']].copy()
    df_nonselected = df_combined[~df_combined['selected_haplotype']].copy()

    # Print epitope distance ranges
    for epitope_col in epitope_cols:
        distance_col_name = epitope_distance_cols[epitope_col]
        print(f"  Selected {distance_col_name}: range {df_selected[distance_col_name].min()}-{df_selected[distance_col_name].max()}")
        if len(df_nonselected) > 0:
            print(f"  Non-selected {distance_col_name}: range {df_nonselected[distance_col_name].min()}-{df_nonselected[distance_col_name].max()}")

    # Reorder columns: strain name, haplotype name, hamming distance, nearest library strain, epitope distances, HA sequence, then others
    # Build priority columns list dynamically to include epitope distance columns
    priority_cols = [
        'representative_strain',
        'derived_haplotype',
        'hamming_distance_nearest_library_strain',
        'nearest_library_strain',
    ]

    # Add epitope distance columns (sorted alphabetically for consistency)
    epitope_distance_col_names = sorted(epitope_distance_cols.values())
    priority_cols.extend(epitope_distance_col_names)

    # Add HA sequence column
    priority_cols.append('representative_strain_ha_sequence')

    # Verify all priority columns exist (they should by definition at this point)
    missing_priority_cols = [col for col in priority_cols if col not in df_combined.columns]
    if missing_priority_cols:
        raise ValueError(
            f"Expected columns are missing from combined dataframe: {missing_priority_cols}\n"
            f"Available columns: {list(df_combined.columns)}"
        )

    # Get remaining columns
    other_cols = [col for col in df_combined.columns if col not in priority_cols]

    # Reorder
    ordered_cols = priority_cols + other_cols

    df_selected = df_selected[ordered_cols]
    df_nonselected = df_nonselected[ordered_cols]

    # Sort the dataframes
    # Library (selected) strains: sorted by hamming distance smallest to largest
    df_selected = df_selected.sort_values('hamming_distance_nearest_library_strain', ascending=True)
    # Non-library (non-selected) strains: sorted by hamming distance largest to smallest
    df_nonselected = df_nonselected.sort_values('hamming_distance_nearest_library_strain', ascending=False)

    # Save outputs
    df_selected.to_csv(selected_output_path, sep='\t', index=False)
    df_nonselected.to_csv(nonselected_output_path, sep='\t', index=False)

    # Write library FASTA file (selected haplotypes only)
    with open(library_fasta_path, 'w') as f:
        for _, row in df_selected.iterrows():
            f.write(f">{row['representative_strain']}\n")
            f.write(f"{row['representative_strain_ha_sequence']}\n")

    print(f"\nOutputs:")
    print(f"  Selected haplotypes: {selected_output_path}")
    print(f"    {len(df_selected)} haplotypes")
    print(f"  Non-selected haplotypes: {nonselected_output_path}")
    print(f"    {len(df_nonselected)} haplotypes")
    print(f"  Library FASTA: {library_fasta_path}")
    print(f"    {len(df_selected)} sequences")

    print(f"\nCuration complete ✓")


if __name__ == '__main__':
    # Snakemake provides these variables
    sys.stdout = sys.stderr = open(snakemake.log[0], "w")
    curate_library_sequences(
        valid_haplotypes_path=snakemake.input.valid_haplotypes,
        additional_haplotypes_paths=snakemake.input.additional_haplotypes,
        site_annotations_path=snakemake.input.site_annotations,
        selected_output_path=snakemake.output.selected_haplotypes,
        nonselected_output_path=snakemake.output.nonselected_haplotypes,
        library_fasta_path=snakemake.output.library_fasta
    )
