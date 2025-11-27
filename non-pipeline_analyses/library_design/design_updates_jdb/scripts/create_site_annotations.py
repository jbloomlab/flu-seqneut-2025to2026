"""
Parse GFF and epitope JSON files to create site annotations.

This script:
1. Parses GFF file to map sequential positions to protein/protein_site
2. Parses epitope JSON files to identify epitope sites
3. Creates a TSV with sequential_site, protein, protein_site, and epitope columns
"""

import sys
import json
import pandas as pd
from pathlib import Path


def parse_gff(gff_path):
    """
    Parse GFF file to extract coordinate mappings.

    Converts nucleotide positions to amino acid positions and assigns
    sequential numbering starting from 1.

    Parameters
    ----------
    gff_path : str
        Path to GFF file

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: sequential_site, protein, protein_site
    """
    # First pass: collect all features and sort by start position
    features = []
    feature_proteins = set()

    with open(gff_path, 'r') as f:
        for line in f:
            # Skip comments and empty lines
            if line.startswith('#') or not line.strip():
                continue

            # GFF format: seqname, source, feature, start, end, score, strand, frame, attribute
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            start_nt = int(fields[3])  # Nucleotide start (1-indexed)
            end_nt = int(fields[4])    # Nucleotide end (1-indexed)
            attributes = fields[8]

            # Parse attributes to get gene/product name
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    # Strip whitespace and quotes from value
                    value = value.strip().strip('"').strip("'")
                    attr_dict[key.strip()] = value

            # Extract protein name (gene_name, gene, or product)
            protein = attr_dict.get('gene_name', attr_dict.get('gene', attr_dict.get('product', feature_type)))

            if protein in feature_proteins:
                raise ValueError(f"Duplicate {protein=} in {gff_path=}")
            feature_proteins.add(protein)

            features.append({
                'protein': protein,
                'start_nt': start_nt,
                'end_nt': end_nt
            })

    # Sort features by start position to ensure consistent sequential numbering
    features.sort(key=lambda x: x['start_nt'])

    # Second pass: create amino acid mappings
    annotations = []
    current_sequential_site = 1

    for feature in features:
        protein = feature['protein']
        start_nt = feature['start_nt']
        end_nt = feature['end_nt']

        # Convert nucleotide range to amino acids
        # Each amino acid is 3 nucleotides (a codon)
        length_nt = end_nt - start_nt + 1

        if length_nt % 3 != 0:
            raise ValueError(
                f"Feature {protein} has nucleotide length {length_nt} which is not divisible by 3!\n"
                f"Nucleotide range: {start_nt}-{end_nt}"
            )

        num_amino_acids = length_nt // 3

        # Create entry for each amino acid in this feature
        for protein_site in range(1, num_amino_acids + 1):
            annotations.append({
                'sequential_site': current_sequential_site,
                'protein': protein,
                'protein_site': protein_site
            })
            current_sequential_site += 1

    df = pd.DataFrame(annotations)

    if len(df) > 0:
        # Check for duplicate sequential_site values (overlapping features)
        duplicates = df[df.duplicated(subset='sequential_site', keep=False)]
        if len(duplicates) > 0:
            raise ValueError(
                f"Found overlapping features in GFF file!\n"
                f"Sequential sites with multiple annotations:\n"
                f"{duplicates.sort_values('sequential_site').head(10)}\n"
                f"Total duplicates: {len(duplicates)}"
            )

        df = df.sort_values('sequential_site').reset_index(drop=True)

        # Check that sequential sites are continuous from 1 to N
        sequential_sites = set(df['sequential_site'])
        expected_sites = set(range(1, len(df) + 1))

        if sequential_sites != expected_sites:
            missing = expected_sites - sequential_sites
            extra = sequential_sites - expected_sites
            error_msg = "Sequential sites are not continuous from 1 to N!\n"
            if missing:
                error_msg += f"  Missing sites: {sorted(missing)[:20]}\n"
            if extra:
                error_msg += f"  Sites outside expected range: {sorted(extra)[:20]}\n"
            error_msg += f"  Expected range: 1 to {len(df)}\n"
            error_msg += f"  Actual range: {min(sequential_sites)} to {max(sequential_sites)}"
            raise ValueError(error_msg)

    return df


def parse_epitope_json(json_path):
    """
    Parse epitope JSON file to extract epitope sites as (protein, protein_site) tuples.

    Parameters
    ----------
    json_path : str
        Path to epitope JSON file

    Returns
    -------
    list of tuples
        List of (protein, protein_site) tuples for epitope sites
    """
    with open(json_path, 'r') as f:
        data = json.load(f)

    epitope_sites = []

    # Expected structure: {"map": {"HA1": {"145": 1, "155": 1, ...}}}
    if 'map' not in data:
        raise ValueError(f"Epitope JSON {json_path} missing 'map' key. Found keys: {list(data.keys())}")

    for protein, positions in data['map'].items():
        if not isinstance(positions, dict):
            raise ValueError(
                f"Expected positions to be a dict for protein {protein}, got {type(positions)}"
            )

        for protein_site_str, value in positions.items():
            # protein_site is a string like "145", convert to int
            protein_site = int(protein_site_str)
            epitope_sites.append((protein, protein_site))

    return epitope_sites


def create_site_annotations(gff_path, epitope_json_paths, output_path, subtype):
    """
    Create site annotations TSV combining GFF and epitope data.

    Parameters
    ----------
    gff_path : str
        Path to GFF file
    epitope_json_paths : list of str
        Paths to epitope JSON files
    output_path : str
        Path to output TSV file
    subtype : str
        Subtype name (e.g., 'H3N2', 'H1N1') to strip from epitope names
    """
    # Parse GFF file
    print(f"Parsing GFF file: {gff_path}")
    df = parse_gff(gff_path)

    if len(df) == 0:
        raise ValueError(f"No annotations found in GFF file: {gff_path}")

    print(f"  Found {len(df)} sequential positions")
    print(f"  Proteins: {df['protein'].unique().tolist()}")

    # Parse each epitope JSON and add as column
    for epitope_json_path in epitope_json_paths:
        # Extract epitope map name from filename
        # Filename format: {subtype}_{epitope_map}.json (e.g., H3N2_Wolf.json)
        epitope_name = Path(epitope_json_path).stem

        # Strip subtype prefix if present (e.g., "H3N2_Wolf" -> "Wolf")
        if epitope_name.startswith(f"{subtype}_"):
            epitope_name = epitope_name[len(subtype) + 1:]  # +1 for the underscore

        column_name = f"{epitope_name}_epitope"

        print(f"Parsing epitope file: {epitope_json_path}")
        epitope_sites = parse_epitope_json(epitope_json_path)
        print(f"  Found {len(epitope_sites)} epitope sites (protein, protein_site pairs)")

        # Create a set of (protein, protein_site) tuples for matching
        epitope_set = set(epitope_sites)

        # Check that all epitope sites are found in the GFF data
        df_sites = set(zip(df['protein'], df['protein_site']))
        missing_sites = epitope_set - df_sites
        if missing_sites:
            raise ValueError(
                f"Epitope sites not found in GFF data for {epitope_name}!\n"
                f"Missing (protein, protein_site) pairs: {sorted(missing_sites)[:10]}\n"
                f"Total missing: {len(missing_sites)}"
            )

        # Add column indicating if (protein, protein_site) is in epitope
        df[column_name] = df.apply(
            lambda row: 1 if (row['protein'], row['protein_site']) in epitope_set else 0,
            axis=1
        )

        print(f"  Matched {df[column_name].sum()} sites")

    # Write output
    df.to_csv(output_path, sep='\t', index=False)
    print(f"\nWrote site annotations to: {output_path}")
    print(f"Columns: {df.columns.tolist()}")
    print(f"Shape: {df.shape}")


if __name__ == '__main__':
    # Snakemake provides these variables
    sys.stdout = sys.stderr = open(snakemake.log[0], "w")
    create_site_annotations(
        gff_path=snakemake.input.gff,
        epitope_json_paths=snakemake.input.epitope_jsons,
        output_path=snakemake.output.tsv,
        subtype=snakemake.wildcards.subtype
    )
