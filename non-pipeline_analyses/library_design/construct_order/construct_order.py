import marimo

__generated_with = "0.17.2"
app = marimo.App(width="medium")


@app.cell
def _(mo):
    mo.md(
        r"""
    # Construct design notebook
    This notebook takes the output of `initial_design` (namely, a list of nucleotide sequences for both H3N2 and H1N1 strains) and generates a list of nucleotide inserts that can be submitted to Twist Biosciences for gene fragment synthesis and cloning. These inserts are designed to fit between the BsmBI-v2 and XbaI cut sites in the Bloom lab vector 2851. They should begin after the 19th codon of WSN upstream signal peptide, and continue all the way through the end of the HA coding region, followed by a double stop codon and a 16-nucleotide barcode.
    """
    )
    return


@app.cell
def _():
    import os
    import marimo as mo
    from pathlib import Path

    import sys
    import yaml
    import pandas as pd
    import Bio
    from Bio import SeqIO
    import random
    random.seed(42) # random number generator seed value
    return Path, mo, os, pd, random, sys, yaml


@app.cell
def _(Path, mo, os):
    # Marimo path to notebook
    notebook_directory: Path = mo.notebook_dir()

    # ID input and output directories
    datadir = notebook_directory / './data/'
    resultsdir = notebook_directory / './results/'
    os.makedirs(datadir, exist_ok=True)
    os.makedirs(resultsdir, exist_ok=True)
    return (notebook_directory,)


@app.cell
def _(
    barcode_index,
    n_barcodes,
    notebook_directory: "Path",
    nucleotides,
    os,
    pd,
    random,
):
    def design_inserts(library, subtype, 
                       insert_filepath, 
                       construct_filepath,
                       ectodomain_start,
                       ectodomain_length,
                       endodomain_sequence,
                       start_codon,
                       special_start_codons = None,
                       append_additional_upstream_sequence = '',
                       append_additional_downstream_sequence = '',
                       virus_id = 1, # start virus naming at this index
                       nucleotides = nucleotides
                      ):

        # Only design if the ordersheet hasn't been generated
        if os.path.exists(notebook_directory / construct_filepath):
            print(f"Already designed '{construct_filepath}', reading that file and NOT regenerating barcodes.")

        elif not os.path.exists(notebook_directory / construct_filepath): 
            print(f'Generating new barcodes at {construct_filepath}...')
            library_strains = pd.read_csv(notebook_directory / insert_filepath, sep='\t') # Input 
            start_codon = start_codon # Define the custom start codon to search for        
            virus_id = virus_id # Define ordersheet name parameters
            inserts = [] # Initialize empty ordersheet to populate with name, sequence
            subtype_specific_library_strains = library_strains.query(f'subtype=="{subtype}"') # Subtype specific strains
            if special_start_codons is None: # Get special start codons if they exist
                special_start_codons = {}

            for index, row in subtype_specific_library_strains.iterrows():  
                # Design barcode 
                i=1 # Initialize barcode counter
                for n in list(range(0,n_barcodes)):
                    for n in list(range(0,100)): # Try 100 times to make a barcode 
                        barcode = ''.join(random.choices(nucleotides, k=16))
                        if barcode[0:2] == 'gg': # Don't use barcodes that start with GG
                            continue
                        if barcode in barcode_index: # Don't use barcodes that have already been used in the library
                            continue
                        if n == 100:
                            print('something really rare happened, try resetting barcode_index')
                        else:
                            barcode_index.append(barcode)
                            break

                    # Get HA sequence
                    record = row['nt_sequence']
                    # Make a strain name with barcode info
                    name = row['strain']
                    name_barcoded = f'{library}_{subtype}_{virus_id}_bc{i}'
                    i+=1
                    # Get Genbank ID (and additional mutations) from FASTA header
                    genbank_id = row['accession_w_aa_muts_added']

                    # Find the position of the first instance of 'ATGAAG' or other custom start
                    codon_to_use = special_start_codons.get(name, start_codon)
                    start_position = record.find(codon_to_use)
                    assert start_position != -1, f"For {name} - no start codon {codon_to_use} found"

                    # Extract the sequence starting from the found position 
                    insert_start = start_position + ectodomain_start 
                    insert_end = start_position + (ectodomain_length) 
                    ectodomain_insert_seq = record[insert_start:insert_end]                
                    # Identify the endodomain region (subtype specific)
                    endodomain = endodomain_sequence

                    # Insert sequence we need to order is just ectodomain, endodomain, and barcode
                    insert_seq = ectodomain_insert_seq + endodomain + barcode

                    # Add upstream and downstream sequences if they exist
                    insert_seq = (append_additional_upstream_sequence or '') + insert_seq + (append_additional_downstream_sequence or '')

                    # Add to inserts list
                    inserts.append([name, genbank_id, name_barcoded, str(insert_seq)])     

                # Add to virus counter
                virus_id+=1

            inserts_df = pd.DataFrame(inserts, columns = ['strain', 'genbank', 'name', 'sequence'])
            inserts_df = inserts_df.to_csv(notebook_directory / construct_filepath, index=False) 

            return inserts
    return (design_inserts,)


@app.cell
def _(notebook_directory: "Path", pd, snakemake, sys, yaml):
    # Load configuration
    def load_config():
        config_path = notebook_directory / "config.yml"

        # If running within Snakemake, it defines a 'snakemake' object
        if "snakemake" in globals():
            config_path = snakemake.configfile
        elif len(sys.argv) > 1:
            # Allow manual override from command-line argument
            config_path = sys.argv[1]
        with open(config_path) as f:
            config = yaml.safe_load(f)
        print(f"Loaded config from: {config_path}")
        return config

    config = load_config()

    # Get barcode configuration
    # Define global variables for insert design from config
    n_barcodes = config['n_barcodes'] # Set the number of barcodes to design for
    nucleotides = config['nucleotides']

    # Load past barcodes from config
    barcode_index = []
    for key in config['past_barcodes_to_avoid']:
        key_barcodes = pd.read_csv(notebook_directory / config['past_barcodes_to_avoid'][key])['barcode'].tolist()
        barcode_index.extend(key_barcodes)
    return barcode_index, config, n_barcodes, nucleotides


@app.cell
def _(config, design_inserts, notebook_directory: "Path", pd):
    # Design all constructs configured in 'orders' key
    order_outputs = []

    for order in config['orders']:

        curr_order = config['orders'][order]

        design = design_inserts(
            library = curr_order['library'],
            subtype = curr_order['subtype'],
            insert_filepath = curr_order['input_file'],
            construct_filepath = curr_order['output_file'],
            ectodomain_start = curr_order['ectodomain_start'],
            ectodomain_length = curr_order['ectodomain_length'],
            endodomain_sequence = curr_order['endodomain_sequence'],
            append_additional_upstream_sequence = curr_order['append_additional_upstream_sequence'],
            start_codon = config['start_codon'],
            special_start_codons = config['special_start_codons']
        )

        print(f'There are inserts to order at {curr_order['output_file']}')

        # Add TSV to list of order outputs
        order_outputs.append(curr_order['output_file'])

    # Add barcode-to-strain output file in top-level directory that aggregatess across all ordersheets in output
    output_dir = notebook_directory / f'./{config['barcode_to_strain']}'
    order_outputs_df = pd.concat([pd.read_csv(f'{notebook_directory}/{f}', sep=',') for f in order_outputs], ignore_index=True)
    order_outputs_df['barcode'] = order_outputs_df['sequence'].str[-16:]
    order_outputs_df.to_csv(output_dir)
    print(f'\nSaved a final barcoce-to-strain map of all designed constructs to {output_dir}, use that for seqneut-pipeline input!')
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
