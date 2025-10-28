import marimo

__generated_with = "0.17.2"
app = marimo.App(width="medium")


@app.cell
def _(mo):
    mo.md(
        r"""
    # Construct design notebook
    This notebook takes the output of `initial_design` -- namely a list of nucleotide sequences for both H3N2 and H1N1 strains -- and generates a list of nucleotide inserts that can be submitted to Twist Biosciences for gene fragment synthesis and cloning. 

    [some more details about the constructs]
    """
    )
    return


@app.cell
def _():
    import marimo as mo
    from pathlib import Path
    import pandas as pd
    import Bio
    from Bio import SeqIO
    import random
    return Path, SeqIO, mo, pd, random


@app.cell
def _(Path, mo, pd):
    ###### Get strains ######
    # marimo path to notebook
    notebook_directory: Path = mo.notebook_dir()
    library_strains_filepath = notebook_directory / "../initial_design/results/aggregated_library_strains/library_strains.tsv"

    # load library strains
    library_strains = pd.read_csv(library_strains_filepath, sep='\t')

    ##### Get past barcodes #####
    barcode = 'NNNNNNNNNNNNNNNN'
    nucleotides = ['a', 'c', 'g', 't']

    # Initialize barcode index (list of barcodes that have been used by a construct already)
    barcode_index = []

    # Add barcodes already used by prior viral libraries
    kikawa2023 = pd.read_csv(notebook_directory / 'data/2023_H3N2_Kikawa.csv')['barcode'].tolist()
    loes2023 = pd.read_csv(notebook_directory / 'data/pdmH1N1_lib2023_loes.csv')['barcode'].tolist()
    kikawa2025 = pd.read_csv(notebook_directory / 'data/flu-seqneut-2025-barcode-to-strain_designed.csv')['barcode'].tolist()
    barcode_index.extend(loes2023)
    barcode_index.extend(kikawa2023)
    barcode_index.extend(kikawa2025)
    return barcode_index, nucleotides


@app.cell
def _():
    ###### Get flanking sequences ######

    # First 19 amino acids of WSN
    wsn_upstream_19aa = 'atgaaggcaaaactactggtcctgttatatgcatttgtagctacagatgcagacaca'

    # First 20 amino acids of WSN
    wsn_upstream_20aa = 'atgaaggcaaaactactggtcctgttatatgcatttgtagctacagatgcagacacaata'

    # H1 endodomain, representing amino acids 521 to 568 (48 amino acids)
    h1_last_46aa_from_WSN = 'aaattggaatcaatgggagtgtatcagattctggcgatatattctacagtggcaagctccttagtactgctagtttctttaggagcgattagcttttggatgtgctccaacggCtcCCtAcaAtgTCgGatTtgTatTTAATAG'

    # H3 endodomain, representing amino acids 521 to 560 (40 amino acids)
    h3_endodomain = 'atcaagggagttgagctgaagtcaggatacaaagattggatcctatggatttcctttgccATGtcTtgCttCCtActGtgCgtAgcACtACtAggCttTatTatgtgggcGtgTcaGaaA'
    # Downstream mutated WSN packaging signal with double stop codon (11 amino acids)
    wsn_downstream = 'ggCtcCCtAcaAtgTCgGatTtgTatTTAATAG'
    return (
        h1_last_46aa_from_WSN,
        h3_endodomain,
        wsn_downstream,
        wsn_upstream_19aa,
        wsn_upstream_20aa,
    )


@app.cell
def _(SeqIO, barcode_index, nucleotides, os, pd, random):
    def design_inserts(subtype, 
                       insert_filepath, 
                       library_nucleotide_sequences,
                       upstream_signalpep,
                       ectodomain_start,
                       ectodomain_length,
                       endodomain_sequence,
                       start_codon = "ATGAAG",
                       append_additional_upstream_sequence = '',
                       append_additional_downstream_sequence = '',
                       virus_id = 1 # start virus naming at this index
                      ):

        # Only design if the ordersheet hasn't been generated
        if os.path.exists(insert_filepath):
            print(f"File '{insert_filepath}' exists, reading that file and NOT regenerating barcodes.")
        else:
            # Input FASTA file of subtype nucleotide sequences
            fasta_file = library_nucleotide_sequences
        
            # Define the custom start codon to search for
            start_codon = start_codon
    
            # Define ordersheet name parameters
            virus_id = virus_id
    
            # Initialize empty ordersheet to populate with name, sequence
            inserts = []
    
            # Set the number of barcodes to design for
            n_barcodes = 2
        
            # Open FASTA file and design constructs for each entry
            with open(fasta_file, "r") as handle:
        
                for record in SeqIO.parse(handle, "fasta"):
                    # Initialize barcode counter
                    i=1
                
                    # Each strain needs barcodes designed 
                    for n in list(range(0,n_barcodes)):
        
                        # Find the position of the first instance of 'ATGAAG'
                        start_position = record.seq.find(start_codon)                
                        assert start_position == -1, f"For {record.id} - no start codon {start_codon} found"

                        # Extract the sequence starting from the found position 
                        insert_start = start_position + ectodomain_start # Insert will start after first 19 amino acids of WSN
                        insert_end = start_position + (ectodomain_length) # Insert stops at amino acid 521
                        ectodomain_insert_seq = record.seq[insert_start:insert_end]
                        ectodomain_insert_translated_seq = ectodomain_insert_seq.translate()                
                        # Identify the endodomain region (subtype specific)
                        endodomain = endodomain_sequence
                                    
                        # Generate a barcode
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
            
                        # Expected sequence, including fixed upstream WSN signal peptide
                        expected_seq = upstream_signalpep + ectodomain_insert_seq + endodomain + barcode
                        # Insert sequence we need to order
                        # Ectodomain, endodomain, and barcode
                        insert_seq = ectodomain_insert_seq + endodomain + barcode

                        # Add upstream sequence
                        if append_additional_upstream_sequence == '':
                            pass
                        else:
                            insert_seq = append_additional_upstream_sequence + insert_seq
                        # Add downstream sequence
                        if append_additional_downstream_sequence == '':
                            pass
                        else:
                            insert_seq = insert_seq + append_additional_downstream_sequence
            
                        # Make a strain name with barcode info
                        name = record.id
                        name_barcoded = f'{subtype}_{virus_id}_bc{i}'
                        i+=1

                        # Get Genbank ID (and additional mutations) from FASTA header
                        genbank_id = record.description[len(record.id):].strip('protein identical to ')
    
                        # Add to inserts list
                        inserts.append([name, genbank_id, name_barcoded, str(insert_seq)])     

                    # Add to virus counter
                    virus_id+=1
    
                inserts_df = pd.DataFrame(inserts, columns = ['strain', 'genbank', 'name', 'sequence'])
                # inserts_df = inserts_df.sort_values(by = 'name').reset_index(drop=True).to_csv(insert_filepath, index=False)
                inserts_df = inserts_df.to_csv(insert_filepath, index=False)
        

    return (design_inserts,)


@app.cell
def _(mo):
    mo.md(
        r"""
    ## Design H1 inserts

    These sequences should be XXXX nucleotides long
    """
    )
    return


@app.cell
def _(
    design_inserts,
    h1_last_46aa_from_WSN,
    ordersheetsdir,
    os,
    pd,
    wsn_upstream_20aa,
):
    design_inserts(
        subtype = 'H1N1',
        insert_filepath = os.path.join(ordersheetsdir, 'h1_inserts.csv'),
        library_nucleotide_sequences = '../results/strains_for_library/h1_nt_seqs_for_library.fasta',
        upstream_signalpep = wsn_upstream_20aa,
        ectodomain_start = 20*3,
        ectodomain_length = 520*3,
        endodomain_sequence = h1_last_46aa_from_WSN,
        append_additional_upstream_sequence = 'catttgtagctacagatgcagacaca' + 'ata', # Overlap with WSN signal peptide
        append_additional_downstream_sequence = 'AGATCGGAAGAGCGTCGTGT', # Overlap with Illumina R1 priming sequence             
    )

    h1_inserts_df = pd.read_csv(os.path.join(ordersheetsdir, 'h1_inserts.csv'))
    h1_inserts_df
    return


@app.cell
def _(mo):
    mo.md(r"""## Design H3 inserts""")
    return


@app.cell
def _(
    design_inserts,
    h3_endodomain,
    ordersheetsdir,
    os,
    pd,
    wsn_downstream,
    wsn_upstream_19aa,
):
    design_inserts(
        subtype = 'H3N2',
        insert_filepath = os.path.join(ordersheetsdir, 'h3_inserts.csv'),
        library_nucleotide_sequences = '../results/strains_for_library/h3_nt_seqs_for_library.fasta',
        upstream_signalpep = wsn_upstream_19aa,
        ectodomain_start = 16*3,
        ectodomain_length = 517*3,
        endodomain_sequence = h3_endodomain + wsn_downstream,
        append_additional_upstream_sequence = 'catttgtagctacagatgcagacaca', # Overlap with WSN signal peptide
        append_additional_downstream_sequence = 'AGATCGGAAGAGCGTCGTGT', # Overlap with Illumina R1 priming sequence                    
    )

    h3_inserts_df = pd.read_csv(os.path.join(ordersheetsdir, 'h3_inserts.csv'))
    h3_inserts_df
    return


if __name__ == "__main__":
    app.run()
