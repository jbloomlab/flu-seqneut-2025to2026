import marimo

__generated_with = "0.17.2"
app = marimo.App(width="medium")


@app.cell
def _():
    import os
    import marimo as mo
    from pathlib import Path

    import sys
    import yaml
    import pandas as pd
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    import pandas as pd
    import glob
    import re
    import numpy as np

    import warnings
    from Bio import BiopythonWarning
    warnings.simplefilter("ignore", BiopythonWarning) # Silence warnings about long plasmid names
    return FeatureLocation, Path, Seq, SeqFeature, SeqIO, SeqRecord, mo, os, pd


@app.cell
def _(Path, mo, os, pd):
    # Marimo path to notebook
    notebook_directory: Path = mo.notebook_dir()

    # ID library name
    library_name = 'flu-seqneut-2025to2026'

    # ID date
    date = '02-FEB-2026'

    # ID input and output directories
    designed_viral_library_df = pd.read_csv(notebook_directory / './flu-seqneut-2025to2026-barcode-to-strain-designed.csv')

    plasmiddir = notebook_directory / './plasmids/'
    os.makedirs(plasmiddir, exist_ok=True)
    return (
        date,
        designed_viral_library_df,
        library_name,
        notebook_directory,
        plasmiddir,
    )


@app.cell
def _(mo):
    mo.md(r"""Now, define a function to make plasmid maps. It should label the H1 and H3 maps slightly differently (given differences in design) and should include pertinent shareable info, like Genbank ID, viral barcode, and library name.""")
    return


@app.cell
def _():
    pHH21_derivative_backbone_upstream = 'actcttcctttttcaatattattgaagcatttatcagggttattgtctcatgagcggatacatatttgaatgtatttagaaaaataaacaaaagagtttgtagaaacgcaaaaaggccatccgtcaggatggccttctgcttaatttgatgcctggcagtttatggcgggcgtcctgcccgccaccctccgggccgttgcttcgcaacgttcaaatccgctcccggcggatttgtcctactcaggagagcgttcaccgacaaacaacagataaaacgaaaggcccagtctttcgactgagcctttcgttttatttgatgcctggcagttccctactctcgcatggggagaccccacactaccatcggcgctacggcgtttcacttctgagttcggcatggggtcaggtgggaccaccgcgctactgccgccaggcaaattctgttttatcagaccgcttctgcgttctgatttaatctgtatcaggctgaaaatttttttGcGGccgccaaaacagccaagctagcggccgatccccaaaaaaaaaaaaaaaaaaagagtccagagtggccccgccgctccgcgccggggggggggggggggggggacactttcggacatctggtcgacctccagcatcgggggaaaaaaaaaaacaaagtgtcgcccggagtactggtcgacctccgaagttgggggggagcaaaagcaggggaaaataaaaacaaccaaa'
    pHH21_derivative_backbone_downstream = 'agatcggaagagcgtcgtgtagggaaagagtgtgcggccgctatctactcaactgtcgccagttcactggtgctttaggtctccctgggggcaatcagtttctggatgtgttctaatgggtctttgcagtgcagaatatgcatctgagattaggatttcagaaatataaggaaaaacacccttgtttctactaataacccggcggcccaaaatgccgactcggagcgaaagatatacctcccccggggccgggaggtcgcgtcaccgaccacgccgccggcccaggcgacgcgcgacacggacacctgtccccaaaaacgccaccatcgcagccacacacggagcgcccggggccctctggtcaaccccaggacacacgcgggagcagcgccgggccggggacgccctcccggccgcccgtgccacacgcagggggccggcccgtgtctccagagcgggagccggaagcattttcggccggcccctcctacgaccgggacacacgagggaccgaaggccggccaggcgcgacctctcgggccgcacgcgcgctcagggagcgctctccgactccgcacggggactcgccagaaaggatcgtgatctgcattaatgaatcaggggataacgcaggaaagaacatgtgagcaaaaggccagcaaaaggccaggaaccgtaaaaaggccgcgttgctggcgtttttccataggctccgcccccctgacgagcatcacaaaaatcgacgctcaagtcagaggtggcgaaacccgacaggactataaagataccaggcgtttccccctggaagctccctcgtgcgctctcctgttccgaccctgccgcttaccggatacctgtccgcctttctcccttcgggaagcgtggcgctttctcatagctcacgctgtaggtatctcagttcggtgtaggtcgttcgctccaagctgggctgtgtgcacgaaccccccgttcagcccgaccgctgcgccttatccggtaactatcgtcttgagtccaacccggtaagacacgacttatcgccactggcagcagccactggtaacaggattagcagagcgaggtatgtaggcggtgctacagagttcttgaagtggtggcctaactacggctacactagaagaacagtatttggtatctgcgctctgctgaagccagttaccttcggaaaaagagttggtagctcttgatccggcaaacaaaccaccgctggtagcggtggtttttttgtttgcaagcagcagattacgcgcagaaaaaaaggatctcaagaagatcctttgatcttttctacggggtctgacgctcagtggaacgaaaactcacgttaagggattttggtcatgagattatcaaaaaggatcttcacctagatccttttaaattaaaaatgaagttttaaatcaatctaaagtatatatgagtaaacttggtctgacagttaccaatgcttaatcagtgaggcacctatctcagcgatctgtctatttcgttcatccatagttgcctgactccccgtcgtgtagataactacgatacgggagggcttaccatctggccccagtgctgcaatgataccgcgagacccacgctcaccggctccagatttatcagcaataaaccagccagccggaagggccgagcgcagaagtggtcctgcaactttatccgcctccatccagtctattaattgttgccgggaagctagagtaagtagttcgccagttaatagtttgcgcaacgttgttgccattgctacaggcatcgtggtgtcacgctcgtcgtttggtatggcttcattcagctccggttcccaacgatcaaggcgagttacatgatcccccatgttgtgcaaaaaagcggttagctccttcggtcctccgatcgttgtcagaagtaagttggccgcagtgttatcactcatggttatggcagcactgcataattctcttactgtcatgccatccgtaagatgcttttctgtgactggtgagtactcaaccaagtcattctgagaatagtgtatgcggcgaccgagttgctcttgcccggcgtcaacacgggataataccgcgccacatagcagaactttaaaagtgctcatcattggaaaacgttcttcggggcgaaaactctcaaggatcttaccgctgttgagatccagttcgatgtaacccactcgtgcacccaactgatcttcagcatcttttactttcaccagcgtttctgggtgagcaaaaacaggaaggcaaaatgccgcaaaaaagggaataagggcgacacggaaatgttgaatactcat'
    return (
        pHH21_derivative_backbone_downstream,
        pHH21_derivative_backbone_upstream,
    )


@app.cell
def _(
    FeatureLocation,
    Seq,
    SeqFeature,
    SeqRecord,
    date,
    library_name,
    pHH21_derivative_backbone_downstream,
    pHH21_derivative_backbone_upstream,
):
    def write_genbank(
        strain_name,
        accession,
        plasmid_name, # the name assigned in the Bloom lab plasmid log
        subtype, 
        ectodomain,
        barcode='NNNNNNNNNNNNNNNN',
        constant_plasmid_sequence_upstream=pHH21_derivative_backbone_upstream,
        constant_plasmid_sequence_downstream=pHH21_derivative_backbone_downstream):


        definition = (
            f"This pHH plasmid contains the HA ectodomain sequence for a {subtype} variant {strain_name}. " +
            f"This plasmid was generated for the {library_name} library. " +
            f"Signal peptide and 3'NCR from WSN, ectodomain from {strain_name} HA with Genbank accession {accession}, and last 46 aa recoded WSN transmembrane and c-terminal domain. " +
            f"With duplicated 5' packaging signals from WSN with a single stop codon in the duplicated packaging signal, with the 16-nucleotide barcode {barcode}. " +
            "This plasmid was cloned and sequence confirmed by Twist, and designed and logged by Caroline Kikawa"
        )

        # Depending on the subtype, the upstream signal peptide and downstream endododomain and CT domain will vary
        if subtype=='H3N2':
            signal_peptide_nc = 'atgaaggcaaaactactggtcctgttatatgcatttgtagctacagatgcagacaca' # first 19 aa
            endodomain_nc = 'atcaagggagttgagctgaagtcaggatacaaagattggatcctatggatttcctttgccATGtcTtgCttCCtActGtgCgtAgcACtACtAggCttTatTatgtgggcGtgTcaGaaA'
            c_term_nc = 'ggCtcCCtAcaAtgTCgGatTtgTatTTAATAG'
        elif subtype=='H1N1':
            signal_peptide_nc = 'atgaaggcaaaactactggtcctgttatatgcatttgtagctacagatgcagacaca'
            endodomain_nc = 'aaattggaatcaatgggagtgtatcagattctggcgatatattctacagtggcaagctccttagtactgctagtttctttaggagcgattagcttttggatgtgctccaacggctccctacaatgtcggatttgtatttaatag'
            c_term_nc = '' # included in endo sequence above

        # Build the expected plasmid sequence
        plasmid_sequence = (
            constant_plasmid_sequence_upstream +
            signal_peptide_nc + 
            ectodomain +
            # endodomain_nc + 
            # c_term_nc + 
            # barcode +
            constant_plasmid_sequence_downstream
        )

        # Write sequence features
        f1 = SeqFeature(FeatureLocation(534, 700, -1), type="terminator", qualifiers = {'label': 'mouse PolI terminator'})
        f2 = SeqFeature(FeatureLocation(700, 712, -1), type="misc_feature", qualifiers = {'label': 'U12'})
        f3 = SeqFeature(FeatureLocation(700, 732, -1), type="misc_feature", qualifiers = {'label': "3' NCR"})

        if subtype=='H3N2':
            f4 = SeqFeature(FeatureLocation(732, 789, +1), type="misc_feature", qualifiers = {'label': 'WSN first 19 aa'})
            f5 = SeqFeature(FeatureLocation(789, 2292, +1), type="misc_feature", qualifiers = {'label': f'HA ectodomain from {strain_name}'})
            f6 = SeqFeature(FeatureLocation(2292, 2412, +1), type="misc_feature", qualifiers = {'label': 'consensus H3 endodomain'})
            f8 = SeqFeature(FeatureLocation(2412, 2445, +1), type="misc_feature", qualifiers = {'label': 'WSN recoded CT'})
            f9 = SeqFeature(FeatureLocation(2445, 2461, +1), type="misc_feature", qualifiers = {'label': 'barcode'})
            f10 = SeqFeature(FeatureLocation(2461, 2494, +1), type="misc_feature", qualifiers = {'label': 'Illumina Read1'})
            f11 = SeqFeature(FeatureLocation(2503, 2608, +1), type="misc_feature", qualifiers = {'label': 'WSN packaging signal'})
            f12 = SeqFeature(FeatureLocation(2608, 2653, -1), type="misc_feature", qualifiers = {'label': "5' NCR"})
            f13 = SeqFeature(FeatureLocation(2641, 2653, -1), type="misc_feature", qualifiers = {'label': 'U13'})
            f14 = SeqFeature(FeatureLocation(2653, 3056, -1), type="misc_feature", qualifiers = {'label': 'Human PolI promoter'})
            f15 = SeqFeature(FeatureLocation(3923, 4781, -1), type="CDS", qualifiers = {'label': 'AmpR'})

        elif subtype=='H1N1':
            f4 = SeqFeature(FeatureLocation(732, 792, +1), type="misc_feature", qualifiers = {'label': 'WSN first 20 aa'})
            f5 = SeqFeature(FeatureLocation(792, 2292, +1), type="misc_feature", qualifiers = {'label': f'HA gene from {strain_name}'})
            f6 = SeqFeature(FeatureLocation(2292, 2436, +1), type="misc_feature", qualifiers = {'label': 'WSN endodomain'})
            f8 = SeqFeature(FeatureLocation(2403, 2436, +1), type="misc_feature", qualifiers = {'label': 'WSN recoded CT'})
            f9 = SeqFeature(FeatureLocation(2436, 2452, +1), type="misc_feature", qualifiers = {'label': 'barcode'})
            f10 = SeqFeature(FeatureLocation(2452, 2485, +1), type="misc_feature", qualifiers = {'label': 'Illumina Read1'})
            f11 = SeqFeature(FeatureLocation(2494, 2599, +1), type="misc_feature", qualifiers = {'label': 'WSN packaging signal'})
            f12 = SeqFeature(FeatureLocation(2599, 2644, -1), type="misc_feature", qualifiers = {'label': "5' NCR"})
            f13 = SeqFeature(FeatureLocation(2632, 2644, -1), type="misc_feature", qualifiers = {'label': 'U13'})
            f14 = SeqFeature(FeatureLocation(2644, 3047, -1), type="misc_feature", qualifiers = {'label': 'Human PolI promoter'})
            f15 = SeqFeature(FeatureLocation(3914, 4772, -1), type="CDS", qualifiers = {'label': 'AmpR'})

        features_list = [f1,f2,f3,f4,f5,f6,f8,f9,f10,f11,f12,f13,f14,f15]

        # Write sequence record and save 
        record = SeqRecord(Seq(plasmid_sequence), 
                           id = '.', 
                           name = plasmid_name, description = definition, 
                           features = features_list, 
                           annotations = {'source': 'synthetic DNA construct',
                                          'organism': 'synthetic DNA construct',
                                          'molecule_type': 'ds-DNA',
                                          'topology': 'circular',
                                          'date': f'{date}',})

        return(plasmid_sequence, definition, record)   
    return (write_genbank,)


@app.cell
def _(designed_viral_library_df):
    # We only need to log the newly ordered strains
    circulating_designed_viral_library_df = designed_viral_library_df.query('strain_type=="circulating_2025to2026"')
    return (circulating_designed_viral_library_df,)


@app.cell
def _(
    SeqIO,
    circulating_designed_viral_library_df,
    notebook_directory: "Path",
    os,
    pd,
    plasmiddir,
    write_genbank,
):
    # Write records for all plasmids

    # Initialize mapping list so we can map strains to Bloom lab plasmid log IDs
    strain_to_plasmid_log_map = []

    # First make an index of plasmid IDs for plate-based logging in Bloom lab
    plasmid_index = []
    # 96-well plate: 8 rows (A-H) × 12 columns (1-12)
    rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    columns = range(1, 13)  # 1 to 12
    plates = range(5, 8)    # plates 5, 6, 7

    # Generate the list
    # First 108 indexes: plate 5 and part of plate 6
    plate = 5
    count = 0
    for p in [5, 6]:
        for col in columns:
            for row in rows:
                plasmid_index.append(f"plate{p}-{row}{col}")
                count += 1
                if count == 108:
                    break
            if count == 108:
                break
        if count == 108:
            break

    # Next 96 indexes: plate 7 (full plate)
    plate = 7
    for col in columns:
        for row in rows:
            plasmid_index.append(f"plate{plate}-{row}{col}")

    # Final 8 wells: plate 8, A1-H1
    plate = 8
    for row in rows:
        plasmid_index.append(f"plate{plate}-{row}1")


    later_names = [
        'flu-seqneut-25to26_H3N2_55_bc1',
        'flu-seqneut-25to26_H3N2_55_bc2',
        'flu-seqneut-25to26_H3N2_56_bc1',
        'flu-seqneut-25to26_H3N2_56_bc2',
        'flu-seqneut-25to26_H3N2_57_bc1',
        'flu-seqneut-25to26_H3N2_57_bc2',
        'flu-seqneut-25to26_H3N2_58_bc1',
        'flu-seqneut-25to26_H3N2_58_bc2',
                  ]

    # Index for loop
    i=0
    # Generate maps
    for name in circulating_designed_viral_library_df.name.unique():

        temp_df = circulating_designed_viral_library_df.query(f'name == "{name}"')

        # Skip 4 later strains that are added to a 4th plate
        if name in later_names:
            pass

        barcode_id = temp_df.name.values[0].split('bc')[1]
        subtype = temp_df.subtype.values[0]
        strain_name = temp_df.strain.values[0].replace(f'_{subtype}', '') 
        plasmid_name = plasmid_index[i] + '_' + strain_name + '_bc' + barcode_id
        plasmid_name = plasmid_name.replace('/', '_')


        map, defintion, record = (write_genbank(
            strain_name = strain_name,
            accession = temp_df.accession.values[0],
            plasmid_name = plasmid_name,
            subtype = subtype,
            barcode = temp_df.barcode.values[0],
            ectodomain = temp_df.nt_sequence_HA_ectodomain.values[0]
            )
                                 )

        outfile = os.path.join(plasmiddir, f'{plasmid_name}.gb')
        with open(outfile, 'w') as f:
            SeqIO.write(record, f, 'genbank')

        strain_to_plasmid_log_map.append([name, f'{plasmid_name}.gb'])

        i+=1

    # Reindex for last 8 plate slots
    n=204

    for name in later_names:

        temp_df = circulating_designed_viral_library_df.query(f'name == "{name}"')

        barcode_id = temp_df.name.values[0].split('bc')[1]
        subtype = temp_df.subtype.values[0]
        strain_name = temp_df.strain.values[0].replace(f'_{subtype}', '') 
        plasmid_name = plasmid_index[n] + '_' + strain_name + '_bc' + barcode_id
        plasmid_name = plasmid_name.replace('/', '_')

        map, defintion, record = (write_genbank(
            strain_name = strain_name,
            accession = temp_df.accession.values[0],
            plasmid_name = plasmid_name,
            subtype = subtype,
            barcode = temp_df.barcode.values[0],
            ectodomain = temp_df.nt_sequence_HA_ectodomain.values[0]
            )
                                 )

        outfile = os.path.join(plasmiddir, f'{plasmid_name}.gb')
        with open(outfile, 'w') as f:
            SeqIO.write(record, f, 'genbank')

        n+=1

    strain_to_plasmid_log_map_df = pd.DataFrame(strain_to_plasmid_log_map, columns=['strain', 'plasmid_log_ID'])
    outfile = notebook_directory / 'flu-seqneut-2025to2026-strain-to-plasmid-log.csv'
    print(f'saving strain to plasmid log ID map to {outfile}...')
    strain_to_plasmid_log_map_df.to_csv(outfile, index=False)
    return


if __name__ == "__main__":
    app.run()
