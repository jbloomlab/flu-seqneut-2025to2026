"""Match proteins to Genbank nucleotide coding sequences.

For command-line options, run::

    python match_prot_to_genbank_nt.py --help

"""


import argparse
import io
import os
import shutil
import subprocess
import textwrap
import zipfile

import Bio.Blast.NCBIXML
import Bio.Data.CodonTable
import Bio.SeqIO

import pandas as pd


CODON_TO_AA = Bio.Data.CodonTable.unambiguous_dna_by_name["Standard"].forward_table

# map each amino acid to its most used codon in human genes
AA_TO_BEST_HUMAN_CODON = {
    'F': 'TTC',  # Phenylalanine
    'L': 'CTG',  # Leucine
    'I': 'ATC',  # Isoleucine
    'M': 'ATG',  # Methionine (only codon)
    'V': 'GTG',  # Valine
    'S': 'AGC',  # Serine
    'P': 'CCC',  # Proline
    'T': 'ACC',  # Threonine
    'A': 'GCC',  # Alanine
    'Y': 'TAC',  # Tyrosine
    'H': 'CAC',  # Histidine
    'Q': 'CAG',  # Glutamine
    'N': 'AAC',  # Asparagine
    'K': 'AAG',  # Lysine
    'D': 'GAC',  # Aspartic acid
    'E': 'GAG',  # Glutamic acid
    'C': 'TGC',  # Cysteine
    'W': 'TGG',  # Tryptophan (only codon)
    'R': 'CGC',  # Arginine
    'G': 'GGC',  # Glycine
    '*': 'TGA'   # Stop codon (TGA slightly more frequent than TAA or TAG)
}


def get_target_fasta(taxon, targets_zip, targets_fasta):
    """Get targets from NCBI datasets, filtering sequences with ambiguous nucleotides."""
    cmd = [
        'datasets', 'download', 'virus', 'genome', 'taxon', taxon,
        '--include', 'genome',
        '--filename', os.path.basename(targets_zip)
    ]

    dirname = os.path.dirname(targets_zip)
    assert os.path.dirname(targets_zip) == os.path.dirname(targets_fasta)

    print(
        f"Running this command in {dirname=} to create {targets_zip=}:\n"
        + " ".join(cmd)
    )
    subprocess.run(cmd, check=True, cwd=dirname)

    # Extract raw FASTA first
    raw_fasta = targets_fasta + ".raw"
    print(f"Extracting genomic.fna from {targets_zip} as {raw_fasta}")
    with zipfile.ZipFile(targets_zip, 'r') as zip_ref:
        for member in zip_ref.namelist():
            if member.endswith('genomic.fna'):
                with zip_ref.open(member) as f_in, open(raw_fasta, 'wb') as f_out:
                    f_out.write(f_in.read())
                break
        else:
            raise ValueError("genomic.fna not found in archive.")

    # Filter with seqkit to remove sequences with ambiguous nucleotides
    print(f"Filtering sequences with ambiguous nucleotides using seqkit to create {targets_fasta=}")
    seqkit_cmd = [
        'seqkit', 'grep',
        '-s', '-r', '-p', '^[ATGCatgc]+$',
        raw_fasta
    ]
    with open(targets_fasta, 'w') as f_out:
        result = subprocess.run(seqkit_cmd, stdout=f_out, stderr=subprocess.PIPE, text=True, check=True)

    # Report filtering stats if available in stderr
    if result.stderr:
        print(result.stderr)

    # Clean up raw file
    os.remove(raw_fasta)
        
        
def make_blastdb(targets_fasta, blastdb):
    """Make BLAST database."""
    cmd = [
        "makeblastdb",
        "-in", targets_fasta,
        "-dbtype", "nucl",
        "-out", blastdb,
        "-parse_seqids",
    ]
    print(f"Running this command to create {blastdb=}:\n{' '.join(cmd)}")
    subprocess.run(cmd, check=True)
        
        
def run_tblastn(query_prots, blastdb, blast_results, evalue, max_target_seqs, threads):
    """Run ``tblastn``"""
    cmd = [
        "tblastn",
        "-query", query_prots,
        "-db", blastdb,
        "-out", blast_results,
        "-outfmt", "5",
        "-evalue", str(evalue),
        "-max_target_seqs", str(max_target_seqs),
        "-num_threads", str(threads),
        "-seg", "no",
    ]
    print(f"Running this command to create {blast_results=}:\n{' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def parse_blast_xml(query_prots, blast_results, blastdb, query_id_name):
    """Parse BLAST results from XML, keeping only full-length gap-free alignments."""
    # parse BLAST result to get highest identity full-length gap-free alignment
    all_hits = {}
    with open(blast_results) as handle:
        for record in Bio.Blast.NCBIXML.parse(handle):
            query_id = record.query
            best_hit = None
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    if (
                        hsp.query_start == 1 and
                        hsp.query_end == record.query_letters and
                        hsp.gaps == 0
                    ):
                        # full-length gap-free alignment
                        prot_identity = hsp.identities / record.query_letters
                        aa_muts = [
                            f"{ha}{r + 1}{qa}"
                            for (r, (qa, ha)) in enumerate(zip(hsp.query, hsp.sbjct))
                            if qa != ha
                        ]
                        if (
                            (best_hit is None)
                            or (prot_identity > best_hit["prot_identity"])
                        ):
                            best_hit = {
                                "accession": alignment.accession,
                                "prot_identity": prot_identity,
                                "aa_muts_accession_to_prot": aa_muts,
                                "nt_start": hsp.sbjct_start,
                                "nt_end": hsp.sbjct_end,
                                "strand": "+" if hsp.frame[1] > 0 else "-",
                            }
                    else:
                        continue
            if best_hit is not None:
                all_hits[query_id] = best_hit
            else:
                all_hits[query_id] = {}

    # make sure we processed all queries and get the nucleotide hit
    extraneous_hits = set(all_hits)
    for query in Bio.SeqIO.parse(query_prots, "fasta"):
        if query.id not in all_hits:
            raise ValueError(f"No results for {query.id=}")
        extraneous_hits.remove(query.id)

        hit_d = all_hits[query.id]

        hit_d["prot_sequence"] = str(query.seq)

        if "accession" in all_hits[query.id]:
            cmds = [
                "blastdbcmd",
                "-db", blastdb,
                "-dbtype", "nucl",
                "-entry", hit_d["accession"],
            ]
            result = subprocess.run(cmds, capture_output=True, text=True, check=True)
            seq_record = Bio.SeqIO.read(io.StringIO(result.stdout), "fasta")
            assert seq_record.id.split(".")[0] == hit_d["accession"]

            nt_start = hit_d["nt_start"]
            nt_end = hit_d["nt_end"]
            strand = hit_d["strand"]
            if (nt_start > nt_end) or (nt_end - nt_start + 1 != 3 * len(query.seq)):
                raise ValueError(f"{nt_start=}, {nt_end=}, {strand=}, {len(query.seq)=}")
            cdna = seq_record[nt_start - 1: nt_end]
            if strand == "-":
                cdna = cdna.reverse_complement()
            elif strand != "+":
                raise ValueError(f"Invalid {strand=}")
            assert len(hit_d["prot_sequence"]) * 3 == len(cdna)

            # add mutations
            cdna = list(str(cdna.seq))
            for mut in hit_d["aa_muts_accession_to_prot"]:
                acc_aa, r, query_aa = mut[0], int(mut[1: -1]), mut[-1]
                assert query_aa == query.seq[r - 1], hit_d
                codon = "".join(cdna[3 * r - 3: 3 * r])
                assert codon in CODON_TO_AA, f"{codon=}, {hit_d=}"
                assert acc_aa == CODON_TO_AA[codon], hit_d
                cdna[3 * r - 3: 3 * r] = AA_TO_BEST_HUMAN_CODON[query_aa]
            cdna = "".join(cdna)

            if query.seq != str(Bio.Seq.Seq(cdna).translate()):
                raise ValueError(f"For {query.id=}, {query.seq=} != {cdna=}")

            hit_d["nt_sequence"] = cdna

    if extraneous_hits:
        raise ValueError(f"{extraneous_hits=}")

    cols_to_keep = [
        "accession_for_nt_sequence",
        "accession_w_aa_muts_added",
        "aa_muts_added_to_accession",
        "prot_sequence",
        "nt_sequence",
    ]
    if query_id_name in cols_to_keep:
        raise ValueError(f"invalid {query_id_name=}, is in {cols_to_keep=}")
    df = (
        pd.DataFrame.from_dict(all_hits, orient="index")
        .assign(
            aa_muts_added_to_accession=lambda x: (
                x["aa_muts_accession_to_prot"].fillna("").map(
                    lambda m: "_".join(m) if m else ""
                )
            ),
            accession_w_aa_muts_added=lambda x: x.apply(
                lambda r: (
                    ""
                    if pd.isnull(r["accession"])
                    else (
                        r["accession"] + "_" + r["aa_muts_added_to_accession"]
                        if r["aa_muts_added_to_accession"]
                        else r["accession"]
                    )
                ),
                axis=1,
            ),
        )
        .rename(columns={"accession": "accession_for_nt_sequence"})
        [cols_to_keep]
        .reset_index(names=query_id_name)
    )

    return df


def main():
    """Main body of script that runs command-line program."""

    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
        pass

    parser = argparse.ArgumentParser(
        description=textwrap.dedent(
            """\
            Match a set of protein sequences against Genbank nucleotide sequences.

            See https://github.com/jbloom/match_prot_to_genbank_nt

            Performs the following steps:
              1. Downloads as search target all NCBI sequences of a given taxon.
              2. Builds a BLAST database.
              3. Queries the proteins against the database with `tblastn`.
              4. Find highest identity match that covers full protein with no indels.
              5. Recover nucleotide sequence of match and add any needed mutations to  make
                 encoded protein fully match query protein (using most common human codon).
              6. Report match results in CSV; any protein missing a match that covers full
                 protein with no gaps will have some empty columns. The CSV is in
                 '<outdir>/match_prot_to_genbank_nt.csv'.

            Example usage:
                
                python match_prot_to_genbank_nt.py \\
                    --query-prots example_prots.fasta \\
                    --taxon "Influenza A virus" \\
                    --outdir example_outdir
            """
        ),
        formatter_class=CustomFormatter,
    )
    parser.add_argument(
        "--query-prots",
        required=True,
        help="FASTA file with query proteins.",
    )
    parser.add_argument(
        "--taxon",
        required=True,
        help='NCBI taxon name (eg, "Influenza A virus") for targets to search'
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="output directory",
    )
    parser.add_argument(
        "--existing-outdir",
        choices={"error", "delete", "keep"},
        default="error",
        help="What if some expected output already exists in outdir?",
    )
    parser.add_argument("--evalue", default=1e-10, help="BLAST E-value")
    parser.add_argument("--max_target_seqs", default=10, help="max BLAST target seqs")
    parser.add_argument("--threads", default=4, help="threads when running BLAST")
    parser.add_argument(
        "--query-id-name", default="strain", help="Name for query id in output CSV"
    )

    args = parser.parse_args()

    outdir = args.outdir
    if os.path.isdir(args.outdir):
        if args.existing_outdir == "error":
            raise IOError(f"{outdir=} already exists")
        elif args.existing_outdir == "delete":
            print(f"Deleting existing {outdir=}")
            shutil.rmtree(outdir)
            os.mkdir(outdir)
        elif args.existing_outdir == "keep":
            print(f"Using existing output in {outdir=} where available")
        else:
            raise ValueError(f"invalid {args.existing_outdir=}")
    else:
        print(f"Creating {outdir=}")
        os.mkdir(outdir)

    targets_zip = os.path.join(outdir, "targets.zip")
    targets_fasta = os.path.join(outdir, "targets.fasta")
    if not os.path.isfile(targets_fasta):
        print(f"Getting {targets_fasta=}")
        get_target_fasta(args.taxon, targets_zip, targets_fasta)
    else:
        print(f"Using existing {targets_fasta=}")

    blastdb = os.path.join(outdir, "blastdb/targets")
    if not os.path.isdir(os.path.dirname(blastdb)):
        print(f"Creating BLAST database {blastdb=} from {targets_fasta=}")
        make_blastdb(targets_fasta, blastdb)
    else:
        print(f"Using existing BLAST database {blastdb=}")

    blast_results = os.path.join(outdir, "blast_results.xml")
    query_prots = args.query_prots
    if not os.path.isfile(blast_results):
        print(f"Querying {query_prots=} against {blastdb=} to create {blast_results=}")
        run_tblastn(
            query_prots,
            blastdb,
            blast_results,
            args.evalue,
            args.max_target_seqs,
            args.threads,
        )
    else:
        print(f"Using existing BLAST results {blast_results=}")

    output_csv = os.path.join(outdir, "match_prot_to_genbank_nt.csv")
    if not os.path.isfile(output_csv):
        print(f"Parsing blast output to create {output_csv=}")
        df = parse_blast_xml(query_prots, blast_results, blastdb, args.query_id_name)
        df.to_csv(output_csv, index=False)
    else:
        print(f"Using existing {output_csv=}")

    print(f"Program complete.")


if __name__ == "__main__":
    main()
