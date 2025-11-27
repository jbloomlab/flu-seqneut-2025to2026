# Match protein sequences to the closest Genbank nucleotide sequences encoding the protein

This repository provides a Python script ([match_prot_to_genbank_nt.py](match_prot_to_genbank_nt.py)) that is designed for the use case when you have a protein sequence and want to find a nucleotide sequence in Genbank encoding that protein.
If there is a Genbank sequence that aligns with high identity to the full-length protein with no indels, it will be found and mutations added if needed to make it match the protein.

This script was written by [Jesse Bloom](https://jbloomlab.org/) for the use case of finding Genbank influenza nucleotide sequences that match influenza protein sequences of interest, but it should work more generally.

To run the script, first build the `conda` environment in [environment.yaml](environment.yaml), which contains the required programs (most prominently [NCBI datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/) and `blast`).
Then activate that conda environment with:

    conda activate match_prot_to_genbank_nt

and run the Python script [match_prot_to_genbank_nt.py](match_prot_to_genbank_nt.py).

To run on the example protein input in [example_prots.fasta](example_prots.fasta), do:

    python match_prot_to_genbank_nt.py --query-prots example_prots.fasta --taxon "Influenza A virus" --outdir example_outdir

After running this, the output directory will contain a variety of files associated with downloaded the NCBI dataset and running BLAST.
The key output will have name `<outdir>/match_prot_to_genbank_nt.csv`; for the example that file is [example_outdir/match_prot_to_genbank_nt.csv](example_outdir/match_prot_to_genbank_nt.csv).

For a detailed usage message with other options, run `python match_prot_to_genbank_nt.py --help`.
That will give the following usage message:

    usage: match_prot_to_genbank_nt.py [-h] --query-prots QUERY_PROTS --taxon TAXON --outdir OUTDIR
                                       [--existing-outdir {error,keep,delete}] [--evalue EVALUE]
                                       [--max_target_seqs MAX_TARGET_SEQS] [--threads THREADS]
                                       [--query-id-name QUERY_ID_NAME]
    
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
    
        python match_prot_to_genbank_nt.py \
            --query-prots example_prots.fasta \
            --taxon "Influenza A virus" \
            --outdir example_outdir
    
    options:
      -h, --help            show this help message and exit
      --query-prots QUERY_PROTS
                            FASTA file with query proteins. (default: None)
      --taxon TAXON         NCBI taxon name (eg, "Influenza A virus") for targets to search (default: None)
      --outdir OUTDIR       output directory (default: None)
      --existing-outdir {error,keep,delete}
                            What if some expected output already exists in outdir? (default: error)
      --evalue EVALUE       BLAST E-value (default: 1e-10)
      --max_target_seqs MAX_TARGET_SEQS
                            max BLAST target seqs (default: 10)
      --threads THREADS     threads when running BLAST (default: 4)
      --query-id-name QUERY_ID_NAME
                            Name for query id in output CSV (default: strain)
