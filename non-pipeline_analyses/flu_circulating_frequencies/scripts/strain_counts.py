"""Get counts of sequences for each strain."""


import collections
import regex
import sys

import Bio.SeqIO

import pandas as pd


sys.stdout = sys.stderr = log = open(snakemake.log[0], "w")

aas = "ACDEFGHIKLMNPQRSTVWY"

trim_start, trim_end = snakemake.params.trim
assert 1 <= trim_start < trim_end

print(f"Analysis for {snakemake.wildcards.protset=}")

print(f"Reading strain_prots from {snakemake.input.strain_prots}")
print(f"Trimming to {trim_start} to {trim_end}")
strain_prots = {}
prot_names = collections.defaultdict(list)
for seq in Bio.SeqIO.parse(snakemake.input.strain_prots, "fasta"):
    name = seq.id
    if name in strain_prots:
        raise ValueError(f"Duplicate {name=} in {snakemake.input.strain_prots}")
    s = str(seq.seq).upper()
    while s.endswith("*"):  # remove trailing stop codons
        s = s[: -1]
    assert regex.fullmatch(f"[{aas}]+", s), f"{name=} has invalid residues:\n{s}\n{seq}"
    assert len(s) >= trim_end, f"{trim_end=}, {len(s)=}"
    trimmed_s = s[trim_start - 1: trim_end]
    if trimmed_s in prot_names:
        print(f"{name} has same trimmed sequence as {prot_names[trimmed_s]}")
    prot_names[trimmed_s].append(name)
    strain_prots[name] = trimmed_s
print(f"Read {len(strain_prots)=} strain proteins\n")
assert len(strain_prots) >= len(set(strain_prots.values()))

# define regexes allowing fuzzy matching but getting best match
maxdiff = int(snakemake.params.maxdiff)
strain_regexes = {
    strain_name: regex.compile(f"(?b)(?:{strain_seq}){{e<={maxdiff}}}")
    for (strain_name, strain_seq) in strain_prots.items()
}

print(f"Reading proteins from {snakemake.input.protset}, matching with {maxdiff=}")
strain_match_records = []
unmatched_seq_counts = collections.defaultdict(int)
unmatched_recent_seq_counts = collections.defaultdict(int)

recent_date = snakemake.params.recent_date
print(f"Defining a recent sequence as one after {recent_date=}")

assert "other" not in strain_prots
accessions = set()
n_incomplete_date = 0
for prot in Bio.SeqIO.parse(snakemake.input.protset, "fasta"):

    p = str(prot.seq).upper()

    try:
        accession, _, _, date = (t.strip() for t in prot.description.split("|"))
    except ValueError:
        raise ValueError(f"Problem parsing header {prot.description=}")

    if not regex.fullmatch(r"\d{4}\-\d{2}\-\d{2}", date):
        if regex.fullmatch(r"\d{4}(?:\-\d{2})?", date):
            n_incomplete_date += 1
            continue
        else:
            raise ValueError(f"invalid {date=} for {prot.description=}")

    assert pd.notnull(accession), prot.description
    if accession in accessions:
        continue  # duplicate accession
    accessions.add(accession)

    matches = collections.defaultdict(list)
    
    for strain_name, strain_regex in strain_regexes.items():
        if m := strain_regex.search(p):
            ndiff = sum(m.fuzzy_counts)
            assert ndiff <= maxdiff
            matches[ndiff].append((strain_name, accession, date))
    if matches:
        # get best match
        matchlist = sorted(matches.items())[0][1]
        weight = 1 / len(matchlist)
        for tup in matchlist:
            strain_match_records.append((*tup, weight))
    else:
        strain_match_records.append(("other", accession, date, 1))
        unmatched_seq_counts[p] += 1 # Get counts of each unique sequence that doesn't appear in library

        # Get counts of each unique sequence that doesn't appear in library AND that is recent
        if date >= recent_date:  # string comparison works for YYYY-MM-DD format?
            unmatched_recent_seq_counts[p] += 1
        

print(f"Dropped {n_incomplete_date} for having an incomplete date")

strain_matches = pd.DataFrame(
    strain_match_records, columns=["variant", "accession", "date", "weight"]
)

print(f"Read {len(strain_matches)} sequences.")
assert strain_matches["accession"].nunique() == strain_matches["weight"].sum()

overall_counts = (
    strain_matches
    .groupby("variant", as_index=False)
    .aggregate(n_sequences=pd.NamedAgg("weight", "sum"))
    .merge(
        pd.Series(strain_prots).rename_axis("variant").rename("sequence").reset_index(),
        on="variant",
        how="outer",
        validate="one_to_one",
    )
    .assign(n_sequences=lambda x: x["n_sequences"].fillna(0).astype(int))
    .sort_values("n_sequences", ascending=False)
)

print("Here are overall number matching each strain:")
print(overall_counts)
print(f"Writing to {snakemake.output.counts_overall}")
overall_counts.to_csv(snakemake.output.counts_overall, index=False)

print(f"\nWriting all strain matches to {snakemake.output.strain_matches}")
strain_matches.to_csv(snakemake.output.strain_matches, index=False)

print(f"\nWriting strain counts by date to {snakemake.output.counts_by_date}")
(
    strain_matches
    .groupby(["variant", "date"], as_index=False)
    .aggregate(sequences=pd.NamedAgg("weight", "sum"))
    .to_csv(snakemake.output.counts_by_date, index=False)
)

print(f"\nWriting counts of sequences with no library matches to {snakemake.output.counts_unmatched_sequences}")
df = pd.DataFrame.from_dict(unmatched_seq_counts, orient="index", columns=["count"])
df.index.name = "sequence"
df.reset_index(inplace=True)
df = df.sort_values("count", ascending=False)
df.to_csv(snakemake.output.counts_unmatched_sequences, index=False)

print(f"\nWriting counts of recent sequences with no library matches to {snakemake.output.counts_recent_unmatched_sequences}")
df = pd.DataFrame.from_dict(unmatched_recent_seq_counts, orient="index", columns=["count"])
df.index.name = "sequence"
df.reset_index(inplace=True)
df = df.sort_values("count", ascending=False)
df.to_csv(snakemake.output.counts_recent_unmatched_sequences, index=False)
