import re

import pandas as pd

import Bio.SeqIO

for subtype, length in [("H3N2", 566), ("H1N1", 566)]:
    haplotypes = {}
    nseq = 0
    for seq in Bio.SeqIO.parse(f"gisaid_seqs/{subtype}.fasta", "fasta"):
        nseq += 1
        prot = str(seq.seq)
        strain = seq.id.split("|")[0]
        date = seq.id.split("|")[-1]
        assert re.fullmatch(r"\d{4}\-\d{2}\-\d{2}", date), f"{seq.id=}, {date=}"
        if prot in haplotypes:
            haplotypes[prot].append((date, strain))
        else:
            haplotypes[prot] = [(date, strain)]
    print(f"For {subtype=}, parsed {nseq=} for {len(haplotypes)=}")
    records = []
    for haplotype, date_strains in haplotypes.items():
        date_strains = sorted(date_strains)
        rep_strain = date_strains[-1][1]
        last_date = date_strains[-1][0]
        records.append((rep_strain, rep_strain, last_date, len(date_strains), haplotype))
    df = (
        pd.DataFrame(
            records,
            columns=["representative_strain", "derived_haplotype", "latest_sequence", "count", "representative_strain_ha_sequence"],
        )
        .sort_values("count", ascending=False)
        .assign(length=lambda x: x["representative_strain_ha_sequence"].map(len))
        .groupby("representative_strain", as_index=False)
        .first()
        .assign(length=lambda x: x["representative_strain_ha_sequence"].map(len))
    )
    print(
        df
        .groupby("length")
        .aggregate(
            n_haplotypes=pd.NamedAgg("count", "count"),
            n_sequences=pd.NamedAgg("count", "sum"),
        )
    )
    print(f"Only keeping {subtype=} of {length=}")
    df = df[df["length"] == length]
    df.to_csv(f"recent_{subtype}_haplotypes.tsv", index=False, sep="\t")

    strain_records = [(haplotype, ", ".join([tup[1] for tup in d])) for haplotype, d in haplotypes.items()]
    (
        pd.DataFrame(strain_records, columns=["haplotype", "strains"])
        .to_csv(f"haplotype_strains_{subtype}.tsv", sep="\t", index=False)
    )
