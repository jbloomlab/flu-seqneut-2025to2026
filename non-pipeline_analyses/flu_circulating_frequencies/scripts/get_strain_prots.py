import Bio.SeqIO


subtype_prots = {subtype: [] for subtype in snakemake.params.subtypes}

for s in Bio.SeqIO.parse(snakemake.input.fasta, "fasta"):
    for subtype in subtype_prots:
        if s.id.endswith(subtype):
            subtype_prots[subtype].append(s)
            break
    else:
        raise ValueError(f"Cannot match subtype of {s.id}")

for subtype, prots in subtype_prots.items():
    outf = snakemake.output[subtype]
    with open(outf, "w") as f:
        Bio.SeqIO.write(prots, f, "fasta")
