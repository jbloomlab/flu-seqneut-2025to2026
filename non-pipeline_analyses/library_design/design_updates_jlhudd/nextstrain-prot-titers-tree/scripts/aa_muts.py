"""Process results of ``treetime ancestral`` to JSON of amino-acid mutations."""

import json
import re
import sys

import Bio.Phylo
import Bio.SeqIO


sys.stderr = sys.stdout = open(snakemake.log[0], "w")

tree = Bio.Phylo.read(snakemake.input.annotated_tree, "nexus")
ancestral_sequences = Bio.SeqIO.parse(snakemake.input.ancestral_sequences, "fasta")
gene = snakemake.params.gene
aa_muts = snakemake.output.aa_muts

mut_pat = re.compile(r'mutations="(?P<mutations>[^"]*)"')

nodes = {}  # JSON nodes for mutations

for clade in tree.find_clades(order="preorder"):
    if (clade.name is None) and (clade.confidence is not None):
        name = str(clade.confidence)  # Nexus wrapper can parse names as confidence
    else:
        name = clade.name
    assert name, f"{clade=} has no name"
    muts = []
    # Biopython puts comment text on .comment for Nexus
    comment = getattr(clade, "comment", None)
    if comment:
        m = mut_pat.search(comment)
        if m:
            assert len(mut_pat.findall(comment)) == 1
            muts = []
            for mut_str in m.group("mutations").split(","):
                assert re.fullmatch(
                    r"[A-Z]\d+[A-Z]", mut_str
                ), f"{mut_str=}\n{comment=}"
                muts.append(mut_str)
            if muts:
                nodes.setdefault(name, {})["aa_muts"] = {gene: muts}

# Get the root (reference) sequence
root = tree.root
assert tree.root == tree.common_ancestor(tree.get_terminals())
root_seqs = [s for s in ancestral_sequences if s.id == root.name]
assert len(root_seqs) == 1, f"{len(root_seqs)=} for {root.name=}"
out = {"nodes": nodes, "reference": {gene: str(root_seqs[0].seq)}}

with open(aa_muts, "w") as f:
    f.write(json.dumps(out, indent=2))
