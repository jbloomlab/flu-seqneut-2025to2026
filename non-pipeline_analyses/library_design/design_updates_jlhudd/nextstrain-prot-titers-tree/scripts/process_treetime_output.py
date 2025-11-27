"""Process results of ``treetime`` output to formats used by ``augur``."""

import copy
import json
import re
import sys

import Bio.Phylo
import Bio.SeqIO

import pandas as pd


sys.stderr = sys.stdout = open(snakemake.log[0], "w")

print("Reading the trees and removing the outgroup.")
trees = {}
root_name = None
for treetype in ["timetree", "divtree"]:
    intree = getattr(snakemake.input, f"{treetype}_w_outgroup")
    outtree = getattr(snakemake.output, treetype)

    print(f"Reading tree from {intree}, pruning outgroup, writing to {outtree}")
    tree_w_outgroup = Bio.Phylo.read(intree, "nexus")

    # Bio.Phylo can incorrectly parse names as confidence; fix that.
    for clade in tree_w_outgroup.find_clades(order="preorder"):
        if (clade.name is None) and (clade.confidence is not None):
            clade.name = clade.confidence
            clade.confidence = None

    # Remove outgroup
    outgroup = [c for c in tree_w_outgroup.find_clades() if c.name == "outgroup"]
    assert len(outgroup) == 1, "not one sequence called outgroup"
    # get common ancestor of all non-outgroup tips
    other_tips = [t for t in tree_w_outgroup.get_terminals() if t.name != "outgroup"]
    assert len(other_tips) + 1 == len(tree_w_outgroup.get_terminals())
    root = tree_w_outgroup.common_ancestor(other_tips)
    root.branch_length = 0
    root.comment = None
    tree = Bio.Phylo.BaseTree.Tree(root=root, rooted=True)
    if any(t.name == "outgroup" for t in tree.get_terminals()):
        raise ValueError("the outgroup is not actually outgroup to all other tips")
    assert root.name, root
    if root_name:
        if root.name != root_name:
            raise ValueError("inconsistent root names between timetree and divtree")
    else:
        root_name = root.name

    # Remove comments from a version of tree to write.
    tree_nocomments = copy.deepcopy(tree)
    for clade in tree_nocomments.find_clades(order="preorder"):
        clade.comment = None

    # write tree without comments to Newick
    Bio.Phylo.write(tree_nocomments, outtree, format="newick")

    # save tree w comments for later use
    trees[treetype] = tree

# Get the root (reference) sequence
root_seqs = [
    s
    for s in Bio.SeqIO.parse(snakemake.input.ancestral_sequences, "fasta")
    if s.id == root_name
]
assert len(root_seqs) == 1, f"{len(root_seqs)=} for {root.name=}"
refseq = str(root_seqs[0].seq)
print(f"Using as the reference sequence {root_name}; {len(refseq)=}")

# Get the site numbering map and do some error checking on it
print(f"Reading site numbering map from {snakemake.input.site_numbering_map}")
site_numbering_map = pd.read_csv(snakemake.input.site_numbering_map, sep="\t")
if len(refseq) != len(site_numbering_map):
    raise ValueError(f"{len(refseq)=}, {len(site_numbering_map)=}")
req_cols = {"sequential_site", "protein", "protein_site"}
if "nuc" in set(site_numbering_map["protein"]):
    raise ValueError("Cannot have a protein called 'nuc' in site_numbering_map")
if not req_cols.issubset(site_numbering_map.columns):
    raise ValueError(f"{site_numbering_map.columns=} lacks {req_cols=}")
for col in ["sequential_site", "protein_site"]:
    if not pd.api.types.is_integer_dtype(site_numbering_map[col]):
        raise ValueError(f"{col=} in site_numbering_map not integer")
if not (
    len(site_numbering_map)
    == site_numbering_map["sequential_site"].nunique()
    == site_numbering_map["sequential_site"].max()
    == len(
        range(
            site_numbering_map["sequential_site"].min(),
            site_numbering_map["sequential_site"].max() + 1,
        )
    )
):
    raise ValueError("sequential_site in site_numbering map not 1, 2, ... integers")
if not (
    len(site_numbering_map[["protein", "protein_site"]].drop_duplicates())
    == len(site_numbering_map)
):
    raise ValueError("each sequential_site not mapped to unique protein, protein_site")
refseq_proteins = {}
protein_coords = {}
for protein, protein_df in site_numbering_map.groupby("protein"):
    protein_start = int(protein_df["sequential_site"].min())
    protein_end = int(protein_df["sequential_site"].max())
    protein_coords[protein] = {"start": 3 * protein_start - 2, "end": protein_end * 3}
    if protein_end - protein_start + 1 != len(protein_df):
        raise ValueError(f"{protein=} in site_numbering_map not sequential consecutive")
    assert protein_start >= 1 and protein_end <= len(refseq)
    refseq_proteins[protein] = refseq[protein_start - 1 : protein_end]
    if not (
        (protein_df["protein_site"].min() == 1)
        and (protein_df["protein_site"].max() == len(protein_df))
        and (protein_df["protein_site"].nunique() == len(protein_df))
    ):
        raise ValueError(
            f"{protein=} does not have sequential consecutive protein_site"
        )
site_numbering_map = (
    site_numbering_map.assign(protein_site=lambda x: x["protein_site"].astype(int))
    .set_index("sequential_site")
    .to_dict(orient="index")
)

# Get amino-acid mutations in JSON node format
aa_muts_nodes = {}

mut_pat = re.compile(r'mutations="(?P<mutations>[^"]*)"')

for clade in trees["divtree"].find_clades(order="preorder"):
    assert clade.name, f"{clade=} has no name"
    # Biopython puts comment text on .comment for Nexus
    comment = getattr(clade, "comment", None)
    if comment:
        m = mut_pat.search(comment)
        if m:
            assert len(mut_pat.findall(comment)) == 1, comment
            muts = {}
            for mut_str in m.group("mutations").split(","):
                m_match = re.fullmatch(
                    r"(?P<parent>[A-Z\-])(?P<site>\d+)(?P<mut>[A-Z\-])", mut_str
                )
                assert m_match, f"{mut_str=}\n{comment=}"
                site = int(m_match.group("site"))
                if site not in site_numbering_map:
                    raise ValueError(
                        f"invalid {site=} in {mut_str=} not in {site_numbering_map=}"
                    )
                prot = site_numbering_map[site]["protein"]
                prot_site = site_numbering_map[site]["protein_site"]
                if prot not in muts:
                    muts[prot] = []
                muts[prot].append(
                    f"{m_match.group('parent')}{prot_site}{m_match.group('mut')}"
                )
            if muts:
                aa_muts_nodes.setdefault(clade.name, {})["aa_muts"] = muts


# Write the amino-acid mutations and reference sequence
print(f"Writing aa mutations to {snakemake.output.aa_muts}")
aa_muts_d = {
    "annotations": {
        "nuc": {"start": 1, "end": len(refseq) * 3, "strand": "+", "type": "source"},
    }
    | {
        prot: coords | {"strand": "+", "type": "gene"}
        for prot, coords in protein_coords.items()
    },
    "nodes": aa_muts_nodes,
    "reference": {"nuc": "N" * len(refseq) * 3} | refseq_proteins,
}
with open(snakemake.output.aa_muts, "w") as f:
    f.write(json.dumps(aa_muts_d, indent=2))

# Get the branch lengths in JSON node format, in units of mutations
brlens_nodes = {}

dates = (
    pd.read_csv(snakemake.input.dates, sep="\t")
    .rename(columns={"numeric date": "num_date"})
    .set_index("#node")
    .to_dict(orient="index")
)
for clade in trees["divtree"].find_clades(order="preorder"):
    assert clade.branch_length is not None, clade
    # Handle '--' values from treetime (dates that couldn't be inferred)
    num_date = dates[clade.name]["num_date"]
    if num_date == "--":
        num_date = None
    else:
        num_date = float(num_date)

    brlens_nodes[clade.name] = {
        "mutation_length": clade.branch_length * len(refseq),
        "date": dates[clade.name]["date"],
        "num_date": num_date,
    }

print(f"Writing branch lengths to {snakemake.output.brlens}")
with open(snakemake.output.brlens, "w") as f:
    f.write(json.dumps({"nodes": brlens_nodes}, indent=2))
