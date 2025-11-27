import pandas as pd

library_strains = pd.read_csv("results/aggregated_library_strains/library_strains.tsv", sep="\t")

existing_library = set(pd.read_csv("../initial_design/results/aggregated_library_strains/library_strains.tsv", sep="\t")["prot_sequence"])

haplotypes = (
    pd.read_csv("data/H3N2_haplotypes.tsv", sep="\t")
    .rename(columns={"representative_strain_ha_sequence": "prot_sequence"})
    [["derived_haplotype", "count", "prot_sequence"]]
    .groupby("prot_sequence", as_index=False)
    .first()
)

library_strains = (
    library_strains
    .query("prot_sequence not in @existing_library")
    .query("subtype == 'H3N2'")
    .merge(haplotypes, how="left", on="prot_sequence", validate="one_to_one")
)

library_strains = library_strains[
    ["derived_haplotype", "count"] + [c for c in library_strains.columns if c not in ["derived_haplotype", "count"]]
]

print(library_strains)

library_strains.to_csv("strains_to_add.tsv", sep="\t", index=False)
