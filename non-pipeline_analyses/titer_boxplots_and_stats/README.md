# titer_boxplots_and_stats

Manual analysis of H3N2 neutralization titers from the 2025–2026 seasonal influenza serology study. Produces boxplots and Mann-Whitney statistics for the 302 'All'-group human sera, used for paper figures.

## Input data

Read from `results/final_titer_data/` in the repo root (relative paths, no download needed):

| File | Contents |
|---|---|
| `human_titers.csv` | Individual titers: one row per serum × virus measurement |
| `human_viruses.csv` | Virus metadata: subclade, derived haplotype, collection date, sequences |

## Notebooks

### `subclade_titers.py`

Boxplots of all titer measurements (serum × virus) for viruses in subclades **K, J.2.4, J.2.3, J.2.2, and J.2**. Tests the hypothesis that K-subclade titers are lower than other subclades.

**Statistics:**
- K vs all others combined: Mann-Whitney U (one-sided, uncorrected — pre-specified directional hypothesis)
- All pairwise subclade comparisons: Mann-Whitney U with Benjamini-Hochberg correction

```bash
marimo run subclade_titers.py
```

---

### `mutation_titers.py`

For K-subclade viruses, compares titers between viruses with and without specific HA mutations. Each panel is one pre-specified comparison (Mann-Whitney U, two-sided, uncorrected).

**Current comparisons:**
- K135E present vs absent
- R189K present vs absent
- Any mutation at sites 96, 207, or 223 vs none

To add a new comparison, extend the `COMPARISONS` list near the top of the notebook with a new dict containing a `filter_fn` lambda — a new panel will be generated automatically.

> **Note:** K135E and R189K each appear in only 1 of 20 K-subclade viruses. Interpret those comparisons cautiously.

```bash
marimo run mutation_titers.py
```

---

To open either notebook in interactive editing mode, use `marimo edit` instead of `marimo run`.