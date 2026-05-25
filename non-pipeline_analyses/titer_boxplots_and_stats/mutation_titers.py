import marimo

__generated_with = "0.17.2"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from scipy import stats
    import re
    from pathlib import Path
    return Path, mo, np, pd, plt, re, stats


@app.cell
def _(mo):
    mo.md(r"""
    # Mutation-stratified titer boxplots

    Compare neutralization titers between viruses **with** and **without**
    specific mutations (or mutations at specific sites) in the HA protein.
    Comparisons can span any combination of subclades and both subtypes.

    Each data point is one (serum × virus) measurement, log₂-transformed.
    Statistics: Mann-Whitney U (two-sided, uncorrected — one pre-specified
    comparison per panel). Fold change = ratio of geometric means.
    """)
    return


@app.cell
def _(Path, pd):
    # ---------------------------------------------------------------------------
    # Paths — relative to this notebook's location
    # ---------------------------------------------------------------------------
    HERE = Path(__file__).parent
    REPO_ROOT = HERE.parent.parent
    DATA_DIR = REPO_ROOT / "results" / "final_titer_data"

    titers_raw = pd.read_csv(DATA_DIR / "human_titers.csv")
    viruses = pd.read_csv(DATA_DIR / "human_viruses.csv")

    results_dir = HERE / "results"
    results_dir.mkdir(exist_ok=True)
    return titers_raw, viruses, results_dir


@app.cell
def _(np, pd, titers_raw, viruses):
    # ---------------------------------------------------------------------------
    # Pre-merge all titer data with virus metadata
    # Covers all subclades/subtypes — comparisons filter by subclade via COMPARISONS
    # ---------------------------------------------------------------------------
    df_all = titers_raw.merge(
        viruses[["virus", "subclade", "derived_haplotype"]],
        on="virus", how="inner",
    )
    df_all["log2_titer"] = np.log2(df_all["titer"])

    haplotype_table = (
        viruses[["virus", "subclade", "derived_haplotype"]]
        .sort_values(["subclade", "derived_haplotype"])
        .reset_index(drop=True)
    )
    return df_all, haplotype_table


@app.cell
def _(mo):
    mo.md("### Virus haplotypes (all subclades)")
    return


@app.cell
def _(haplotype_table):
    haplotype_table.sort_values("derived_haplotype").reset_index(drop=True)
    return


@app.cell
def _(re):
    # ---------------------------------------------------------------------------
    # Helper functions for mutation parsing
    # ---------------------------------------------------------------------------

    def parse_mutations(haplotype: str) -> list[str]:
        """Return list of extra mutation strings, e.g. ['A272T', 'S279T']."""
        h = str(haplotype)
        if ":" not in h:
            return []
        return h.split(":")[1].split(",")

    def has_mutation(haplotype: str, mutation: str) -> bool:
        """Check if haplotype contains a specific mutation (exact match)."""
        return mutation in parse_mutations(haplotype)

    def has_mutation_at_sites(haplotype: str, sites: list[int]) -> bool:
        """Check if haplotype has any mutation at any of the given numeric sites."""
        muts = parse_mutations(haplotype)
        for m in muts:
            match = re.search(r"\d+", m)
            if match and int(match.group()) in sites:
                return True
        return False
    return has_mutation, has_mutation_at_sites


@app.cell
def _(mo):
    mo.md("""
    ---
    ## Comparison definitions

    Edit the list below to add or remove comparisons.
    Each entry is a dict with:
    - `title`: panel title
    - `filter_fn`: function(haplotype) → bool, True = "with mutation" group
    - `subclades`: **list** of subclades to draw viruses from (can mix subtypes)
    - `with_label` / `without_label`: group labels for the two boxes
    """)
    return


@app.cell
def _(has_mutation, has_mutation_at_sites):
    # ---------------------------------------------------------------------------
    # Comparison definitions — add more here as needed
    # ---------------------------------------------------------------------------

    # H3N2 — K subclade
    H3N2_K_COMPARISONS = [
        {
            "title": "K135E",
            "subclades": ["K"],
            "with_label": "K:K135E",
            "without_label": "K: no K135E",
            "filter_fn": lambda h: has_mutation(h, "K135E"),
        },
        {
            "title": "R189K",
            "subclades": ["K"],
            "with_label": "K:R189K",
            "without_label": "K: no R189K",
            "filter_fn": lambda h: has_mutation(h, "R189K"),
        },
        {
            "title": "Any mutation at sites 96, 207, 223, or 261",
            "subclades": ["K"],
            "with_label": "K: mut@96/207/223/261",
            "without_label": "K: no mut@96/207/223/261",
            "filter_fn": lambda h: has_mutation_at_sites(h, [96, 207, 223, 261]),
        },
    ]

    # H3N2 — J.2.4 subclade
    H3N2_J24_COMPARISONS = [
        {
            "title": "K135N",
            "subclades": ["J.2.4"],
            "with_label": "J.2.4:K135N",
            "without_label": "J.2.4: no K135N",
            "filter_fn": lambda h: has_mutation(h, "K135N"),
        },
    ]

    # H1N1 — D.3.1 + D.3.1.1 (pooled so mutations can span both subclades)
    H1N1_COMPARISONS = [
        {
            "title": "Mut at sites 155 or 157",
            "subclades": ["D.3.1", "D.3.1.1"],
            "with_label": "mut@155/157",
            "without_label": "no mut@155/157",
            "filter_fn": lambda h: has_mutation_at_sites(h, [155, 157]),
        },
    ]

    return H1N1_COMPARISONS, H3N2_J24_COMPARISONS, H3N2_K_COMPARISONS


@app.cell
def _(H1N1_COMPARISONS, H3N2_J24_COMPARISONS, H3N2_K_COMPARISONS, df_all, mo, np, pd, plt, results_dir, stats, viruses):
    # ---------------------------------------------------------------------------
    # Build figures — one per comparison group, original layout preserved
    # ---------------------------------------------------------------------------
    COLOR_WITH    = "#0072B2"
    COLOR_WITHOUT = "#E69F00"

    def _sig_label(p):
        if p < 0.001: return "***"
        elif p < 0.01: return "**"
        elif p < 0.05: return "*"
        else: return "ns"

    def _make_figure(comparisons, suptitle=None):
        n_panels = len(comparisons)
        fig, axes = plt.subplots(1, n_panels, figsize=(4 * n_panels, 5.5), sharey=True)
        if n_panels == 1:
            axes = [axes]

        stat_rows = []

        tick_values = [40, 80, 160, 320, 640, 1280, 2560, 5120, 10240, 16384]
        log2_ticks  = np.log2(tick_values)

        for ax, comp in zip(axes, comparisons):
            sc_viruses = viruses[viruses["subclade"].isin(comp["subclades"])][
                ["virus", "derived_haplotype"]
            ]

            with_viruses = sc_viruses[
                sc_viruses["derived_haplotype"].apply(comp["filter_fn"])
            ]["virus"].tolist()
            without_viruses = sc_viruses[
                ~sc_viruses["derived_haplotype"].apply(comp["filter_fn"])
            ]["virus"].tolist()

            df_with    = df_all[df_all["virus"].isin(with_viruses)]["log2_titer"].values
            df_without = df_all[df_all["virus"].isin(without_viruses)]["log2_titer"].values

            n_with_v    = len(with_viruses)
            n_without_v = len(without_viruses)

            if len(df_with) > 0 and len(df_without) > 0:
                _, p_val    = stats.mannwhitneyu(df_with, df_without, alternative="two-sided")
                fold_change = 2 ** (np.mean(df_with) - np.mean(df_without))
            else:
                p_val       = float("nan")
                fold_change = float("nan")

            stat_rows.append({
                "comparison": comp["title"],
                "subclades": ", ".join(comp["subclades"]),
                "n_viruses_with": n_with_v,
                "n_viruses_without": n_without_v,
                "n_obs_with": len(df_with),
                "n_obs_without": len(df_without),
                "fold_change_with_vs_without": fold_change,
                "p_mannwhitney": p_val,
                "significance": _sig_label(p_val),
            })

            data_groups = [df_without, df_with]
            colors      = [COLOR_WITHOUT, COLOR_WITH]
            positions   = [0, 1]

            bp = ax.boxplot(
                data_groups,
                patch_artist=True,
                positions=positions,
                widths=0.5,
                medianprops=dict(color="black", linewidth=2),
                whiskerprops=dict(linewidth=1.2),
                capprops=dict(linewidth=1.2),
                flierprops=dict(marker="o", markersize=3, alpha=0.4, linestyle="none"),
                boxprops=dict(linewidth=1.2),
            )
            for patch, color in zip(bp["boxes"], colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.75)

            rng = np.random.default_rng(42)
            for pos, vals, color in zip(positions, data_groups, colors):
                jitter = rng.uniform(-0.18, 0.18, size=len(vals))
                ax.scatter(pos + jitter, vals, color=color, s=3, alpha=0.25, zorder=3)

            ax.set_xticks([0, 1])
            ax.set_xticklabels(
                [f"{comp['without_label']}\n({n_without_v} viruses)",
                 f"{comp['with_label']}\n({n_with_v} virus{'es' if n_with_v != 1 else ''})"],
                fontsize=8.5
            )

            y_top = max(
                np.max(df_with)    if len(df_with)    > 0 else 0,
                np.max(df_without) if len(df_without) > 0 else 0,
            )
            bar_y = y_top + 0.3
            fc_str = f"{fold_change:.2f}x" if not np.isnan(fold_change) else ""
            ax.plot([0, 0, 1, 1], [bar_y, bar_y + 0.1, bar_y + 0.1, bar_y],
                    lw=1.5, color="black")
            ax.text(0.5, bar_y + 0.15,
                    f"{_sig_label(p_val)}, {fc_str}\np={p_val:.2e}",
                    ha="center", va="bottom", fontsize=9)

            ax.set_title(comp["title"], fontsize=10, pad=8)
            ax.spines[["top", "right"]].set_visible(False)

            if ax == axes[0]:
                ax.set_yticks(log2_ticks)
                ax.set_yticklabels([str(v) for v in tick_values], fontsize=9)
                ax.set_ylabel("neutralization titer", fontsize=11)

        if suptitle:
            fig.suptitle(suptitle, fontsize=11, y=1.01)
        plt.tight_layout()
        return fig, pd.DataFrame(stat_rows)

    fig_h3n2_k,   stats_h3n2_k   = _make_figure(H3N2_K_COMPARISONS)
    fig_h3n2_j24, stats_h3n2_j24 = _make_figure(H3N2_J24_COMPARISONS)
    fig_h1n1,     stats_h1n1     = _make_figure(H1N1_COMPARISONS)

    fig_h3n2_k.savefig(results_dir / "mutation_titers_H3N2_K.svg", bbox_inches="tight")
    fig_h3n2_j24.savefig(results_dir / "mutation_titers_H3N2_J24.svg", bbox_inches="tight")
    fig_h1n1.savefig(results_dir / "mutation_titers_H1N1.svg", bbox_inches="tight")

    mo.mpl.interactive(fig_h3n2_k)
    return fig_h3n2_j24, fig_h1n1, stats_h3n2_k, stats_h3n2_j24, stats_h1n1


@app.cell
def _(fig_h3n2_j24, mo):
    mo.mpl.interactive(fig_h3n2_j24)
    return


@app.cell
def _(fig_h1n1, mo):
    mo.mpl.interactive(fig_h1n1)
    return


@app.cell
def _(mo):
    mo.md("### Statistical summary")
    return


@app.cell
def _(mo, stats_h3n2_k, stats_h3n2_j24, stats_h1n1):
    mo.vstack([
        mo.md("**K subclade**"),
        stats_h3n2_k.style.format({"fold_change_with_vs_without": "{:.2f}", "p_mannwhitney": "{:.3e}"}),
        mo.md("**J.2.4 subclade**"),
        stats_h3n2_j24.style.format({"fold_change_with_vs_without": "{:.2f}", "p_mannwhitney": "{:.3e}"}),
        mo.md("**H1N1**"),
        stats_h1n1.style.format({"fold_change_with_vs_without": "{:.2f}", "p_mannwhitney": "{:.3e}"}),
    ])
    return


@app.cell
def _(mo):
    mo.md("""
    ---
    **Notes:**
    - K135E and R189K each appear in only **1 of 20 K-subclade viruses** — interpret these comparisons cautiously.
    - Sites 96/207/223 appear in **3 of 20 K-subclade viruses**.
    - All p-values are uncorrected (one pre-specified comparison per panel).
    - Titer = 40 (minimum) may be upper-bound censored (`titer_bound == "upper"`).
    """)
    return


if __name__ == "__main__":
    app.run()
