import marimo

__generated_with = "0.17.2"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from scipy import stats
    from statsmodels.stats.multitest import multipletests
    import itertools
    import re
    from pathlib import Path
    return Path, itertools, mo, multipletests, np, pd, plt, stats


@app.cell
def _(mo):
    mo.md(r"""
    # Subclade titer boxplots

    Neutralization titers for all 302 'All'-group human sera against viruses
    from subclades of interest, separately for H3N2 and H1N1.

    **Statistics (both subtypes):**
    - Primary test: hypothesis subclade vs all others combined (Mann-Whitney U, one-sided, uncorrected)
    - Pairwise tests: all subclade pairs (Mann-Whitney U, BH-corrected)

    Each data point is one (serum × virus) measurement, log₂-transformed.
    Fold changes are ratios of geometric means.
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

    titers_path = DATA_DIR / "human_titers.csv"
    viruses_path = DATA_DIR / "human_viruses.csv"

    titers_raw = pd.read_csv(titers_path)
    viruses = pd.read_csv(viruses_path)

    results_dir = HERE / "results"
    results_dir.mkdir(exist_ok=True)
    return titers_raw, viruses, results_dir


@app.cell
def _(mo, titers_raw, viruses):
    mo.md(f"""
    **Data loaded:**
    - {len(titers_raw):,} titer measurements ({titers_raw['serum'].nunique()} sera × {titers_raw['virus'].nunique()} viruses)
    - {len(viruses)} viruses with metadata
    """)
    return


@app.cell
def _(np, titers_raw, viruses):
    # ---------------------------------------------------------------------------
    # H3N2 — filter to subclades of interest
    # ---------------------------------------------------------------------------
    H3N2_SUBCLADE_ORDER = ["K", "J.2.4", "J.2.3", "J.2.2", "J.2"]

    # Wong (2011) colorblind-safe palette
    H3N2_SUBCLADE_COLORS = {
        "K":     "#0072B2",
        "J.2.4": "#009E73",
        "J.2.3": "#D55E00",
        "J.2.2": "#CC79A7",
        "J.2":   "#56B4E9",
    }

    virus_meta_h3n2 = viruses[viruses["subclade"].isin(H3N2_SUBCLADE_ORDER)][
        ["virus", "subclade"]
    ].copy()

    df_h3n2 = titers_raw.merge(virus_meta_h3n2, on="virus", how="inner")
    df_h3n2["log2_titer"] = np.log2(df_h3n2["titer"])

    counts_h3n2 = (
        df_h3n2.groupby("subclade")
        .agg(n_viruses=("virus", "nunique"), n_sera=("serum", "nunique"), n_obs=("titer", "count"))
        .reindex(H3N2_SUBCLADE_ORDER)
    )
    return H3N2_SUBCLADE_COLORS, H3N2_SUBCLADE_ORDER, counts_h3n2, df_h3n2


@app.cell
def _(mo):
    mo.md("## H3N2\n### Observations per subclade")
    return


@app.cell
def _(counts_h3n2):
    counts_h3n2
    return


@app.cell
def _(df_h3n2, H3N2_SUBCLADE_ORDER, itertools, mo, multipletests, np, pd, stats):
    # ---------------------------------------------------------------------------
    # H3N2 statistics
    # ---------------------------------------------------------------------------

    # 1. Primary: K vs all others combined (uncorrected, one-sided)
    k_titers = df_h3n2.loc[df_h3n2["subclade"] == "K", "log2_titer"].values
    other_titers_h3n2 = df_h3n2.loc[df_h3n2["subclade"] != "K", "log2_titer"].values
    _, p_k_vs_all = stats.mannwhitneyu(k_titers, other_titers_h3n2, alternative="less")
    fold_k_vs_all = 2 ** (np.mean(k_titers) - np.mean(other_titers_h3n2))

    # 2. All pairwise Mann-Whitney (BH corrected)
    subclades_present_h3n2 = df_h3n2["subclade"].unique()
    pairs_h3n2 = list(itertools.combinations(
        [s for s in H3N2_SUBCLADE_ORDER if s in subclades_present_h3n2], 2
    ))
    pairwise_raw_h3n2 = []
    for s1, s2 in pairs_h3n2:
        g1 = df_h3n2.loc[df_h3n2["subclade"] == s1, "log2_titer"].values
        g2 = df_h3n2.loc[df_h3n2["subclade"] == s2, "log2_titer"].values
        _, p = stats.mannwhitneyu(g1, g2, alternative="two-sided")
        pairwise_raw_h3n2.append(p)

    reject_h3n2, pvals_bh_h3n2, _, _ = multipletests(pairwise_raw_h3n2, method="fdr_bh")

    fold_changes_h3n2 = [
        2 ** (np.mean(df_h3n2.loc[df_h3n2["subclade"] == s1, "log2_titer"].values) -
              np.mean(df_h3n2.loc[df_h3n2["subclade"] == s2, "log2_titer"].values))
        for s1, s2 in pairs_h3n2
    ]

    pairwise_df_h3n2 = pd.DataFrame({
        "subclade_1": [p[0] for p in pairs_h3n2],
        "subclade_2": [p[1] for p in pairs_h3n2],
        "fold_change_1_vs_2": fold_changes_h3n2,
        "p_raw": pairwise_raw_h3n2,
        "p_bh": pvals_bh_h3n2,
        "significant": reject_h3n2,
    }).sort_values("p_bh")

    mo.md(f"""
    ### Primary test: K vs all others combined
    Mann-Whitney U (one-sided, K < others): **{fold_k_vs_all:.2f}×, p = {p_k_vs_all:.2e}**
    *(uncorrected — pre-specified directional hypothesis)*
    """)
    return fold_k_vs_all, p_k_vs_all, pairwise_df_h3n2


@app.cell
def _(mo):
    mo.md("### Pairwise tests (BH-corrected)")
    return


@app.cell
def _(pairwise_df_h3n2):
    pairwise_df_h3n2.style.format({"fold_change_1_vs_2": "{:.2f}", "p_raw": "{:.3e}", "p_bh": "{:.3e}"})
    return


@app.cell
def _(
    H3N2_SUBCLADE_COLORS,
    H3N2_SUBCLADE_ORDER,
    df_h3n2,
    fold_k_vs_all,
    mo,
    np,
    p_k_vs_all,
    pairwise_df_h3n2,
    plt,
    results_dir,
    stats,
):
    # ---------------------------------------------------------------------------
    # Shared figure helper — reused for H3N2 and H1N1
    # ---------------------------------------------------------------------------
    def _sig_label(p):
        if p < 0.001:   return "***"
        elif p < 0.01:  return "**"
        elif p < 0.05:  return "*"
        else:           return "ns"

    def make_subclade_figure(
        df, subclade_order, subclade_colors,
        primary_subclade, pairwise_df, p_primary, fold_primary,
        title=None,
    ):
        """
        Draw a subclade boxplot figure.

        Parameters
        ----------
        df               : tidy dataframe with 'subclade' and 'log2_titer' columns
        subclade_order   : list of subclades in display order
        subclade_colors  : dict mapping subclade → hex color
        primary_subclade : subclade for the primary hypothesis test (shown as spanning bar)
        pairwise_df      : DataFrame with cols subclade_1, subclade_2, fold_change_1_vs_2,
                           p_bh, significant
        p_primary        : p-value for primary test (uncorrected)
        fold_primary     : fold change for primary test
        title            : optional figure title string
        """
        tick_values = [40, 80, 160, 320, 640, 1280, 2560, 5120, 10240, 16384]
        log2_ticks  = np.log2(tick_values)
        y_data_min  = log2_ticks[0]
        y_data_max  = log2_ticks[-1]

        subclades_ordered = [s for s in subclade_order if s in df["subclade"].unique()]
        grouped = [df.loc[df["subclade"] == s, "log2_titer"].values for s in subclades_ordered]

        fig, ax = plt.subplots(figsize=(7, 5))

        bp = ax.boxplot(
            grouped,
            patch_artist=True,
            positions=range(len(subclades_ordered)),
            widths=0.55,
            medianprops=dict(color="black", linewidth=2),
            whiskerprops=dict(linewidth=1.2),
            capprops=dict(linewidth=1.2),
            flierprops=dict(marker="o", markersize=2.5, alpha=0.4, linestyle="none"),
            boxprops=dict(linewidth=1.2),
        )
        for patch, sc in zip(bp["boxes"], subclades_ordered):
            patch.set_facecolor(subclade_colors[sc])
            patch.set_alpha(0.75)

        rng = np.random.default_rng(42)
        for i, (sc, vals) in enumerate(zip(subclades_ordered, grouped)):
            jitter = rng.uniform(-0.22, 0.22, size=len(vals))
            ax.scatter(i + jitter, vals, color=subclade_colors[sc], s=3, alpha=0.25, zorder=3)

        ax.set_yticks(log2_ticks)
        ax.set_yticklabels([str(v) for v in tick_values], fontsize=9)
        ax.set_ylabel("neutralization titer", fontsize=11)

        ax.set_xticks(range(len(subclades_ordered)))
        n_labels = [
            f"{sc}\n(n={df[df['subclade']==sc]['virus'].nunique()} viruses)"
            for sc in subclades_ordered
        ]
        ax.set_xticklabels(n_labels, fontsize=10)
        ax.set_xlabel("subclade", fontsize=11)

        if title:
            ax.set_title(title, fontsize=11)

        # Primary hypothesis annotation
        p_idx     = subclades_ordered.index(primary_subclade)
        last_idx  = len(subclades_ordered) - 1
        y_top     = ax.get_ylim()[1]
        bar_y     = y_top + 0.15
        ax.annotate("", xy=(last_idx, bar_y), xytext=(p_idx, bar_y),
                    arrowprops=dict(arrowstyle="-", lw=1.5, color="black"))
        ax.text(
            (p_idx + last_idx) / 2, bar_y + 0.05,
            f"{primary_subclade} vs others: {_sig_label(p_primary)}, {fold_primary:.2f}× "
            f"(p={p_primary:.2e}, uncorrected)",
            ha="center", va="bottom", fontsize=8.5,
        )

        # BH-corrected pairwise annotations (significant pairs only)
        sig_pairs   = pairwise_df[pairwise_df["significant"]]
        step        = 0.55
        annotation_y = bar_y + 0.55
        for _, row in sig_pairs.iterrows():
            if row["subclade_1"] not in subclades_ordered or row["subclade_2"] not in subclades_ordered:
                continue
            x1 = subclades_ordered.index(row["subclade_1"])
            x2 = subclades_ordered.index(row["subclade_2"])
            ax.annotate("", xy=(x2, annotation_y), xytext=(x1, annotation_y),
                        arrowprops=dict(arrowstyle="-", lw=1.2, color="#444444"))
            ax.text(
                (x1 + x2) / 2, annotation_y + 0.04,
                f"{_sig_label(row['p_bh'])}, {row['fold_change_1_vs_2']:.2f}× "
                f"(BH p={row['p_bh']:.2e})",
                ha="center", va="bottom", fontsize=7.5, color="#444444",
            )
            annotation_y += step

        ax.set_ylim(ax.get_ylim()[0], annotation_y + 0.3)
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines["left"].set_bounds(y_data_min, y_data_max)
        y_label_center = (y_data_min + y_data_max) / 2
        ax.yaxis.set_label_coords(
            -0.08,
            (y_label_center - ax.get_ylim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0]),
        )
        plt.tight_layout()
        return fig

    # ---------------------------------------------------------------------------
    # H3N2 figure
    # ---------------------------------------------------------------------------
    fig_h3n2 = make_subclade_figure(
        df=df_h3n2,
        subclade_order=H3N2_SUBCLADE_ORDER,
        subclade_colors=H3N2_SUBCLADE_COLORS,
        primary_subclade="K",
        pairwise_df=pairwise_df_h3n2,
        p_primary=p_k_vs_all,
        fold_primary=fold_k_vs_all,
    )

    fig_h3n2.savefig(results_dir / "subclade_titers_H3N2.svg", bbox_inches="tight")
    mo.mpl.interactive(fig_h3n2)
    return make_subclade_figure,


@app.cell
def _(mo):
    mo.md("---\n## H1N1\n### Observations per subclade")
    return


@app.cell
def _(np, titers_raw, viruses):
    # ---------------------------------------------------------------------------
    # H1N1 — D.3.1 and D.3.1.1 only (other subclades are singletons)
    # ---------------------------------------------------------------------------
    H1N1_SUBCLADE_ORDER = ["D.3.1.1", "D.3.1"]

    # Restart color scale — remaining Wong colors not used in H3N2
    H1N1_SUBCLADE_COLORS = {
        "D.3.1.1": "#E69F00",   # orange
        "D.3.1":   "#000000",   # black
    }

    virus_meta_h1n1 = viruses[viruses["subclade"].isin(H1N1_SUBCLADE_ORDER)][
        ["virus", "subclade"]
    ].copy()

    df_h1n1 = titers_raw.merge(virus_meta_h1n1, on="virus", how="inner")
    df_h1n1["log2_titer"] = np.log2(df_h1n1["titer"])

    counts_h1n1 = (
        df_h1n1.groupby("subclade")
        .agg(n_viruses=("virus", "nunique"), n_sera=("serum", "nunique"), n_obs=("titer", "count"))
        .reindex(H1N1_SUBCLADE_ORDER)
    )
    return H1N1_SUBCLADE_COLORS, H1N1_SUBCLADE_ORDER, counts_h1n1, df_h1n1


@app.cell
def _(counts_h1n1):
    counts_h1n1


@app.cell
def _(df_h1n1, H1N1_SUBCLADE_ORDER, mo, multipletests, np, pd, stats):
    # ---------------------------------------------------------------------------
    # H1N1 statistics
    # ---------------------------------------------------------------------------

    # 1. Primary: D.3.1.1 vs D.3.1 (one-sided, uncorrected)
    d311_titers = df_h1n1.loc[df_h1n1["subclade"] == "D.3.1.1", "log2_titer"].values
    d31_titers  = df_h1n1.loc[df_h1n1["subclade"] == "D.3.1",   "log2_titer"].values
    _, p_d311_vs_d31 = stats.mannwhitneyu(d311_titers, d31_titers, alternative="less")
    fold_d311_vs_d31 = 2 ** (np.mean(d311_titers) - np.mean(d31_titers))

    # 2. Pairwise (only one pair here, so BH correction is trivial but kept for consistency)
    reject_h1n1, pvals_bh_h1n1, _, _ = multipletests([p_d311_vs_d31], method="fdr_bh")

    pairwise_df_h1n1 = pd.DataFrame({
        "subclade_1": ["D.3.1.1"],
        "subclade_2": ["D.3.1"],
        "fold_change_1_vs_2": [fold_d311_vs_d31],
        "p_raw": [p_d311_vs_d31],
        "p_bh": pvals_bh_h1n1,
        "significant": reject_h1n1,
    })

    mo.md(f"""
    ### Primary test: D.3.1.1 vs D.3.1
    Mann-Whitney U (one-sided, D.3.1.1 < D.3.1): **{fold_d311_vs_d31:.2f}×, p = {p_d311_vs_d31:.2e}**
    *(uncorrected — pre-specified directional hypothesis)*
    """)
    return fold_d311_vs_d31, p_d311_vs_d31, pairwise_df_h1n1


@app.cell
def _(
    H1N1_SUBCLADE_COLORS,
    H1N1_SUBCLADE_ORDER,
    df_h1n1,
    fold_d311_vs_d31,
    make_subclade_figure,
    mo,
    p_d311_vs_d31,
    pairwise_df_h1n1,
    results_dir,
):
    fig_h1n1 = make_subclade_figure(
        df=df_h1n1,
        subclade_order=H1N1_SUBCLADE_ORDER,
        subclade_colors=H1N1_SUBCLADE_COLORS,
        primary_subclade="D.3.1.1",
        pairwise_df=pairwise_df_h1n1,
        p_primary=p_d311_vs_d31,
        fold_primary=fold_d311_vs_d31,
    )
    fig_h1n1.savefig(results_dir / "subclade_titers_H1N1.svg", bbox_inches="tight")
    mo.mpl.interactive(fig_h1n1)
    return


@app.cell
def _(mo):
    mo.md("""
    ---
    *Note: Each data point is one serum–virus measurement.
    Titer values at the minimum (40) may represent upper-bound censored observations (`titer_bound == "upper"`).
    Statistical tests use log₂-transformed titers.*
    """)
    return


if __name__ == "__main__":
    app.run()
