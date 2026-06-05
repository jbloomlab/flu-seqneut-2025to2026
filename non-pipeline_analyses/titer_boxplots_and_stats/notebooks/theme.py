"""
Matplotlib theme for titer boxplot figures.
Translated from the lab's Altair theme (ha-preference-shifts/notebooks/theme.py).

Usage
-----
    import theme          # applying the rcParams on import is sufficient
    from theme import format_pvalue
"""

import matplotlib as mpl
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Typography
# ---------------------------------------------------------------------------
FONT       = "Helvetica"
FONTSIZE   = 11        # base size — axes labels, tick labels, annotations
TITLE_SIZE = 11

# ---------------------------------------------------------------------------
# Apply rcParams globally
# ---------------------------------------------------------------------------
mpl.rcParams.update({
    # Font
    "font.family":          "sans-serif",
    "font.sans-serif":      ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size":            FONTSIZE,
    # Axes
    "axes.spines.top":      False,
    "axes.spines.right":    False,
    "axes.linewidth":       0.8,
    "axes.labelsize":       FONTSIZE,
    "axes.titlesize":       TITLE_SIZE,
    "axes.titleweight":     "normal",
    "axes.titlepad":        6,
    # Ticks
    "xtick.labelsize":      FONTSIZE,
    "ytick.labelsize":      FONTSIZE,
    "xtick.major.size":     3,
    "ytick.major.size":     3,
    "xtick.major.width":    0.8,
    "ytick.major.width":    0.8,
    "xtick.direction":      "out",
    "ytick.direction":      "out",
    # Grid — off by default (matching reference)
    "axes.grid":            False,
    # Figure
    "figure.dpi":           150,
    "savefig.dpi":          300,
    "savefig.bbox":         "tight",
    # Lines / patches
    "lines.linewidth":      1.0,
    "patch.linewidth":      0.8,
})


# ---------------------------------------------------------------------------
# p-value formatter — matches reference notebook style
# ---------------------------------------------------------------------------
def format_pvalue(p: float) -> str:
    """
    Format a p-value for figure annotation.

    Examples
    --------
    >>> format_pvalue(1.6e-37)
    'p < 1 × 10⁻⁵'
    >>> format_pvalue(3.2e-4)
    'p = 3 × 10⁻⁴'
    >>> format_pvalue(0.032)
    'p = 0.032'
    >>> format_pvalue(0.41)
    'p = 0.41'
    """
    superscripts = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")
    if p < 1e-5:
        return "< 1 × 10⁻⁵"
    elif p < 1e-3:
        base, exp = f"{p:.0e}".split("e")
        exp_str = str(int(exp)).translate(superscripts)
        return f"{base} × 10{exp_str}"
    else:
        return f"p = {p:.2g}"


def sig_stars(p: float) -> str:
    """Return significance stars for a p-value."""
    if p < 0.001: return "***"
    elif p < 0.01: return "**"
    elif p < 0.05: return "*"
    else: return "ns"


def make_stats_table(df, path, col_widths=None):
    """
    Save a styled matplotlib table as SVG.

    Parameters
    ----------
    df         : DataFrame to render — columns are used as headers
    path       : Path to save the SVG
    col_widths : list of column widths in inches (defaults to equal widths of 1.1")
    """
    import matplotlib.pyplot as plt

    n_rows = len(df)
    n_cols = len(df.columns)

    if col_widths is None:
        col_widths = [1.1] * n_cols

    fig_w = sum(col_widths) + 0.3
    fig_h = (n_rows + 1) * 0.38 + 0.2

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off")

    tbl = ax.table(
        cellText=df.values,
        colLabels=df.columns,
        cellLoc="center",
        loc="center",
        colWidths=[w / fig_w for w in col_widths],
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)

    for j in range(n_cols):
        cell = tbl[0, j]
        cell.set_facecolor("#2c2c2c")
        cell.set_text_props(color="white", fontweight="bold")
        cell.set_edgecolor("white")
        cell.set_height(0.13)

    sig_col = df.columns.get_loc("Sig.") if "Sig." in df.columns else None
    for i in range(n_rows):
        sig = sig_col is not None and df.iloc[i, sig_col] != "ns"
        for j in range(n_cols):
            cell = tbl[i + 1, j]
            cell.set_edgecolor("#cccccc")
            cell.set_height(0.10)
            if sig:
                cell.set_facecolor("#e8f4e8")
            elif i % 2 == 0:
                cell.set_facecolor("#f7f7f7")
            else:
                cell.set_facecolor("white")

    plt.tight_layout(pad=0.2)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
