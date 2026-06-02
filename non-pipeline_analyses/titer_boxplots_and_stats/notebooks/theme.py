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
        return "p < 1 × 10⁻⁵"
    elif p < 1e-3:
        base, exp = f"{p:.0e}".split("e")
        exp_str = str(int(exp)).translate(superscripts)
        return f"p = {base} × 10{exp_str}"
    else:
        return f"p = {p:.2g}"
