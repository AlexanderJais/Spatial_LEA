"""Shared figure style, so every panel in the manuscript matches.

Nature-style: white ground, dark grey text, thin consistent rules, no gridlines,
no titles inside axes.  Colours are Okabe-Ito, which is colour-blind safe and
survives greyscale conversion by lightness.  Each colour carries one biological
meaning throughout the manuscript and is never reused for anything else.
"""

from __future__ import annotations

import glob

import matplotlib as mpl
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt

# Age is an ordinal variable, so it is encoded sequentially in one hue: dark for
# aged, light for adult.  This is colour-blind safe by construction, survives
# greyscale as a lightness difference, and -- unlike a blue/pink pair -- cannot
# be misread as the sex of the cohort, which is confounded with panel here.
AGED = "#08519C"
ADULT = "#9ECAE1"

# Cell identity. The population this manuscript is about; grey is everything
# else, so the eye is never asked to compare two saturated colours.
POP = "#4B2E83"
TISSUE = "#DDDDDD"      # cells not central to the message
OTHER_CELLS = "#BFBFBF"  # typed cells in a zoomed field, other than POP
INK = "#1A1A1A"
MUTED = "#666666"
RULE = "#4D4D4D"

MM = 1 / 25.4
SINGLE, ONE_HALF, FULL = 89 * MM, 120 * MM, 183 * MM   # Nature column widths


def use_style() -> None:
    for path in glob.glob("/usr/share/fonts/**/NimbusSans*.otf", recursive=True):
        try:
            fm.fontManager.addfont(path)
        except Exception:
            pass
    family = ("Nimbus Sans" if any(f.name == "Nimbus Sans" for f in fm.fontManager.ttflist)
              else "DejaVu Sans")
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [family, "Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 6.5,
        "axes.labelsize": 6.5, "axes.titlesize": 6.5,
        "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
        "text.color": INK, "axes.labelcolor": INK,
        "xtick.color": RULE, "ytick.color": RULE,
        "axes.edgecolor": RULE, "axes.linewidth": .5,
        "xtick.major.width": .5, "ytick.major.width": .5,
        "xtick.major.size": 2.2, "ytick.major.size": 2.2,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": False,
        "legend.frameon": False,
        "figure.facecolor": "white", "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "figure.dpi": 400, "savefig.dpi": 400, "savefig.bbox": "tight",
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "lines.linewidth": .8,
    })


def panel(ax, letter: str, dx: float = -0.24, dy: float = 1.10) -> None:
    """Bold lower-case panel letter, positioned in axes fraction."""
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=8,
            fontweight="bold", va="top", ha="left", color=INK)


def scalebar(ax, length_um: float, label: str, frac=(0.06, 0.06), lw=1.2) -> None:
    """Horizontal scale bar in data units, placed by axes fraction."""
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
    x = x0 + (x1 - x0) * frac[0]
    y = y0 + (y1 - y0) * frac[1]
    ax.plot([x, x + length_um], [y, y], lw=lw, color=INK, solid_capstyle="butt",
            zorder=10)
    ax.text(x + length_um / 2, y + (y1 - y0) * 0.025, label, ha="center",
            va="bottom", fontsize=5.5, color=INK, zorder=10)


def bare(ax) -> None:
    """Spatial panels carry a scale bar, not axes."""
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
