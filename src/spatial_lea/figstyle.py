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
# The receptor this study is about, wherever it is shown as a measurement rather
# than as a group difference: two tints, so a cell carrying one transcript and a
# cell carrying several are not drawn as the same observation.
GALR1 = "#B4531A"
GALR1_LOW = "#E3B48F"
TISSUE = "#DDDDDD"      # cells not central to the message
OTHER_CELLS = "#BFBFBF"  # typed cells in a zoomed field, other than POP
INK = "#1A1A1A"
MUTED = "#1A1A1A"      # panel sub-headings are set in ink, not grey
RULE = "#4D4D4D"

# Registered hypothalamic subregions.  Muted, distinguishable, and always
# accompanied by a direct label on the map, so colour is never the only channel
# carrying the identity.  Territory that is not a named nucleus stays grey.
NUCLEUS_COLOUR = {
    "ARC": "#9C6114",
    "ME_3V": "#C9A227",
    "VMH": "#2C6FA8",
    "DMH": "#6A3D9A",
    "LHA": "#2E7D5B",
    "ZI": "#A34E6B",
    "DHA_PH": "#5B6B7A",
    "TUseg": "#CFCFCF",
    "fibre": "#E2E2E2",
    "edge": "#EFEFEF",
}
NUCLEUS_LABEL = {
    "ARC": "ARC", "ME_3V": "ME / 3V", "VMH": "VMH", "DMH": "DMH",
    "LHA": "LHA", "ZI": "ZI", "DHA_PH": "DHA / PH",
}

MM = 1 / 25.4
SINGLE, ONE_HALF, FULL = 89 * MM, 120 * MM, 183 * MM   # Nature column widths


def use_style() -> None:
    for pattern in ("/usr/share/fonts/**/NimbusSans*.otf",
                    "/usr/share/fonts/**/NimbusSans-*.otf"):
        for path in glob.glob(pattern, recursive=True):
            try:
                fm.fontManager.addfont(path)
            except Exception:
                pass
    # Nimbus Sans is the Helvetica metric clone; Liberation Sans is the Arial
    # one and is the next best thing.  DejaVu is matplotlib's default and looks
    # nothing like either, so it is the last resort rather than a silent one.
    have = {f.name for f in fm.fontManager.ttflist}
    family = next((f for f in ("Nimbus Sans", "Helvetica", "Liberation Sans",
                               "Arial", "DejaVu Sans") if f in have), "DejaVu Sans")
    if family != "Nimbus Sans":
        print(f"[figstyle] Nimbus Sans not installed; using {family}. "
              "apt install fonts-urw-base35")
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
    """Bold lower-case panel letter, positioned in axes fraction.

    Prefer ``panel_letters`` for a multi-panel figure: positioning each letter
    against its own axes puts them at ragged heights and offsets, because axes
    differ in how much room their tick labels take and an equal-aspect panel
    shrinks inside the cell it was given.
    """
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=8,
            fontweight="bold", va="top", ha="left", color=INK)


def panel_letters(fig, pairs, dx: float = 0.0, dy: float = 0.012) -> None:
    """Place panel letters on the grid rather than on the axes.

    Each letter goes at the top-left of the cell its panel was allotted, so
    letters in a row share a baseline and letters in a column share a margin,
    whatever the axes inside those cells do.
    """
    for letter, ax in pairs:
        cell = ax.get_subplotspec().get_position(fig)
        fig.text(cell.x0 + dx, cell.y1 + dy, letter, fontsize=8,
                 fontweight="bold", va="bottom", ha="left", color=INK)


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


def wash(colour: str, amount: float = .55) -> tuple:
    """Mix a colour towards white.

    Subregions are the ground of a spatial panel, not its subject, so they are
    drawn as washes and full saturation is left to the cells plotted on them.
    """
    rgb = [int(colour[i:i + 2], 16) / 255 for i in (1, 3, 5)]
    return tuple(c + (1 - c) * amount for c in rgb)
