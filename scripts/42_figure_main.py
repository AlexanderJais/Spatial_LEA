#!/usr/bin/env python3
"""Figure 1: which MBH cells change Galr1 between adult and aged mice.

    Galr1 is expressed across the mediobasal hypothalamus.  In which cells does
    it change with age?

    One: a glutamatergic Otp+/Cbln1+/Prdm8+ cell type of the dorsal
    hypothalamus.  What that population is, and what else changes in it, is
    Figure 2 (45_figure_population.py).

Panels:
  a  the registered subregions, every cell from all eight animals pooled in the
     shared anatomical frame, with the population that changes on top
  b  Galr1 in adult against aged, in every cell type the design can test.  Group
     means joined, one small point per animal; this is the comparison the study
     is about, so it is shown for all cell types at once and not summarised to
     a fold change
  c  the same data per animal: Galr1 z-scored across the eight animals, adult
     block then aged block, cell types clustered by their pattern.  A cell type
     that separates the groups shows as a left-right split in its row
  d  Galr1 in the cell type that changes, one point per animal, paired within
     block

Choices made by rule rather than by eye:
  * a cell type enters the screen when every animal contributes at least
    MIN_CELLS_PER_ANIMAL of it, and when at least MIN_PCT_POS of its cells carry
    Galr1 at all.  The first is the blocked design: a type one animal lacks
    cannot be tested.  The second is detection: a fold change in a gene found in
    a few per cent of cells measures the detection process, not regulation.
    Both are applied before any statistics are looked at, and what they exclude
    is printed.
  * "changes" means all four blocks agreeing in sign, aged and adult ranges not
    overlapping, and a permutation P <= 0.05.  0.0286 is the smallest P eight
    animals can produce, so this is the strictest the design allows.
  * panel a paints the subregions largest first.  Pooled over eight sections,
    cells share pixels and whichever region is drawn last takes them; in the
    order the colour table lists them the DMH was painted over by LHA, ZI and
    DHA_PH, its three larger neighbours.
  * panels b, c, e are computed here, so the figure regenerates from the data.

Adult is plotted first throughout, as the reference condition.  Per-animal
panels are drawn on a linear CPM axis from zero: a log axis cropped to the data
makes any difference look like whatever the crop chose.
"""

from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.cluster.hierarchy import leaves_list, linkage, optimal_leaf_ordering
from scipy.spatial.distance import pdist

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.figstyle import (  # noqa: E402
    ADULT, AGED as C_AGED, FULL, GALR1, INK, NUCLEUS_COLOUR, NUCLEUS_LABEL, POP,
    TISSUE, bare, panel, scalebar, use_style,
)
from spatial_lea.io import (  # noqa: E402
    ADULT as A_ADULT, AGED as A_AGED, BLOCKS, ONE_PER_MOUSE, counts_matrix,
)

PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "figures"
SRC = OUT / "source_data"
POPNAME = "Glut Prdm8/Cbln1"
ANIMALS = list(A_ADULT) + list(A_AGED)
# Entry criteria for the screen, both applied before any statistic is looked at.
MIN_CELLS_PER_ANIMAL = 30
MIN_PCT_POS = 10.0
# What "changes" means: every block agreeing, the two groups not overlapping,
# and the smallest P eight animals can give.
MAX_P = 0.05
# Diverging, colour-blind safe, and mapped through zero so no shift of colour
# happens at a value the data did not cross.
LFC_CMAP = LinearSegmentedColormap.from_list("lfc", ["#2C6FA8", "#F4F4F4", GALR1])


def wash(colour: str, amount: float = .55) -> tuple:
    """Mix a colour towards white.

    Subregions are the ground of panel a, not its subject, so they are drawn as
    washes and full saturation is left to the cells plotted on them.  Without
    this the DMH's colour and the population's are both purple at full strength
    and the population disappears into the region it sits in.
    """
    rgb = np.array([int(colour[i:i + 2], 16) / 255 for i in (1, 3, 5)])
    return tuple(rgb + (1 - rgb) * amount)


def blocked_stats(cpm: pd.DataFrame) -> pd.DataFrame:
    lfc = pd.DataFrame({b: cpm.loc[x] - cpm.loc[y] for b, (x, y) in BLOCKS.items()}).T
    mean = lfc.mean()
    agree = (np.sign(lfc) == np.sign(mean)).sum()
    vals = cpm.loc[ANIMALS].to_numpy()
    n_ad = len(A_ADULT)
    idx = list(range(len(ANIMALS)))
    obs = vals[n_ad:].mean(axis=0) - vals[:n_ad].mean(axis=0)
    null = np.array([vals[list(c)].mean(axis=0)
                     - vals[[i for i in idx if i not in c]].mean(axis=0)
                     for c in combinations(idx, n_ad)])
    p = (np.abs(null) >= np.abs(obs)).mean(axis=0)
    aged_v, adult_v = vals[n_ad:], vals[:n_ad]
    margin = np.maximum(aged_v.min(axis=0) - adult_v.max(axis=0),
                        adult_v.min(axis=0) - aged_v.max(axis=0))
    out = pd.DataFrame({"lfc": mean, "blocks": agree, "p": p, "margin": margin})
    for b in BLOCKS:
        out[f"lfc_{b}"] = lfc.loc[b]
    return out


def galr1_screen(counts: np.ndarray, var: np.ndarray, cell_types: np.ndarray,
                 animals: np.ndarray) -> tuple[pd.DataFrame, pd.DataFrame, list]:
    """Galr1 in every cell type the design can test, by one rule throughout.

    Returns the per-type statistics, the per-animal Galr1 level behind them, and
    the types the entry criteria excluded -- the last so it is reported rather
    than left implicit.
    """
    gene = int(np.flatnonzero(var == "Galr1")[0])
    rows, per_animal_cpm, excluded = [], {}, []
    for name in pd.unique(cell_types):
        if str(name).startswith(("unlabelled", "unresolved")):
            continue
        sel = cell_types == name
        n_per_animal = [int((sel & (animals == a)).sum()) for a in ANIMALS]
        pct_pos = float((counts[sel][:, gene] > 0).mean() * 100)
        if min(n_per_animal) < MIN_CELLS_PER_ANIMAL:
            excluded.append((name, f"{min(n_per_animal)} cells in the thinnest animal"))
            continue
        if pct_pos < MIN_PCT_POS:
            excluded.append((name, f"{pct_pos:.1f}% of cells Galr1+"))
            continue
        mat = pd.DataFrame({a: counts[sel & (animals == a)].sum(axis=0)
                            for a in ANIMALS}, index=var).T
        linear = mat.div(mat.sum(axis=1), axis=0) * 1e6
        r = blocked_stats(np.log2(linear + 1)).loc["Galr1"]
        per_animal_cpm[name] = linear["Galr1"]
        rows.append({"cell_type": name, "n_cells": int(sel.sum()),
                     "min_cells": min(n_per_animal), "pct_pos": round(pct_pos, 1),
                     "cpm": round(float(counts[sel][:, gene].sum()
                                        / counts[sel].sum() * 1e6), 1),
                     "lfc": round(float(r.lfc), 3), "fold": round(float(2 ** r.lfc), 2),
                     "blocks": int(r.blocks), "p": round(float(r.p), 4),
                     "margin": round(float(r.margin), 3)})
    out = pd.DataFrame(rows)
    out["significant"] = ((out.blocks == 4) & (out.margin > 0) & (out.p <= MAX_P))
    out = out.sort_values("lfc", ascending=False)
    return out, pd.DataFrame(per_animal_cpm).T.loc[out.cell_type], excluded


def paired_panel(ax, values: dict, ylabel: str, title: str) -> None:
    """One point per animal, adult joined to its block partner.

    Linear CPM from zero.  The effect is what it is; an axis cropped to the
    points would set its apparent size by the crop.
    """
    for a_aged, a_adult in BLOCKS.values():
        ax.plot([0, 1], [values[a_adult], values[a_aged]], color="#C9C9C9",
                lw=.6, zorder=1)
    for j, (members, colour) in enumerate(((A_ADULT, ADULT), (A_AGED, C_AGED))):
        ys = [values[a] for a in members]
        ax.scatter([j] * len(ys), ys, s=15, color=colour, zorder=3, linewidths=0)
        ax.plot([j - .22, j + .22], [np.mean(ys)] * 2, color=colour, lw=1.4,
                solid_capstyle="butt", zorder=2)
    ax.set_xlim(-.45, 1.45); ax.set_xticks([0, 1])
    ax.set_xticklabels(["adult", "aged"])
    ax.set_ylim(0, max(values.values()) * 1.16)
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", pad=3)


def main() -> int:
    use_style()
    OUT.mkdir(parents=True, exist_ok=True); SRC.mkdir(parents=True, exist_ok=True)

    win = sc.read_h5ad(PROC / "hypothalamus_nuclei.h5ad")
    win = win[win.obs["section"].astype(str).isin(ONE_PER_MOUSE)].copy()
    cell_types = win.obs["cell_type"].astype(str).to_numpy()
    if POPNAME not in set(cell_types):
        raise SystemExit(f"'{POPNAME}' is not among this run's cell types; see "
                         "results/annotation/cluster_label_matching.csv")
    is_pop = cell_types == POPNAME
    animals = win.obs["animal"].astype(str).to_numpy()
    ml, dv = win.obs["ml"].to_numpy(), win.obs["dv"].to_numpy()
    nuc = win.obs["nucleus_ext"].astype(str).to_numpy()
    counts = counts_matrix(win)
    var = win.var_names.to_numpy()

    screen, per_animal, excluded = galr1_screen(counts, var, cell_types, animals)
    screen.to_csv(SRC / "galr1_by_celltype.csv", index=False)
    per_animal.to_csv(SRC / "galr1_by_celltype_per_animal.csv")

    mat = pd.DataFrame({a: counts[is_pop & (animals == a)].sum(axis=0)
                        for a in ANIMALS}, index=var).T
    pop_cpm = mat.div(mat.sum(axis=1), axis=0) * 1e6
    pop_stats = blocked_stats(np.log2(pop_cpm + 1))[mat.sum(axis=0) >= 200]
    pop_stats.to_csv(SRC / "fig_main_stats.csv")
    pd.DataFrame({"abundance": {a: (is_pop & (animals == a)).sum()
                                / (animals == a).sum() * 100 for a in ANIMALS}}
                 ).loc[ANIMALS].to_csv(SRC / "fig_main_abundance.csv")

    fig = plt.figure(figsize=(FULL, 3.05))
    gs = fig.add_gridspec(1, 4, width_ratios=[1.55, 1.15, 1.25, .68], wspace=.92)

    # (a) the registered subregions, pooled, with the population on them
    ax = fig.add_subplot(gs[0]); panel(ax, "a", dx=-0.04, dy=1.10)
    present = [k for k in NUCLEUS_LABEL if (nuc == k).any()]
    for key in ("edge", "fibre", "TUseg"):
        m = nuc == key
        if m.any():
            ax.scatter(ml[m], dv[m], s=.45, c=TISSUE, linewidths=0, rasterized=True)
    for key in sorted(present, key=lambda k: int((nuc == k).sum()), reverse=True):
        m = nuc == key
        ax.scatter(ml[m], dv[m], s=.45, c=[wash(NUCLEUS_COLOUR[key])],
                   linewidths=0, rasterized=True)
    ax.scatter(ml[is_pop], dv[is_pop], s=4.0, c=POP, linewidths=.2,
               edgecolors="white", rasterized=True)
    for key in present:
        m = nuc == key
        side = ml[m] > 0 if (ml[m] > 0).sum() > 30 else ml[m] < 0
        ax.annotate(NUCLEUS_LABEL[key],
                    (np.median(ml[m][side]), np.median(dv[m][side])),
                    fontsize=5.4, color=INK, ha="center", va="center",
                    fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.10", fc="white", ec="none",
                              alpha=.72))
    ax.set_xlim(-1500, 1500); ax.set_ylim(-100, 1800)
    ax.set_aspect("equal"); bare(ax)
    scalebar(ax, 500, "500 µm")
    ax.legend(handles=[Line2D([], [], marker="o", linestyle="none", markersize=2.4,
                              color=POP, label=POPNAME)],
              loc="upper left", bbox_to_anchor=(0, .04), frameon=False,
              fontsize=5.4, handletextpad=.3, borderpad=0, borderaxespad=0)
    ax.set_title("subregions, 8 animals pooled", loc="left", pad=2, x=.04)

    # (b) adult against aged, in every cell type tested
    ax = fig.add_subplot(gs[1]); panel(ax, "b", dx=-0.62, dy=1.10)
    order = screen.sort_values("lfc")
    ys = np.arange(len(order))
    for y, name in zip(ys, order.cell_type):
        vals = per_animal.loc[name]
        a_mean, g_mean = vals[list(A_ADULT)].mean(), vals[list(A_AGED)].mean()
        ax.plot([a_mean, g_mean], [y, y], color="#C9C9C9", lw=.9, zorder=1)
        ax.scatter(vals[list(A_ADULT)], [y] * len(A_ADULT), s=3.5, color=ADULT,
                   linewidths=0, zorder=2, alpha=.75)
        ax.scatter(vals[list(A_AGED)], [y] * len(A_AGED), s=3.5, color=C_AGED,
                   linewidths=0, zorder=2, alpha=.75)
        ax.scatter([a_mean], [y], s=17, color=ADULT, linewidths=0, zorder=3)
        ax.scatter([g_mean], [y], s=17, color=C_AGED, linewidths=0, zorder=3)
    ax.set_yticks(ys)
    ax.set_yticklabels([f"{n} *" if s else n
                        for n, s in zip(order.cell_type, order.significant)],
                       fontsize=5.2)
    for tick, sig in zip(ax.get_yticklabels(), order.significant):
        tick.set_fontweight("bold" if sig else "normal")
        tick.set_color(POP if sig else INK)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-.7, len(order) - .3)
    ax.set_xlim(0, float(per_animal.to_numpy().max()) * 1.06)
    ax.set_xlabel("$\\it{Galr1}$ (CPM)")
    # Named in the corner in their own colours rather than in a legend box:
    # every row carries points, so a box placed anywhere inside sits on data.
    # Placed against the rows whose points sit at low CPM, so the words land on
    # empty axis rather than on data or on panel c's row labels.
    for y, label, colour in ((.46, "aged", C_AGED), (.39, "adult", ADULT)):
        ax.annotate(label, xy=(.97, y), xycoords="axes fraction", fontsize=5.8,
                    color=colour, fontweight="bold", ha="right", va="top")
    ax.set_title("adult vs aged, per cell type", loc="left", pad=2)

    # (c) the same data, per animal
    ax = fig.add_subplot(gs[2]); panel(ax, "c", dx=-0.60, dy=1.10)
    z = per_animal.sub(per_animal.mean(axis=1), axis=0).div(
        per_animal.std(axis=1).replace(0, np.nan), axis=0)
    link = optimal_leaf_ordering(linkage(z.to_numpy(), "average"), pdist(z.to_numpy()))
    z = z.iloc[leaves_list(link)]
    lim = float(np.nanmax(np.abs(z.to_numpy())))
    im = ax.imshow(z.to_numpy(), aspect="auto", cmap=LFC_CMAP,
                   norm=TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim))
    ax.axvline(len(A_ADULT) - .5, color="white", lw=1.4)
    ax.set_xticks(range(len(ANIMALS)))
    ax.set_xticklabels(ANIMALS, fontsize=4.8, rotation=90)
    for x, label in ((len(A_ADULT) / 2 - .5, "adult"),
                     (len(A_ADULT) + len(A_AGED) / 2 - .5, "aged")):
        ax.annotate(label, xy=(x, -.85), xycoords=("data", "data"), fontsize=5.6,
                    ha="center", va="bottom", color=INK, annotation_clip=False)
    sig = screen.set_index("cell_type").significant.reindex(z.index)
    ax.set_yticks(range(len(z)))
    ax.set_yticklabels([f"{n} *" if s else n for n, s in zip(z.index, sig)],
                       fontsize=5.2)
    for tick, s_ in zip(ax.get_yticklabels(), sig):
        tick.set_fontweight("bold" if s_ else "normal")
        tick.set_color(POP if s_ else INK)
    ax.tick_params(axis="both", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=.05, pad=.03)
    cb.ax.tick_params(labelsize=5.2, length=1.5)
    cb.set_label("$\\it{Galr1}$ (z across animals)", fontsize=5.6)
    cb.outline.set_visible(False)
    # Pad clears the adult/aged group labels drawn above the columns.
    ax.set_title("* changes with age", loc="left", pad=17)

    # (d) the cell type that changes
    ax = fig.add_subplot(gs[3]); panel(ax, "d", dx=-0.70, dy=1.10)
    r = pop_stats.loc["Galr1"]
    paired_panel(ax, pop_cpm["Galr1"].to_dict(), "$\\it{Galr1}$ (CPM)", POPNAME)
    ax.title.set_color(POP)
    ax.title.set_fontweight("bold")
    ax.annotate(f"{2 ** r.lfc:.2f}×\n$P$ = {r.p:.3f}", xy=(.04, .98),
                xycoords="axes fraction", fontsize=6, color=INK, ha="left",
                va="top", linespacing=1.4)

    fig.savefig(OUT / "figure1_screen.pdf")
    fig.savefig(OUT / "figure1_screen.png")
    plt.close(fig)

    print(f"Galr1 screen: {len(screen)} cell types tested "
          f"(>= {MIN_CELLS_PER_ANIMAL} cells in every animal, "
          f">= {MIN_PCT_POS:.0f}% of cells positive)")
    with pd.option_context("display.width", 200):
        print(screen[["cell_type", "n_cells", "pct_pos", "cpm", "fold", "blocks",
                      "p", "margin", "significant"]].to_string(index=False))
    print(f"\nchanges with age: {list(screen[screen.significant].cell_type)}")
    print(f"\nexcluded before testing ({len(excluded)}):")
    for name, why in excluded:
        print(f"   {name:28s} {why}")
    print(f"\nwrote {OUT / 'figure1_screen.pdf'} (+ .png)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
