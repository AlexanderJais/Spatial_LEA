#!/usr/bin/env python3
"""Figure 1: which MBH cells change Galr1 between adult and aged mice.

    Galr1 is expressed across the mediobasal hypothalamus.  In which cells does
    it change with age?

    One: a glutamatergic Otp+/Cbln1+/Prdm8+ cell type of the dorsal
    hypothalamus, which also changes 21 other genes.

Panels:
  a  the delineated subregions, every cell from all eight animals pooled in the
     shared anatomical frame, with what was analysed
  b  where the cell type that changes sits, in the same frame
  c  Galr1 per animal, z-scored across the eight, adult block then aged block,
     for every cell type the design can test.  Ordered by the age test, so what
     changes is at the top
  d  Galr1 in that cell type, one point per animal, paired within block
  e  the same cells in expression space, carrying their HypoMap identity
  f  how well that identity holds, against the next five candidates
  g  every other gene that changes in this population, per animal

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
  * every panel is computed here, so the figure regenerates from the data.

Effect sizes and P values belong in the legend, not on the panels.  They are
printed at the end of this script and written to results/figures/source_data/.
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
HYPOMAP = REPO / "results" / "hypomap" / "hypomap_C185_named_all_celltypes.csv"
POPNAME = "Glut Prdm8/Cbln1"
ANIMALS = list(A_ADULT) + list(A_AGED)
# The manifest keys an animal "G_073"; the animal is G073.
ANIMAL_LABEL = {a: a.replace("_", "") for a in ANIMALS}
HYPOMAP_SHOWN = 6
MIN_CELLS_PER_ANIMAL = 30
MIN_PCT_POS = 10.0
MIN_TOTAL_COUNTS = 200
MAX_P = 0.05
LFC_CMAP = LinearSegmentedColormap.from_list("lfc", ["#2C6FA8", "#F4F4F4", GALR1])


def blocked_stats(cpm: pd.DataFrame) -> pd.DataFrame:
    """Within-block log2 change, block agreement, permutation P, and the gap
    between the two groups."""
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
    return pd.DataFrame({"lfc": mean, "blocks": agree, "p": p, "margin": margin})


def galr1_screen(counts, var, cell_types, animals):
    """Galr1 in every cell type the design can test, by one rule throughout."""
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
    out["significant"] = (out.blocks == 4) & (out.margin > 0) & (out.p <= MAX_P)
    # Ordered by the test the panel is about: what changes first, then by how
    # cleanly the groups separate, then by effect.  Not by clustering -- a
    # reader looking for the result should not have to hunt for it.
    out = out.sort_values(["significant", "margin", "lfc"],
                          ascending=[False, False, False])
    return out, pd.DataFrame(per_animal_cpm).T.loc[out.cell_type], excluded


def zscore(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.sub(frame.mean(axis=1), axis=0).div(
        frame.std(axis=1).replace(0, np.nan), axis=0)


def animal_heatmap(fig, ax, z: pd.DataFrame, cbar_label: str, ylabels,
                   bold=None) -> None:
    """Rows against the eight animals, adult block then aged block."""
    lim = float(np.nanmax(np.abs(z.to_numpy())))
    im = ax.imshow(z.to_numpy(), aspect="auto", cmap=LFC_CMAP,
                   norm=TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim))
    ax.axvline(len(A_ADULT) - .5, color="white", lw=1.4)
    ax.set_xticks(range(len(ANIMALS)))
    ax.set_xticklabels([ANIMAL_LABEL[a] for a in ANIMALS], fontsize=4.8,
                       rotation=90)
    for x, label in ((len(A_ADULT) / 2 - .5, "adult"),
                     (len(A_ADULT) + len(A_AGED) / 2 - .5, "aged")):
        ax.annotate(label, xy=(x, -.8), fontsize=5.8, ha="center", va="bottom",
                    color=INK, annotation_clip=False)
    ax.set_yticks(range(len(z)))
    ax.set_yticklabels(ylabels, fontsize=5.2)
    if bold is not None:
        for tick, flag in zip(ax.get_yticklabels(), bold):
            tick.set_fontweight("bold" if flag else "normal")
            tick.set_color(POP if flag else INK)
    ax.tick_params(axis="both", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=.05, pad=.03)
    cb.ax.tick_params(labelsize=5.2, length=1.5)
    cb.set_label(cbar_label, fontsize=5.6)
    cb.outline.set_visible(False)


def frame_panel(ax, ml, dv, nuc, present, coloured: bool):
    """The shared anatomical frame; subregions coloured, or all tissue grey."""
    ax.scatter(ml, dv, s=.45, c=TISSUE, linewidths=0, rasterized=True)
    if coloured:
        for key in sorted(present, key=lambda k: int((nuc == k).sum()),
                          reverse=True):
            m = nuc == key
            ax.scatter(ml[m], dv[m], s=.45, c=NUCLEUS_COLOUR[key], linewidths=0,
                       rasterized=True)
    ax.set_xlim(-1500, 1500); ax.set_ylim(-100, 1800)
    ax.set_aspect("equal"); bare(ax)
    scalebar(ax, 500, "500 µm")


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
    testable = mat.sum(axis=0) >= MIN_TOTAL_COUNTS
    pop_cpm = mat.div(mat.sum(axis=1), axis=0) * 1e6
    pop_stats = blocked_stats(np.log2(pop_cpm + 1))[testable]
    pop_stats["fold"] = 2 ** pop_stats.lfc
    pop_stats["changes"] = ((pop_stats.blocks == 4) & (pop_stats.margin > 0)
                            & (pop_stats.p <= MAX_P))
    pop_stats.sort_values("lfc", ascending=False).to_csv(
        SRC / "population_age_genes.csv")
    movers = pop_stats[pop_stats.changes].sort_values("lfc", ascending=False)

    fig = plt.figure(figsize=(FULL, 5.7))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.0], hspace=.45)
    top = outer[0].subgridspec(1, 3, width_ratios=[1.15, 1.15, 1.40], wspace=.78)
    bot = outer[1].subgridspec(1, 4, width_ratios=[.62, 1.05, 1.15, 1.25],
                               wspace=.80)
    present = [k for k in NUCLEUS_LABEL if (nuc == k).any()]

    # (a) what was analysed
    ax = fig.add_subplot(top[0]); panel(ax, "a", dx=-0.04, dy=1.12)
    frame_panel(ax, ml, dv, nuc, present, coloured=True)
    for key in present:
        m = nuc == key
        side = ml[m] > 0 if (ml[m] > 0).sum() > 30 else ml[m] < 0
        ax.annotate(NUCLEUS_LABEL[key],
                    (np.median(ml[m][side]), np.median(dv[m][side])),
                    fontsize=5.4, color=INK, ha="center", va="center",
                    fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.10", fc="white", ec="none",
                              alpha=.72))
    ax.annotate(f"{win.n_obs:,} cells", xy=(0, -0.02), xycoords="axes fraction",
                fontsize=5.6, color=INK, ha="left", va="top")

    # (b) where the cell type that changes sits
    ax = fig.add_subplot(top[1]); panel(ax, "b", dx=-0.04, dy=1.12)
    frame_panel(ax, ml, dv, nuc, present, coloured=False)
    ax.scatter(ml[is_pop], dv[is_pop], s=4.2, c=POP, linewidths=.2,
               edgecolors="white", rasterized=True)
    ax.annotate(POPNAME, xy=(0, -0.02), xycoords="axes fraction", fontsize=5.6,
                color=POP, ha="left", va="top", fontweight="bold")

    # (c) Galr1 per animal, every testable cell type
    ax = fig.add_subplot(top[2]); panel(ax, "c", dx=-0.60, dy=1.12)
    z = zscore(per_animal)
    animal_heatmap(fig, ax, z, "$\\it{Galr1}$ (z across animals)",
                   list(z.index), bold=list(screen.significant))

    # (d) Galr1 in that cell type
    ax = fig.add_subplot(bot[0]); panel(ax, "d", dx=-0.78, dy=1.12)
    values = pop_cpm["Galr1"].to_dict()
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
    ax.set_ylim(0, max(values.values()) * 1.12)
    ax.set_ylabel("$\\it{Galr1}$ (CPM)")
    ax.set_title(POPNAME, loc="left", pad=3, color=POP, fontweight="bold")

    # (e) the same cells in expression space, with the HypoMap identity
    ax = fig.add_subplot(bot[1]); panel(ax, "e", dx=-0.10, dy=1.12)
    ref = sc.read_h5ad(PROC / "mbh_roi_annotated.h5ad")
    umap = ref.obsm["X_umap"]
    pop_u = ref.obs["cell_type"].astype(str).to_numpy() == POPNAME
    ax.scatter(umap[~pop_u, 0], umap[~pop_u, 1], s=.6, c=TISSUE, linewidths=0,
               rasterized=True)
    ax.scatter(umap[pop_u, 0], umap[pop_u, 1], s=2.6, c=POP, linewidths=0,
               rasterized=True)
    hmap = pd.read_csv(HYPOMAP, index_col=0)
    best = hmap[POPNAME].sort_values(ascending=False)
    ax.annotate(POPNAME, xy=(np.median(umap[pop_u, 0]), np.median(umap[pop_u, 1])),
                xytext=(0, 14), textcoords="offset points", fontsize=5.6,
                color=POP, fontweight="bold", ha="center", va="bottom")
    ax.set_aspect("equal"); bare(ax)
    ax.set_title("Xenium UMAP", loc="left", pad=2, x=.06)
    del ref

    # (f) how well that identity holds
    ax = fig.add_subplot(bot[2]); panel(ax, "f", dx=-0.52, dy=1.12)
    shown = best.head(HYPOMAP_SHOWN)[::-1]
    ax.barh(range(len(shown)), shown.values, height=.70,
            color=["#C4C4C4"] * (len(shown) - 1) + [POP])
    ax.set_yticks(range(len(shown)))
    ax.set_yticklabels([i.split(": ", 1)[-1] for i in shown.index], fontsize=5.4)
    ax.get_yticklabels()[-1].set_fontweight("bold")
    ax.get_yticklabels()[-1].set_color(POP)
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(.6, float(shown.max()) * 1.05)
    ax.set_xlabel("Spearman ρ")

    # (g) what else changes in this population
    ax = fig.add_subplot(bot[3]); panel(ax, "g", dx=-0.46, dy=1.12)
    lv = np.log2(pop_cpm[movers.index] + 1).T
    zg = zscore(lv)
    if len(zg) > 2:
        link = optimal_leaf_ordering(linkage(zg.to_numpy(), "average"),
                                     pdist(zg.to_numpy()))
        zg = zg.iloc[leaves_list(link)]
    animal_heatmap(fig, ax, zg, "z across animals",
                   [f"$\\it{{{g}}}$" for g in zg.index])

    fig.savefig(OUT / "figure1.pdf")
    fig.savefig(OUT / "figure1.png")
    plt.close(fig)

    print(f"{win.n_obs:,} cells, {len(ANIMALS)} mice, 1 section each")
    print(f"{int(is_pop.sum()):,} {POPNAME} cells\n")
    print(f"Galr1 screen: {len(screen)} cell types tested "
          f"(>= {MIN_CELLS_PER_ANIMAL} cells in every animal, "
          f">= {MIN_PCT_POS:.0f}% positive)")
    with pd.option_context("display.width", 200):
        print(screen[["cell_type", "n_cells", "pct_pos", "cpm", "fold", "blocks",
                      "p", "margin", "significant"]].to_string(index=False))
    print(f"\nexcluded before testing ({len(excluded)}):")
    for name, why in excluded:
        print(f"   {name:28s} {why}")
    print(f"\n=== for the legend: {POPNAME} ===")
    print(f"HypoMap C185 {best.index[0].split(': ', 1)[-1]}, "
          f"Spearman ρ = {best.iloc[0]:.3f} "
          f"(next {best.index[1].split(': ', 1)[-1]}, {best.iloc[1]:.3f})")
    print(f"{len(movers)} of {int(testable.sum())} testable genes change:")
    print(movers[["fold", "lfc", "blocks", "p", "margin"]].round(3).to_string())
    print(f"\nwrote {OUT / 'figure1.pdf'} (+ .png)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
