#!/usr/bin/env python3
"""Figure 2: what the cell type that changes is, and what else changes in it.

Figure 1 asks which MBH cells change Galr1 between adult and aged mice and finds
one.  This figure is about that population:

  a  where it sits in expression space -- the reference clustering's UMAP
  b  what it corresponds to in HypoMap, at the C185 level
  c  what else changes in it with age: every expressed gene, effect against how
     cleanly the two groups separate.  Genes that meet the same criterion Galr1
     met are named
  d  those genes per animal, z-scored, adult block then aged block

Same rule as Figure 1 for what counts as a change: all four blocks agreeing in
sign, aged and adult ranges not overlapping, permutation P <= 0.05.  0.0286 is
the smallest P eight animals can give.

Genes are tested only if they carry enough signal to be tested at all
(MIN_TOTAL_COUNTS across the population), which is stated here and printed at
the end with what it excluded.
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
    FULL, GALR1, INK, POP, TISSUE, bare, panel, use_style,
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
HYPOMAP_SHOWN = 6
MIN_TOTAL_COUNTS = 200
LABELLED = 7                     # named in panel c; panel d names them all
MAX_P = 0.05
LFC_CMAP = LinearSegmentedColormap.from_list("lfc", ["#2C6FA8", "#F4F4F4", GALR1])


def blocked_stats(cpm: pd.DataFrame) -> pd.DataFrame:
    """Within-block log2 change, block agreement, permutation P, and the gap
    between the groups.  Identical to Figure 1's test."""
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
    counts = counts_matrix(win)
    var = win.var_names.to_numpy()

    mat = pd.DataFrame({a: counts[is_pop & (animals == a)].sum(axis=0)
                        for a in ANIMALS}, index=var).T
    testable = mat.sum(axis=0) >= MIN_TOTAL_COUNTS
    cpm = mat.div(mat.sum(axis=1), axis=0) * 1e6
    stats = blocked_stats(np.log2(cpm + 1))[testable]
    stats["fold"] = 2 ** stats.lfc
    stats["changes"] = (stats.blocks == 4) & (stats.margin > 0) & (stats.p <= MAX_P)
    stats.sort_values("lfc", ascending=False).to_csv(SRC / "population_age_genes.csv")
    movers = stats[stats.changes].sort_values("lfc", ascending=False)

    fig = plt.figure(figsize=(FULL, 3.05))
    gs = fig.add_gridspec(1, 4, width_ratios=[1.0, 1.35, 1.25, 1.05], wspace=.80)

    # (a) the population in expression space
    ax = fig.add_subplot(gs[0]); panel(ax, "a", dx=-0.06, dy=1.10)
    ref = sc.read_h5ad(PROC / "mbh_roi_annotated.h5ad")
    umap = ref.obsm["X_umap"]
    pop_u = ref.obs["cell_type"].astype(str).to_numpy() == POPNAME
    ax.scatter(umap[~pop_u, 0], umap[~pop_u, 1], s=.6, c=TISSUE, linewidths=0,
               rasterized=True)
    ax.scatter(umap[pop_u, 0], umap[pop_u, 1], s=2.6, c=POP, linewidths=0,
               rasterized=True)
    ax.annotate(POPNAME, (np.median(umap[pop_u, 0]), np.median(umap[pop_u, 1])),
                xytext=(0, 20), textcoords="offset points", fontsize=5.6,
                color=POP, fontweight="bold", ha="center",
                arrowprops=dict(arrowstyle="-", color=POP, lw=.5))
    ax.annotate(f"{ref.n_obs:,} cells, {ref.obs['section'].nunique()} ROI sections",
                xy=(.02, .02), xycoords="axes fraction", fontsize=5.2, color=INK,
                ha="left", va="bottom")
    ax.set_aspect("equal"); bare(ax)
    ax.set_title("UMAP", loc="left", pad=2, x=.10)
    del ref

    # (b) what it corresponds to in HypoMap
    ax = fig.add_subplot(gs[1]); panel(ax, "b", dx=-0.52, dy=1.10)
    hmap = pd.read_csv(HYPOMAP, index_col=0)
    best = hmap[POPNAME].sort_values(ascending=False).head(HYPOMAP_SHOWN)[::-1]
    ax.barh(range(len(best)), best.values, height=.70,
            color=["#C4C4C4"] * (len(best) - 1) + [POP])
    ax.set_yticks(range(len(best)))
    ax.set_yticklabels([i.split(": ", 1)[-1] for i in best.index], fontsize=5.4)
    ax.get_yticklabels()[-1].set_fontweight("bold")
    ax.get_yticklabels()[-1].set_color(POP)
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(.6, float(best.max()) * 1.05)
    ax.set_xlabel("Spearman ρ to HypoMap C185")
    ax.set_title("HypoMap identity", loc="left", pad=2)

    # (c) what else changes with age in this population
    ax = fig.add_subplot(gs[2]); panel(ax, "c", dx=-0.42, dy=1.10)
    rest = stats[~stats.changes]
    ax.axhline(0, color=INK, lw=.5)
    ax.axvline(0, color=INK, lw=.5)
    ax.scatter(rest.lfc, rest.margin, s=7, c="#C4C4C4", linewidths=0)
    ax.scatter(movers.lfc, movers.margin, s=20, c=POP, linewidths=0, zorder=3)
    # The movers sit almost on top of each other just above zero, so labelling
    # all of them in place is unreadable.  A few are named on a ladder with
    # leader lines; panel d names every one of them.
    # Galr1 is named whatever its rank: it is the gene the study is about, and
    # its position among the others is the point of the panel.
    pick = list(movers.lfc.abs().sort_values(ascending=False).index[:LABELLED])
    if "Galr1" in movers.index and "Galr1" not in pick:
        pick = pick[:-1] + ["Galr1"]
    named = movers.reindex(pick).sort_values("lfc")
    top = float(movers.margin.max())
    ladder = np.linspace(top + .25, top + .25 + .16 * (len(named) - 1), len(named))
    for (gene, row), y in zip(named.iterrows(), ladder):
        ax.annotate(f"$\\it{{{gene}}}$", (row.lfc, row.margin), xytext=(row.lfc, y),
                    textcoords="data", fontsize=5.4, color=POP, fontweight="bold",
                    ha="center", va="bottom",
                    arrowprops=dict(arrowstyle="-", color=POP, lw=.35,
                                    shrinkA=1, shrinkB=2))
    ax.set_ylim(float(movers.margin.min()) - .1 if len(rest) == 0
                else float(stats.margin.min()) - .15, ladder[-1] + .3)
    ax.set_xlabel("aged / adult (log$_2$)")
    ax.set_ylabel("gap between the groups (log$_2$)")
    ax.set_title(f"{len(movers)} of {len(stats)} genes change", loc="left", pad=2)

    # (d) those genes, per animal
    ax = fig.add_subplot(gs[3]); panel(ax, "d", dx=-0.46, dy=1.10)
    lv = np.log2(cpm[movers.index] + 1).T
    z = lv.sub(lv.mean(axis=1), axis=0).div(lv.std(axis=1).replace(0, np.nan), axis=0)
    if len(z) > 2:
        link = optimal_leaf_ordering(linkage(z.to_numpy(), "average"),
                                     pdist(z.to_numpy()))
        z = z.iloc[leaves_list(link)]
    lim = float(np.nanmax(np.abs(z.to_numpy())))
    im = ax.imshow(z.to_numpy(), aspect="auto", cmap=LFC_CMAP,
                   norm=TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim))
    ax.axvline(len(A_ADULT) - .5, color="white", lw=1.4)
    ax.set_xticks(range(len(ANIMALS)))
    ax.set_xticklabels(ANIMALS, fontsize=4.8, rotation=90)
    for x, label in ((len(A_ADULT) / 2 - .5, "adult"),
                     (len(A_ADULT) + len(A_AGED) / 2 - .5, "aged")):
        ax.annotate(label, xy=(x, -.85), fontsize=5.6, ha="center", va="bottom",
                    color=INK, annotation_clip=False)
    ax.set_yticks(range(len(z)))
    ax.set_yticklabels([f"$\\it{{{g}}}$" for g in z.index], fontsize=5.4)
    ax.tick_params(axis="both", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=.055, pad=.04)
    cb.ax.tick_params(labelsize=5.2, length=1.5)
    cb.set_label("z across animals", fontsize=5.6)
    cb.outline.set_visible(False)
    ax.set_title("per animal", loc="left", pad=17)

    fig.savefig(OUT / "figure2_population.pdf")
    fig.savefig(OUT / "figure2_population.png")
    plt.close(fig)

    print(f"{int(is_pop.sum())} {POPNAME} cells across 8 animals")
    print(f"{int(testable.sum())} of {len(var)} genes carry >= {MIN_TOTAL_COUNTS} "
          "counts in this population and were tested\n")
    print(f"=== changes with age ({len(movers)}) ===")
    print(movers[["fold", "lfc", "blocks", "p", "margin"]].round(3).to_string())
    print("\n=== nearest misses (4 blocks, groups overlap) ===")
    near = stats[(stats.blocks == 4) & (~stats.changes)].sort_values(
        "margin", ascending=False).head(8)
    print(near[["fold", "lfc", "blocks", "p", "margin"]].round(3).to_string())
    print(f"\nwrote {OUT / 'figure2_population.pdf'} (+ .png)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
