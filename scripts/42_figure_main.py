#!/usr/bin/env python3
"""Figure 1: which MBH cells upregulate Galr1 with age.

    Galr1 is expressed across the mediobasal hypothalamus.  In which cells does
    it change between adult and aged mice?

    One population: a glutamatergic Otp+/Cbln1+/Prdm8+ cell type of the dorsal
    hypothalamus, which gains Galr1 and Ghsr with age.

Panels:
  a  the registered subregions, every cell from all eight animals pooled in the
     shared anatomical frame, with this population on top
  b  Galr1 expression across the cell types the design can test
  c  the age effect, per block: a clustered heatmap of Galr1 aged-minus-adult in
     each of the four animal pairs, so consistency is visible rather than
     summarised away.  Cell types that change significantly are marked
  d  the same cells in expression space (UMAP), with the population that moved
  e  its HypoMap C185 correspondence
  f  Galr1 in that population, one point per animal, paired within block
  g  Ghsr, the second receptor these neurons carry, behaves the same way

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
HYPOMAP = REPO / "results" / "hypomap" / "hypomap_C185_named_all_celltypes.csv"
POPNAME = "Glut Prdm8/Cbln1"
ANIMALS = list(A_ADULT) + list(A_AGED)
KEY_GENES = ["Galr1", "Ghrh", "Ghsr", "Galr3"]
MAIN_GENES = ["Galr1", "Ghsr"]
HYPOMAP_SHOWN = 6
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
                 animals: np.ndarray) -> tuple[pd.DataFrame, list]:
    """Galr1 in every cell type the design can test, by one rule throughout.

    Returns the table and the types the entry criteria excluded, so the second
    is reported rather than left implicit.
    """
    gene = int(np.flatnonzero(var == "Galr1")[0])
    rows, excluded = [], []
    for name in pd.unique(cell_types):
        if str(name).startswith(("unlabelled", "unresolved")):
            continue
        sel = cell_types == name
        per_animal = [int((sel & (animals == a)).sum()) for a in ANIMALS]
        pct_pos = float((counts[sel][:, gene] > 0).mean() * 100)
        if min(per_animal) < MIN_CELLS_PER_ANIMAL:
            excluded.append((name, f"{min(per_animal)} cells in the thinnest animal"))
            continue
        if pct_pos < MIN_PCT_POS:
            excluded.append((name, f"{pct_pos:.1f}% of cells $Galr1$+"))
            continue
        mat = pd.DataFrame({a: counts[sel & (animals == a)].sum(axis=0)
                            for a in ANIMALS}, index=var).T
        cpm = np.log2(mat.div(mat.sum(axis=1), axis=0) * 1e6 + 1)
        r = blocked_stats(cpm).loc["Galr1"]
        row = {"cell_type": name, "n_cells": int(sel.sum()),
               "min_cells": min(per_animal), "pct_pos": round(pct_pos, 1),
               "cpm": round(float(np.log2(counts[sel][:, gene].sum()
                                          / counts[sel].sum() * 1e6 + 1)), 2),
               "lfc": round(float(r.lfc), 3), "fold": round(float(2 ** r.lfc), 2),
               "blocks": int(r.blocks), "p": round(float(r.p), 4),
               "margin": round(float(r.margin), 3)}
        for b in BLOCKS:
            row[f"lfc_{b}"] = round(float(r[f"lfc_{b}"]), 3)
        rows.append(row)
    out = pd.DataFrame(rows)
    out["significant"] = ((out.blocks == 4) & (out.margin > 0) & (out.p <= MAX_P))
    return out.sort_values("lfc", ascending=False), excluded


def population_markers(counts: np.ndarray, var: np.ndarray,
                       is_pop: np.ndarray) -> pd.Series:
    inside = counts[is_pop].sum(axis=0) / counts[is_pop].sum() * 1e6
    outside = counts[~is_pop].sum(axis=0) / counts[~is_pop].sum() * 1e6
    return pd.Series(np.log2((inside + 1) / (outside + 1)), index=var)


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
        ax.scatter([j] * len(ys), ys, s=13, color=colour, zorder=3, linewidths=0)
        ax.plot([j - .22, j + .22], [np.mean(ys)] * 2, color=colour, lw=1.3,
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

    mat = pd.DataFrame({a: counts[is_pop & (animals == a)].sum(axis=0)
                        for a in ANIMALS}, index=var).T
    expressed = mat.sum(axis=0) >= 200
    frac = mat.div(mat.sum(axis=1), axis=0) * 1e6
    cpm = np.log2(frac + 1)
    stats = blocked_stats(cpm)[expressed]
    stats.to_csv(SRC / "fig_main_stats.csv")
    cpm.loc[ANIMALS, [g for g in ("Fos",) if g in cpm.columns]].to_csv(
        SRC / "fig_main_fos_per_animal.csv")
    pct = {a: (is_pop & (animals == a)).sum() / (animals == a).sum() * 100
           for a in ANIMALS}
    pd.DataFrame({"abundance": pct}).loc[ANIMALS].to_csv(SRC / "fig_main_abundance.csv")

    screen, excluded = galr1_screen(counts, var, cell_types, animals)
    screen.to_csv(SRC / "galr1_by_celltype.csv", index=False)
    markers = population_markers(counts, var, is_pop)
    markers.sort_values(ascending=False).to_csv(SRC / "population_markers.csv")

    fig = plt.figure(figsize=(FULL, 5.4))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.0], hspace=.55)
    top = outer[0].subgridspec(1, 3, width_ratios=[1.45, 1.15, 1.30], wspace=.60)
    bot = outer[1].subgridspec(1, 4, width_ratios=[1.15, 1.30, .85, .85], wspace=1.02)

    # (a) the registered subregions, pooled, with the population on them
    ax = fig.add_subplot(top[0]); panel(ax, "a", dx=-0.04, dy=1.13)
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
              loc="upper left", bbox_to_anchor=(0, .03), frameon=False,
              fontsize=5.4, handletextpad=.3, borderpad=0, borderaxespad=0)
    ax.set_title("subregions, 8 animals pooled", loc="left", pad=2, x=.04)

    # (b) Galr1 expression across the testable cell types
    ax = fig.add_subplot(top[1]); panel(ax, "b", dx=-0.52, dy=1.13)
    expr = screen.sort_values("cpm")
    ys = np.arange(len(expr))
    ax.barh(ys, expr.cpm, height=.70,
            color=[POP if s else "#C4C4C4" for s in expr.significant])
    ax.set_yticks(ys); ax.set_yticklabels(expr.cell_type, fontsize=5.2)
    for tick, sig in zip(ax.get_yticklabels(), expr.significant):
        tick.set_fontweight("bold" if sig else "normal")
        tick.set_color(POP if sig else INK)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("$\\it{Galr1}$ (log$_2$ CPM)")
    ax.set_title("expression by cell type", loc="left", pad=2)

    # (c) the age effect per block, clustered
    ax = fig.add_subplot(top[2]); panel(ax, "c", dx=-0.46, dy=1.13)
    block_cols = [f"lfc_{b}" for b in BLOCKS]
    hm = screen.set_index("cell_type")[block_cols]
    link = optimal_leaf_ordering(linkage(hm.to_numpy(), "average"),
                                 pdist(hm.to_numpy()))
    hm = hm.iloc[leaves_list(link)]
    lim = float(np.abs(hm.to_numpy()).max())
    im = ax.imshow(hm.to_numpy(), aspect="auto", cmap=LFC_CMAP,
                   norm=TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim))
    ax.set_xticks(range(len(block_cols)))
    ax.set_xticklabels([b for b in BLOCKS], fontsize=5.6)
    ax.set_xlabel("animal pair (block)")
    sig = screen.set_index("cell_type").significant.reindex(hm.index)
    ax.set_yticks(range(len(hm)))
    ax.set_yticklabels([f"{n} *" if s else n for n, s in zip(hm.index, sig)],
                       fontsize=5.2)
    for tick, s in zip(ax.get_yticklabels(), sig):
        tick.set_fontweight("bold" if s else "normal")
        tick.set_color(POP if s else INK)
    ax.tick_params(axis="both", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=.055, pad=.04)
    cb.ax.tick_params(labelsize=5.2, length=1.5)
    cb.set_label("$\\it{Galr1}$ aged − adult (log$_2$)", fontsize=5.6)
    cb.outline.set_visible(False)
    ax.set_title("* all 4 blocks, $P$ ≤ 0.05", loc="left", pad=2)

    # (d) the same cells in expression space
    ax = fig.add_subplot(bot[0]); panel(ax, "d", dx=-0.10)
    ref = sc.read_h5ad(PROC / "mbh_roi_annotated.h5ad")
    umap = ref.obsm["X_umap"]
    ref_types = ref.obs["cell_type"].astype(str).to_numpy()
    pop_u = ref_types == POPNAME
    ax.scatter(umap[~pop_u, 0], umap[~pop_u, 1], s=.6, c=TISSUE, linewidths=0,
               rasterized=True)
    ax.scatter(umap[pop_u, 0], umap[pop_u, 1], s=2.4, c=POP, linewidths=0,
               rasterized=True)
    ax.annotate(POPNAME, (np.median(umap[pop_u, 0]), np.median(umap[pop_u, 1])),
                xytext=(0, 16), textcoords="offset points", fontsize=5.6,
                color=POP, fontweight="bold", ha="center",
                arrowprops=dict(arrowstyle="-", color=POP, lw=.5))
    ax.set_aspect("equal"); bare(ax)
    # Said plainly on the panel: this is the reference clustering, which is the
    # four hand-drawn ROI sections, not the eight the rest of the figure uses.
    ax.annotate(f"{ref.n_obs:,} cells, {ref.obs['section'].nunique()} ROI sections",
                xy=(.02, .02), xycoords="axes fraction",
                fontsize=5.4, color=INK, ha="left", va="bottom")
    ax.set_title("UMAP", loc="left", pad=2)
    del ref

    # (e) HypoMap correspondence
    ax = fig.add_subplot(bot[1]); panel(ax, "e", dx=-0.44)
    hmap = pd.read_csv(HYPOMAP, index_col=0)
    best = hmap[POPNAME].sort_values(ascending=False).head(HYPOMAP_SHOWN)[::-1]
    names = [i.split(": ", 1)[-1] for i in best.index]
    ax.barh(range(len(best)), best.values, height=.70,
            color=["#C4C4C4"] * (len(best) - 1) + [POP])
    ax.set_yticks(range(len(best))); ax.set_yticklabels(names, fontsize=5.4)
    ax.get_yticklabels()[-1].set_fontweight("bold")
    ax.get_yticklabels()[-1].set_color(POP)
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(.6, float(best.max()) * 1.06)
    ax.set_xlabel("Spearman ρ to HypoMap C185")
    ax.set_title("HypoMap identity", loc="left", pad=2)

    # (f, g) the receptors in that population
    for j, gene in enumerate(MAIN_GENES):
        ax = fig.add_subplot(bot[j + 2]); panel(ax, "fg"[j], dx=-0.46)
        r = stats.loc[gene]
        paired_panel(ax, frac[gene].to_dict(), f"$\\it{{{gene}}}$ (CPM)",
                     f"{2 ** r.lfc:.2f}×  $P$ = {r.p:.3f}")

    fig.savefig(OUT / "figure_main.pdf")
    fig.savefig(OUT / "figure_main.png")
    plt.close(fig)

    print(f"{int(is_pop.sum())} {POPNAME} neurons across 8 animals")
    for g in KEY_GENES:
        if g in stats.index:
            r = stats.loc[g]
            print(f"  {g:6s} {2 ** r.lfc:.2f}x  blocks {int(r.blocks)}/4  P {r.p:.4f}")

    print(f"\nGalr1 screen: {len(screen)} cell types tested "
          f"(>= {MIN_CELLS_PER_ANIMAL} cells in every animal, "
          f">= {MIN_PCT_POS:.0f}% of cells positive)")
    with pd.option_context("display.width", 200):
        print(screen[["cell_type", "n_cells", "pct_pos", "cpm", "fold", "blocks",
                      "p", "margin", "significant"]].to_string(index=False))
    print(f"\nsignificant: {list(screen[screen.significant].cell_type)}")
    print(f"\nexcluded before testing ({len(excluded)}):")
    for name, why in excluded:
        print(f"   {name:28s} {why}")
    print(f"\nwrote {OUT / 'figure_main.pdf'} (+ .png)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
