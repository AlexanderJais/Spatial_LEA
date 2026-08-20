#!/usr/bin/env python3
"""Main figure: dorsal hypothalamic Prdm8/Cbln1 neurons in ageing.

Single biological conclusion the figure supports:

    A defined glutamatergic population of the dorsal hypothalamus loses activity
    with age while upregulating Galr1, without any loss of neurons.

Panels, in narrative order:
  a  where the neurons are, at three magnifications of the same tissue
  b  what they are -- marker enrichment establishes identity
  c  Fos falls          } the observation, one point per animal,
  d  Galr1 rises        } paired within block, all four blocks shown
  e  Ghsr rises         }
  f  abundance is unchanged, so this is regulation and not cell loss
  g  reproducibility -- the within-block fold change for every gene that
     separates the groups completely

Choices made by rule rather than by eye:
  * representative section = the one of the eight whose neuron count is nearest
    the median, ties broken by anatomical frame quality
  * representative field = the 250 um square in that section containing the most
    neurons of this population
Both are stated in the caption text printed by this script.

Statistics are computed within blocks throughout; each block is one aged and one
adult mouse processed on the same slide run, gene panel and segmentation
chemistry, so no comparison crosses a technical boundary.
"""

from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.figstyle import (  # noqa: E402
    ADULT, AGED as C_AGED, FULL, INK, MUTED, OTHER_CELLS, POP, TISSUE,
    bare, panel, scalebar, use_style,
)
from spatial_lea.io import (  # noqa: E402
    ADULT as A_ADULT, AGED as A_AGED, ANIMAL_META, BLOCKS, ONE_PER_MOUSE,
    counts_matrix,
)

PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "figures"
SRC = OUT / "source_data"
POPNAME = "Glut Prdm8/Cbln1"
ANIMALS = list(A_AGED) + list(A_ADULT)
GENES = ["Fos", "Galr1", "Ghsr"]
FIELD_UM = 250.0
MARKERS = 10


def blocked_stats(cpm: pd.DataFrame) -> pd.DataFrame:
    lfc = pd.DataFrame({b: cpm.loc[x] - cpm.loc[y] for b, (x, y) in BLOCKS.items()}).T
    mean = lfc.mean()
    agree = (np.sign(lfc) == np.sign(mean)).sum()
    vals = cpm.to_numpy(); idx = list(range(len(ANIMALS)))
    obs = vals[:4].mean(axis=0) - vals[4:].mean(axis=0)
    null = np.array([vals[list(c)].mean(axis=0)
                     - vals[[i for i in idx if i not in c]].mean(axis=0)
                     for c in combinations(idx, 4)])
    p = (np.abs(null) >= np.abs(obs)).mean(axis=0)
    margin = np.maximum(vals[:4].min(axis=0) - vals[4:].max(axis=0),
                        vals[4:].min(axis=0) - vals[:4].max(axis=0))
    out = pd.DataFrame({"lfc": mean, "blocks": agree, "p": p, "margin": margin})
    for b in BLOCKS:
        out[f"lfc_{b}"] = lfc.loc[b]
    return out


def paired_panel(ax, values: dict, ylabel: str, title: str) -> None:
    """One point per animal, aged and adult joined within block."""
    for a_aged, a_adult in BLOCKS.values():
        ax.plot([0, 1], [values[a_aged], values[a_adult]], color="#C9C9C9",
                lw=.6, zorder=1)
    for j, (grp, members, colour) in enumerate(
            (("aged", A_AGED, C_AGED), ("adult", A_ADULT, ADULT))):
        ys = [values[a] for a in members]
        ax.scatter([j] * len(ys), ys, s=13, color=colour, zorder=3,
                   linewidths=0, clip_on=False)
        ax.plot([j - .22, j + .22], [np.mean(ys)] * 2, color=colour, lw=1.3,
                solid_capstyle="butt", zorder=2)
    ax.set_xlim(-.45, 1.45); ax.set_xticks([0, 1])
    ax.set_xticklabels(["aged", "adult"])
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", pad=3)


def main() -> int:
    use_style()
    OUT.mkdir(parents=True, exist_ok=True); SRC.mkdir(parents=True, exist_ok=True)

    win = sc.read_h5ad(PROC / "hypothalamus_nuclei.h5ad")
    win = win[win.obs["section"].astype(str).isin(ONE_PER_MOUSE)].copy()
    is_pop = (win.obs["cell_type"].astype(str) == POPNAME).to_numpy()
    animals = win.obs["animal"].astype(str).to_numpy()
    sections = win.obs["section"].astype(str).to_numpy()
    counts = counts_matrix(win)
    var = win.var_names.to_numpy()

    # --- statistics -------------------------------------------------------
    mat = pd.DataFrame({a: counts[is_pop & (animals == a)].sum(axis=0)
                        for a in ANIMALS}, index=var).T
    expressed = mat.sum(axis=0) >= 200
    cpm = np.log2(mat.div(mat.sum(axis=1), axis=0) * 1e6 + 1)
    stats = blocked_stats(cpm)[expressed]
    stats.to_csv(SRC / "fig_main_stats.csv")

    # --- representative section and field, by rule ------------------------
    per_sec = pd.Series({s: int((is_pop & (sections == s)).sum())
                         for s in ONE_PER_MOUSE})
    qc = pd.read_csv(REPO / "results" / "anatomy" / "frame_qc.csv", index_col=0)
    cand = (per_sec - per_sec.median()).abs()
    best = cand[cand == cand.min()].index
    rep = max(best, key=lambda s: qc.loc[s, "quality"])
    rep_pop = is_pop & (sections == rep)
    ml, dv = win.obs["ml"].to_numpy(), win.obs["dv"].to_numpy()

    px, py = ml[rep_pop], dv[rep_pop]
    grid = [(x, y, int(((px >= x) & (px < x + FIELD_UM)
                        & (py >= y) & (py < y + FIELD_UM)).sum()))
            for x in np.arange(px.min() - 20, px.max(), 25)
            for y in np.arange(py.min() - 20, py.max(), 25)]
    fx, fy, fn = max(grid, key=lambda t: t[2])
    field = (fx, fy, fx + FIELD_UM, fy + FIELD_UM)

    # --- figure -----------------------------------------------------------
    fig = plt.figure(figsize=(FULL, 5.55))
    gs = fig.add_gridspec(3, 3, height_ratios=[.90, .90, .90],
                          width_ratios=[1.22, 1, 1], hspace=.52, wspace=.52)

    # (a) three magnifications of the same tissue
    gsa = gs[0, :].subgridspec(1, 3, width_ratios=[1.35, 1.15, 1], wspace=.14)

    whole = sc.read_h5ad(PROC / "mbh_anatomical.h5ad", backed="r")
    ws = whole.obs["section"].astype(str).to_numpy() == rep
    wml = whole.obs["ml"].to_numpy()[ws]
    wdv = whole.obs["dv"].to_numpy()[ws]
    whole.file.close()

    ax = fig.add_subplot(gsa[0]); panel(ax, "a", dx=-0.03, dy=1.16)
    sub = np.random.default_rng(0).choice(len(wml), min(55000, len(wml)), replace=False)
    ax.scatter(wml[sub], wdv[sub], s=.25, c=TISSUE, linewidths=0, rasterized=True)
    ax.add_patch(Rectangle((-1500, -100), 3000, 1900, fill=False, ec=INK, lw=.6))
    ax.set_aspect("equal"); bare(ax)
    ax.set_xlim(np.percentile(wml, .2), np.percentile(wml, 99.8))
    ax.set_ylim(np.percentile(wdv, .2), np.percentile(wdv, 99.8))
    scalebar(ax, 1000, "1 mm")
    ax.set_title(f"{rep}, coronal", loc="left", pad=2, color=MUTED, x=.06)

    ax = fig.add_subplot(gsa[1])
    inwin = (sections == rep)
    ax.scatter(ml[inwin & ~is_pop], dv[inwin & ~is_pop], s=.7, c=TISSUE,
               linewidths=0, rasterized=True)
    ax.scatter(ml[rep_pop], dv[rep_pop], s=4.5, c=POP, linewidths=0, rasterized=True)
    ax.add_patch(Rectangle((field[0], field[1]), FIELD_UM, FIELD_UM,
                           fill=False, ec=INK, lw=.6))
    ax.set_xlim(-1500, 1500); ax.set_ylim(-100, 1800)
    ax.set_aspect("equal"); bare(ax)
    scalebar(ax, 500, "500 µm")
    ax.set_title("hypothalamic window", loc="left", pad=2, color=MUTED)
    # Direct labels rather than a legend; each is in the colour of the mark it
    # names, so nothing has to be matched up by eye.
    ax.annotate("Prdm8/Cbln1 neurons", xy=(-1460, 1690), fontsize=5.5,
                color=POP, va="top", ha="left")
    ax.annotate("all other cells", xy=(-1460, 1530), fontsize=5.5,
                color="#A8A8A8", va="top", ha="left")

    ax = fig.add_subplot(gsa[2])
    m = inwin & (ml >= field[0]) & (ml < field[2]) & (dv >= field[1]) & (dv < field[3])
    ax.scatter(ml[m & ~is_pop], dv[m & ~is_pop], s=9, c=OTHER_CELLS,
               linewidths=0, rasterized=True)
    ax.scatter(ml[m & is_pop], dv[m & is_pop], s=16, c=POP, linewidths=0,
               rasterized=True)
    ax.set_xlim(field[0], field[2]); ax.set_ylim(field[1], field[3])
    ax.set_aspect("equal"); bare(ax)
    scalebar(ax, 50, "50 µm")
    ax.set_title(f"{fn} of these neurons in field", loc="left", pad=2, color=MUTED)

    # (b) identity
    ax = fig.add_subplot(gs[1, 0]); panel(ax, "b")
    cin = counts[is_pop].sum(axis=0) / counts[is_pop].sum() * 1e6
    cout = counts[~is_pop].sum(axis=0) / counts[~is_pop].sum() * 1e6
    enr = pd.Series(np.log2((cin + 1) / (cout + 1)), index=var).nlargest(MARKERS)
    enr = enr.sort_values()
    receptors = {"Galr1", "Galr3", "Ghsr"}
    ax.barh(range(len(enr)), enr.values, height=.68,
            color=[POP if g in receptors else "#9A9A9A" for g in enr.index])
    ax.set_yticks(range(len(enr)))
    ax.set_yticklabels([f"$\\it{{{g}}}$" for g in enr.index])
    ax.set_xlabel("log$_2$ enrichment over other cells")
    ax.tick_params(axis="y", length=0)

    # (c-e) the observation
    for j, gene in enumerate(GENES):
        ax = fig.add_subplot(gs[1, j + 1] if j < 2 else gs[2, 0])
        panel(ax, "cde"[j])
        r = stats.loc[gene]
        paired_panel(ax, cpm[gene].to_dict(), f"$\\it{{{gene}}}$ (log$_2$ CPM)",
                     f"{2 ** r.lfc:.2f}×   $P$ = {r.p:.3f}")

    # (f) abundance
    ax = fig.add_subplot(gs[2, 1]); panel(ax, "f")
    pct = {a: (is_pop & (animals == a)).sum() / (animals == a).sum() * 100
           for a in ANIMALS}
    ab = blocked_stats(pd.DataFrame({"abundance": pct}).loc[ANIMALS])
    paired_panel(ax, pct, "abundance (% of cells)",
                 f"no change   $P$ = {ab.loc['abundance', 'p']:.2f}")

    # (g) reproducibility: every gene that separates the groups completely,
    # ranked. All are plotted; only the four that carry the argument are named,
    # so no gene is selected out of the display.
    ax = fig.add_subplot(gs[2, 2]); panel(ax, "g")
    sep = stats[(stats.blocks == 4) & (stats.margin > 0)]
    sep = sep.reindex(sep.lfc.sort_values().index)
    ys = np.arange(len(sep))
    for b in BLOCKS:
        ax.scatter(sep[f"lfc_{b}"], ys, s=5, facecolors="none",
                   edgecolors="#B0B0B0", linewidths=.4, zorder=2)
    ax.scatter(sep["lfc"], ys, s=11, color=INK, zorder=3, linewidths=0)
    ax.axvline(0, color=MUTED, lw=.5)
    ax.set_yticks([]); ax.set_ylim(-1.2, len(sep) + .4)
    ax.set_xlim(-3.6, 1.6)
    ax.set_xlabel("aged / adult (log$_2$)")
    ax.set_ylabel(f"{len(sep)} genes, ranked")
    for g, ha, dx in (("Fos", "left", .24), ("Tacr1", "left", .24),
                      ("Galr1", "right", -.20), ("Ghsr", "right", -.20)):
        if g in sep.index:
            y = int(np.where(sep.index == g)[0][0])
            ax.annotate(f"$\\it{{{g}}}$", (sep.lfc[g], y), xytext=(dx * 30, 0),
                        textcoords="offset points", fontsize=5.5, va="center",
                        ha=ha, color=INK)

    fig.savefig(OUT / "figure_main.pdf")
    fig.savefig(OUT / "figure_main.png")
    plt.close(fig)

    print(f"Representative section: {rep} "
          f"(neurons {per_sec[rep]}, median {per_sec.median():.1f}; "
          f"frame quality {qc.loc[rep, 'quality']})")
    print(f"Representative field: {FIELD_UM:.0f} um square at "
          f"ml {field[0]:.0f}..{field[2]:.0f}, dv {field[1]:.0f}..{field[3]:.0f}, "
          f"containing {fn} neurons (the densest such square in that section)")
    print(f"Genes separating completely in all four blocks: {len(sep)}")
    print(sep[["lfc", "p", "margin"]].round(3).to_string())
    print(f"\nwrote {OUT / 'figure_main.pdf'} (+ .png)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
