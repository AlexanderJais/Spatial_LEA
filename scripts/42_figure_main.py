#!/usr/bin/env python3
"""Main figure: Galr1 upregulation in dorsal hypothalamic Ghrh neurons with age.

Single biological conclusion the figure supports:

    A defined glutamatergic population of the dorsal hypothalamus -- Galr1+,
    Ghrh+, Ghsr+, Galr3+ -- upregulates Galr1 with age, without loss of neurons.

Panels:
  a  where the neurons are: one intact section, then every neuron of this type
     from all eight animals in the shared anatomical frame
  c  Galr1 tested in every cell type with enough cells: of twenty, one moves
     consistently and separates the animals, which is what makes the result a
     statement about this population rather than about ageing hypothalamus
  d  Galr1 in that population, one point per animal, paired within block
  e  Ghsr, the second receptor these neurons carry, behaves the same way

Choices made by rule rather than by eye:
  * representative section = the section with the highest whole-section bilateral
    balance, i.e. the most intact tissue.  Balance is the ratio of cells left
    and right of the midline; a torn or folded section scores low.  Restricting
    that test to the analysed window is not enough -- a section can be perfectly
    symmetric inside the window and badly damaged outside it.

Adult is plotted first throughout, as the reference condition.  Fos falls
steeply in these neurons in aged mice, but these animals are untreated, so
immediate-early gene expression is hard to interpret without a stimulus; it is
reported in Extended Data rather than as a main-figure claim.
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

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.figstyle import (  # noqa: E402
    ADULT, AGED as C_AGED, FULL, INK, NUCLEUS_COLOUR, NUCLEUS_LABEL, POP,
    TISSUE, bare, panel, scalebar, use_style,
)
from spatial_lea.io import (  # noqa: E402
    ADULT as A_ADULT, AGED as A_AGED, BLOCKS, ONE_PER_MOUSE, counts_matrix,
)

PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "figures"
SRC = OUT / "source_data"
POPNAME = "Glut Prdm8/Cbln1"
# Adult first: the reference condition.
ANIMALS = list(A_ADULT) + list(A_AGED)
# Priority order for this study; Galr1 leads, then the peptide that names the
# population, then the other two receptors it carries.
KEY_GENES = ["Galr1", "Ghrh", "Ghsr", "Galr3"]
MAIN_GENES = ["Galr1", "Ghsr"]
CONTEXT_MARKERS = 8


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


def paired_panel(ax, values: dict, ylabel: str, title: str) -> None:
    """One point per animal; adult on the left, joined to its block partner."""
    for a_aged, a_adult in BLOCKS.values():
        ax.plot([0, 1], [values[a_adult], values[a_aged]], color="#C9C9C9",
                lw=.6, zorder=1)
    for j, (members, colour) in enumerate(((A_ADULT, ADULT), (A_AGED, C_AGED))):
        ys = [values[a] for a in members]
        ax.scatter([j] * len(ys), ys, s=13, color=colour, zorder=3,
                   linewidths=0, clip_on=False)
        ax.plot([j - .22, j + .22], [np.mean(ys)] * 2, color=colour, lw=1.3,
                solid_capstyle="butt", zorder=2)
    ax.set_xlim(-.45, 1.45); ax.set_xticks([0, 1])
    ax.set_xticklabels(["adult", "aged"])
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
    ml, dv = win.obs["ml"].to_numpy(), win.obs["dv"].to_numpy()
    counts = counts_matrix(win)
    var = win.var_names.to_numpy()

    mat = pd.DataFrame({a: counts[is_pop & (animals == a)].sum(axis=0)
                        for a in ANIMALS}, index=var).T
    expressed = mat.sum(axis=0) >= 200
    cpm = np.log2(mat.div(mat.sum(axis=1), axis=0) * 1e6 + 1)
    stats = blocked_stats(cpm)[expressed]
    stats.to_csv(SRC / "fig_main_stats.csv")
    # Fos is reported in Extended Data, so its per-animal values travel too.
    cpm.loc[ANIMALS, [g for g in ("Fos",) if g in cpm.columns]].to_csv(
        SRC / "fig_main_fos_per_animal.csv")

    # Abundance is reported in Extended Data; its per-animal values travel too.
    pct = {a: (is_pop & (animals == a)).sum() / (animals == a).sum() * 100
           for a in ANIMALS}
    pd.DataFrame({"abundance": pct}).loc[ANIMALS].to_csv(SRC / "fig_main_abundance.csv")

    # --- representative section: the most intact tissue, whole section --------
    whole = sc.read_h5ad(PROC / "mbh_anatomical.h5ad", backed="r")
    wsec = whole.obs["section"].astype(str).to_numpy()
    wml = whole.obs["ml"].to_numpy()
    balance = {}
    for s in ONE_PER_MOUSE:
        x = wml[wsec == s]
        left, right = (x < -200).sum(), (x > 200).sum()
        balance[s] = min(left, right) / max(left, right)
    rep = max(balance, key=balance.get)
    whole.file.close()

    fig = plt.figure(figsize=(FULL, 4.6))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.12, .88], hspace=.52)
    top = outer[0].subgridspec(1, 3, width_ratios=[1.30, 1.05, .95], wspace=.36)
    bot = outer[1].subgridspec(1, 3, width_ratios=[1, 1, 1.32], wspace=.60)

    # (a) the registered subregions this panel resolves
    ax = fig.add_subplot(top[0]); panel(ax, "a", dx=-0.05, dy=1.13)
    nuc = win.obs["nucleus_ext"].astype(str).to_numpy()
    for key in ("edge", "fibre", "TUseg"):
        m = nuc == key
        if m.any():
            ax.scatter(ml[m], dv[m], s=.45, c=TISSUE, linewidths=0, rasterized=True)
    for key in NUCLEUS_LABEL:
        m = nuc == key
        if not m.any():
            continue
        ax.scatter(ml[m], dv[m], s=.45, c=NUCLEUS_COLOUR[key], linewidths=0,
                   rasterized=True)
        # Direct label on the right-hand side of the bilateral structure, so the
        # map is readable without a colour key.
        side = ml[m] > 0 if (ml[m] > 0).sum() > 30 else ml[m] < 0
        ax.annotate(NUCLEUS_LABEL[key],
                    (np.median(ml[m][side]), np.median(dv[m][side])),
                    fontsize=5.6, color=INK, ha="center", va="center",
                    fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none",
                              alpha=.72))
    ax.set_xlim(-1500, 1500); ax.set_ylim(-100, 1800)
    ax.set_aspect("equal"); bare(ax)
    scalebar(ax, 500, "500 µm")
    ax.set_title("registered subregions, 8 animals pooled", loc="left", pad=2)

    # (b) representative image
    ax = fig.add_subplot(top[1]); panel(ax, "b", dx=-0.09, dy=1.13)
    m = sections == rep
    ax.scatter(ml[m & ~is_pop], dv[m & ~is_pop], s=.8, c=TISSUE, linewidths=0,
               rasterized=True)
    ax.scatter(ml[m & is_pop], dv[m & is_pop], s=5.5, c=POP, linewidths=0,
               rasterized=True)
    ax.set_xlim(-1500, 1500); ax.set_ylim(-100, 1800)
    ax.set_aspect("equal"); bare(ax)
    scalebar(ax, 500, "500 µm")
    ax.annotate("Otp$^+$ Cbln1$^+$ neurons", xy=(-1460, 1740), fontsize=5.8,
                color=POP, va="top", ha="left")
    ax.set_title("representative image", loc="left", pad=2, x=.03)

    # (c) is the increase specific?  Galr1 in every cell type that has enough
    # cells in all eight animals.  Nothing is filtered on significance.
    ax = fig.add_subplot(top[2]); panel(ax, "c", dx=-0.30, dy=1.13)
    by_ct = pd.read_csv(SRC / "galr1_by_celltype.csv")
    by_ct = by_ct.sort_values("lfc")
    ys = np.arange(len(by_ct))
    passes = ((by_ct.blocks == 4) & (by_ct.margin > 0)).to_numpy()
    ax.axvline(0, color=INK, lw=.5)
    ax.scatter(by_ct.lfc[~passes], ys[~passes], s=11, color="#B0B0B0",
               linewidths=0, zorder=3)
    ax.scatter(by_ct.lfc[passes], ys[passes], s=26, color=POP, linewidths=0,
               zorder=4)
    for y, row in zip(ys[passes], by_ct[passes].itertuples()):
        ax.annotate(row.cell_type, (row.lfc, y), xytext=(7, 0),
                    textcoords="offset points", fontsize=5.4, va="center",
                    ha="left", color=POP, fontweight="bold")
    ax.set_yticks([]); ax.set_ylim(-1, len(by_ct))
    ax.set_xlim(-0.75, 1.75)
    ax.set_xlabel("$\\it{Galr1}$, aged / adult (log$_2$)")
    ax.set_ylabel(f"{len(by_ct)} cell types")
    ax.set_title(f"{int(passes.sum())} of {len(by_ct)} cell types change",
                 loc="left", pad=2)

    # (d, e) the receptors in that population
    for j, gene in enumerate(MAIN_GENES):
        ax = fig.add_subplot(bot[j]); panel(ax, "de"[j], dx=-0.34)
        r = stats.loc[gene]
        paired_panel(ax, cpm[gene].to_dict(), f"$\\it{{{gene}}}$ (log$_2$ CPM)",
                     f"{2 ** r.lfc:.2f}×  $P$ = {r.p:.3f}")

    fig.savefig(OUT / "figure_main.pdf")
    fig.savefig(OUT / "figure_main.png")
    plt.close(fig)

    print("whole-section bilateral balance (1.0 = symmetric):")
    for s, v in sorted(balance.items(), key=lambda kv: -kv[1]):
        print(f"   {s:8s} {v:.3f}{'   <- shown' if s == rep else ''}")
    print(f"\n{int(is_pop.sum())} neurons across 8 animals")
    for g in KEY_GENES:
        if g in stats.index:
            r = stats.loc[g]
            print(f"  {g:6s} {2 ** r.lfc:.2f}x  blocks {int(r.blocks)}/4  P {r.p:.4f}")
    print(f"\nwrote {OUT / 'figure_main.pdf'} (+ .png)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
