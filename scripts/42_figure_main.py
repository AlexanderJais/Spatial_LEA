#!/usr/bin/env python3
"""Main figure: Galr1 upregulation in dorsal hypothalamic Ghrh neurons with age.

Single biological conclusion the figure supports:

    A defined glutamatergic population of the dorsal hypothalamus -- Galr1+,
    Ghrh+, Ghsr+, Galr3+ -- upregulates Galr1 with age, without loss of neurons.

Panels:
  a  where the neurons are: one intact section, then every neuron of this type
     from all eight animals in the shared anatomical frame
  b  what they are, with the four genes that define them for this study
  c  Galr1 rises      } one point per animal, adult and aged joined within block
  d  Ghsr rises       }
  e  abundance is unchanged, so this is regulation and not cell loss
  f  every gene that separates the groups completely, for context

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
MAIN_GENES = ["Galr1", "Gal", "Ghsr"]
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
    outer = fig.add_gridspec(2, 1, height_ratios=[1.06, .94], hspace=.50)
    top = outer[0].subgridspec(1, 3, width_ratios=[1.30, 1.05, .95], wspace=.36)
    bot = outer[1].subgridspec(1, 4, width_ratios=[1, 1, 1, 1.05], wspace=.62)

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

    # (c) what defines them
    ax = fig.add_subplot(top[2]); panel(ax, "c", dx=-0.40, dy=1.13)
    frac = pd.Series((counts[is_pop] > 0).mean(axis=0) * 100, index=var)
    show = ["Cbln1", "Slc17a6", "Otp", "Prdm8", "Bdnf", "Galr1", "Ghsr",
            "Galr3", "Gal", "Ghrh"]
    show = [g for g in show if g in frac.index][::-1]
    ax.barh(range(len(show)), frac[show].values, height=.7,
            color=[POP if g in ("Galr1", "Gal") else "#B9B9B9" for g in show])
    ax.set_yticks(range(len(show)))
    ax.set_yticklabels([f"$\\it{{{g}}}$" for g in show])
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("% of these neurons positive")
    ax.set_xlim(0, 105)

    # (d-f) the galanin and ghrelin receptors
    for j, gene in enumerate(MAIN_GENES):
        ax = fig.add_subplot(bot[j]); panel(ax, "def"[j], dx=-0.48)
        r = stats.loc[gene]
        sig = "$P$ = %.3f" % r.p
        paired_panel(ax, cpm[gene].to_dict(), f"$\\it{{{gene}}}$ (log$_2$ CPM)",
                     f"{2 ** r.lfc:.2f}×  {sig}")

    # (g) every gene that separates the groups completely
    ax = fig.add_subplot(bot[3]); panel(ax, "g", dx=-0.34)
    sep = stats[(stats.blocks == 4) & (stats.margin > 0)]
    sep = sep.reindex(sep.lfc.sort_values().index)
    ys = np.arange(len(sep))
    for b in BLOCKS:
        ax.scatter(sep[f"lfc_{b}"], ys, s=4, facecolors="none",
                   edgecolors="#B8B8B8", linewidths=.4, zorder=2)
    ax.scatter(sep["lfc"], ys, s=9, color=INK, zorder=3, linewidths=0)
    ax.axvline(0, color=INK, lw=.5)
    ax.set_yticks([]); ax.set_ylim(-1.2, len(sep) + .4)
    ax.set_xlabel("aged / adult (log$_2$)")
    ax.set_ylabel(f"{len(sep)} genes, ranked")
    for g, dyy in (("Galr1", 6), ("Ghsr", -6)):
        if g in sep.index:
            y = int(np.where(sep.index == g)[0][0])
            ax.annotate(f"$\\it{{{g}}}$", (sep.lfc[g], y), xytext=(7, dyy),
                        textcoords="offset points", fontsize=5.5, va="center",
                        ha="left", color=INK)
    ax.set_xlim(sep[[f"lfc_{b}" for b in BLOCKS]].min().min() - .4,
                sep[[f"lfc_{b}" for b in BLOCKS]].max().max() + 1.4)

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
