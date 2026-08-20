#!/usr/bin/env python3
"""Main figure: Galr1 upregulation in dorsal hypothalamic Ghrh neurons with age.

Single biological conclusion the figure supports:

    A defined glutamatergic population of the dorsal hypothalamus -- Galr1+,
    Ghrh+, Ghsr+, Galr3+ -- upregulates Galr1 with age, without loss of neurons.

Panels:
  a  where the neurons can be placed: the registered subregions, every cell
     from all eight animals pooled in the shared anatomical frame.  This is not
     a section -- panel b is the section
  b  one representative section, with every neuron of this type on it
  c  which subregion they belong to: enrichment over the whole window, DMH
     highest
  d  what they are, matched against HypoMap at the C185 level
  e  Galr1 tested in every cell type with enough cells: of twenty, one moves
     consistently and separates the animals, which is what makes the result a
     statement about this population rather than about ageing hypothalamus
  f  Galr1 in that population, one point per animal, paired within block
  g  Ghsr, the second receptor these neurons carry, behaves the same way

Choices made by rule rather than by eye:
  * representative section = the section with the highest whole-section bilateral
    balance, i.e. the most intact tissue.  Balance is the ratio of cells left
    and right of the midline; a torn or folded section scores low.  Restricting
    that test to the analysed window is not enough -- a section can be perfectly
    symmetric inside the window and badly damaged outside it.
  * panel a paints the subregions largest first.  Eight sections are pooled into
    one frame, so cells from different animals share pixels and whichever region
    is drawn last takes them.  In the order the colour table happens to list
    them, the DMH (9 707 cells in the window) was painted over by its three
    larger neighbours -- LHA 28 530, ZI 13 532, DHA_PH 10 648, all drawn after
    it -- so it lost pixels at exactly the borders the parcellation is least
    sure of, and the nucleus panel c makes a claim about read as smaller than
    the assignment behind that claim.  Size order is a rule that applies to
    every region rather than a thumb on the scale for this one.

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

    fig = plt.figure(figsize=(FULL, 5.0))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.05, .95], hspace=.60)
    top = outer[0].subgridspec(1, 3, width_ratios=[1.30, 1.05, .95], wspace=.36)
    bot = outer[1].subgridspec(1, 4, width_ratios=[1.35, 1.15, .95, .95], wspace=.72)

    # (a) the registered subregions this panel resolves
    ax = fig.add_subplot(top[0]); panel(ax, "a", dx=-0.05, dy=1.13)
    nuc = win.obs["nucleus_ext"].astype(str).to_numpy()
    for key in ("edge", "fibre", "TUseg"):
        m = nuc == key
        if m.any():
            ax.scatter(ml[m], dv[m], s=.45, c=TISSUE, linewidths=0, rasterized=True)
    present = [k for k in NUCLEUS_LABEL if (nuc == k).any()]
    # Largest first, so no region is hidden by a bigger one drawn after it.
    for key in sorted(present, key=lambda k: int((nuc == k).sum()), reverse=True):
        m = nuc == key
        ax.scatter(ml[m], dv[m], s=.45, c=NUCLEUS_COLOUR[key], linewidths=0,
                   rasterized=True)
    for key in present:
        m = nuc == key
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

    # (c) which subregion
    ax = fig.add_subplot(top[2]); panel(ax, "c", dx=-0.42, dy=1.13)
    enr = pd.read_csv(SRC / "population_regional_enrichment.csv")
    enr = enr.sort_values("enrichment")
    ys = np.arange(len(enr))
    ax.barh(ys, enr.enrichment, height=.68,
            color=[POP if n == "DMH" else "#C4C4C4" for n in enr.nucleus])
    ax.axvline(1, color=INK, lw=.5)
    ax.set_yticks(ys)
    ax.set_yticklabels([NUCLEUS_LABEL.get(n, n) for n in enr.nucleus])
    for tick, n in zip(ax.get_yticklabels(), enr.nucleus):
        tick.set_fontweight("bold" if n == "DMH" else "normal")
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("enrichment over the whole window")
    ax.set_title("most enriched in DMH", loc="left", pad=2)

    # (d) HypoMap identity, C185 level
    ax = fig.add_subplot(bot[0]); panel(ax, "d", dx=-0.40)
    hm = pd.read_csv(REPO / "results" / "hypomap" /
                     "hypomap_C185_named_all_celltypes.csv", index_col=0)
    top_hm = hm[POPNAME].sort_values(ascending=False).head(6)[::-1]
    names = [i.split(": ", 1)[-1] for i in top_hm.index]
    # top_hm is reversed for plotting, so the best match is the LAST bar.
    ax.barh(range(len(top_hm)), top_hm.values, height=.68,
            color=["#C4C4C4"] * (len(top_hm) - 1) + [POP])
    ax.set_yticks(range(len(top_hm)))
    ax.set_yticklabels(names, fontsize=5.4)
    for tick, keep in zip(ax.get_yticklabels(),
                          [False] * (len(top_hm) - 1) + [True]):
        tick.set_fontweight("bold" if keep else "normal")
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(.6, .95)
    ax.set_xlabel("Spearman ρ to HypoMap C185")
    ax.set_title("HypoMap identity", loc="left", pad=2)

    # (e) is the increase specific?  Galr1 in every cell type with enough cells.
    ax = fig.add_subplot(bot[1]); panel(ax, "e", dx=-0.34)
    by_ct = pd.read_csv(SRC / "galr1_by_celltype.csv").sort_values("lfc")
    ys = np.arange(len(by_ct))
    passes = ((by_ct.blocks == 4) & (by_ct.margin > 0)).to_numpy()
    ax.axvline(0, color=INK, lw=.5)
    ax.scatter(by_ct.lfc[~passes], ys[~passes], s=10, color="#B0B0B0",
               linewidths=0, zorder=3)
    ax.scatter(by_ct.lfc[passes], ys[passes], s=24, color=POP, linewidths=0,
               zorder=4)
    for y, row in zip(ys[passes], by_ct[passes].itertuples()):
        ax.annotate("this population", (row.lfc, y), xytext=(6, 0),
                    textcoords="offset points", fontsize=5.4, va="center",
                    ha="left", color=POP, fontweight="bold")
    ax.set_yticks([]); ax.set_ylim(-1, len(by_ct))
    ax.set_xlim(-0.75, 1.9)
    ax.set_xlabel("$\\it{Galr1}$, aged / adult (log$_2$)")
    ax.set_ylabel(f"{len(by_ct)} cell types")
    ax.set_title(f"{int(passes.sum())} of {len(by_ct)} change", loc="left", pad=2)

    # (f, g) the receptors in that population
    for j, gene in enumerate(MAIN_GENES):
        ax = fig.add_subplot(bot[j + 2]); panel(ax, "fg"[j], dx=-0.40)
        r = stats.loc[gene]
        paired_panel(ax, cpm[gene].to_dict(), f"$\\it{{{gene}}}$ (log$_2$ CPM)",
                     f"{2 ** r.lfc:.2f}×  $P$ = {r.p:.3f}")

    fig.savefig(OUT / "figure_main.pdf")
    fig.savefig(OUT / "figure_main.png")
    plt.close(fig)

    print("whole-section bilateral balance (1.0 = symmetric):")
    for s, v in sorted(balance.items(), key=lambda kv: -kv[1]):
        print(f"   {s:8s} {v:.3f}{'   <- shown' if s == rep else ''}")
    print("\nsubregions drawn in panel a (pooled over the eight sections):")
    for key in sorted(present, key=lambda k: int((nuc == k).sum()), reverse=True):
        m = nuc == key
        a_ml, d = np.abs(ml[m]), dv[m]
        print(f"   {key:7s} {int(m.sum()):6d} cells   "
              f"|ml| {np.percentile(a_ml, 5):4.0f}-{np.percentile(a_ml, 95):4.0f}"
              f"   dv {np.percentile(d, 5):4.0f}-{np.percentile(d, 95):4.0f} um")

    print(f"\n{int(is_pop.sum())} neurons across 8 animals")
    for g in KEY_GENES:
        if g in stats.index:
            r = stats.loc[g]
            print(f"  {g:6s} {2 ** r.lfc:.2f}x  blocks {int(r.blocks)}/4  P {r.p:.4f}")
    print(f"\nwrote {OUT / 'figure_main.pdf'} (+ .png)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
