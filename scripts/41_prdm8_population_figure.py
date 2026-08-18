#!/usr/bin/env python3
"""Where the Glut Prdm8/Cbln1 neurons are, and how they change with age.

Note on the location.  The tissue-derived domain that most of these cells fall
into is labelled ME_3V, but that domain is a broad periventricular band (its
mediolateral extent runs to 1341 um), so the label overstates how medial they
are.  Measured directly they sit in two bilateral columns: median |ml| 591 um,
only 18% within 250 um of the midline and 58% beyond 500 um, at dv ~1000 um.
That is dorsal hypothalamic level, not median eminence, and the figure says so
rather than repeating the domain name.

One slide per mouse, four aged against four adult, in four blocks.  Every
statistic is computed within blocks, so no comparison crosses the panel or
segmentation-chemistry boundary between cohorts.

  a  where they sit, in the anatomical frame, pooled over the eight sections
  b  which nuclei they fall in
  c  what they are -- enrichment over the rest of the hypothalamic window
  d  the three headline genes, one point per animal, joined within block
  e  every gene, effect against consistency
  f  abundance is unchanged, so this is expression and not loss of neurons
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

from spatial_lea.io import (  # noqa: E402
    ADULT, AGED, ANIMAL_META, BLOCKS, ONE_PER_MOUSE, counts_matrix,
)

PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "figures"
SRC = OUT / "source_data"
POP = "Glut Prdm8/Cbln1"
ANIMALS = list(AGED) + list(ADULT)
GRP = {"aged": "#B4531A", "adult": "#0E5A61"}
HEADLINE = ["Fos", "Galr1", "Ghsr"]
MARKERS = 12

plt.rcParams.update({
    "font.size": 7.5, "axes.titlesize": 8, "axes.labelsize": 7.5,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
    "font.family": "sans-serif",
})


def label(ax, letter, dx=-0.18, dy=1.14):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="top", ha="left")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True); SRC.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(PROC / "hypothalamus_nuclei.h5ad")
    adata = adata[adata.obs["section"].astype(str).isin(ONE_PER_MOUSE)].copy()
    sel = (adata.obs["cell_type"].astype(str) == POP).to_numpy()
    animals = adata.obs["animal"].astype(str).to_numpy()
    counts = counts_matrix(adata)
    var = adata.var_names.to_numpy()

    # Per-animal pseudobulk, blocked statistics.
    mat = pd.DataFrame({a: counts[sel & (animals == a)].sum(axis=0) for a in ANIMALS},
                       index=var).T
    expressed = mat.sum(axis=0) >= 200
    cpm = np.log2(mat.div(mat.sum(axis=1), axis=0) * 1e6 + 1)
    lfc = pd.DataFrame({b: cpm.loc[x] - cpm.loc[y] for b, (x, y) in BLOCKS.items()}).T
    mean = lfc.mean()
    agree = (np.sign(lfc) == np.sign(mean)).sum()
    vals = cpm.to_numpy(); idx = list(range(8))
    obs = vals[:4].mean(axis=0) - vals[4:].mean(axis=0)
    null = np.array([vals[list(c)].mean(axis=0)
                     - vals[[i for i in idx if i not in c]].mean(axis=0)
                     for c in combinations(idx, 4)])
    p_exact = pd.Series((np.abs(null) >= np.abs(obs)).mean(axis=0), index=var)
    stats = pd.DataFrame({"lfc": mean, "blocks": agree, "p": p_exact})[expressed]
    stats.to_csv(SRC / "prdm8_age_stats.csv")

    fig = plt.figure(figsize=(7.4, 8.2))
    gs = fig.add_gridspec(3, 3, height_ratios=[1.15, .88, .88],
                          width_ratios=[1.25, 1, 1], hspace=.78, wspace=.62)

    # (a) where they are
    ax = fig.add_subplot(gs[0, :2]); label(ax, "a", dx=-0.09)
    ml, dv = adata.obs["ml"].to_numpy(), adata.obs["dv"].to_numpy()
    bg = np.random.default_rng(0).choice(np.where(~sel)[0], 40000, replace=False)
    ax.scatter(ml[bg], dv[bg], s=.6, c="#E4E9E9", rasterized=True, linewidths=0)
    ax.scatter(ml[sel], dv[sel], s=7, c="#6A3D9A", rasterized=True, linewidths=0)
    ax.set_xlim(-1450, 1450); ax.set_ylim(-100, 1800); ax.set_aspect("equal")
    ax.set_xlabel("mediolateral from 3rd ventricle (µm)")
    ax.set_ylabel("dorsal from ventral surface (µm)")
    ax.set_title(f"{int(sel.sum())} {POP} neurons, 8 sections pooled", loc="left")
    absml = np.abs(ml[sel])
    ax.annotate(f"two bilateral columns\nmedian |ml| {np.median(absml):.0f} µm, "
                f"dv {np.median(dv[sel]):.0f} µm\n"
                f"{(absml > 500).mean()*100:.0f}% lie beyond 500 µm from the midline",
                xy=(.02, .97), xycoords="axes fraction", fontsize=6,
                color="#4A5656", va="top")

    # (b) nuclei
    ax = fig.add_subplot(gs[0, 2]); label(ax, "b", dx=-0.34)
    nx = adata.obs["nucleus_ext"].astype(str)
    share = (nx[sel].value_counts(normalize=True) * 100).head(6).sort_values()
    ax.barh(range(len(share)), share.values, color="#6A3D9A")
    ax.set_yticks(range(len(share))); ax.set_yticklabels(share.index, fontsize=6.5)
    ax.set_xlabel("% of these neurons")
    ax.set_title("Domains they fall in\n(these are broad bands)", loc="left")

    # (c) identity
    ax = fig.add_subplot(gs[1, 0]); label(ax, "c", dx=-0.30)
    cin = counts[sel].sum(axis=0) / counts[sel].sum() * 1e6
    cout = counts[~sel].sum(axis=0) / counts[~sel].sum() * 1e6
    enr = pd.Series(np.log2((cin + 1) / (cout + 1)), index=var).nlargest(MARKERS)
    enr = enr.sort_values()
    cols = ["#B4531A" if g in ("Galr1", "Galr3", "Ghsr", "Ghrh") else "#9AA8A8"
            for g in enr.index]
    ax.barh(range(len(enr)), enr.values, color=cols)
    ax.set_yticks(range(len(enr)))
    ax.set_yticklabels(enr.index, fontsize=6, style="italic")
    ax.set_xlabel("log2 vs rest of window")
    ax.set_title("Identity", loc="left")

    # (d) headline genes, paired within block
    for j, gene in enumerate(HEADLINE):
        ax = fig.add_subplot(gs[1, j + 1] if j < 2 else gs[2, 0])
        label(ax, "def"[j], dx=-0.34)
        for b, (a_aged, a_adult) in BLOCKS.items():
            ya, yd = cpm.loc[a_aged, gene], cpm.loc[a_adult, gene]
            ax.plot([0, 1], [ya, yd], color="#C8D0D0", lw=1.0, zorder=1)
            ax.scatter([0], [ya], s=34, color=GRP["aged"], zorder=3)
            ax.scatter([1], [yd], s=34, color=GRP["adult"], zorder=3)
        for k, grp in enumerate(("aged", "adult")):
            members = AGED if grp == "aged" else ADULT
            ax.plot([k - .2, k + .2], [cpm.loc[list(members), gene].mean()] * 2,
                    color=GRP[grp], lw=2.4, solid_capstyle="butt", zorder=2)
        ax.set_xlim(-.5, 1.5); ax.set_xticks([0, 1])
        ax.set_xticklabels(["aged", "adult"], fontsize=6.5)
        ax.set_ylabel(f"{gene} (log2 CPM)")
        r = stats.loc[gene]
        ax.set_title(f"$\\it{{{gene}}}$  {2**r.lfc:.2f}×\n{int(r.blocks)}/4 blocks · "
                     f"p {r.p:.3f}", loc="left")

    # (e) every gene
    ax = fig.add_subplot(gs[2, 1]); label(ax, "g", dx=-0.32)
    four = stats.blocks == 4
    ax.scatter(stats.lfc[~four], -np.log10(stats.p[~four] + 1e-3), s=8,
               c="#D3DADA", linewidths=0)
    ax.scatter(stats.lfc[four], -np.log10(stats.p[four] + 1e-3), s=11,
               c="#0E5A61", linewidths=0)
    # Fos and Tacr1 sit left, Galr1 and Ghsr right and nearly on top of each
    # other, so the two right-hand labels are stacked rather than overlaid.
    nudge = {"Fos": (-4, 5, "right"), "Tacr1": (-4, 5, "right"),
             "Galr1": (5, 6, "left"), "Ghsr": (5, -8, "left")}
    for g in HEADLINE + ["Tacr1"]:
        if g not in stats.index:
            continue
        x, y = stats.lfc[g], -np.log10(stats.p[g] + 1e-3)
        ax.scatter([x], [y], s=26, c="#B4531A", zorder=4, linewidths=0)
        dx, dy, ha = nudge.get(g, (4, 3, "left"))
        ax.annotate(g, (x, y), xytext=(dx, dy), textcoords="offset points",
                    fontsize=6, style="italic", color="#B4531A", ha=ha)
    ax.axvline(0, color="#999", lw=.6)
    ax.set_xlabel("log2 aged / adult"); ax.set_ylabel("−log10 exact p")
    ax.set_title(f"{int(four.sum())}/{len(stats)} genes agree\nin all 4 blocks",
                 loc="left")

    # (f) abundance
    ax = fig.add_subplot(gs[2, 2]); label(ax, "h", dx=-0.34)
    pct = {a: (sel & (animals == a)).sum() / (animals == a).sum() * 100
           for a in ANIMALS}
    for b, (a_aged, a_adult) in BLOCKS.items():
        ax.plot([0, 1], [pct[a_aged], pct[a_adult]], color="#C8D0D0", lw=1.0, zorder=1)
        ax.scatter([0], [pct[a_aged]], s=34, color=GRP["aged"], zorder=3)
        ax.scatter([1], [pct[a_adult]], s=34, color=GRP["adult"], zorder=3)
    for k, grp in enumerate(("aged", "adult")):
        members = AGED if grp == "aged" else ADULT
        ax.plot([k - .2, k + .2], [np.mean([pct[a] for a in members])] * 2,
                color=GRP[grp], lw=2.4, solid_capstyle="butt", zorder=2)
    ax.set_xlim(-.5, 1.5); ax.set_xticks([0, 1])
    ax.set_xticklabels(["aged", "adult"], fontsize=6.5)
    ax.set_ylabel("% of window cells")
    ax.set_title("Abundance unchanged —\nexpression, not cell loss", loc="left")

    fig.text(.02, .995, "Dorsal hypothalamic Prdm8/Cbln1 neurons: less active, more receptor",
             fontsize=10.5, fontweight="bold", va="top")
    fig.text(.02, .967, "One slide per mouse · 4 aged vs 4 adult in 4 blocks · every "
             "statistic computed within blocks", fontsize=7, color="#4A5656", va="top")
    fig.savefig(OUT / "prdm8_cbln1_population.png")
    fig.savefig(OUT / "prdm8_cbln1_population.pdf")
    plt.close(fig)
    print(f"wrote {OUT / 'prdm8_cbln1_population.png'} (+ .pdf)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
