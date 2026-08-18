#!/usr/bin/env python3
"""Per-section and per-animal expression for one cell population, any genes.

Layout follows figure 5a: one point per section, positioned by animal, with the
section number as the marker so any section can be identified and gone back to.
The bar is the animal mean.

The right-hand column is the animal level, which is where inference happens --
sections from one mouse are pseudoreplicates.  Animals are paired by block, and
every statistic is computed within blocks, because the two cohorts ran different
panels and different segmentation chemistry and both are aligned with sex.  A
within-block fold change never crosses that boundary.

    python scripts/40_population_expression.py --cell-type "Glut Prdm8/Cbln1" \\
        --genes Galr1 Gal
"""

from __future__ import annotations

import argparse
import sys
from itertools import combinations, product
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
    ADULT, AGED, ANIMAL_META, BLOCKS, ONE_PER_MOUSE, SECTION_ANIMAL, counts_matrix,
)

PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "figures"
SRC = OUT / "source_data"
GRP_COLOUR = {"aged": "#B4531A", "adult": "#0E5A61"}
ANIMALS = list(AGED) + list(ADULT)
MIN_CELLS = 15

plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
    "font.family": "sans-serif",
})


def per_section(adata, cell_type: str, genes: list[str]) -> pd.DataFrame:
    is_type = (adata.obs["cell_type"].astype(str) == cell_type).to_numpy()
    sections = adata.obs["section"].astype(str).to_numpy()
    counts = counts_matrix(adata)
    var = adata.var_names.to_numpy()
    gi = {g: int(np.where(var == g)[0][0]) for g in genes if g in var}
    rows = []
    for section in sorted(SECTION_ANIMAL):
        sel = is_type & (sections == section)
        if sel.sum() < MIN_CELLS:
            continue
        c = counts[sel]
        animal = SECTION_ANIMAL[section]
        row = {"section": section, "section_no": int(section.rsplit("_", 1)[1]),
               "animal": animal, "group": ANIMAL_META[animal]["age_group"],
               "n_cells": int(sel.sum())}
        for g, i in gi.items():
            row[f"{g}_cpm"] = float(c[:, i].sum() / c.sum() * 1e6)
            row[f"{g}_detect"] = float((c[:, i] > 0).mean() * 100)
        rows.append(row)
    return pd.DataFrame(rows)


def stats_for(adata, cell_type: str, gene: str) -> dict:
    """Blocked and exact tests on per-animal pseudobulk."""
    is_type = (adata.obs["cell_type"].astype(str) == cell_type).to_numpy()
    animals = adata.obs["animal"].astype(str).to_numpy()
    counts = counts_matrix(adata)
    gi = int(np.where(adata.var_names.to_numpy() == gene)[0][0])
    cpm = {}
    for a in ANIMALS:
        c = counts[is_type & (animals == a)]
        cpm[a] = float(np.log2(c[:, gi].sum() / c.sum() * 1e6 + 1))
    lfc = {b: cpm[x] - cpm[y] for b, (x, y) in BLOCKS.items()}
    mean = float(np.mean(list(lfc.values())))
    flips = np.array(list(product([1, -1], repeat=len(BLOCKS))))
    null = flips @ np.array(list(lfc.values())) / len(BLOCKS)
    p_blocked = float((np.abs(null) >= abs(mean)).mean())

    vals = np.array([cpm[a] for a in ANIMALS])
    obs = vals[:len(AGED)].mean() - vals[len(AGED):].mean()
    idx = list(range(len(ANIMALS)))
    null2 = np.array([vals[list(c)].mean() - vals[[i for i in idx if i not in c]].mean()
                      for c in combinations(idx, len(AGED))])
    p_exact = float((np.abs(null2) >= abs(obs)).mean())
    return {"gene": gene, "cpm": cpm, "lfc": lfc, "mean_lfc": mean,
            "p_blocked": p_blocked, "p_exact": p_exact,
            "agree": int(sum(np.sign(v) == np.sign(mean) for v in lfc.values()))}


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cell-type", default="Glut Prdm8/Cbln1")
    p.add_argument("--genes", nargs="+", default=["Galr1", "Gal"])
    p.add_argument("--one-per-mouse", action="store_true")
    p.add_argument("--suffix", default="")
    args = p.parse_args()
    OUT.mkdir(parents=True, exist_ok=True); SRC.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(PROC / "hypothalamus_nuclei.h5ad")
    if args.one_per_mouse:
        keep_s = adata.obs["section"].astype(str).isin(ONE_PER_MOUSE).to_numpy()
        adata = adata[keep_s].copy()
        print("One slide per mouse: " + ", ".join(sorted(ONE_PER_MOUSE)))
    per = per_section(adata, args.cell_type, args.genes)
    stats = {g: stats_for(adata, args.cell_type, g) for g in args.genes}
    tag = args.cell_type.lower().replace(" ", "_").replace("/", "_") + args.suffix
    per.to_csv(SRC / f"{tag}_per_section.csv", index=False)

    n = len(args.genes)
    fig, axes = plt.subplots(n, 2, figsize=(8.4, 3.1 * n),
                             gridspec_kw={"width_ratios": [2.05, 1]})
    axes = np.atleast_2d(axes)
    fig.subplots_adjust(top=.80, hspace=.55, wspace=.30)
    # Sections nudged apart on x so labels never collide; the offset means nothing.
    OFF = {1: -.22, 2: .0, 3: .22}

    for r, gene in enumerate(args.genes):
        st = stats[gene]
        col = f"{gene}_cpm"

        ax = axes[r, 0]
        for i, animal in enumerate(ANIMALS):
            sub = per[per.animal == animal]
            if not len(sub):
                continue
            colour = GRP_COLOUR[ANIMAL_META[animal]["age_group"]]
            ax.plot([i - .32, i + .32], [sub[col].mean()] * 2, color=colour, lw=2.4,
                    solid_capstyle="butt", zorder=2)
            for _, row in sub.iterrows():
                ax.text(i + OFF[row.section_no], row[col], str(row.section_no),
                        color=colour, fontsize=8.5, fontweight="bold",
                        ha="center", va="center", zorder=3)
        lo, hi = per[col].min(), per[col].max()
        pad = (hi - lo) * .16
        ax.set_ylim(lo - pad, hi + pad)
        ax.set_xlim(-.6, len(ANIMALS) - .4)
        ax.set_xticks(range(len(ANIMALS)))
        ax.set_xticklabels([f"{a}\n{ANIMAL_META[a]['age_group']}\n{ANIMAL_META[a]['batch']}"
                            for a in ANIMALS], fontsize=6)
        ax.set_ylabel(f"{gene} (CPM)")
        ax.set_title(f"{gene} — one label per section", loc="left")
        ax.text(-0.10, 1.12, "ab"[r] if r < 2 else chr(97 + r), transform=ax.transAxes,
                fontsize=11, fontweight="bold", va="top")

        # Animal level, paired by block.
        ax = axes[r, 1]
        adult_y = {b: st["cpm"][d] for b, (_, d) in BLOCKS.items()}
        span = max(adult_y.values()) - min(adult_y.values()) or 1.0
        # Nudge block labels apart where two adults sit at nearly the same value,
        # otherwise they print on top of each other.
        placed, offsets = [], {}
        for b, y in sorted(adult_y.items(), key=lambda kv: kv[1]):
            shift = 0.0
            while any(abs((y + shift) - q) < span * .055 for q in placed):
                shift += span * .055
            offsets[b] = shift
            placed.append(y + shift)
        for b, (a_aged, a_adult) in BLOCKS.items():
            ya, yd = st["cpm"][a_aged], st["cpm"][a_adult]
            ax.plot([0, 1], [ya, yd], color="#B8C2C2", lw=1.0, zorder=1)
            ax.scatter([0], [ya], s=46, color=GRP_COLOUR["aged"], zorder=3)
            ax.scatter([1], [yd], s=46, color=GRP_COLOUR["adult"], zorder=3)
            ax.annotate(b, (1, yd + offsets[b]), xytext=(8, 0),
                        textcoords="offset points", fontsize=6.5, va="center",
                        color="#5A6A6A")
        for j, grp in enumerate(("aged", "adult")):
            members = AGED if grp == "aged" else ADULT
            m = np.mean([st["cpm"][a] for a in members])
            ax.plot([j - .18, j + .18], [m] * 2, color=GRP_COLOUR[grp], lw=2.4,
                    solid_capstyle="butt", zorder=2)
        ax.set_xlim(-.45, 1.6); ax.set_xticks([0, 1])
        ax.set_xticklabels(["aged", "adult"])
        ax.set_ylabel(f"{gene} (log2 CPM), animal")
        ax.set_title(f"{st['agree']}/4 blocks agree · {2**st['mean_lfc']:.2f}×\n"
                     f"p(blocked) {st['p_blocked']:.3f} · p(exact) {st['p_exact']:.3f}",
                     loc="left", fontsize=8)
        ax.text(-0.28, 1.12, "cd"[r] if r < 2 else chr(99 + r), transform=ax.transAxes,
                fontsize=11, fontweight="bold", va="top")

    fig.text(.02, .99, f"{args.cell_type}", fontsize=11, fontweight="bold", va="top")
    fig.text(.02, .935, f"{len(per)} sections, {per.animal.nunique()} animals "
             f"({per[per.group=='aged'].animal.nunique()} aged vs "
             f"{per[per.group=='adult'].animal.nunique()} adult), "
             f"{int(per.n_cells.min())}–{int(per.n_cells.max())} cells per section. "
             "Left: number = section, bar = animal mean.", fontsize=7, color="#4A5656",
             va="top")
    fig.text(.02, .885, "Right: one point per animal, joined within block — each block "
             "is one aged and one adult mouse on the same slide run, panel and "
             "segmentation chemistry.", fontsize=7, color="#4A5656", va="top")

    fig.savefig(OUT / f"{tag}_expression.png")
    fig.savefig(OUT / f"{tag}_expression.pdf")
    plt.close(fig)

    print(f"=== {args.cell_type} ===")
    print(per.to_string(index=False))
    for g, st in stats.items():
        print(f"\n{g}: mean blocked lfc {st['mean_lfc']:+.3f} "
              f"({2**st['mean_lfc']:.2f}x), {st['agree']}/4 blocks agree, "
              f"p_blocked {st['p_blocked']:.4f}, p_exact {st['p_exact']:.4f}")
        for b, (a, d) in BLOCKS.items():
            print(f"   {b}  {a:6s} {st['cpm'][a]:6.2f}  vs  {d:6s} "
                  f"{st['cpm'][d]:6.2f}   lfc {st['lfc'][b]:+.2f}")
    print(f"\nwrote {OUT / f'{tag}_expression.png'} (+ .pdf)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
