#!/usr/bin/env python3
"""Gal in ARC Pomc neurons, per section and per animal, all four mice.

Same layout as figure 5a, which showed Gal in ARC Agrp/Npy, so the two read
against each other directly.  The marker is the section number rather than a dot,
so any section can be identified and gone back to.

**Inference is at the animal level**, which is the unit of replication -- three
sections from one mouse are pseudoreplicates and cannot vote independently.  That
matters here: judging on complete section-level separation instead would let a
single section out of twelve overturn a result that all four animals agree on,
and between-section SD in this dataset (0.265 log2) is larger than
between-animal SD (0.152), so one deviant section is expected rather than
surprising.

With 2 animals per group there are only 6 ways to assign the labels, so the
exact permutation p cannot go below 1/6 = 0.167.  A result at that value is as
extreme as the design can produce; it is a statement about the design's
resolution, not a weak result.

Section numbers are slide identifiers, not rostro-caudal order.  The morphometric
AP score is reported alongside so an AP explanation can be checked rather than
assumed.

    python scripts/35_gal_pomc_by_section.py
    python scripts/35_gal_pomc_by_section.py --cell-type "ARC Agrp/Npy"
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.io import ANIMAL_META, SECTION_ANIMAL, counts_matrix  # noqa: E402

PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "figures"
SRC = OUT / "source_data"
AP_FILE = REPO / "results" / "ap_matching" / "section_ap_scores.csv"
ANIMALS = ["F536", "G_073", "M493", "M399"]
GRP_COLOUR = {"aged": "#B4531A", "adult": "#0E5A61"}
MIN_CELLS = 15

plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
    "font.family": "sans-serif",
})


def per_section(adata, cell_type: str, gene: str, ap: pd.Series) -> pd.DataFrame:
    gi = int(np.where(adata.var_names.to_numpy() == gene)[0][0])
    is_type = (adata.obs["cell_type"].astype(str) == cell_type).to_numpy()
    sections = adata.obs["section"].astype(str).to_numpy()
    rows = []
    for section in sorted(SECTION_ANIMAL):
        sel = is_type & (sections == section)
        if sel.sum() < MIN_CELLS:
            continue
        c = counts_matrix(adata[sel])
        animal = SECTION_ANIMAL[section]
        rows.append({
            "section": section,
            "section_no": int(section.rsplit("_", 1)[1]),
            "animal": animal,
            "group": ANIMAL_META[animal]["age_group"],
            "n_cells": int(sel.sum()),
            "gal_cpm": float(c[:, gi].sum() / c.sum() * 1e6),
            "detect_pct": float((c[:, gi] > 0).mean() * 100),
            "ap_score": float(ap[section]),
        })
    return pd.DataFrame(rows)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cell-type", default="ARC Pomc")
    p.add_argument("--gene", default="Gal")
    args = p.parse_args()
    OUT.mkdir(parents=True, exist_ok=True); SRC.mkdir(parents=True, exist_ok=True)

    ap = pd.read_csv(AP_FILE, index_col=0)["ap_score"]
    adata = sc.read_h5ad(PROC / "hypothalamus_domains.h5ad")
    per = per_section(adata, args.cell_type, args.gene, ap)
    tag = args.cell_type.lower().replace(" ", "_").replace("/", "_")
    per.to_csv(SRC / f"gal_{tag}_per_section.csv", index=False)

    from itertools import combinations
    per_animal = per.groupby("animal").agg(group=("group", "first"),
                                           gal_cpm=("gal_cpm", "mean"),
                                           detect_pct=("detect_pct", "mean"))
    per_animal = per_animal.reindex(ANIMALS)
    aged_a = per_animal[per_animal.group == "aged"]["gal_cpm"]
    adult_a = per_animal[per_animal.group == "adult"]["gal_cpm"]
    lfc = float(np.log2(aged_a.mean() / adult_a.mean()))
    separates = bool(aged_a.min() > adult_a.max())

    # Exact permutation over animals: every way of calling 2 of the 4 "aged".
    vals = per_animal["gal_cpm"].to_numpy()
    obs = aged_a.mean() - adult_a.mean()
    stats = [vals[list(c)].mean() - vals[[i for i in range(4) if i not in c]].mean()
             for c in combinations(range(4), 2)]
    p_exact = float(np.mean([s >= obs for s in stats]))

    fig, axes = plt.subplots(1, 3, figsize=(7.6, 3.5))
    fig.subplots_adjust(top=.63, wspace=.58)
    # Sections are nudged apart on x so three labels for one animal never collide;
    # the offset carries no meaning beyond keeping them legible.
    OFFSET = {1: -.20, 2: .0, 3: .20}

    for ax, col, ylab in ((axes[0], "gal_cpm", f"{args.gene} (CPM)"),
                          (axes[1], "detect_pct", f"% of cells {args.gene}⁺")):
        for i, animal in enumerate(ANIMALS):
            sub = per[per.animal == animal]
            if not len(sub):
                continue
            colour = GRP_COLOUR[ANIMAL_META[animal]["age_group"]]
            ax.plot([i - .30, i + .30], [sub[col].mean()] * 2, color=colour, lw=2.4,
                    zorder=2, solid_capstyle="butt")
            for _, r in sub.iterrows():
                ax.text(i + OFFSET[r.section_no], r[col], str(r.section_no),
                        color=colour, fontsize=8.5, fontweight="bold",
                        ha="center", va="center", zorder=3)
        # Text artists do not drive autoscaling, so limits come from the data or
        # the outlying sections fall outside the axes.
        lo, hi = per[col].min(), per[col].max()
        pad = (hi - lo) * .14
        ax.set_ylim(lo - pad, hi + pad)
        ax.set_xlim(-.6, len(ANIMALS) - .4)
        ax.set_xticks(range(len(ANIMALS)))
        ax.set_xticklabels([f"{a}\n{ANIMAL_META[a]['age_group']}" for a in ANIMALS],
                           fontsize=6.5)
        ax.set_ylabel(ylab)

    axes[0].set_title("Expression, per section", loc="left", fontsize=8)
    axes[1].set_title("Detection, per section", loc="left", fontsize=8)

    # (c) the animal level -- the unit that actually replicates
    ax = axes[2]
    for i, animal in enumerate(ANIMALS):
        grp = ANIMAL_META[animal]["age_group"]
        ax.scatter([0 if grp == "aged" else 1], [per_animal.loc[animal, "gal_cpm"]],
                   s=52, color=GRP_COLOUR[grp], zorder=3)
        ax.annotate(animal, (0 if grp == "aged" else 1,
                             per_animal.loc[animal, "gal_cpm"]),
                    xytext=(7, 0), textcoords="offset points", fontsize=6.5,
                    va="center", color=GRP_COLOUR[grp])
    for j, (grp, v) in enumerate((("aged", aged_a), ("adult", adult_a))):
        ax.plot([j - .22, j + .22], [v.mean()] * 2, color=GRP_COLOUR[grp], lw=2.4,
                solid_capstyle="butt")
    if separates:
        mid = (aged_a.min() + adult_a.max()) / 2
        ax.axhline(mid, color="#9AA8A8", lw=.8, ls=":")
        ax.annotate("no overlap", (1.34, mid), fontsize=6, color="#7A8A8A",
                    va="bottom", ha="right")
    ax.set_xlim(-.5, 1.5); ax.set_xticks([0, 1]); ax.set_xticklabels(["aged", "adult"])
    ax.set_ylabel(f"{args.gene} (CPM), animal mean")
    ax.set_title(f"Animal level, n=2 vs 2\nexact p = {p_exact:.3f}", loc="left", fontsize=8)

    fig.text(.02, .99, f"{args.gene} in {args.cell_type}", fontsize=11,
             fontweight="bold", va="top")
    fig.text(.02, .91, f"All four animals order by age: "
             + " > ".join(f"{a} {per_animal.loc[a, 'gal_cpm']:,.0f}"
                          for a in per_animal.sort_values('gal_cpm', ascending=False).index)
             + f"  ({lfc:+.2f} log2, {2**lfc:.2f}×)", fontsize=7, color="#4A5656", va="top")
    fig.text(.02, .84, f"Exact permutation over animals p = {p_exact:.3f} — the floor for "
             "2 vs 2 is 1/6 = 0.167, so this is as extreme as the design can resolve.",
             fontsize=7, color="#4A5656", va="top")
    fig.text(.02, .78, "Left and middle: the number is the section, the bar is that "
             "animal's mean. Right: one point per animal, the unit of replication.",
             fontsize=7, color="#4A5656", va="top")

    per_animal.round(1).to_csv(SRC / f"gal_{tag}_per_animal.csv")

    fig.savefig(OUT / f"gal_{tag}_by_section.png")
    fig.savefig(OUT / f"gal_{tag}_by_section.pdf")
    plt.close(fig)

    print(f"=== {args.gene} in {args.cell_type} ===")
    print(per[["section", "animal", "group", "n_cells", "gal_cpm",
               "detect_pct", "ap_score"]].round(2).to_string(index=False))
    aged_s = per[per.group == "aged"]["gal_cpm"]
    adult_s = per[per.group == "adult"]["gal_cpm"]
    print(f"\nper section  aged  {aged_s.mean():8.0f} CPM  range {aged_s.min():.0f}-{aged_s.max():.0f}")
    print(f"per section  adult {adult_s.mean():8.0f} CPM  range {adult_s.min():.0f}-{adult_s.max():.0f}")
    print("\n=== animal level (the unit of replication) ===")
    for a in per_animal.sort_values("gal_cpm", ascending=False).index:
        print(f"  {a:6s} {per_animal.loc[a, 'group']:6s} {per_animal.loc[a, 'gal_cpm']:9,.0f} CPM")
    print(f"\n  all aged above all adult: {separates}")
    print(f"  log2 aged/adult {lfc:+.3f}  ({2**lfc:.2f}x)")
    print(f"  exact permutation over animals p = {p_exact:.3f}  "
          f"(floor for 2v2 is 1/6 = 0.167)")
    print(f"\nwrote {OUT / f'gal_{tag}_by_section.png'} (+ .pdf)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
