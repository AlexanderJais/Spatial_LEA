#!/usr/bin/env python3
"""Gal in ARC Pomc neurons, one label per section, all four animals.

Same layout as figure 5a, which showed Gal in ARC Agrp/Npy, so the two read
against each other directly.  The marker is the section number rather than a
dot, so a point that sits apart from its animal's other two can be identified
and gone back to.

Section numbers are slide identifiers, not rostro-caudal order -- the
morphometric AP score for each is printed to the console and written to the
source-data file, since that is what actually orders them.

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

    fig, axes = plt.subplots(1, 2, figsize=(6.0, 3.3))
    fig.subplots_adjust(top=.70, wspace=.50)
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
            ax.plot([i - .30, i + .30], [sub[col].mean()] * 2, color=colour, lw=2,
                    zorder=2, solid_capstyle="butt")
            for _, r in sub.iterrows():
                ax.text(i + OFFSET[r.section_no], r[col], str(r.section_no),
                        color=colour, fontsize=8.5, fontweight="bold",
                        ha="center", va="center", zorder=3)
        # Text artists do not drive autoscaling, so the limits are set from the
        # data or the outlying sections fall outside the axes.
        lo, hi = per[col].min(), per[col].max()
        pad = (hi - lo) * .14
        ax.set_ylim(lo - pad, hi + pad)
        ax.set_xlim(-.6, len(ANIMALS) - .4)
        ax.set_xticks(range(len(ANIMALS)))
        ax.set_xticklabels([f"{a}\n{ANIMAL_META[a]['age_group']}" for a in ANIMALS],
                           fontsize=6.5)
        ax.set_ylabel(ylab)

    axes[0].set_title("Expression", loc="left", fontsize=8)
    axes[1].set_title("Detection", loc="left", fontsize=8)

    aged = per[per.group == "aged"]["gal_cpm"]
    adult = per[per.group == "adult"]["gal_cpm"]
    lfc = np.log2(aged.mean() / adult.mean())
    overlap = aged.min() <= adult.max()
    fig.text(.02, .99, f"{args.gene} in {args.cell_type}", fontsize=11,
             fontweight="bold", va="top")
    fig.text(.02, .90, f"{len(per)} sections, 4 animals. The number is the section; "
             "the bar is the animal mean.", fontsize=7, color="#4A5656", va="top")
    fig.text(.02, .82, f"aged {aged.mean():,.0f} vs adult {adult.mean():,.0f} CPM "
             f"({lfc:+.2f} log2, {2**lfc:.2f}×). Section ranges "
             f"{aged.min():,.0f}–{aged.max():,.0f} and {adult.min():,.0f}–{adult.max():,.0f} "
             + ("overlap." if overlap else "do not overlap."),
             fontsize=7, color="#4A5656", va="top")

    fig.savefig(OUT / f"gal_{tag}_by_section.png")
    fig.savefig(OUT / f"gal_{tag}_by_section.pdf")
    plt.close(fig)

    print(f"=== {args.gene} in {args.cell_type} ===")
    print(per[["section", "animal", "group", "n_cells", "gal_cpm",
               "detect_pct", "ap_score"]].round(2).to_string(index=False))
    print(f"\naged   mean {aged.mean():8.1f} CPM   range {aged.min():.0f}–{aged.max():.0f}")
    print(f"adult  mean {adult.mean():8.1f} CPM   range {adult.min():.0f}–{adult.max():.0f}")
    print(f"log2 aged/adult {lfc:+.3f}")
    per_animal = per.groupby("animal")["gal_cpm"].mean()
    print("\nper-animal means:")
    for a in ANIMALS:
        print(f"  {a:6s} {ANIMAL_META[a]['age_group']:6s} {per_animal[a]:8.1f}")
    print(f"\nwrote {OUT / f'gal_{tag}_by_section.png'} (+ .pdf)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
