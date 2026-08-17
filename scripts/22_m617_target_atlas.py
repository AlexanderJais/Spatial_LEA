#!/usr/bin/env python3
"""Galr1 target atlas — what M617 can act on, and where.

M617 is a Galr1-selective agonist.  For interpreting a treatment experiment the
question is not "what changed with age" but "which cells carry the receptor,
what are they, and where do they sit".  That is descriptive: 165k cells, no age
contrast, none of the n=2 limitation.

Four things a treatment study needs and this data can give:

1. **Where the receptor is.** Not just which cell types are Galr1-positive, but
   what share of the total Galr1 signal each population carries -- a rare type
   with high expression and an abundant type with low expression are very
   different drug targets.
2. **What those cells do.** Galr1 is Gi/Go-coupled, so an agonist inhibits the
   cell that carries it. Whether Galr1 sits on GABAergic or glutamatergic
   neurons decides whether M617 inhibits a circuit or disinhibits it.
3. **Off-target risk.** Galr3 is also Gi-coupled and also in the panel;
   populations co-expressing both are where Galr1 selectivity matters least.
4. **Endogenous tone.** Where Gal is expressed relative to the receptor sets the
   baseline the agonist adds to -- a receptor already saturated by local ligand
   should respond differently from an unoccupied one.

Adult vs aged is reported throughout as a descriptive contrast of the target
population, per animal, with no p-value claimed.
"""

from __future__ import annotations

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

from spatial_lea.io import ANIMAL_META, counts_matrix  # noqa: E402

PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "m617_atlas"
ANIMALS = ["F536", "G_073", "M493", "M399"]
NUCLEI = ("ARC", "ME_3V", "VMH", "DMH", "LHA", "ZI", "DHA_PH")
GABA = ["Gad1", "Gad2", "Slc32a1"]
GLUT = ["Slc17a6", "Slc17a7"]
MIN_CELLS = 150


def transmitter_class(counts: pd.DataFrame) -> pd.Series:
    """GABAergic, glutamatergic, both or neither, per cell."""
    gaba = counts[[g for g in GABA if g in counts]].sum(axis=1) > 0
    glut = counts[[g for g in GLUT if g in counts]].sum(axis=1) > 0
    out = pd.Series("non-neuronal", index=counts.index)
    out[gaba & ~glut] = "GABAergic"
    out[glut & ~gaba] = "glutamatergic"
    out[gaba & glut] = "mixed"
    return out


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(PROC / "hypothalamus_nuclei.h5ad")
    adata = adata[adata.obs["nucleus_ext"].isin(NUCLEI)].copy()
    genes = adata.var_names.to_numpy()

    need = ["Galr1", "Galr3", "Gal"] + GABA + GLUT
    present = [g for g in need if g in genes]
    mat = pd.DataFrame(counts_matrix(adata[:, present]), columns=present,
                       index=adata.obs_names)
    total = pd.Series(counts_matrix(adata).sum(axis=1), index=adata.obs_names)

    obs = adata.obs.copy()
    obs["galr1"] = mat["Galr1"].to_numpy()
    obs["galr3"] = mat["Galr3"].to_numpy()
    obs["gal"] = mat["Gal"].to_numpy()
    obs["nt"] = transmitter_class(mat).to_numpy()
    obs["total"] = total.to_numpy()
    obs["pop"] = obs["nucleus_ext"].astype(str) + " | " + obs["cell_type"].astype(str)

    galr1_total = obs["galr1"].sum()
    print(f"{adata.n_obs:,} cells across {len(NUCLEI)} nuclei; "
          f"{int((obs.galr1 > 0).sum()):,} are Galr1-positive "
          f"({(obs.galr1 > 0).mean()*100:.1f}%)\n")

    # ---- 1. where the receptor is -----------------------------------------
    rows = []
    for pop, sub in obs.groupby("pop", observed=True):
        if len(sub) < MIN_CELLS:
            continue
        pos = sub["galr1"] > 0
        per_animal = sub.groupby("animal", observed=True).apply(
            lambda x: (x["galr1"] > 0).mean() * 100, include_groups=False)
        aged = np.mean([per_animal.get(a, np.nan) for a in ["F536", "G_073"]])
        adult = np.mean([per_animal.get(a, np.nan) for a in ["M493", "M399"]])
        rows.append({
            "population": pop,
            "nucleus": sub["nucleus_ext"].iloc[0],
            "cell_type": sub["cell_type"].iloc[0],
            "transmitter": sub["nt"].mode().iloc[0],
            "n_cells": len(sub),
            "cells_per_section": round(len(sub) / 12, 1),
            "pct_galr1_pos": round(pos.mean() * 100, 1),
            "mean_galr1_in_pos": round(sub.loc[pos, "galr1"].mean(), 2) if pos.any() else 0,
            # The drug-relevant quantity: how much of all Galr1 in the region
            # this population accounts for.
            "share_of_galr1_pct": round(sub["galr1"].sum() / galr1_total * 100, 2),
            "pct_galr3_pos": round((sub["galr3"] > 0).mean() * 100, 1),
            "pct_gal_pos": round((sub["gal"] > 0).mean() * 100, 1),
            "pct_galr1_and_galr3": round(((sub["galr1"] > 0) & (sub["galr3"] > 0)).mean() * 100, 1),
            "aged_pct_pos": round(aged, 1),
            "adult_pct_pos": round(adult, 1),
        })
    atlas = pd.DataFrame(rows).sort_values("share_of_galr1_pct", ascending=False)
    atlas.to_csv(OUT / "galr1_target_atlas.csv", index=False)

    print("=== Where M617 can act: populations ranked by share of total Galr1 ===")
    top = atlas.head(18)
    with pd.option_context("display.width", 250, "display.max_colwidth", 34):
        print(top[["population", "transmitter", "cells_per_section", "pct_galr1_pos",
                   "share_of_galr1_pct", "pct_galr3_pos", "aged_pct_pos",
                   "adult_pct_pos"]].to_string(index=False))
    print(f"\n  top 5 populations carry {atlas.head(5).share_of_galr1_pct.sum():.0f}% "
          f"of all Galr1 signal; top 10 carry {atlas.head(10).share_of_galr1_pct.sum():.0f}%")

    # ---- 2. what the drug will do -----------------------------------------
    print("\n=== Galr1 is Gi-coupled, so who carries it decides the net effect ===")
    nt = obs.groupby("nt", observed=True).apply(
        lambda x: pd.Series({
            "cells": len(x),
            "pct_galr1_pos": (x["galr1"] > 0).mean() * 100,
            "share_of_galr1": x["galr1"].sum() / galr1_total * 100,
        }), include_groups=False).round(1)
    print(nt.to_string())
    print("  An agonist inhibits the cells carrying the receptor, so signal on")
    print("  GABAergic neurons predicts net disinhibition of their targets.")

    # ---- 3. per nucleus ----------------------------------------------------
    print("\n=== Galr1 by nucleus (target density, and adult vs aged) ===")
    nuc_rows = []
    for nucleus, sub in obs.groupby("nucleus_ext", observed=True):
        per_animal = sub.groupby("animal", observed=True).apply(
            lambda x: (x["galr1"] > 0).mean() * 100, include_groups=False)
        nuc_rows.append({
            "nucleus": nucleus, "cells": len(sub),
            "galr1_pos_per_section": round((sub["galr1"] > 0).sum() / 12),
            "pct_galr1_pos": round((sub["galr1"] > 0).mean() * 100, 1),
            "share_of_galr1_pct": round(sub["galr1"].sum() / galr1_total * 100, 1),
            **{a: round(per_animal.get(a, np.nan), 1) for a in ANIMALS},
        })
    nuc = pd.DataFrame(nuc_rows).sort_values("share_of_galr1_pct", ascending=False)
    nuc.to_csv(OUT / "galr1_by_nucleus.csv", index=False)
    print(nuc.to_string(index=False))
    print("  (last four columns: % Galr1-positive per animal — F536/G_073 aged, "
          "M493/M399 adult)")

    # ---- 4. selectivity and endogenous tone --------------------------------
    print("\n=== Galr1/Galr3 overlap — where M617 selectivity matters least ===")
    overlap = atlas[atlas.pct_galr1_pos >= 10].nlargest(8, "pct_galr1_and_galr3")
    print(overlap[["population", "pct_galr1_pos", "pct_galr3_pos",
                   "pct_galr1_and_galr3"]].to_string(index=False))

    print("\n=== Endogenous Gal tone at the top targets ===")
    print(atlas.head(8)[["population", "pct_galr1_pos", "pct_gal_pos",
                         "share_of_galr1_pct"]].to_string(index=False))

    _plot(atlas, nuc, obs, galr1_total)
    print(f"\nWrote {OUT}")
    return 0


def _plot(atlas, nuc, obs, galr1_total) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(17, 6))
    nt_colour = {"GABAergic": "#B4531A", "glutamatergic": "#0E5A61",
                 "mixed": "#6A3D9A", "non-neuronal": "#9AA8A8"}

    ax = axes[0]
    top = atlas.head(16).iloc[::-1]
    ax.barh(range(len(top)), top["share_of_galr1_pct"],
            color=[nt_colour.get(t, "#999") for t in top["transmitter"]])
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels([p[:40] for p in top["population"]], fontsize=7)
    ax.set_xlabel("% of all Galr1 signal in the region")
    ax.set_title("Where M617 can act\n(colour = transmitter class)")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in nt_colour.values()]
    ax.legend(handles, nt_colour.keys(), fontsize=7, loc="lower right")

    ax = axes[1]
    sub = atlas[atlas.n_cells >= 150]
    ax.scatter(sub["pct_galr1_pos"], sub["cells_per_section"],
               s=sub["share_of_galr1_pct"] * 25 + 8,
               c=[nt_colour.get(t, "#999") for t in sub["transmitter"]], alpha=.75)
    for _, r in sub.nlargest(7, "share_of_galr1_pct").iterrows():
        ax.annotate(r["cell_type"][:22], (r["pct_galr1_pos"], r["cells_per_section"]),
                    fontsize=6.5, textcoords="offset points", xytext=(5, 4))
    ax.set_xlabel("% of cells Galr1-positive")
    ax.set_ylabel("cells per section")
    ax.set_yscale("log")
    ax.set_title("Target abundance vs receptor density\n(size = share of total Galr1)")

    ax = axes[2]
    order = nuc.sort_values("share_of_galr1_pct")
    y = np.arange(len(order))
    for i, animal in enumerate(ANIMALS):
        grp = ANIMAL_META[animal]["age_group"]
        ax.scatter(order[animal], y + (i - 1.5) * .16, s=34,
                   color="#B4531A" if grp == "aged" else "#0E5A61")
    ax.set_yticks(y); ax.set_yticklabels(order["nucleus"])
    ax.set_xlabel("% of cells Galr1-positive")
    ax.set_title("Galr1-positive fraction per nucleus\none point per animal "
                 "(orange aged, teal adult)")
    ax.grid(axis="x", alpha=.25)

    fig.tight_layout()
    fig.savefig(OUT / "m617_target_atlas.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
