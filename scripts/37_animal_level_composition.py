#!/usr/bin/env python3
"""Cell-population abundance with the animal as the unit of replication.

The animal-level counterpart of script 26, which screened on complete separation
of the six aged sections from the six adult sections and found nothing in the
evenly-sampled nuclei.  Same reasoning as script 36: sections from one mouse are
pseudoreplicates, and a rule one section can veto is measuring section variance.

Abundance is a share of its **own nucleus**, so ROI geometry cannot leak in.
Separation at the animal level means both aged animals fall one side of both
adults -- exact one-sided p = 1/6 = 0.167, the floor for 2 versus 2.

Results stay split by how evenly each nucleus is sampled between animals, because
that has not changed: ARC, VMH and ME/3V run 10-15% between-animal CV, while DMH,
LHA, ZI and DHA/PH run 30-49% because their cross-sectional area changes steeply
with rostro-caudal level.  In those, an abundance difference is confounded with
section plane before any biology is involved -- and unlike expression, which is
normalised within the population, abundance cannot be protected from it.

    python scripts/37_animal_level_composition.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.io import ANIMAL_META  # noqa: E402

PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "animal_level"
ANIMALS = ["F536", "G_073", "M493", "M399"]
AGED, ADULT = ["F536", "G_073"], ["M493", "M399"]
EVEN = ("ARC", "VMH", "ME_3V")
AP_LIMITED = ("DMH", "LHA", "ZI", "DHA_PH")
GLIA = ("Astrocyte", "Microglia", "Oligodendrocyte", "OPC", "Tanycyte", "Ependymal")
MIN_MEAN_PCT, MIN_CELLS = 0.5, 300


def screen(obs: pd.DataFrame, within: bool) -> pd.DataFrame:
    groups = ([(n, obs[obs.nucleus_ext == n]) for n in obs.nucleus_ext.unique()]
              if within else [("all nuclei", obs)])
    rows = []
    for nucleus, sub in groups:
        tab = pd.crosstab(sub["animal"], sub["cell_type"]).reindex(ANIMALS, fill_value=0)
        if (tab.sum(axis=1) < 200).any():
            continue
        pct = tab.div(tab.sum(axis=1), axis=0) * 100
        keep = (pct.mean() >= MIN_MEAN_PCT) & (tab.sum() >= MIN_CELLS)
        if not keep.any():
            continue
        pct = pct.loc[:, keep]
        a, d = pct.loc[AGED], pct.loc[ADULT]
        sep = (a.min() > d.max()) | (a.max() < d.min())
        lfc = np.log2((a.mean() + .01) / (d.mean() + .01))
        rows.append(pd.DataFrame({
            "nucleus": nucleus, "cell_type": pct.columns,
            "aged_pct": a.mean().round(2).to_numpy(),
            "adult_pct": d.mean().round(2).to_numpy(),
            "lfc": np.round(lfc.to_numpy(), 2),
            "separates": sep.to_numpy(),
            "n_cells": tab.loc[:, keep].sum().to_numpy(),
            **{f"pct_{an}": pct.loc[an].round(2).to_numpy() for an in ANIMALS},
        }))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(PROC / "hypothalamus_nuclei.h5ad")
    obs = adata.obs[adata.obs["nucleus_ext"].isin(EVEN + AP_LIMITED)].copy()
    obs["animal"] = obs["animal"].astype(str)
    obs["cell_type"] = obs["cell_type"].astype(str)
    obs = obs[~obs["cell_type"].str.startswith("unlabelled")]

    both = pd.concat([screen(obs, True), screen(obs, False)], ignore_index=True)
    both.to_csv(OUT / "composition_animal_level.csv", index=False)
    print(f"Tested {len(both)} populations\n")

    print("=== Calibration: how often does the criterion fire? ===")
    print("   (chance for exchangeable populations is 1/6 = 16.7%)")
    for tag, sel in (("evenly sampled (ARC, VMH, ME/3V)", both.nucleus.isin(EVEN)),
                     ("AP-limited (DMH, LHA, ZI, DHA/PH)", both.nucleus.isin(AP_LIMITED)),
                     ("region-wide", both.nucleus == "all nuclei")):
        sub = both[sel]
        if len(sub):
            print(f"  {tag:36s} {int(sub.separates.sum()):2d}/{len(sub):3d} "
                  f"({sub.separates.mean()*100:4.1f}%)")

    print("\n=== Separating populations in evenly-sampled nuclei and region-wide ===")
    hits = both[both.separates & both.nucleus.isin(list(EVEN) + ["all nuclei"])]
    hits = hits.reindex(hits["lfc"].abs().sort_values(ascending=False).index)
    cols = ["nucleus", "cell_type", "n_cells", "aged_pct", "adult_pct", "lfc",
            "pct_F536", "pct_G_073", "pct_M493", "pct_M399"]
    print(hits[cols].to_string(index=False) if len(hits) else "  none")

    print("\n=== Same, in AP-limited nuclei (abundance is confounded with section plane) ===")
    risky = both[both.separates & both.nucleus.isin(AP_LIMITED)]
    risky = risky.reindex(risky["lfc"].abs().sort_values(ascending=False).index)
    print(risky[cols].head(15).to_string(index=False) if len(risky) else "  none")

    print("\n=== Glia and immune populations, evenly-sampled nuclei ===")
    glia = both[both.cell_type.isin(GLIA) & both.nucleus.isin(list(EVEN) + ["all nuclei"])]
    print(glia.sort_values(["cell_type", "nucleus"])[
        ["nucleus", "cell_type", "n_cells", "aged_pct", "adult_pct", "lfc",
         "separates"]].to_string(index=False))
    print(f"\nWrote {OUT / 'composition_animal_level.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
