#!/usr/bin/env python3
"""Is any cell population changing in abundance with age?

Composition, not expression.  Same discipline as the expression sweep:

* abundance is a share of its **own nucleus**, so ROI geometry cannot leak in
* the criterion is **complete separation of the six aged from the six adult
  sections**, consistent inside both blocks -- no model, no covariate
* results are split by whether the nucleus is evenly sampled between animals.
  ARC, VMH and ME/3V have 10-15% between-animal CV; DMH, LHA, ZI and DHA/PH
  have 30-49%, because their cross-sectional area changes steeply with
  rostro-caudal level.  A composition result in those is confounded with
  section plane before any biology is involved.
* calibration is the rate across all populations tested: if a third of them
  separate, the criterion is measuring anatomy rather than ageing.

Glial and immune populations are also reported explicitly, since those are the
a priori expected ageing changes and deserve to be seen whether or not they
clear the bar.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.io import ANIMAL_META, SECTION_ANIMAL  # noqa: E402

PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "composition"
BLOCKS = {"B1": ("F536", "M493"), "B2": ("G_073", "M399")}
EVEN = ("ARC", "VMH", "ME_3V")        # between-animal CV 10-15%
AP_LIMITED = ("DMH", "LHA", "ZI", "DHA_PH")   # CV 30-49%
GLIA = ("Astrocyte", "Microglia", "Oligodendrocyte", "OPC", "Tanycyte", "Ependymal")
MIN_MEAN_PCT = 0.5      # ignore populations under 0.5% of their nucleus
MIN_CELLS = 300


def separation(pct: pd.DataFrame) -> pd.DataFrame:
    """pct: sections x populations, % of the section's nucleus."""
    animals = pd.Series([SECTION_ANIMAL[s] for s in pct.index], index=pct.index)
    grp = animals.map(lambda a: ANIMAL_META[a]["age_group"])
    aged, adult = pct[grp == "aged"], pct[grp == "adult"]
    up = aged.min() > adult.max()
    dn = aged.max() < adult.min()

    block_ok = pd.Series(True, index=pct.columns)
    for a, d in BLOCKS.values():
        va, vd = pct[animals == a], pct[animals == d]
        block_ok &= ((va.min() > vd.max()) & up) | ((va.max() < vd.min()) & dn)

    lfc = np.mean([np.log2((pct[animals == a].mean() + .01)
                           / (pct[animals == d].mean() + .01))
                   for a, d in BLOCKS.values()], axis=0)
    return pd.DataFrame({
        "aged_pct": aged.mean().round(2), "adult_pct": adult.mean().round(2),
        "lfc": np.round(lfc, 2), "separates": (up | dn) & block_ok,
        "mean_pct": pct.mean().round(2),
    })


def run(obs: pd.DataFrame, within: str | None) -> pd.DataFrame:
    """Composition within each nucleus, or across the whole region."""
    rows = []
    groups = [(n, obs[obs.nucleus_ext == n]) for n in obs.nucleus_ext.unique()] \
        if within == "nucleus" else [("all nuclei", obs)]
    for nucleus, sub in groups:
        tab = pd.crosstab(sub["section"], sub["cell_type"])
        tab = tab.reindex(sorted(SECTION_ANIMAL), fill_value=0)
        if (tab.sum(axis=1) < 200).any():
            continue
        pct = tab.div(tab.sum(axis=1), axis=0) * 100
        keep = (pct.mean() >= MIN_MEAN_PCT) & (tab.sum() >= MIN_CELLS)
        if not keep.any():
            continue
        res = separation(pct.loc[:, keep])
        res["nucleus"] = nucleus
        res["n_cells"] = tab.loc[:, keep].sum()
        rows.append(res.reset_index().rename(columns={"index": "cell_type"}))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(PROC / "hypothalamus_nuclei.h5ad")
    obs = adata.obs[adata.obs["nucleus_ext"].isin(EVEN + AP_LIMITED)].copy()
    obs["section"] = obs["section"].astype(str)
    obs["cell_type"] = obs["cell_type"].astype(str)
    obs = obs[~obs["cell_type"].str.startswith("unlabelled")]

    per_nuc = run(obs, "nucleus")
    region = run(obs, None)
    both = pd.concat([per_nuc, region], ignore_index=True)
    both.to_csv(OUT / "composition_age.csv", index=False)

    print(f"Tested {len(both)} populations "
          f"({per_nuc.nucleus.nunique()} nuclei + region-wide)\n")

    print("=== Calibration: how often does the criterion fire? ===")
    for tag, sel in (("evenly sampled (ARC, VMH, ME/3V)", both.nucleus.isin(EVEN)),
                     ("AP-limited (DMH, LHA, ZI, DHA/PH)", both.nucleus.isin(AP_LIMITED)),
                     ("region-wide", both.nucleus == "all nuclei")):
        sub = both[sel]
        if not len(sub):
            continue
        print(f"  {tag:36s} {int(sub.separates.sum()):2d}/{len(sub):3d} "
              f"({sub.separates.mean()*100:4.1f}%)")

    print("\n=== Populations that fully separate, in evenly-sampled nuclei ===")
    hits = both[both.separates & both.nucleus.isin(list(EVEN) + ["all nuclei"])]
    hits = hits.sort_values("lfc", key=abs, ascending=False)
    if len(hits):
        print(hits[["nucleus", "cell_type", "n_cells", "aged_pct",
                    "adult_pct", "lfc"]].to_string(index=False))
    else:
        print("  none")

    print("\n=== Same, in AP-limited nuclei (interpret with caution) ===")
    risky = both[both.separates & both.nucleus.isin(AP_LIMITED)]
    if len(risky):
        print(risky.sort_values("lfc", key=abs, ascending=False)[
            ["nucleus", "cell_type", "n_cells", "aged_pct", "adult_pct",
             "lfc"]].to_string(index=False))
    else:
        print("  none")

    print("\n=== Glia and immune populations, whether or not they separate ===")
    glia = both[both.cell_type.isin(GLIA) & both.nucleus.isin(list(EVEN) + ["all nuclei"])]
    print(glia.sort_values(["cell_type", "nucleus"])[
        ["nucleus", "cell_type", "n_cells", "aged_pct", "adult_pct",
         "lfc", "separates"]].to_string(index=False))
    print(f"\nWrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
