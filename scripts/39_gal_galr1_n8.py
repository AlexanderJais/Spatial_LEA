#!/usr/bin/env python3
"""Cell types where Gal or Galr1 changes with age, across all eight animals.

Four aged against four adult, in four blocks.  Only the galanin ligand and its
two panel receptors are reported; every other gene appears solely as calibration.

**Everything is computed within blocks.**  Each block is one aged and one adult
animal run on the same slide, panel and segmentation chemistry, so a within-block
fold change never crosses a technical boundary.  This matters because the two
cohorts used different panels (315 of 325 targets shared) and different
segmentation chemistry, and both are perfectly aligned with sex -- so raw values
are not comparable between cohorts, while within-block fold changes are.

Two statistics, reported side by side because they fail differently:

  blocked    sign of the fold change in each of the four blocks.  All four
             agreeing is exact one-sided p = 1/16 = 0.0625 under sign-flipping,
             and it is immune to the panel and chemistry difference because no
             comparison ever leaves a block.
  separation both aged animals... all four aged animals on one side of all four
             adults.  Exact one-sided p = 1/70 = 0.014 over the C(8,4) label
             assignments -- but it compares animals across cohorts, so it is the
             more powerful and the more fragile of the two.

A gene clearing both has p = 0.014 by the stronger test and 0.0625 by the test
that cannot be blamed on the panel.

Calibration is the same pair of tests applied to all shared panel genes in the
same population, so every hit carries the local false-positive rate.

    python scripts/39_gal_galr1_n8.py
    python scripts/39_gal_galr1_n8.py --min-cells 25
"""

from __future__ import annotations

import argparse
import sys
from itertools import combinations, product
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.io import AGED, ADULT, ANIMAL_META, BLOCKS, counts_matrix  # noqa: E402

PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "gal_n8"
GENES = ["Gal", "Galr1", "Galr3"]
ANIMALS = list(ANIMAL_META)
MIN_COUNTS = 200
MIN_MARGIN = 0.25


def blocked_stats(cpm: pd.DataFrame) -> pd.DataFrame:
    """cpm: animals x genes, log2 CPM.  Fold change inside each block."""
    lfc = pd.DataFrame({b: cpm.loc[a] - cpm.loc[d] for b, (a, d) in BLOCKS.items()}).T
    mean_lfc = lfc.mean()
    agree = (np.sign(lfc) == np.sign(mean_lfc)).sum()

    # Exact sign-flip permutation within blocks: 2^4 = 16 arrangements.
    flips = np.array(list(product([1, -1], repeat=len(BLOCKS))))
    null = flips @ lfc.to_numpy() / len(BLOCKS)
    p_blocked = (np.abs(null) >= np.abs(mean_lfc.to_numpy())).mean(axis=0)

    aged_v = cpm.loc[list(AGED)].to_numpy()
    adult_v = cpm.loc[list(ADULT)].to_numpy()
    margin = np.maximum(aged_v.min(axis=0) - adult_v.max(axis=0),
                        adult_v.min(axis=0) - aged_v.max(axis=0))

    # Exact permutation over the C(8,4) = 70 ways of splitting eight animals.
    vals = cpm.to_numpy()
    obs = aged_v.mean(axis=0) - adult_v.mean(axis=0)
    idx = list(range(len(ANIMALS)))
    null2 = np.array([vals[list(c)].mean(axis=0)
                      - vals[[i for i in idx if i not in c]].mean(axis=0)
                      for c in combinations(idx, 4)])
    p_exact = (np.abs(null2) >= np.abs(obs)).mean(axis=0)

    return pd.DataFrame({
        "mean_lfc": mean_lfc, "blocks_agree": agree, "p_blocked": p_blocked,
        "margin": margin, "separates": margin > 0, "p_exact": p_exact,
        **{f"lfc_{b}": lfc.loc[b] for b in BLOCKS},
        **{f"cpm_{a}": cpm.loc[a] for a in ANIMALS},
    })


def screen(adata, key: str, min_cells: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    counts = counts_matrix(adata)
    genes = adata.var_names.to_numpy()
    animals = adata.obs["animal"].astype(str).to_numpy()
    labels = adata.obs[key].astype(str).to_numpy() if key != "cell_type_only" \
        else adata.obs["cell_type"].astype(str).to_numpy()

    rows, cal = [], []
    for pop in pd.unique(labels):
        if "unlabelled" in pop:
            continue
        sel = labels == pop
        n = {a: int((sel & (animals == a)).sum()) for a in ANIMALS}
        if min(n.values()) < min_cells:
            continue
        mat = pd.DataFrame({a: counts[sel & (animals == a)].sum(axis=0) for a in ANIMALS},
                           index=genes).T
        expressed = mat.sum(axis=0) >= MIN_COUNTS
        if expressed.sum() < 30:
            continue
        cpm = np.log2(mat.div(mat.sum(axis=1), axis=0) * 1e6 + 1)
        res = blocked_stats(cpm)
        # Two bars.  The blocked one is primary: it never compares across
        # cohorts, so it cannot be produced by the panel or chemistry
        # difference.  Complete separation is stronger but does compare across
        # them, so it is reported as a bonus rather than a requirement.
        consistent = res.blocks_agree == 4
        strict = consistent & res.separates & (res.margin >= MIN_MARGIN)

        cal.append({"population": pop, "min_cells": min(n.values()),
                    "genes_tested": int(expressed.sum()),
                    "pct_blocks_agree_4": round(float((res.blocks_agree[expressed] == 4).mean() * 100), 1),
                    "pct_separating": round(float(res.separates[expressed].mean() * 100), 1),
                    "pct_strict": round(float(strict[expressed].mean() * 100), 1)})
        rank = res.loc[expressed, "mean_lfc"].abs().rank(pct=True) * 100
        for g in GENES:
            if g not in res.index or not expressed[g]:
                continue
            r = res.loc[g]
            rows.append({
                "population": pop, "gene": g, "min_cells": min(n.values()),
                "mean_lfc": round(float(r.mean_lfc), 3),
                "fold": round(float(2 ** r.mean_lfc), 2),
                "blocks_agree": f"{int(r.blocks_agree)}/4",
                "p_blocked": round(float(r.p_blocked), 4),
                "separates": bool(r.separates), "margin": round(float(r.margin), 3),
                "p_exact": round(float(r.p_exact), 4),
                "lfc_pct_rank": round(float(rank.get(g, np.nan)), 1),
                "all_blocks_agree": bool(consistent[g]),
                "passes_strict": bool(strict[g]),
                "panel_strict_pct": cal[-1]["pct_strict"],
                **{f"lfc_{b}": round(float(r[f"lfc_{b}"]), 2) for b in BLOCKS},
            })
    return pd.DataFrame(rows), pd.DataFrame(cal)


def report(res: pd.DataFrame, cal: pd.DataFrame, title: str) -> None:
    print(f"\n{'='*94}\n{title}\n{'='*94}")
    print(f"{cal.shape[0]} populations tested\n")
    print("Calibration across all shared panel genes:")
    print(f"  all 4 blocks agree      {cal.pct_blocks_agree_4.mean():5.1f}% of genes  "
          "(chance 2/16 = 12.5%)")
    print(f"  aged separate from adult{cal.pct_separating.mean():5.1f}% of genes  "
          "(chance 2/70 = 2.9%)")
    print(f"  both, margin >= 0.25    {cal.pct_strict.mean():5.1f}% of genes  "
          "<- the rate any hit has to beat")

    cols = ["population", "min_cells", "mean_lfc", "fold", "p_blocked", "p_exact",
            "lfc_pct_rank", "separates", "lfc_B1", "lfc_B2", "lfc_B3", "lfc_B4"]
    for gene in GENES:
        sub = res[res.gene == gene]
        if not len(sub):
            continue
        hits = sub[sub.all_blocks_agree]
        hits = hits.reindex(hits["mean_lfc"].abs().sort_values(ascending=False).index)
        print(f"\n--- {gene}: {len(hits)} of {len(sub)} populations change in the same "
              f"direction in all four blocks ---")
        if not len(hits):
            print("    none")
        else:
            with pd.option_context("display.width", 230, "display.max_colwidth", 34):
                print(hits[cols].to_string(index=False))
            n_strict = int(hits.passes_strict.sum())
            print(f"    of these, {n_strict} also separate the four aged animals "
                  "completely from the four adults")
        near = sub[~sub.all_blocks_agree & (sub.blocks_agree == "3/4")]
        near = near.reindex(near["mean_lfc"].abs().sort_values(ascending=False).index)
        if len(near):
            print(f"    3 of 4 blocks agreeing ({len(near)}), largest first:")
            with pd.option_context("display.width", 230, "display.max_colwidth", 34):
                print(near[cols].head(5).to_string(index=False))


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--min-cells", type=int, default=30)
    args = p.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(PROC / "hypothalamus_nuclei.h5ad")
    adata.obs["pop"] = (adata.obs["nucleus_ext"].astype(str) + " | "
                        + adata.obs["cell_type"].astype(str))

    res_ct, cal_ct = screen(adata, "cell_type_only", args.min_cells)
    report(res_ct, cal_ct, "BY CELL TYPE  (pooled across the whole hypothalamic window)")
    res_ct.to_csv(OUT / "gal_galr1_by_cell_type.csv", index=False)
    cal_ct.to_csv(OUT / "calibration_by_cell_type.csv", index=False)

    keep = adata.obs["nucleus_ext"].isin(
        ["ARC", "ME_3V", "VMH", "DMH", "LHA", "ZI", "DHA_PH"]).to_numpy()
    res_np, cal_np = screen(adata[keep], "pop", args.min_cells)
    report(res_np, cal_np, "BY NUCLEUS x CELL TYPE")
    res_np.to_csv(OUT / "gal_galr1_by_nucleus_cell_type.csv", index=False)
    cal_np.to_csv(OUT / "calibration_by_nucleus_cell_type.csv", index=False)

    print(f"\nWrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
