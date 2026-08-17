#!/usr/bin/env python3
"""Gal / Galr1 analysis on the single best AP-matched section per mouse.

Selection: for each animal keep the one section that minimises the spread of the
morphometric AP score across the four chosen sections; discard the other two.
All 3^4 = 81 combinations are enumerated, so the chosen set is optimal, not
merely reasonable.

With one section per animal the section-level separation criterion is gone --
there are four numbers, one per mouse. What remains is the blocked fold change
in each of the two batches, its consistency, and the panel-wide rate of genes
clearing the same bar in the same population, which is the empirical
false-positive rate for this design.
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.io import ANIMAL_META, SECTION_ANIMAL, counts_matrix  # noqa: E402

PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "matched_gal"
AP_FILE = REPO / "results" / "ap_matching" / "section_ap_scores.csv"
BLOCKS = {"B1": ("F536", "M493"), "B2": ("G_073", "M399")}
ANIMALS = ["F536", "G_073", "M493", "M399"]
NUCLEI = ("ARC", "ME_3V", "VMH", "DMH", "LHA", "ZI", "DHA_PH")
GENES = ["Gal", "Galr1", "Galr3"]
MIN_CELLS = 40
MIN_COUNTS = 150
LFC_THRESHOLD = 0.58


def choose(ap: pd.Series) -> tuple[list[str], float]:
    by_animal = {a: [s for s in ap.index if SECTION_ANIMAL[s] == a] for a in ANIMALS}
    best = None
    for combo in itertools.product(*[by_animal[a] for a in ANIMALS]):
        spread = float(ap[list(combo)].max() - ap[list(combo)].min())
        if best is None or spread < best[0]:
            best = (spread, list(combo))
    return best[1], best[0]


def population_table(adata, sections) -> dict:
    """Per-animal summed counts per gene, for each population."""
    keep = adata.obs["section"].astype(str).isin(sections).to_numpy()
    sub = adata[keep]
    counts = counts_matrix(sub)
    pops = (sub.obs["nucleus_ext"].astype(str) + " | "
            + sub.obs["cell_type"].astype(str)).to_numpy()
    animals = sub.obs["animal"].astype(str).to_numpy()
    out = {}
    for pop in pd.unique(pops):
        sel = pops == pop
        if not set(animals[sel]) == set(ANIMALS):
            continue
        mat, n = {}, {}
        for a in ANIMALS:
            m = sel & (animals == a)
            n[a] = int(m.sum())
            mat[a] = counts[m].sum(axis=0)
        if min(n.values()) < MIN_CELLS:
            continue
        out[pop] = (pd.DataFrame(mat, index=adata.var_names).T, n)
    return out


def blocked(mat: pd.DataFrame) -> pd.DataFrame:
    cpm = mat.div(mat.sum(axis=1), axis=0) * 1e6
    lfc = {b: np.log2((cpm.loc[a] + 1) / (cpm.loc[d] + 1)) for b, (a, d) in BLOCKS.items()}
    b1, b2 = lfc["B1"].to_numpy(), lfc["B2"].to_numpy()
    return pd.DataFrame({
        "lfc_B1": b1, "lfc_B2": b2, "mean_lfc": (b1 + b2) / 2,
        "consistent": np.sign(b1) == np.sign(b2),
        "min_abs": np.minimum(np.abs(b1), np.abs(b2)),
        "counts": mat.sum(axis=0).to_numpy(),
    }, index=mat.columns)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    ap = pd.read_csv(AP_FILE, index_col=0)["ap_score"]
    chosen, spread = choose(ap)

    print("=== Selected sections (1 per mouse, AP spread minimised over 81 sets) ===")
    for s in sorted(chosen, key=lambda x: ap[x]):
        a = SECTION_ANIMAL[s]
        print(f"  {s:8s} {a:6s} {ANIMAL_META[a]['age_group']:6s} AP {ap[s]:+.2f}")
    dropped = [s for s in ap.index if s not in chosen]
    print(f"  dropped: {', '.join(sorted(dropped))}")
    print(f"  AP spread of the chosen set {spread:.2f}  (all 12 sections span "
          f"{ap.max() - ap.min():.2f})")
    nearest = sorted(ap[[s for s in chosen if SECTION_ANIMAL[s] != 'F536']])
    print(f"  F536_2 sits {ap['F536_2'] - nearest[-1]:+.2f} from the next closest "
          f"section, so the set is matched only as well as F536 allows.")

    adata = sc.read_h5ad(PROC / "hypothalamus_nuclei.h5ad")
    adata = adata[adata.obs["nucleus_ext"].isin(NUCLEI)].copy()
    pops = population_table(adata, chosen)
    print(f"\n{len(pops)} populations with >={MIN_CELLS} cells in every mouse\n")

    rows = []
    for pop, (mat, n) in pops.items():
        res = blocked(mat)
        expressed = res["counts"] >= MIN_COUNTS
        passes = res["consistent"] & (res["min_abs"] >= LFC_THRESHOLD) & expressed
        rate = float(passes.sum() / max(expressed.sum(), 1))
        for gene in GENES:
            if gene not in res.index or not expressed[gene]:
                continue
            r = res.loc[gene]
            rows.append({
                "population": pop, "gene": gene, "min_cells": min(n.values()),
                "lfc_B1": round(r.lfc_B1, 2), "lfc_B2": round(r.lfc_B2, 2),
                "mean_lfc": round(r.mean_lfc, 2),
                "consistent": bool(r.consistent),
                "passes": bool(passes[gene]),
                "panel_pass_rate_pct": round(rate * 100, 1),
                "panel_genes_tested": int(expressed.sum()),
            })
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "matched_gal_results.csv", index=False)

    print("=== Gal / Galr1 clearing consistency + >=1.5x, on the matched sections ===")
    hits = res[res["passes"]].sort_values("panel_pass_rate_pct")
    if len(hits):
        with pd.option_context("display.width", 220, "display.max_colwidth", 34):
            print(hits[["population", "gene", "min_cells", "lfc_B1", "lfc_B2",
                        "mean_lfc", "panel_pass_rate_pct"]].to_string(index=False))
        print("\n  panel_pass_rate_pct = % of expressed panel genes clearing the same")
        print("  bar in that population, i.e. the false-positive rate to beat.")
    else:
        print("  none")

    print(f"\n=== Overall ===")
    for gene in GENES:
        sub = res[res.gene == gene]
        if len(sub):
            print(f"  {gene:6s} tested in {len(sub):2d} populations, "
                  f"{int(sub.consistent.sum()):2d} consistent, "
                  f"{int(sub.passes.sum()):2d} clearing the bar "
                  f"(panel rate {sub.panel_pass_rate_pct.mean():.0f}%)")

    print("\n=== Does matching change the Gal/AgRP result? ===")
    for pop in ("ARC | ARC Agrp/Npy", "ME_3V | ARC Agrp/Npy"):
        r = res[(res.population == pop) & (res.gene == "Gal")]
        if len(r):
            r = r.iloc[0]
            print(f"  {pop:22s} matched  B1 {r.lfc_B1:+.2f}  B2 {r.lfc_B2:+.2f}  "
                  f"mean {r.mean_lfc:+.2f}  passes={r.passes}")
    print("  (all 12 sections gave ARC +0.79/+0.66 and ME_3V +1.18/+1.23)")
    print(f"\nWrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
