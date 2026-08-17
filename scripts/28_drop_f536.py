#!/usr/bin/env python3
"""G073 (aged) against M399 and M493 (adult), with F536 dropped.

Dropping F536 buys AP comparability: the remaining three animals overlap in
morphometric AP score, which no four-animal set does.  It costs the blocked
design and leaves one aged animal, so **any difference is confounded with
G073's individual identity** -- this cannot be an age test.

What makes it worth running anyway is that the two adults provide their own
control.  M399 versus M493 measures how much two same-age, AP-comparable animals
differ.  A gene only counts here if G073 departs from *both* adults in the same
direction *and* by more than the adults depart from each other.

Two readouts:
  1. animal level -- G073 vs each adult, against the adult-adult difference
  2. section level -- do G073's three sections separate from the six adult ones
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.io import SECTION_ANIMAL, counts_matrix  # noqa: E402

PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "drop_f536"
AP_FILE = REPO / "results" / "ap_matching" / "section_ap_scores.csv"
AGED, ADULTS = "G_073", ["M399", "M493"]
KEEP = [AGED] + ADULTS
NUCLEI = ("ARC", "ME_3V", "VMH", "DMH", "LHA", "ZI", "DHA_PH")
GENES = ["Gal", "Galr1", "Galr3"]
MIN_CELLS, MIN_PER_SECTION, MIN_COUNTS = 60, 15, 200


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    ap = pd.read_csv(AP_FILE, index_col=0)["ap_score"]
    kept = [s for s in ap.index if SECTION_ANIMAL[s] in KEEP]

    print("=== AP comparability after dropping F536 ===")
    for a in KEEP:
        v = ap[[s for s in kept if SECTION_ANIMAL[s] == a]]
        print(f"  {a:6s} {'aged ' if a == AGED else 'adult'}  AP {v.min():+.2f} .. {v.max():+.2f}")
    span = ap[kept].max() - ap[kept].min()
    print(f"  spread across the three animals {span:.2f}  "
          f"(was {ap.max() - ap.min():.2f} with F536)")

    adata = sc.read_h5ad(PROC / "hypothalamus_nuclei.h5ad")
    adata = adata[adata.obs["nucleus_ext"].isin(NUCLEI)
                  & adata.obs["animal"].astype(str).isin(KEEP)].copy()
    counts = counts_matrix(adata)
    pops = (adata.obs["nucleus_ext"].astype(str) + " | "
            + adata.obs["cell_type"].astype(str)).to_numpy()
    animals = adata.obs["animal"].astype(str).to_numpy()
    sections = adata.obs["section"].astype(str).to_numpy()
    genes = adata.var_names.to_numpy()

    rows = []
    for pop in pd.unique(pops):
        sel = pops == pop
        n = {a: int((sel & (animals == a)).sum()) for a in KEEP}
        if min(n.values()) < MIN_CELLS:
            continue
        mat = pd.DataFrame({a: counts[sel & (animals == a)].sum(axis=0) for a in KEEP},
                           index=genes).T
        expressed = mat.sum(axis=0) >= MIN_COUNTS
        cpm = np.log2(mat.div(mat.sum(axis=1), axis=0) * 1e6 + 1)

        d1 = cpm.loc[AGED] - cpm.loc["M399"]
        d2 = cpm.loc[AGED] - cpm.loc["M493"]
        adult_gap = (cpm.loc["M399"] - cpm.loc["M493"]).abs()
        # Must exceed the adult-adult difference in both comparisons, and agree
        # in sign, or it is within the range two same-age animals already show.
        beats = (np.sign(d1) == np.sign(d2)) & (d1.abs() > adult_gap) & (d2.abs() > adult_gap)

        # Section level: G073's three sections vs the six adult sections.
        ssel = [s for s in pd.unique(sections[sel]) if (sections[sel] == s).sum() >= MIN_PER_SECTION]
        sep = pd.Series(False, index=genes)
        if len(ssel) >= 9:
            smat = np.vstack([counts[sel & (sections == s)].sum(axis=0) for s in ssel])
            scpm = np.log2(smat / smat.sum(axis=1, keepdims=True) * 1e6 + 1)
            is_aged = np.array([SECTION_ANIMAL[s] == AGED for s in ssel])
            a_, d_ = scpm[is_aged], scpm[~is_aged]
            sep = pd.Series((a_.min(axis=0) > d_.max(axis=0))
                            | (a_.max(axis=0) < d_.min(axis=0)), index=genes)

        strict = beats & sep & expressed
        rate = float(strict.sum() / max(expressed.sum(), 1))
        for gene in GENES:
            if gene not in genes or not bool(expressed[gene]):
                continue
            rows.append({
                "population": pop, "gene": gene, "min_cells": min(n.values()),
                "vs_M399": round(float(d1[gene]), 2), "vs_M493": round(float(d2[gene]), 2),
                "adult_gap": round(float(adult_gap[gene]), 2),
                "beats_adult_gap": bool(beats[gene]),
                "sections_separate": bool(sep[gene]),
                "passes": bool(strict[gene]),
                "panel_pass_rate_pct": round(rate * 100, 1),
                "panel_tested": int(expressed.sum()),
            })

    res = pd.DataFrame(rows)
    res.to_csv(OUT / "drop_f536_results.csv", index=False)
    print(f"\n{res['population'].nunique()} populations with >={MIN_CELLS} cells "
          f"in all three animals\n")

    print("=== Gal / Galr1: G073 vs both adults, beating the adult-adult gap ===")
    hits = res[res["passes"]].sort_values("panel_pass_rate_pct")
    if len(hits):
        with pd.option_context("display.width", 230, "display.max_colwidth", 32):
            print(hits[["population", "gene", "min_cells", "vs_M399", "vs_M493",
                        "adult_gap", "panel_pass_rate_pct"]].to_string(index=False))
    else:
        print("  none")

    print("\n=== Counts ===")
    for gene in GENES:
        s = res[res.gene == gene]
        if len(s):
            print(f"  {gene:6s} {len(s):2d} populations | "
                  f"{int(s.beats_adult_gap.sum()):2d} beat the adult gap | "
                  f"{int(s.sections_separate.sum()):2d} separate at section level | "
                  f"{int(s.passes.sum()):2d} both  (panel rate "
                  f"{s.panel_pass_rate_pct.mean():.1f}%)")

    print("\n=== The AgRP populations specifically ===")
    agrp = res[res.population.str.contains("Agrp") & (res.gene == "Gal")]
    if len(agrp):
        print(agrp[["population", "min_cells", "vs_M399", "vs_M493", "adult_gap",
                    "beats_adult_gap", "sections_separate"]].to_string(index=False))
    else:
        print("  below the cell threshold in this subset")
    print(f"\nWrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
