#!/usr/bin/env python3
"""Rank the animal-level hits by whether rostro-caudal position explains them.

Script 36 recovers 357 strong hits the section-level rule had discarded, and
script 37 finds that half to two thirds of abundance comparisons separate.  Both
need a way to rank, because separation alone fires at 41% and 61% against a 16.7%
chance rate.

The confound worth ranking on is real and measurable: the aged sections average
+1.19 on the morphometric AP score and the adult sections -1.19, a 2.38-unit
offset.  Anything that varies with rostro-caudal level will therefore look like
an age effect.

**The earlier AP test failed audit and is not reused here.**  It compared
sections across animals, where the AP gap is perfectly confounded with block
(AUC 1.00), and the regression correction it motivated extrapolated 5.4x.  This
is a different test: each animal contributes three sections spanning its own AP
range, so the slope is fitted **within** an animal, where age, batch and animal
identity are all constant by construction.  Nothing is corrected or normalised --
the within-animal slope is reported next to the observed age difference so the
two can be compared directly.

For each hit:

  slope        change per AP unit, averaged over the four within-animal fits
  sign_agree   how many of the four animals agree on that direction (0-4)
  ap_predicted slope x 2.38, the difference AP alone would produce
  ratio        observed age difference / ap_predicted

A hit with |ratio| well above 1 is larger than AP can account for.  A hit with
|ratio| near or below 1 and four animals agreeing on the slope is the pattern a
sampling difference produces.  Sign agreement of 2/4 means the within-animal
slopes disagree, so AP is not driving that population and the ratio is noise.

    python scripts/38_within_animal_ap.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.io import ANIMAL_META, SECTION_ANIMAL, counts_matrix  # noqa: E402

PROC = REPO / "data" / "processed"
RES = REPO / "results" / "animal_level"
AP_FILE = REPO / "results" / "ap_matching" / "section_ap_scores.csv"
ANIMALS = ["F536", "G_073", "M493", "M399"]
AGED, ADULT = ["F536", "G_073"], ["M493", "M399"]
MIN_PER_SECTION = 15
TOP_N = 30


def within_animal_slope(values: pd.Series, ap: pd.Series) -> tuple[float, int, list]:
    """Mean slope of value on AP, fitted separately inside each animal."""
    slopes = []
    for animal in ANIMALS:
        secs = [s for s in values.index if SECTION_ANIMAL[s] == animal]
        if len(secs) < 3:
            continue
        x, y = ap[secs].to_numpy(), values[secs].to_numpy()
        if np.ptp(x) < 1e-6:
            continue
        slopes.append(float(np.polyfit(x, y, 1)[0]))
    if not slopes:
        return np.nan, 0, []
    mean = float(np.mean(slopes))
    agree = int(max(sum(s > 0 for s in slopes), sum(s < 0 for s in slopes)))
    return mean, agree, [round(s, 3) for s in slopes]


def main() -> int:
    ap = pd.read_csv(AP_FILE, index_col=0)["ap_score"]
    ap_gap = (ap[[s for s in ap.index if SECTION_ANIMAL[s] in AGED]].mean()
              - ap[[s for s in ap.index if SECTION_ANIMAL[s] in ADULT]].mean())
    print(f"Aged sections average AP {ap[[s for s in ap.index if SECTION_ANIMAL[s] in AGED]].mean():+.2f}, "
          f"adult {ap[[s for s in ap.index if SECTION_ANIMAL[s] in ADULT]].mean():+.2f}  "
          f"→ offset {ap_gap:+.2f} AP units\n")

    adata = sc.read_h5ad(PROC / "hypothalamus_nuclei.h5ad")
    counts = counts_matrix(adata)
    genes = adata.var_names.to_numpy()
    sections = adata.obs["section"].astype(str).to_numpy()
    pops = (adata.obs["nucleus_ext"].astype(str) + " | "
            + adata.obs["cell_type"].astype(str)).to_numpy()

    sweep = pd.read_csv(RES / "animal_level_sweep.csv")
    strong = sweep[sweep.animal_separates & (sweep.lfc_pct_rank >= 95)]
    strong = strong.reindex(strong["lfc"].abs().sort_values(ascending=False).index)
    strong = strong.head(TOP_N)

    rows = []
    for _, h in strong.iterrows():
        sel = pops == h.population
        gi = int(np.where(genes == h.gene)[0][0])
        vals = {}
        for s in pd.unique(sections[sel]):
            m = sel & (sections == s)
            if m.sum() < MIN_PER_SECTION:
                continue
            c = counts[m]
            vals[s] = float(np.log2(c[:, gi].sum() / c.sum() * 1e6 + 1))
        v = pd.Series(vals)
        if len(v) < 12:
            continue
        slope, agree, per = within_animal_slope(v, ap)
        pred = slope * ap_gap
        rows.append({
            "population": h.population, "gene": h.gene, "observed_lfc": h.lfc,
            "slope_per_ap": round(slope, 3), "sign_agree": f"{agree}/4",
            "ap_predicted": round(pred, 2),
            "ratio": round(abs(h.lfc) / abs(pred), 2) if abs(pred) > 1e-6 else np.inf,
            "verdict": ("AP explains it" if agree >= 4 and abs(pred) >= abs(h.lfc) * .7
                        else "larger than AP" if abs(pred) < abs(h.lfc) * .5
                        else "partly AP"),
        })

    out = pd.DataFrame(rows)
    out.to_csv(RES / "within_animal_ap_check.csv", index=False)
    with pd.option_context("display.width", 220, "display.max_colwidth", 30):
        print("=== Top animal-level expression hits, checked within animals ===")
        print(out.to_string(index=False))

    print("\n=== Summary ===")
    for v, n in out.verdict.value_counts().items():
        print(f"  {v:18s} {n:3d}")

    print("\n=== The galanin system, same check ===")
    gal = sweep[sweep.gene.isin(["Gal", "Galr1", "Galr3"]) & sweep.animal_separates]
    gal = gal.reindex(gal["lfc"].abs().sort_values(ascending=False).index).head(12)
    rows = []
    for _, h in gal.iterrows():
        sel = pops == h.population
        gi = int(np.where(genes == h.gene)[0][0])
        vals = {}
        for s in pd.unique(sections[sel]):
            m = sel & (sections == s)
            if m.sum() < MIN_PER_SECTION:
                continue
            c = counts[m]
            vals[s] = float(np.log2(c[:, gi].sum() / c.sum() * 1e6 + 1))
        v = pd.Series(vals)
        if len(v) < 12:
            continue
        slope, agree, per = within_animal_slope(v, ap)
        pred = slope * ap_gap
        rows.append({"population": h.population, "gene": h.gene,
                     "observed_lfc": h.lfc, "slope_per_ap": round(slope, 3),
                     "sign_agree": f"{agree}/4", "ap_predicted": round(pred, 2),
                     "ratio": round(abs(h.lfc) / abs(pred), 2) if abs(pred) > 1e-6 else np.inf})
    galout = pd.DataFrame(rows)
    galout.to_csv(RES / "within_animal_ap_galanin.csv", index=False)
    with pd.option_context("display.width", 220, "display.max_colwidth", 30):
        print(galout.to_string(index=False))
    print(f"\nWrote {RES}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
