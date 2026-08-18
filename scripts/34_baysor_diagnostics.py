#!/usr/bin/env python3
"""Two ways the Baysor comparison could be flattering itself, tested.

**1. Depth.**  Baysor cells carry 87 median transcripts against the vendor's 218.
A cell with fewer transcripts has fewer chances to show both a GABAergic and a
glutamatergic marker, so part of the drop in mixed identity could be arithmetic
rather than better segmentation.  The test: compare the two within bins of equal
transcript count.  If Baysor still wins at matched depth, the improvement is
real; if the curves lie on top of each other, the whole effect was depth.

**2. Fragments.**  Baysor returns 1.85 cells per DAPI nucleus.  Some of that is
expected -- a cell can cross the section plane while its nucleus does not, which
for a 10.7 um cell, a 7.12 um nucleus and a 5.54 um section predicts 1.28 -- but
1.85 is more than geometry accounts for.  The excess is either real cytoplasm
whose nucleus is out of plane, or fragments of neuropil.  Fragments would show up
as small, low-count, and biased toward the cell types with the most diffuse
processes, which would make the astrocyte gain in the composition an artefact.

Both matter for how much of the comparison to believe, so both are printed
whether or not they are flattering.
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

SEG = REPO / "data" / "segmentation"
PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "baysor"
GABA = ["Gad1", "Gad2", "Slc32a1"]
GLUT = ["Slc17a6", "Slc17a7"]
BINS = [10, 30, 60, 100, 150, 220, 320, 500, 10**6]


def mixed_by_depth(counts: pd.DataFrame) -> pd.DataFrame:
    tot = counts.sum(axis=1)
    gaba = counts[[g for g in GABA if g in counts]].sum(axis=1)
    glut = counts[[g for g in GLUT if g in counts]].sum(axis=1)
    neuronal = (gaba >= 2) | (glut >= 2)
    mixed = (gaba >= 2) & (glut >= 2)
    b = pd.cut(tot, BINS, right=False)
    return pd.DataFrame({
        "n_cells": tot.groupby(b, observed=True).size(),
        "n_marker_pos": neuronal.groupby(b, observed=True).sum(),
        "pct_mixed": (mixed.groupby(b, observed=True).sum()
                      / neuronal.groupby(b, observed=True).sum() * 100).round(1),
    })


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    ref = sc.read_h5ad(PROC / "hypothalamus_nuclei.h5ad")

    for section in sys.argv[1:]:
        print(f"\n{'='*76}\n{section}\n{'='*76}")
        bay = pd.read_csv(OUT / f"{section}_baysor_cells.csv", index_col=0)
        counts = pd.read_csv(SEG / f"{section}_baysor" / "segmentation_counts.tsv",
                             sep="\t", index_col=0).T
        counts = counts.loc[counts.index.isin(bay.index)]
        bay = bay.loc[counts.index]

        ven = ref[ref.obs["section"].astype(str) == section]
        ven_counts = pd.DataFrame(counts_matrix(ven), columns=ven.var_names)

        # --- 1. mixed identity at matched depth ---------------------------
        b = mixed_by_depth(counts).rename(columns=lambda c: "baysor_" + c)
        v = mixed_by_depth(ven_counts).rename(columns=lambda c: "vendor_" + c)
        tab = b.join(v, how="outer")
        tab = tab[(tab["baysor_n_marker_pos"] >= 100)
                  & (tab["vendor_n_marker_pos"] >= 100)]
        print("1. mixed GABA/glut rate within bins of equal transcript count")
        print(f"   {'transcripts/cell':>18s} {'vendor':>18s} {'Baysor':>18s}")
        for idx, r in tab.iterrows():
            print(f"   {str(idx):>18s} "
                  f"{r.vendor_pct_mixed:6.1f}% (n={int(r.vendor_n_marker_pos):5d}) "
                  f"{r.baysor_pct_mixed:6.1f}% (n={int(r.baysor_n_marker_pos):5d})")
        wins = int((tab["baysor_pct_mixed"] < tab["vendor_pct_mixed"]).sum())
        print(f"   Baysor lower in {wins} of {len(tab)} comparable bins")
        tab.to_csv(OUT / f"{section}_mixed_by_depth.csv")

        # --- 2. are the nucleus-free cells fragments? ---------------------
        has_nuc = bay["prior_cell"] > 0
        print("\n2. Baysor cells with and without a DAPI nucleus in the plane")
        print(f"   {'':22s} {'with nucleus':>14s} {'without':>14s}")
        for name, series in [("n cells", pd.Series(1, index=bay.index)),
                             ("median transcripts", counts.sum(axis=1)),
                             ("median area um2", bay["area"]),
                             ("median diameter um", 2 * np.sqrt(bay["area"] / np.pi))]:
            f = (lambda s: s.sum()) if name == "n cells" else (lambda s: s.median())
            print(f"   {name:22s} {f(series[has_nuc]):>14.1f} {f(series[~has_nuc]):>14.1f}")

        share = (bay.loc[~has_nuc, "cell_type"].value_counts(normalize=True) * 100)
        share_n = (bay.loc[has_nuc, "cell_type"].value_counts(normalize=True) * 100)
        comp = pd.DataFrame({"with_nucleus_pct": share_n.round(1),
                             "without_pct": share.round(1)}).fillna(0)
        comp["enrichment"] = (comp["without_pct"]
                              / comp["with_nucleus_pct"].clip(lower=.1)).round(2)
        print("\n   cell types over-represented among the nucleus-free cells:")
        for t, r in comp.nlargest(5, "enrichment").iterrows():
            print(f"     {t:34s} {r.with_nucleus_pct:5.1f}% -> {r.without_pct:5.1f}%"
                  f"   ({r.enrichment}x)")
        comp.to_csv(OUT / f"{section}_nucleus_free.csv")

        # --- 3. like for like: only the Baysor cells that have a nucleus --
        # Every vendor cell has a nucleus by construction, so this is the only
        # comparison of matched objects.  It is also matched in depth: these
        # cells carry a median transcript count close to the vendor's.
        def mixed_of(df):
            g = df[[x for x in GABA if x in df]].sum(axis=1)
            l = df[[x for x in GLUT if x in df]].sum(axis=1)
            n = (g >= 2) | (l >= 2)
            return ((g >= 2) & (l >= 2)).sum() / max(n.sum(), 1) * 100

        sub = counts[has_nuc.to_numpy()]
        print(f"\n3. like for like -- the {int(has_nuc.sum()):,} Baysor cells that sit on a"
              f" measured nucleus\n   (median {int(sub.sum(axis=1).median())} transcripts"
              f" against the vendor's {int(ven_counts.sum(axis=1).median())})")
        print(f"   mixed rate   vendor {mixed_of(ven_counts):.1f}%"
              f"   ->  Baysor {mixed_of(sub):.1f}%")
        print(f"   (using all Baysor cells instead gives {mixed_of(counts):.1f}%, but those"
              f" include\n    {int((~has_nuc).sum()):,} nucleus-free objects at a median of"
              f" {int(counts[~has_nuc.to_numpy()].sum(axis=1).median())} transcripts, and the"
              f" mixed rate\n    falls with depth, so that number flatters Baysor.)")

        pd.DataFrame([{
            "section": section,
            "vendor": round(mixed_of(ven_counts), 1),
            "baysor_matched": round(mixed_of(sub), 1),
            "baysor_all": round(mixed_of(counts), 1),
            "n_matched": int(has_nuc.sum()),
            "n_nucleus_free": int((~has_nuc).sum()),
            "median_counts_matched": int(sub.sum(axis=1).median()),
            "median_counts_free": int(counts[~has_nuc.to_numpy()].sum(axis=1).median()),
            "median_counts_vendor": int(ven_counts.sum(axis=1).median()),
        }]).to_csv(OUT / f"{section}_mixed_matched.csv", index=False)

        # Composition on the same matched set: does the astrocyte gain survive?
        from importlib.util import module_from_spec, spec_from_file_location
        spec = spec_from_file_location("cmp", REPO / "scripts" / "31_compare_segmentation.py")
        cmp = module_from_spec(spec); spec.loader.exec_module(cmp)
        prof = cmp.reference_profiles(ref)
        comp2 = pd.DataFrame({
            "vendor_pct": cmp.transfer(ven_counts, prof).value_counts(normalize=True) * 100,
            "baysor_all_pct": cmp.transfer(counts, prof).value_counts(normalize=True) * 100,
            "baysor_nuc_pct": cmp.transfer(sub, prof).value_counts(normalize=True) * 100,
        }).fillna(0).round(2)
        comp2["delta_matched"] = (comp2["baysor_nuc_pct"] - comp2["vendor_pct"]).round(2)
        comp2.to_csv(OUT / f"{section}_composition_matched.csv")
        print("\n   composition, largest shifts on the matched set:")
        print(f"   {'':34s} {'vendor':>8s} {'Baysor':>8s} {'Baysor':>8s}")
        print(f"   {'':34s} {'':>8s} {'(all)':>8s} {'(w/ nuc)':>8s}")
        for ct, r in comp2.reindex(
                comp2["delta_matched"].abs().sort_values(ascending=False).index).head(6).iterrows():
            print(f"     {ct:32s} {r.vendor_pct:7.2f}% {r.baysor_all_pct:7.2f}%"
                  f" {r.baysor_nuc_pct:7.2f}%  ({r.delta_matched:+.2f})")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
