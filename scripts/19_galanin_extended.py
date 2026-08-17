#!/usr/bin/env python3
"""Galanin system across the extended parcellation, with the AP test applied.

ARC, VMH, DMH, LHA, ZI, DHA/PH and the ME/3V wall.  LHA and ZI were outside the
previous MBH definition entirely, so they are new territory.

Every candidate is put through the test that separated signal from artifact
earlier: correlate the effect size against the rostro-caudal mismatch of the
section pair it was computed from.  An effect that grows with mismatch is
anatomy, not age.

Sampling balance is reported per nucleus first, because it decides which nuclei
can support a cross-animal comparison at all: ARC and VMH are sampled evenly
(CV ~10%), while DMH, LHA, ZI and DHA vary 30-50% between animals.
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats as st

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.io import ANIMAL_META, SECTION_ANIMAL, counts_matrix  # noqa: E402

PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "galanin_extended"
AP_FILE = REPO / "results" / "ap_matching" / "section_ap_scores.csv"
GENES = ["Gal", "Galr1", "Galr3"]
NUCLEI = ("ARC", "ME_3V", "VMH", "DMH", "LHA", "ZI", "DHA_PH")
BLOCKS = {"B1": ("F536", "M493"), "B2": ("G_073", "M399")}
MIN_CELLS, MIN_POS, MIN_SECTION = 40, 15, 15
LFC_THRESHOLD = 0.58


def per_animal(adata, gi: int, mask: np.ndarray) -> pd.DataFrame:
    counts = counts_matrix(adata[mask])
    obs = adata.obs[mask]
    df = pd.DataFrame({"animal": obs["animal"].astype(str).to_numpy(),
                       "g": counts[:, gi], "tot": counts.sum(axis=1)})
    g = df.groupby("animal")
    # Aggregated with named columns rather than groupby.apply, which pandas
    # deprecates and which emitted a warning per call across ~500 calls here.
    agg = g.agg(gg=("g", "sum"), tt=("tot", "sum"), n=("g", "size"),
                n_pos=("g", lambda x: int((x > 0).sum())))
    return pd.DataFrame({"cpm": agg["gg"] / agg["tt"] * 1e6,
                         "n": agg["n"], "n_pos": agg["n_pos"]})


def per_section(adata, gi: int, mask: np.ndarray) -> pd.Series:
    counts = counts_matrix(adata[mask])
    obs = adata.obs[mask]
    df = pd.DataFrame({"section": obs["section"].astype(str).to_numpy(),
                       "g": counts[:, gi], "tot": counts.sum(axis=1)})
    g = df.groupby("section").agg(gg=("g", "sum"), tt=("tot", "sum"), n=("g", "size"))
    g = g[g["n"] >= MIN_SECTION]
    return g["gg"] / g["tt"] * 1e6


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(PROC / "hypothalamus_nuclei.h5ad")
    ap = pd.read_csv(AP_FILE, index_col=0)
    genes = adata.var_names.to_numpy()

    tab = pd.crosstab(adata.obs["animal"], adata.obs["nucleus_ext"])
    keep = [n for n in NUCLEI if n in tab.columns]
    cv = (tab[keep].std() / tab[keep].mean() * 100).round(1)
    print("=== Sampling balance between animals (CV of cells per nucleus) ===")
    for n in keep:
        verdict = "usable" if cv[n] < 20 else "AP-limited"
        print(f"  {n:8s} CV {cv[n]:5.1f}%   {verdict}")
    print("  Nuclei above ~20% are unevenly sampled between animals, so a "
          "cross-animal\n  comparison in them is confounded with section level.")

    rows = []
    for nucleus in keep:
        in_nuc = (adata.obs["nucleus_ext"] == nucleus).to_numpy()
        types = adata.obs.loc[in_nuc, "cell_type"].value_counts()
        targets = ["<all cells>"] + [t for t in types.index if types[t] >= 240][:10]
        for cell_type in targets:
            mask = in_nuc if cell_type == "<all cells>" else (
                in_nuc & (adata.obs["cell_type"].astype(str) == cell_type).to_numpy())
            if mask.sum() < MIN_CELLS * 4:
                continue
            for gene in GENES:
                gi = int(np.where(genes == gene)[0][0])
                pa = per_animal(adata, gi, mask)
                if len(pa) < 4 or pa["n"].min() < MIN_CELLS or pa["n_pos"].min() < MIN_POS:
                    continue
                lfc = {b: float(np.log2((pa.loc[a, "cpm"] + 1) / (pa.loc[d, "cpm"] + 1)))
                       for b, (a, d) in BLOCKS.items()}
                consistent = np.sign(lfc["B1"]) == np.sign(lfc["B2"])
                min_abs = min(abs(lfc["B1"]), abs(lfc["B2"]))
                if not (consistent and min_abs >= LFC_THRESHOLD):
                    continue

                ps = per_section(adata, gi, mask)
                pairs = []
                for block, (aged, adult) in BLOCKS.items():
                    a_s = [s for s in ps.index if SECTION_ANIMAL[s] == aged]
                    d_s = [s for s in ps.index if SECTION_ANIMAL[s] == adult]
                    for si, sj in itertools.product(a_s, d_s):
                        pairs.append((abs(ap.loc[si, "ap_score"] - ap.loc[sj, "ap_score"]),
                                      float(np.log2((ps[si] + 1) / (ps[sj] + 1)))))
                if len(pairs) < 6:
                    continue
                gaps, lfcs = np.array([p[0] for p in pairs]), np.array([p[1] for p in pairs])
                r, p = st.pearsonr(gaps, lfcs)
                best = lfcs[np.argmin(gaps)]
                rows.append({
                    "nucleus": nucleus, "cell_type": cell_type, "gene": gene,
                    "n_min": int(pa["n"].min()),
                    "lfc_B1": round(lfc["B1"], 2), "lfc_B2": round(lfc["B2"], 2),
                    "mean_lfc": round((lfc["B1"] + lfc["B2"]) / 2, 2),
                    "n_pairings": len(pairs),
                    "same_sign_pairings": int((np.sign(lfcs) == np.sign(np.mean(lfcs))).sum()),
                    "best_matched_lfc": round(float(best), 2),
                    "r_vs_ap_gap": round(float(r), 2), "p_vs_ap_gap": round(float(p), 3),
                    "ap_independent": bool(p > 0.05),
                    "sampling_cv": float(cv[nucleus]),
                })

    res = pd.DataFrame(rows)
    if res.empty:
        print("\nNo candidate cleared the blocked filters.")
        return 0
    res = res.sort_values(["ap_independent", "mean_lfc"], ascending=[False, False])
    res.to_csv(OUT / "extended_candidates.csv", index=False)

    print(f"\n=== Candidates clearing consistency + >=1.5x, across {len(keep)} nuclei ===")
    with pd.option_context("display.width", 250, "display.max_colwidth", 26):
        print(res[["nucleus", "cell_type", "gene", "n_min", "lfc_B1", "lfc_B2",
                   "best_matched_lfc", "same_sign_pairings", "n_pairings",
                   "r_vs_ap_gap", "p_vs_ap_gap", "ap_independent",
                   "sampling_cv"]].to_string(index=False))

    survivors = res[res["ap_independent"] & (res["sampling_cv"] < 20)]
    print(f"\n=== Surviving both the AP test and even sampling: {len(survivors)} ===")
    if len(survivors):
        print(survivors[["nucleus", "cell_type", "gene", "lfc_B1", "lfc_B2",
                         "r_vs_ap_gap", "p_vs_ap_gap"]].to_string(index=False))
    else:
        print("  (none)")
    ap_dep = res[~res["ap_independent"]]
    if len(ap_dep):
        print(f"\n  {len(ap_dep)} candidate(s) track AP mismatch and are therefore anatomy:")
        for _, r in ap_dep.iterrows():
            print(f"    {r.nucleus:7s} {r.cell_type:24s} {r.gene:6s} "
                  f"r={r.r_vs_ap_gap:+.2f} p={r.p_vs_ap_gap:.3f}")
    print(f"\nWrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
