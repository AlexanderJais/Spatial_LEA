#!/usr/bin/env python3
"""Map spatial domains onto named nuclei, and extend beyond ARC/VMH/DMH.

The Gaussian model delineated three nuclei from marker anchors and under-called
the DMH.  The 14 tissue-derived domains cover the whole hypothalamic window, so
they can carry a wider and better-supported parcellation -- provided each domain
is assigned on evidence rather than by eye.

Each domain is scored against nucleus-diagnostic marker sets and placed by its
position in the anatomical frame.  Both are printed, so a disputed assignment
can be argued from the table.

One domain is split rather than assigned whole: domain 8 spans |ml| 460-1011 um
and is Hcrt-dominated laterally (64.5% Hcrt+ beyond 500 um) but Grp-enriched
medially (18.6% Grp+ within 500 um).  It is one niche by composition but two
structures by anatomy -- the medial fringe of the DMH and the LHA.
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

from spatial_lea.io import counts_matrix  # noqa: E402

PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "extended_nuclei"
SPLIT_ML_UM = 500.0   # domain 8: medial fringe vs lateral hypothalamus

# Nucleus definitions.  Each is a marker set plus the position it occupies in
# the intrinsic anatomical frame, in micrometres.  Domains are matched to these
# by evidence at run time.
#
# This replaces a hardcoded {KMeans domain id -> nucleus} table.  That table was
# written against the domains derived from the first cohort; KMeans labels are
# arbitrary integers, so when the domains were recomputed on the full set the
# ids no longer meant what they had, and the old table was applied to the wrong
# clusters.  Only the ARC happened to survive.  A mapping keyed on cluster
# number cannot be reused across runs, so the assignment is now derived.
NUCLEUS_DEF = {
    "ARC":    {"markers": ["Agrp", "Pomc", "Ghrh", "Slc6a3"], "ml": 200, "dv": 150},
    "ME_3V":  {"markers": ["Gpr50", "Spag16"],                "ml": 60,  "dv": 400},
    "VMH":    {"markers": ["Adcyap1", "Calb1", "Rasgrf2"],    "ml": 380, "dv": 560},
    "DMH":    {"markers": ["Grp", "Ppp1r17"],                 "ml": 260, "dv": 1050},
    "LHA":    {"markers": ["Hcrt"],                           "ml": 900, "dv": 900},
    "ZI":     {"markers": ["Cacna2d2", "Pvalb", "Slc32a1"],   "ml": 700, "dv": 1450},
    "DHA_PH": {"markers": ["Slc17a6", "Otp"],                 "ml": 500, "dv": 1330},
}
# How far a domain may sit from a nucleus centre before position stops
# supporting the call.  Generous: it separates nuclei, it does not define them.
POS_SCALE_UM = 700.0
# A domain has to be positively supported to be given a nucleus name.  Below
# this it is left as unassigned tuberal territory rather than absorbed into the
# nearest nucleus, which is what inflated the VMH on the first pass.
MIN_SCORE = 0.45


def domain_profile(adata, obs: pd.DataFrame) -> tuple:
    """Per-domain marker z-scores, position, and composition."""
    doms = sorted(obs["domain"].astype(str).unique(), key=int)
    counts = counts_matrix(adata)
    genes = adata.var_names.to_numpy()
    dom = obs["domain"].astype(str).to_numpy()
    cpm = {}
    for d in doms:
        m = dom == d
        c = counts[m]
        cpm[d] = np.log2(c.sum(axis=0) / max(c.sum(), 1) * 1e6 + 1)
    cpm = pd.DataFrame(cpm, index=genes).T
    z = (cpm - cpm.mean()) / cpm.std().replace(0, np.nan)
    pos = obs.assign(dom=dom).groupby("dom").agg(
        absml=("ml", lambda x: float(np.median(np.abs(x)))),
        dv=("dv", "median"), cells=("dv", "size"))
    return z, pos.loc[doms]


def score_domains(z: pd.DataFrame, pos: pd.DataFrame) -> pd.DataFrame:
    """Evidence for each nucleus in each domain: markers plus position."""
    rows = {}
    for nucleus, spec in NUCLEUS_DEF.items():
        present = [g for g in spec["markers"] if g in z.columns]
        marker = z[present].mean(axis=1) if present else pd.Series(0.0, index=z.index)
        dist = np.hypot(pos["absml"] - spec["ml"], pos["dv"] - spec["dv"])
        rows[nucleus] = marker - dist / POS_SCALE_UM
    return pd.DataFrame(rows)


def assign_by_evidence(scores: pd.DataFrame, obs: pd.DataFrame,
                       composition: pd.DataFrame) -> dict:
    """Best-supported nucleus per domain; glial-dominated domains stay unnamed."""
    out = {}
    for d in scores.index:
        top = composition.loc[d].idxmax() if d in composition.index else ""
        frac = composition.loc[d].max() if d in composition.index else 0.0
        if top in ("ME / meningeal fibroblast", "VLMC / fibroblast",
                   "Arachnoid fibroblast") and frac > .35:
            out[d] = ("edge", f"{top} {frac*100:.0f}% — tissue edge")
            continue
        if top == "Oligodendrocyte" and frac > .35:
            out[d] = ("fibre", f"oligodendrocyte {frac*100:.0f}% — fibre tract")
            continue
        best = scores.loc[d].idxmax()
        value = scores.loc[d].max()
        if value < MIN_SCORE:
            out[d] = ("TUseg", f"no nucleus supported (best {best} {value:+.2f})")
        else:
            out[d] = (best, f"{best} score {value:+.2f}")
    return out


ANALYSIS_NUCLEI = ("ARC", "ME_3V", "VMH", "DMH", "LHA", "ZI", "DHA_PH")


def assign(obs: pd.DataFrame, mapping: dict) -> pd.Series:
    dom = obs["domain"].astype(str)
    return pd.Series([mapping.get(d, ("unassigned", ""))[0] for d in dom],
                     index=obs.index)


def extent(obs: pd.DataFrame, nucleus: str) -> dict:
    sub = obs[obs["nucleus_ext"] == nucleus]
    if not len(sub):
        return {}
    ml, dv = sub["ml"].abs(), sub["dv"]
    return {
        "nucleus": nucleus, "cells": len(sub), "per_section": round(len(sub) / 12),
        "ml_p5": round(ml.quantile(.05)), "ml_p95": round(ml.quantile(.95)),
        "dv_p5": round(dv.quantile(.05)), "dv_p95": round(dv.quantile(.95)),
        "width_um": round(ml.quantile(.95) * 2),
        "height_um": round(dv.quantile(.95) - dv.quantile(.05)),
    }


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(PROC / "hypothalamus_domains.h5ad")
    obs = adata.obs
    z, pos = domain_profile(adata, obs)
    scores = score_domains(z, pos)
    comp = pd.crosstab(obs["domain"].astype(str), obs["cell_type"].astype(str),
                       normalize="index")
    mapping = assign_by_evidence(scores, obs, comp)
    adata.obs["nucleus_ext"] = assign(obs, mapping).to_numpy()

    print("=== Domain -> nucleus assignment, derived from the data ===")
    print(f"  {'dom':>3s} {'cells':>6s} {'|ml|':>6s} {'dv':>6s}  {'call':8s} evidence")
    for d in sorted(mapping, key=int):
        nucleus, why = mapping[d]
        print(f"  {d:>3s} {int(pos.loc[d, 'cells']):6d} "
              f"{pos.loc[d, 'absml']:6.0f} {pos.loc[d, 'dv']:6.0f}  {nucleus:8s} {why}")
    print("\n=== Score matrix (marker z minus distance penalty) ===")
    print(scores.round(2).to_string())

    ext = pd.DataFrame([extent(adata.obs, n) for n in ANALYSIS_NUCLEI]).dropna(how="all")
    ext.to_csv(OUT / "nucleus_extents.csv", index=False)
    print("\n=== Extent of each nucleus (p5-p95, anatomical frame) ===")
    print(ext.to_string(index=False))

    if "DMH" in set(ext.nucleus):
        dmh = ext[ext.nucleus == "DMH"].iloc[0]
        print("\n=== DMH against published dimensions ===")
        print(f"  derived here        : {dmh.cells:,} cells, "
              f"{dmh.width_um} x {dmh.height_um} um")
        print(f"  published mouse DMH : ~800-1200 x 500-700 um")

    print("\n=== Cells per nucleus per animal (balance across the design) ===")
    tab = pd.crosstab(adata.obs["animal"], adata.obs["nucleus_ext"])
    keep = [n for n in ANALYSIS_NUCLEI if n in tab.columns]
    print(tab[keep].to_string())
    cv = (tab[keep].std() / tab[keep].mean() * 100).round(1)
    print("\n  between-animal CV (%):", cv.to_dict())

    adata.write_h5ad(PROC / "hypothalamus_nuclei.h5ad")
    tab.to_csv(OUT / "cells_per_nucleus_per_animal.csv")
    _plot(adata)
    print(f"\nWrote {PROC / 'hypothalamus_nuclei.h5ad'} and {OUT}")
    return 0


def _plot(adata) -> None:
    colours = {"ARC": "#B4531A", "ME_3V": "#E0A468", "VMH": "#0E5A61", "DMH": "#6A3D9A",
               "LHA": "#2E7D32", "ZI": "#C2185B", "DHA_PH": "#5D6D7E", "TUseg": "#B0BEC5",
               "fibre": "#ECEFF1", "edge": "#F5F5F5", "unassigned": "#FAFAFA"}
    obs = adata.obs
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.5))

    ax = axes[0]
    sample = obs.sample(min(120000, len(obs)), random_state=0)
    ax.scatter(sample["ml"], sample["dv"], s=.8, alpha=.6,
               c=[colours.get(n, "#EEE") for n in sample["nucleus_ext"]], rasterized=True)
    for nucleus in ("ARC", "ME_3V", "VMH", "DMH", "LHA", "ZI", "DHA_PH"):
        sub = obs[obs["nucleus_ext"] == nucleus]
        if len(sub) < 100:
            continue
        ax.annotate(nucleus, (sub["ml"].abs().median(), sub["dv"].median()),
                    fontsize=11, fontweight="bold", ha="center", color="#111",
                    bbox=dict(boxstyle="round,pad=.18", fc="white", ec="none", alpha=.8))
    ax.set_xlim(-1500, 1500); ax.set_ylim(-150, 1800); ax.set_aspect("equal")
    ax.set_xlabel("ml (µm)"); ax.set_ylabel("dv above ventral surface (µm)")
    ax.set_title("Extended parcellation, all 12 sections pooled")

    ax = axes[1]
    d8 = obs[obs["domain"].astype(str) == "8"]
    ax.scatter(d8["ml"], d8["dv"], s=1.2, alpha=.5,
               c=["#6A3D9A" if abs(m) < SPLIT_ML_UM else "#2E7D32" for m in d8["ml"]],
               rasterized=True)
    for x in (-SPLIT_ML_UM, SPLIT_ML_UM):
        ax.axvline(x, color="k", lw=1, ls="--")
    ax.set_xlim(-1500, 1500); ax.set_ylim(600, 1600); ax.set_aspect("equal")
    ax.set_xlabel("ml (µm)"); ax.set_ylabel("dv (µm)")
    ax.set_title("Domain 8 split: medial fringe (DMH, purple)\nvs lateral hypothalamus (green)")

    fig.tight_layout()
    fig.savefig(OUT / "extended_nuclei.png", dpi=145)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
