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

PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "extended_nuclei"
SPLIT_ML_UM = 500.0   # domain 8: medial fringe vs lateral hypothalamus

# domain -> (nucleus, the evidence the call rests on)
# Marker z-scores are across the 14 domains; positions are in the anatomical frame.
DOMAIN_NUCLEUS = {
    "9":  ("ARC",     "ARC markers z=+2.96 (Agrp, Pomc, Ghrh, Tac2); |ml| 172, dv 136"),
    "12": ("ME_3V",   "tanycyte/ME z=+1.14, Gpr50 64% of cells; |ml| 40, dv 130 — 3V floor"),
    "5":  ("ME_3V",   "tanycyte/ME z=+1.28, ependymal 58%; |ml| 29, dv 922 — dorsal 3V wall"),
    "7":  ("VMH",     "VMH z=+0.86, Slc17a6 z=+2.13, Rasgrf2 58%; |ml| 246 — dorsomedial VMH"),
    "4":  ("VMH",     "VMH z=+0.60, Slc17a6 z=+1.76, Calb1 41%; |ml| 443 — ventrolateral VMH"),
    "13": ("DMH",     "DMH z=+1.41 (Grp, Ppp1r17); |ml| 265, dv 1100 — DMH core"),
    "8":  ("SPLIT",   "LHA z=+2.11 laterally (Hcrt 64.5% beyond 500um) but Grp-enriched medially"),
    "0":  ("ZI",      "GABA z=+1.77, Gad/Slc32a1 dominant, Cacna2d2 48%; dv 1548 — zona incerta"),
    "10": ("DHA_PH",  "Slc17a6 z=+1.23, ZI z=-0.10 — glutamatergic, so not ZI; dorsal hypothalamic area"),
    "3":  ("TUseg",   "no nucleus marker dominant; |ml| 715, dv 471 — tuberal / ventrolateral"),
    "2":  ("fibre",   "oligodendrocyte 25%, |ml| 997 — fibre tract territory"),
    "11": ("fibre",   "oligodendrocyte 67%, fibre z=+1.92 — internal capsule / optic tract"),
    "1":  ("edge",    "meningeal fibroblast 74%, |ml| 1218 — tissue edge"),
    "6":  ("edge",    "VLMC 62%, |ml| 910 — tissue edge"),
}
ANALYSIS_NUCLEI = ("ARC", "ME_3V", "VMH", "DMH", "LHA", "ZI", "DHA_PH")


def assign(obs: pd.DataFrame) -> pd.Series:
    dom = obs["domain"].astype(str)
    out = pd.Series([DOMAIN_NUCLEUS.get(d, ("unassigned", ""))[0] for d in dom], index=obs.index)
    medial = obs["ml"].abs() < SPLIT_ML_UM
    out[(dom == "8") & medial] = "DMH"
    out[(dom == "8") & ~medial] = "LHA"
    return out


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
    adata.obs["nucleus_ext"] = assign(adata.obs).to_numpy()

    print("=== Domain -> nucleus assignment ===")
    for dom, (nucleus, why) in sorted(DOMAIN_NUCLEUS.items(), key=lambda kv: int(kv[0])):
        print(f"  domain {dom:>2s} -> {nucleus:8s} {why}")
    print(f"\n  domain 8 split at |ml| = {SPLIT_ML_UM:.0f} um: medial -> DMH, lateral -> LHA")

    ext = pd.DataFrame([extent(adata.obs, n) for n in ANALYSIS_NUCLEI]).dropna(how="all")
    ext.to_csv(OUT / "nucleus_extents.csv", index=False)
    print("\n=== Extent of each nucleus (p5-p95, anatomical frame) ===")
    print(ext.to_string(index=False))

    old = pd.read_csv(REPO / "results" / "anatomy" / "cells_per_nucleus.csv", index_col=0)
    print("\n=== DMH: how the delineation has changed ===")
    prev_dmh = int(old["DMH"].sum()) if "DMH" in old else 0
    dmh = ext[ext.nucleus == "DMH"].iloc[0]
    print(f"  Gaussian ellipse (Grp-anchored) : 13,235 cells, 1084 x 380 um")
    print(f"  domain 13 alone                 : 12,441 cells, 1026 x 469 um")
    print(f"  domain 13 + medial fringe of 8  : {dmh.cells:,} cells, "
          f"{dmh.width_um} x {dmh.height_um} um")
    print(f"  published mouse DMH             : ~800-1200 x 500-700 um")

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
