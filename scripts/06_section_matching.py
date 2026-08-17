#!/usr/bin/env python3
"""Are the four sections cut at comparable rostro-caudal levels?

This matters more than it looks.  Cell-type abundance in the hypothalamus
changes steeply along the AP axis, so if the aged and adult sections sit at
different levels, a "cell type X is more abundant in aged" result can be pure
sampling geometry, and any nucleus-specific claim inherits the same problem.

Two readouts:

* Composition correlation between sections -- how similar the sections are as a
  whole, before any group comparison.
* Per-cell-type coefficient of variation across animals, split by whether the
  variation follows age or not.  A population that varies 3-15x between animals
  without regard to age is telling you about the cut, not about biology.
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
OUT = REPO / "results" / "section_matching"
BLOCKS = {"B1": ("F536", "M493"), "B2": ("G_073", "M399")}

# Populations with a known AP bias in the mouse hypothalamus; their relative
# abundance is a coarse proxy for where each section was cut.
AP_INFORMATIVE = {
    "Oxt/Avp magnocellular": "rostral (PVN/SON)",
    "ARC Agrp/Npy": "mid ARC",
    "ARC Pomc": "mid ARC",
    "VMH-like Glut Rasgrf2": "mid VMH",
    "DMH Grp/Ppp1r17": "caudal (DMH)",
    "LHA Hcrt (orexin)": "caudal-lateral (LHA)",
}


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(PROC / "mbh_roi_annotated.h5ad")

    frac = pd.crosstab(adata.obs["cell_type"], adata.obs["animal"], normalize="columns") * 100
    frac = frac.loc[frac.sum(axis=1) > 0]
    order = ["F536", "G_073", "M493", "M399"]  # aged, aged, adult, adult
    frac = frac[order]

    corr = frac.corr(method="spearman")
    corr.to_csv(OUT / "composition_correlation.csv")
    print("=== Composition correlation between sections (Spearman, cell-type %) ===")
    print(corr.round(3).to_string())
    aged, adult = ["F536", "G_073"], ["M493", "M399"]
    within = np.mean([corr.loc[aged[0], aged[1]], corr.loc[adult[0], adult[1]]])
    between = np.mean([corr.loc[a, b] for a in aged for b in adult])
    print(f"\n  mean within-group r = {within:.3f} | mean between-group r = {between:.3f}")
    print("  (if within >> between the groups differ compositionally; if similar, "
          "section-to-section variation dominates)")

    cv = pd.DataFrame({
        "mean_pct": frac.mean(axis=1).round(3),
        "cv_pct": (frac.std(axis=1) / frac.mean(axis=1) * 100).round(1),
        "max_fold_range": (frac.max(axis=1) / frac.replace(0, np.nan).min(axis=1)).round(1),
    })
    for block, (a, b) in BLOCKS.items():
        cv[f"lfc_{block}"] = np.log2((frac[a] + .01) / (frac[b] + .01)).round(2)
    cv["age_consistent"] = np.sign(cv["lfc_B1"]) == np.sign(cv["lfc_B2"])
    cv["ap_hint"] = [AP_INFORMATIVE.get(t, "") for t in cv.index]
    cv = cv.sort_values("max_fold_range", ascending=False)
    cv.to_csv(OUT / "composition_variability.csv")

    print("\n=== Cell-type abundance variability across the four sections ===")
    with pd.option_context("display.width", 220):
        print(cv.head(14).to_string())

    unstable = cv[(cv["max_fold_range"] >= 2.0) & (~cv["age_consistent"])]
    print(f"\n{len(unstable)} cell types vary >=2x across sections WITHOUT following age.")
    print("These are sampling/section-plane effects, and any age result in them is unsafe:")
    print("  " + ", ".join(unstable.index[:12]))

    print("\n=== AP-informative populations (% of each section's ROI cells) ===")
    ap = frac.loc[[t for t in AP_INFORMATIVE if t in frac.index]].round(3)
    ap["hint"] = [AP_INFORMATIVE[t] for t in ap.index]
    print(ap.to_string())

    _plot(frac, cv, order)
    print(f"\nWrote tables and figures to {OUT}")
    return 0


def _plot(frac, cv, order) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(17, 6.5))
    colors = {"F536": "#B4531A", "G_073": "#8C3F12", "M493": "#0E5A61", "M399": "#0A4046"}

    ax = axes[0]
    top = frac.loc[frac.mean(axis=1).sort_values(ascending=False).index[:18]]
    y = np.arange(len(top))
    for i, animal in enumerate(order):
        ax.barh(y + (i - 1.5) * .2, top[animal], height=.2, color=colors[animal], label=animal)
    ax.set_yticks(y); ax.set_yticklabels(top.index, fontsize=8); ax.invert_yaxis()
    ax.set_xlabel("% of ROI cells"); ax.set_title("Cell-type composition per section")
    ax.legend(fontsize=8)

    ax = axes[1]
    sub = cv[cv["mean_pct"] > 0.2].sort_values("max_fold_range")
    ax.barh(np.arange(len(sub)), sub["max_fold_range"],
            color=["#B4531A" if c else "#9AA8A8" for c in sub["age_consistent"]])
    ax.axvline(2, color="k", ls="--", lw=1)
    ax.set_yticks(np.arange(len(sub))); ax.set_yticklabels(sub.index, fontsize=7)
    ax.set_xlabel("max / min abundance across sections")
    ax.set_title("Abundance variability\n(orange = direction follows age)")

    ax = axes[2]
    ax.scatter(cv["max_fold_range"], cv[["lfc_B1", "lfc_B2"]].abs().min(axis=1),
               c=["#B4531A" if c else "#9AA8A8" for c in cv["age_consistent"]], s=34)
    ax.set_xscale("log")
    ax.set_xlabel("abundance variability across sections (fold)")
    ax.set_ylabel("min |log2 FC| between blocks")
    ax.set_title("Age effect vs sampling noise")
    ax.grid(alpha=.25)

    fig.tight_layout()
    fig.savefig(OUT / "section_matching.png", dpi=155)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
