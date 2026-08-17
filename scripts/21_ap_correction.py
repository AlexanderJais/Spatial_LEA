#!/usr/bin/env python3
"""Replace the AP diagnostic with a within-animal AP correction.

The audit (20_audit_ap_test.py) found the pooled "effect size vs AP mismatch"
correlation unusable: the AP gaps available in block B1 (2.36-4.94) and block B2
(0.04-2.00) do not overlap at all, so the correlation tests a block difference;
the 18 pairings share sections, making the naive p-value anticonservative by
roughly 40x; and 38% of panel genes fire at p<0.05, against the 5% a calibrated
test would give.

This uses the one AP relationship that cannot be confounded with age: the slope
**within** an animal.  Age is constant inside an animal, so any dependence of
expression on rostro-caudal level measured there is anatomy by construction.

  1. fit  log2 CPM ~ AP + animal  on the 12 sections, taking only the common
     within-animal slope beta
  2. correct every section to a common AP level:  y* = y - beta * (ap - ap_mean)
  3. recompute the blocked age effect on the corrected values

If an apparent age effect is anatomical it shrinks toward zero; if it is real it
survives.  Everything is calibrated against the whole panel so the size of the
change has a reference distribution.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc
import statsmodels.formula.api as smf

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.io import ANIMAL_META, SECTION_ANIMAL, counts_matrix  # noqa: E402

PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "ap_correction"
AP_FILE = REPO / "results" / "ap_matching" / "section_ap_scores.csv"
BLOCKS = {"B1": ("F536", "M493"), "B2": ("G_073", "M399")}
CASES = [("Galr1", "GABA Gal/Galr1"), ("Gal", "ARC Agrp/Npy"), ("Gal", "GABA Cacna2d2")]
MIN_CELLS = 15
MIN_COUNTS = 200


def section_table(adata, cell_type: str):
    sel = (adata.obs["cell_type"].astype(str) == cell_type).to_numpy()
    sub = adata[sel]
    counts = counts_matrix(sub)
    sections = sub.obs["section"].astype(str).to_numpy()
    keep = [s for s in pd.unique(sections) if (sections == s).sum() >= MIN_CELLS]
    mat = np.vstack([counts[sections == s].sum(axis=0) for s in keep])
    return pd.DataFrame(mat, index=keep, columns=adata.var_names)


def fit_and_correct(cpm_log: pd.Series, ap: pd.Series):
    """Common within-animal slope on AP, and the age effect before/after."""
    df = pd.DataFrame({"y": cpm_log, "ap": ap[cpm_log.index],
                       "animal": [SECTION_ANIMAL[s] for s in cpm_log.index]})
    if df["animal"].nunique() < 4:
        return None
    # animal as a fixed effect absorbs every between-animal difference, including
    # age, so the AP slope is estimated purely from variation inside animals.
    model = smf.ols("y ~ ap + C(animal)", df).fit()
    beta = float(model.params["ap"])
    p_beta = float(model.pvalues["ap"])

    df["y_corr"] = df["y"] - beta * (df["ap"] - df["ap"].mean())
    out = {"beta": beta, "p_beta": p_beta}
    for label, col in (("raw", "y"), ("corr", "y_corr")):
        per_animal = df.groupby("animal")[col].mean()
        for block, (aged, adult) in BLOCKS.items():
            out[f"{label}_{block}"] = float(per_animal[aged] - per_animal[adult])
        out[f"{label}_mean"] = np.mean([out[f"{label}_{b}"] for b in BLOCKS])
    out["shrinkage"] = 1 - abs(out["corr_mean"]) / max(abs(out["raw_mean"]), 1e-9)
    return out


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(PROC / "hypothalamus_nuclei.h5ad")
    ap = pd.read_csv(AP_FILE, index_col=0)["ap_score"]

    summary = []
    panel_ref = {}
    for gene, cell_type in CASES:
        mat = section_table(adata, cell_type)
        if len(mat) < 10:
            continue
        totals = mat.sum(axis=1)

        # Panel reference: the same correction applied to every gene here, so
        # the shrinkage of the gene of interest has a distribution to sit in.
        shrinks, betas = [], []
        for g in mat.columns:
            if mat[g].sum() < MIN_COUNTS:
                continue
            y = np.log2(mat[g] / totals * 1e6 + 1)
            if y.std() < 1e-9:
                continue
            res = fit_and_correct(y, ap)
            if res is None:
                continue
            shrinks.append({"gene": g, **res})
        panel = pd.DataFrame(shrinks).set_index("gene")
        panel_ref[(gene, cell_type)] = panel
        panel.to_csv(OUT / f"panel_{cell_type.replace('/', '_').replace(' ', '_')}.csv")

        row = panel.loc[gene]
        pct = float((panel["shrinkage"] >= row["shrinkage"]).mean() * 100)
        summary.append({
            "gene": gene, "cell_type": cell_type,
            "ap_slope": round(row["beta"], 3), "p_slope": round(row["p_beta"], 4),
            "raw_B1": round(row["raw_B1"], 2), "raw_B2": round(row["raw_B2"], 2),
            "corr_B1": round(row["corr_B1"], 2), "corr_B2": round(row["corr_B2"], 2),
            "raw_mean": round(row["raw_mean"], 2), "corr_mean": round(row["corr_mean"], 2),
            "shrinkage_pct": round(row["shrinkage"] * 100, 1),
            "panel_rank_pct": round(pct, 1),
            "sign_kept": bool(np.sign(row["corr_B1"]) == np.sign(row["corr_B2"])
                              == np.sign(row["raw_mean"])),
        })

    res = pd.DataFrame(summary)
    res.to_csv(OUT / "ap_corrected_candidates.csv", index=False)
    print("=== Age effect before and after within-animal AP correction ===")
    print("   ap_slope is the anatomical gradient, estimated where age is constant.")
    print("   shrinkage_pct is how much of the effect the correction removes.")
    print("   panel_rank_pct: % of panel genes shrinking at least as much (100 = typical).\n")
    with pd.option_context("display.width", 240):
        print(res.to_string(index=False))

    print("\n=== Reading ===")
    for _, r in res.iterrows():
        if r.shrinkage_pct > 60:
            verdict = "largely anatomical"
        elif r.shrinkage_pct > 30:
            verdict = "partly anatomical"
        else:
            verdict = "survives correction"
        print(f"  {r.gene:6s} in {r.cell_type:16s}  {r.raw_mean:+.2f} -> {r.corr_mean:+.2f} "
              f"({r.shrinkage_pct:5.1f}% removed, panel rank {r.panel_rank_pct:5.1f}%)  {verdict}")

    _plot(res, panel_ref)
    print(f"\nWrote {OUT}")
    return 0


def _plot(res, panel_ref) -> None:
    fig, axes = plt.subplots(1, len(panel_ref), figsize=(4.6 * len(panel_ref), 4.2))
    axes = np.atleast_1d(axes)
    for ax, ((gene, cell_type), panel) in zip(axes, panel_ref.items()):
        ax.hist(panel["shrinkage"] * 100, bins=40, color="#C3CCCC")
        ax.axvline(panel.loc[gene, "shrinkage"] * 100, color="#B4531A", lw=1.8)
        ax.annotate(gene, (panel.loc[gene, "shrinkage"] * 100, ax.get_ylim()[1] * .85),
                    color="#B4531A", fontsize=9, style="italic", ha="left")
        ax.set_xlabel("% of age effect removed by AP correction")
        ax.set_ylabel("panel genes")
        ax.set_title(f"{gene} in {cell_type}", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "ap_correction.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
