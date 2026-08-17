#!/usr/bin/env python3
"""Estimate each section's rostro-caudal level, and pick an AP-matched core set.

The rostro-caudal gradient is what dismantled the Galr1 lead: one animal's own
sections span more Galr1 than the age groups differ by.  The fix is to compare
sections cut at the same level.

The AP score is built from **morphometry** -- physical measurements of the
section in its own anatomical frame -- not from cell-type composition.  That
distinction matters.  Matching on composition would risk regressing out the very
age differences the study is looking for, since composition is exactly what an
age effect would change.  Ventricle length and tissue width are structural, and
far less plausible as age-regulated quantities.

Features (all in micrometres, all in the anatomical frame):
  * dorsoventral extent of the third-ventricle lining
  * tissue half-width at four dorsoventral levels
  * dorsoventral extent of tissue on the midline
  * dorsoventral position and mediolateral width of the VMH cell mass
  * cell density in the window

Validation asks whether the score orders each animal's sections the way they
were cut, which is the only ground truth available.
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats as st

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.io import ANIMAL_META, SECTION_ANIMAL  # noqa: E402

PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "ap_matching"
ANIMALS = ["F536", "G_073", "M493", "M399"]
DV_LEVELS = (200, 600, 1000, 1400)
LINING_TYPES = ("Tanycyte", "Ependymal")
VMH_TYPES = ("VMH-like Glut Rasgrf2", "VMH-like Glut Calb1", "VMH-like Glut Tac1")


def morphometry(obs: pd.DataFrame) -> dict:
    ml, dv = obs["ml"].to_numpy(), obs["dv"].to_numpy()
    ct = obs["cell_type"].astype(str).to_numpy()

    feats = {}
    lining = np.isin(ct, LINING_TYPES) & (np.abs(ml) < 200)
    if lining.sum() > 30:
        feats["v3_length"] = float(np.quantile(dv[lining], .95) - np.quantile(dv[lining], .05))
        feats["v3_top"] = float(np.quantile(dv[lining], .95))
    else:
        feats["v3_length"] = feats["v3_top"] = np.nan

    for level in DV_LEVELS:
        band = np.abs(dv - level) < 75
        feats[f"halfwidth_dv{level}"] = (
            float(np.quantile(np.abs(ml[band]), .95)) if band.sum() > 50 else np.nan
        )

    midline = np.abs(ml) < 150
    feats["midline_dv_extent"] = (
        float(np.quantile(dv[midline], .98)) if midline.sum() > 50 else np.nan
    )

    vmh = np.isin(ct, VMH_TYPES)
    if vmh.sum() > 50:
        feats["vmh_dv_centre"] = float(np.median(dv[vmh]))
        feats["vmh_ml_width"] = float(np.quantile(np.abs(ml[vmh]), .90))
        feats["vmh_dv_extent"] = float(np.quantile(dv[vmh], .95) - np.quantile(dv[vmh], .05))
    else:
        feats["vmh_dv_centre"] = feats["vmh_ml_width"] = feats["vmh_dv_extent"] = np.nan

    feats["n_cells_window"] = float(len(obs))
    return feats


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(PROC / "hypothalamus_domains.h5ad")
    obs = adata.obs

    rows = {}
    for section, sub in obs.groupby("section", observed=True):
        rows[str(section)] = morphometry(sub)
    feats = pd.DataFrame(rows).T
    feats["animal"] = [SECTION_ANIMAL[s] for s in feats.index]
    feats["idx"] = [int(s.rsplit("_", 1)[1]) for s in feats.index]
    feats["group"] = feats["animal"].map(lambda a: ANIMAL_META[a]["age_group"])

    num = feats.drop(columns=["animal", "idx", "group"]).astype(float)
    num = num.loc[:, num.notna().all()]
    z = (num - num.mean()) / num.std(ddof=0)

    u, s, vt = np.linalg.svd(z.to_numpy(), full_matrices=False)
    ap = pd.Series(u[:, 0] * s[0], index=z.index)
    # Orient so that increasing score means increasing section index on average,
    # giving the axis a consistent direction across runs.
    if np.corrcoef(ap, feats["idx"])[0, 1] < 0:
        ap = -ap
    feats["ap_score"] = ap
    var = (s**2 / (s**2).sum())[0]

    print("=== Morphometric AP score ===")
    print(f"  built from {z.shape[1]} features; PC1 explains {var*100:.0f}% of variance")
    print(f"  loadings: " + ", ".join(
        f"{c} {vt[0][i]:+.2f}" for i, c in enumerate(z.columns)))
    print()
    show = feats[["animal", "group", "idx", "ap_score", "v3_length",
                  "vmh_dv_centre", "vmh_ml_width", "midline_dv_extent"]]
    print(show.sort_values("ap_score").round(1).to_string())

    print("\n=== Validation: does the score order sections as they were cut? ===")
    agree = 0
    for animal, sub in feats.groupby("animal", observed=True):
        r = st.spearmanr(sub["idx"], sub["ap_score"]).statistic
        agree += r > 0
        print(f"  {animal:6s} idx {list(sub.sort_values('idx')['idx'])} -> "
              f"ap {[round(v,1) for v in sub.sort_values('idx')['ap_score']]}   rho = {r:+.1f}")
    print(f"  {agree}/4 animals ordered consistently with cut order")

    print("\n=== Is the AP score confounded with age? ===")
    for grp, sub in feats.groupby("group", observed=True):
        print(f"  {grp:6s} ap_score {sub.ap_score.min():+.2f} .. {sub.ap_score.max():+.2f}")
    aged, adult = feats[feats.group == "aged"], feats[feats.group == "adult"]
    overlap = min(aged.ap_score.max(), adult.ap_score.max()) - max(aged.ap_score.min(), adult.ap_score.min())
    print(f"  ranges overlap: {overlap > 0}  (overlap width {overlap:+.2f})")

    matched = _select_matched(feats)
    feats.to_csv(OUT / "section_ap_scores.csv")
    matched.to_csv(OUT / "matched_sets.csv", index=False)
    _plot(feats, matched)
    print(f"\nWrote {OUT}")
    return 0


def _select_matched(feats: pd.DataFrame) -> pd.DataFrame:
    """Exhaustively pick the tightest set with a fixed number of sections/animal."""
    by_animal = {a: list(sub.index) for a, sub in feats.groupby("animal", observed=True)}
    rows = []
    for per_animal in (1, 2):
        combos = itertools.product(*[
            list(itertools.combinations(by_animal[a], per_animal)) for a in ANIMALS
        ])
        best = None
        for combo in combos:
            picked = [s for group in combo for s in group]
            vals = feats.loc[picked, "ap_score"]
            spread = float(vals.max() - vals.min())
            if best is None or spread < best[0]:
                best = (spread, picked)
        spread, picked = best
        rows.append({
            "per_animal": per_animal, "n_sections": len(picked),
            "ap_spread": round(spread, 2),
            "sections": ", ".join(sorted(picked)),
        })
        print(f"\n=== Tightest set with {per_animal} section(s) per animal ===")
        print(f"  sections : {', '.join(sorted(picked))}")
        print(f"  AP spread: {spread:.2f}  (full 12-section spread "
              f"{feats.ap_score.max()-feats.ap_score.min():.2f})")
        sub = feats.loc[picked, ["animal", "group", "idx", "ap_score"]].sort_values("ap_score")
        print(sub.round(2).to_string())
    return pd.DataFrame(rows)


def _plot(feats, matched) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    colours = {"aged": "#B4531A", "adult": "#0E5A61"}

    ax = axes[0]
    for animal, sub in feats.groupby("animal", observed=True):
        sub = sub.sort_values("idx")
        c = colours[sub["group"].iloc[0]]
        ax.plot(sub["idx"], sub["ap_score"], "-o", color=c, label=f"{animal} ({sub.group.iloc[0]})")
        for _, r in sub.iterrows():
            ax.annotate(r.name.split("_")[-1], (r["idx"], r["ap_score"]),
                        textcoords="offset points", xytext=(6, 0), fontsize=8)
    ax.set_xlabel("section index (cut order)"); ax.set_ylabel("morphometric AP score")
    ax.set_title("AP score vs cut order\n(monotone within animal = the score tracks position)")
    ax.legend(fontsize=8); ax.grid(alpha=.25)

    ax = axes[1]
    best = matched[matched.per_animal == 1]["sections"].iloc[0].split(", ")
    for _, r in feats.iterrows():
        chosen = r.name in best
        ax.scatter(r["ap_score"], r["animal"], s=150 if chosen else 55,
                   color=colours[r["group"]], edgecolor="k" if chosen else "none",
                   linewidth=1.6 if chosen else 0, zorder=3 if chosen else 2)
    lo, hi = feats.loc[best, "ap_score"].min(), feats.loc[best, "ap_score"].max()
    ax.axvspan(lo, hi, color="#9AA8A8", alpha=.22, zorder=1)
    ax.set_xlabel("morphometric AP score"); ax.set_ylabel("")
    ax.set_title("AP-matched core set (ringed)\nshaded band = span of the matched set")
    ax.grid(axis="x", alpha=.25)

    fig.tight_layout()
    fig.savefig(OUT / "ap_matching.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
