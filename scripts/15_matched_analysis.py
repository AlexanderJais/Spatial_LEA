#!/usr/bin/env python3
"""Re-test the galanin candidates on AP-matched sections, block by block.

Matching is done *within block* rather than across all four animals.  The reason
is in the AP scores: F536's sections sit at +2.2 to +4.2 while the other three
animals span -2.1 to -0.1, so no four-animal set can be well matched.  Within
block B2 (G_073 vs M399) a near-perfect pair exists; within block B1 (F536 vs
M493) the best available pair is still far apart, and that residual is reported
rather than hidden.

For each candidate this compares three things:
  * the effect using all sections
  * the effect using the best-matched pair in each block
  * the residual AP difference of that pair, so the reader can judge it
"""

from __future__ import annotations

import itertools
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

from spatial_lea.io import ANIMAL_META, SECTION_ANIMAL, counts_matrix  # noqa: E402

PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "matched"
AP = REPO / "results" / "ap_matching" / "section_ap_scores.csv"
BLOCKS = {"B1": ("F536", "M493"), "B2": ("G_073", "M399")}

CANDIDATES = [
    ("DMH Galr1 lead", "GABA Gal/Galr1", "Galr1", "13"),
    ("best Gal candidate", "ARC Agrp/Npy", "Gal", "9"),
    ("Gal, DMH GABA Cacna2d2", "GABA Cacna2d2", "Gal", None),
]


def best_pairs(ap: pd.DataFrame) -> dict:
    out = {}
    for block, (aged, adult) in BLOCKS.items():
        a = ap[ap.animal == aged]
        d = ap[ap.animal == adult]
        pairs = [
            (abs(a.loc[i, "ap_score"] - d.loc[j, "ap_score"]), i, j)
            for i, j in itertools.product(a.index, d.index)
        ]
        gap, si, sj = min(pairs)
        out[block] = {"aged_section": si, "adult_section": sj, "ap_gap": round(gap, 2)}
    return out


def measure(adata, cell_type, gene, sections=None) -> pd.Series:
    sel = (adata.obs["cell_type"].astype(str) == cell_type).to_numpy()
    if sections is not None:
        sel &= adata.obs["section"].astype(str).isin(sections).to_numpy()
    sub = adata[sel]
    if sub.n_obs == 0:
        return pd.Series(dtype=float)
    gi = int(np.where(adata.var_names.to_numpy() == gene)[0][0])
    counts = counts_matrix(sub)
    df = pd.DataFrame({
        "animal": sub.obs["animal"].astype(str).to_numpy(),
        "g": counts[:, gi],
        "tot": counts.sum(axis=1),
    })
    g = df.groupby("animal")
    return pd.Series({
        **{f"cpm_{a}": v for a, v in (g["g"].sum() / g["tot"].sum() * 1e6).items()},
        **{f"det_{a}": v for a, v in (g.apply(lambda x: (x["g"] > 0).mean() * 100)).items()},
        **{f"n_{a}": v for a, v in g.size().items()},
    })


def lfc(vals: pd.Series, prefix: str) -> dict:
    out = {}
    for block, (aged, adult) in BLOCKS.items():
        ka, kd = f"{prefix}_{aged}", f"{prefix}_{adult}"
        if ka in vals and kd in vals:
            out[block] = float(np.log2((vals[ka] + 1) / (vals[kd] + 1)))
    return out


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    ap = pd.read_csv(AP, index_col=0)
    adata = sc.read_h5ad(PROC / "hypothalamus_domains.h5ad")

    pairs = best_pairs(ap)
    print("=== Best AP-matched pair within each block ===")
    for block, spec in pairs.items():
        aged, adult = BLOCKS[block]
        print(f"  {block}: {spec['aged_section']} ({aged}, aged, ap "
              f"{ap.loc[spec['aged_section'],'ap_score']:+.2f})  vs  "
              f"{spec['adult_section']} ({adult}, adult, ap "
              f"{ap.loc[spec['adult_section'],'ap_score']:+.2f})   gap = {spec['ap_gap']}")
    print(f"\n  For scale: the full 12-section AP spread is "
          f"{ap.ap_score.max()-ap.ap_score.min():.2f}, and a typical adjacent pair "
          f"within one animal differs by "
          f"{np.median([abs(np.diff(sorted(s.ap_score))).mean() for _, s in ap.groupby('animal')]):.2f}.")

    rows = []
    for label, cell_type, gene, _ in CANDIDATES:
        allsec = measure(adata, cell_type, gene)
        if allsec.empty:
            print(f"\n{label}: cell type '{cell_type}' not present, skipped")
            continue
        matched_sections = [s for spec in pairs.values()
                            for s in (spec["aged_section"], spec["adult_section"])]
        match = measure(adata, cell_type, gene, matched_sections)

        for metric, prefix in (("counts per million", "cpm"), ("detection %", "det")):
            l_all, l_match = lfc(allsec, prefix), lfc(match, prefix)
            rows.append({
                "candidate": label, "cell_type": cell_type, "gene": gene, "metric": metric,
                "all_B1": round(l_all.get("B1", np.nan), 2),
                "all_B2": round(l_all.get("B2", np.nan), 2),
                "matched_B1": round(l_match.get("B1", np.nan), 2),
                "matched_B2": round(l_match.get("B2", np.nan), 2),
                "all_consistent": np.sign(l_all.get("B1", 0)) == np.sign(l_all.get("B2", 0)),
                "matched_consistent": np.sign(l_match.get("B1", 0)) == np.sign(l_match.get("B2", 0)),
            })

    res = pd.DataFrame(rows)
    res.to_csv(OUT / "matched_vs_all.csv", index=False)
    print("\n=== Effect with all sections vs AP-matched pairs (log2 aged/adult) ===")
    with pd.option_context("display.width", 240):
        print(res.to_string(index=False))

    print("\n=== Reading ===")
    for _, r in res[res.metric == "counts per million"].iterrows():
        b1_keeps = np.sign(r.all_B1) == np.sign(r.matched_B1)
        b2_keeps = np.sign(r.all_B2) == np.sign(r.matched_B2)
        verdict = ("holds in both blocks" if b1_keeps and b2_keeps else
                   "holds in B2 only" if b2_keeps else
                   "holds in B1 only" if b1_keeps else "reverses in both")
        print(f"  {r.candidate:24s} {r.gene:6s}  all ({r.all_B1:+.2f}, {r.all_B2:+.2f}) -> "
              f"matched ({r.matched_B1:+.2f}, {r.matched_B2:+.2f})   {verdict}")
    print("\n  B2 (G_073 vs M399) is matched to within 0.05 AP units and is the one to trust.")
    print("  B1 (F536 vs M493) cannot be matched: F536's block was cut at a different level.")

    _plot(ap, pairs)
    print(f"\nWrote {OUT}")
    return 0


def _plot(ap, pairs) -> None:
    fig, ax = plt.subplots(figsize=(9, 4.5))
    colours = {"aged": "#B4531A", "adult": "#0E5A61"}
    order = ["F536", "M493", "G_073", "M399"]
    chosen = {s for spec in pairs.values() for s in (spec["aged_section"], spec["adult_section"])}
    for y, animal in enumerate(order):
        sub = ap[ap.animal == animal]
        for sec, r in sub.iterrows():
            picked = sec in chosen
            ax.scatter(r["ap_score"], y, s=190 if picked else 60,
                       color=colours[r["group"]], edgecolor="k" if picked else "none",
                       linewidth=1.8 if picked else 0, zorder=3 if picked else 2)
            if picked:
                ax.annotate(sec, (r["ap_score"], y), textcoords="offset points",
                            xytext=(0, 13), ha="center", fontsize=8)
    for block, spec in pairs.items():
        xa = ap.loc[spec["aged_section"], "ap_score"]
        xd = ap.loc[spec["adult_section"], "ap_score"]
        ya = order.index(BLOCKS[block][0]); yd = order.index(BLOCKS[block][1])
        ax.plot([xa, xd], [ya, yd], color="#555", lw=1.6, ls="--", zorder=1)
        ax.annotate(f"{block}: gap {spec['ap_gap']}", ((xa + xd) / 2, (ya + yd) / 2),
                    fontsize=9, ha="center", va="bottom",
                    bbox=dict(boxstyle="round,pad=.2", fc="white", ec="none", alpha=.8))
    ax.set_yticks(range(4)); ax.set_yticklabels(order)
    ax.set_xlabel("morphometric AP score"); ax.grid(axis="x", alpha=.25)
    ax.set_title("AP-matched pairing within each block\n"
                 "B2 pairs almost exactly; B1 cannot be matched")
    fig.tight_layout()
    fig.savefig(OUT / "matched_pairs.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
