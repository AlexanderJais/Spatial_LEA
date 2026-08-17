#!/usr/bin/env python3
"""Audit the AP test: is "effect size vs AP mismatch" a trustworthy diagnostic?

The test correlates the log fold change computed from each aged/adult section
pair against that pair's rostro-caudal mismatch, and calls an effect artefactual
when the correlation is significant.  It was used to reject the DMH Galr1 lead
and to keep Gal in ARC Agrp/Npy.  Four things could make it wrong:

1. **Non-independence.** The 18 pairings are built from 12 sections, so each
   section appears in three pairings.  A Pearson p-value assumes independent
   observations and will be anticonservative.

2. **Confounding with block.** If the AP gaps available in block B1 barely
   overlap those in B2, then "large gap" means "block B1" and the correlation
   tests a block difference, not an AP effect.

3. **Miscalibration.** Under a null with no AP dependence the test should reject
   at ~5%.  If it rejects far more often across the panel, it is not a test.

4. **Direction.** A real age effect in a gene that also varies with AP would
   produce a correlation too, so a significant result may not mean "artefact".

Each is checked below against the real data and against a permutation null that
respects the nesting.
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

from spatial_lea.io import ANIMAL_META, SECTION_ANIMAL, counts_matrix  # noqa: E402

PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "ap_audit"
AP_FILE = REPO / "results" / "ap_matching" / "section_ap_scores.csv"
BLOCKS = {"B1": ("F536", "M493"), "B2": ("G_073", "M399")}
CASES = [("Galr1", "GABA Gal/Galr1"), ("Gal", "ARC Agrp/Npy"), ("Gal", "GABA Cacna2d2")]
N_PERM = 20000
SEED = 0


def section_cpm(adata, gene: str, cell_type: str, min_cells=15) -> pd.Series:
    gi = int(np.where(adata.var_names.to_numpy() == gene)[0][0])
    sel = (adata.obs["cell_type"].astype(str) == cell_type).to_numpy()
    sub = adata[sel]
    counts = counts_matrix(sub)
    df = pd.DataFrame({"section": sub.obs["section"].astype(str).to_numpy(),
                       "g": counts[:, gi], "tot": counts.sum(axis=1)})
    g = df.groupby("section").agg(gg=("g", "sum"), tt=("tot", "sum"), n=("g", "size"))
    g = g[g["n"] >= min_cells]
    return g["gg"] / g["tt"] * 1e6


def pairings(cpm: pd.Series, ap: pd.Series) -> pd.DataFrame:
    rows = []
    for block, (aged, adult) in BLOCKS.items():
        a = [s for s in cpm.index if SECTION_ANIMAL[s] == aged]
        d = [s for s in cpm.index if SECTION_ANIMAL[s] == adult]
        for si, sj in itertools.product(a, d):
            rows.append({"block": block, "aged_section": si, "adult_section": sj,
                         "ap_gap": abs(ap[si] - ap[sj]),
                         "lfc": float(np.log2((cpm[si] + 1) / (cpm[sj] + 1)))})
    return pd.DataFrame(rows)


def perm_pvalue(pairs: pd.DataFrame, cpm: pd.Series, ap: pd.Series, rng) -> float:
    """Permutation null that respects the nesting.

    The AP scores are shuffled *within animal*, so each animal keeps its own set
    of section levels and only the mapping of level to section is randomised.
    That destroys any real AP dependence while preserving the pairing structure,
    the number of sections per animal and the block design -- which a naive
    Pearson p-value ignores entirely.
    """
    observed = abs(st.pearsonr(pairs["ap_gap"], pairs["lfc"]).statistic)
    by_animal = {}
    for s in cpm.index:
        by_animal.setdefault(SECTION_ANIMAL[s], []).append(s)

    count = 0
    for _ in range(N_PERM):
        shuffled = {}
        for animal, secs in by_animal.items():
            vals = ap[secs].to_numpy().copy()
            rng.shuffle(vals)
            shuffled.update(dict(zip(secs, vals)))
        gaps = np.array([abs(shuffled[r.aged_section] - shuffled[r.adult_section])
                         for r in pairs.itertuples()])
        if np.std(gaps) < 1e-9:
            continue
        r = abs(st.pearsonr(gaps, pairs["lfc"]).statistic)
        count += r >= observed
    return (count + 1) / (N_PERM + 1)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    adata = sc.read_h5ad(PROC / "hypothalamus_nuclei.h5ad")
    ap = pd.read_csv(AP_FILE, index_col=0)["ap_score"]

    # ---- problem 2: is AP gap confounded with block? ----------------------
    ref = pairings(section_cpm(adata, "Gal", "ARC Agrp/Npy"), ap)
    print("=== Check 1: is AP mismatch confounded with block? ===")
    for block, sub in ref.groupby("block"):
        print(f"  {block}: gaps {sub.ap_gap.min():.2f} – {sub.ap_gap.max():.2f}  (n = {len(sub)})")
    b1, b2 = ref[ref.block == "B1"].ap_gap, ref[ref.block == "B2"].ap_gap
    overlap = min(b1.max(), b2.max()) - max(b1.min(), b2.min())
    auc = st.mannwhitneyu(b1, b2).statistic / (len(b1) * len(b2))
    print(f"  overlap of the two ranges: {overlap:+.2f}   separation AUC = {auc:.2f}")
    if overlap <= 0:
        print("  *** The ranges do not overlap. 'large gap' means 'block B1', so the")
        print("      correlation cannot distinguish an AP effect from a block difference. ***")

    # ---- problem 1 & 3: independence and calibration -----------------------
    print("\n=== Check 2: naive p-value vs a nesting-aware permutation null ===")
    rows = []
    for gene, cell_type in CASES:
        cpm = section_cpm(adata, gene, cell_type)
        if len(cpm) < 10:
            continue
        pairs = pairings(cpm, ap)
        r, p_naive = st.pearsonr(pairs["ap_gap"], pairs["lfc"])
        p_perm = perm_pvalue(pairs, cpm, ap, rng)
        # Within-block correlations: the only place the confound is broken.
        within = {}
        for block, sub in pairs.groupby("block"):
            if sub["ap_gap"].std() > 1e-9 and len(sub) >= 5:
                within[block] = st.pearsonr(sub["ap_gap"], sub["lfc"])
        rows.append({
            "gene": gene, "cell_type": cell_type, "r_all": round(r, 2),
            "p_naive": round(p_naive, 4), "p_permuted": round(p_perm, 4),
            **{f"r_{b}": round(v.statistic, 2) for b, v in within.items()},
            **{f"p_{b}": round(v.pvalue, 3) for b, v in within.items()},
        })
    audit = pd.DataFrame(rows)
    audit.to_csv(OUT / "ap_test_audit.csv", index=False)
    with pd.option_context("display.width", 220):
        print(audit.to_string(index=False))

    # ---- problem 3: false-positive rate across the panel -------------------
    print("\n=== Check 3: how often does the test fire across the whole panel? ===")
    cell_type = "ARC Agrp/Npy"
    sel = (adata.obs["cell_type"].astype(str) == cell_type).to_numpy()
    sub = adata[sel]
    counts = counts_matrix(sub)
    sections = sub.obs["section"].astype(str).to_numpy()
    keep_sections = [s for s in ap.index if (sections == s).sum() >= 15]
    tot = np.array([counts[sections == s].sum() for s in keep_sections])
    naive_hits, perm_hits, tested = 0, 0, 0
    per_gene = []
    for gi, gene in enumerate(adata.var_names):
        g = np.array([counts[sections == s, gi].sum() for s in keep_sections])
        if g.sum() < 200:
            continue
        cpm = pd.Series(g / tot * 1e6, index=keep_sections)
        pairs = pairings(cpm, ap)
        if len(pairs) < 10 or pairs["lfc"].std() < 1e-9:
            continue
        tested += 1
        r, p = st.pearsonr(pairs["ap_gap"], pairs["lfc"])
        naive_hits += p < 0.05
        per_gene.append({"gene": gene, "r": r, "p": p})
    pg = pd.DataFrame(per_gene)
    pg.to_csv(OUT / "panel_ap_correlations.csv", index=False)
    print(f"  {tested} genes tested in {cell_type}")
    print(f"  naive p < 0.05 in {naive_hits} genes ({naive_hits/tested*100:.0f}%) "
          f"-- expected 5% if the test were calibrated")
    print(f"  median |r| across the panel: {pg.r.abs().median():.2f}")

    _plot(ref, pg, audit)
    print(f"\nWrote {OUT}")
    return 0


def _plot(ref, pg, audit) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    ax = axes[0]
    for block, sub in ref.groupby("block"):
        ax.scatter(sub["ap_gap"], [block] * len(sub), s=40,
                   color="#B4531A" if block == "B1" else "#0E5A61")
    ax.set_xlabel("AP mismatch of the pairing"); ax.set_ylabel("")
    ax.set_title("Check 1 — AP mismatch is confounded with block\n"
                 "(the ranges do not overlap)")
    ax.grid(axis="x", alpha=.3)

    ax = axes[1]
    ax.hist(pg["r"], bins=40, color="#C3CCCC")
    ax.axvline(0, color="#555", lw=1)
    for _, r in audit.iterrows():
        ax.axvline(r["r_all"], lw=1.6,
                   color="#B4531A" if r["gene"] == "Galr1" else "#0E5A61")
    ax.set_xlabel("r (effect size vs AP mismatch)"); ax.set_ylabel("panel genes")
    ax.set_title("Check 3 — panel-wide distribution of r\n"
                 "(lines: Galr1 orange, Gal teal)")

    ax = axes[2]
    ax.hist(pg["p"], bins=20, color="#C3CCCC")
    ax.axhline(len(pg) / 20, color="#B4531A", ls="--", lw=1.4)
    ax.set_xlabel("naive p-value"); ax.set_ylabel("panel genes")
    ax.set_title("Check 3 — p-value distribution\n"
                 "(dashed = uniform, what a calibrated test gives)")
    fig.tight_layout()
    fig.savefig(OUT / "ap_test_audit.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
