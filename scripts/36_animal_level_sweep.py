#!/usr/bin/env python3
"""Re-run the expression screen with the animal as the unit of replication.

Scripts 24, 26, 27 and 28 all screened on **complete separation of the six aged
sections from the six adult sections**.  That criterion is inconsistent with the
statistical position taken everywhere else here: three sections from one mouse
are pseudoreplicates, and between-section SD in this dataset (0.265 log2) is
larger than between-animal SD (0.152).  A rule that lets one section in twelve
veto a result all four animals agree on is measuring section variance.

Here the criterion is separation at the **animal** level: pseudobulk each
animal's cells in a population, and ask whether both aged animals fall on one
side of both adults.  With 2 versus 2 there are six ways to assign the labels, so
a separating gene has exact one-sided p = 1/6 = 0.167 -- the floor.  Nothing in
this design can do better, and that is a fact about the design, not the gene.

Because 1/6 of exchangeable genes separate by chance, separation alone is close
to meaningless.  Two things make it informative:

  * **panel calibration** -- the same test is applied to all 297 genes in the
    same population, so every hit carries the observed rate.  A population where
    40% of the panel separates is telling you about the population.
  * **effect size rank** -- a gene's |log2 FC| percentile among the expressed
    panel in that population.  Separating *and* being in the top few percent of
    effects is a much rarer joint event than separating alone.

Both criteria are computed on the same populations, so the last table answers
the question that prompted this: what did the section-level rule throw away?

    python scripts/36_animal_level_sweep.py
    python scripts/36_animal_level_sweep.py --min-cells 25
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.io import ANIMAL_META, SECTION_ANIMAL, counts_matrix  # noqa: E402

PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "animal_level"
NUCLEI = ("ARC", "ME_3V", "VMH", "DMH", "LHA", "ZI", "DHA_PH")
ANIMALS = ["F536", "G_073", "M493", "M399"]
AGED = ["F536", "G_073"]
ADULT = ["M493", "M399"]
BLOCKS = {"B1": ("F536", "M493"), "B2": ("G_073", "M399")}
FOCUS = ["Gal", "Galr1", "Galr3"]
MIN_COUNTS = 200          # gene must be expressed in the population
MIN_PER_SECTION = 15      # for the section-level comparison only
TOP_PCT = 95.0            # effect-size percentile that makes a separator notable


def pseudobulk(counts: np.ndarray, groups: np.ndarray, order: list) -> np.ndarray:
    """log2 CPM per group, rows in `order`."""
    mat = np.vstack([counts[groups == g].sum(axis=0) for g in order])
    return np.log2(mat / mat.sum(axis=1, keepdims=True) * 1e6 + 1)


def separates(values: np.ndarray, is_aged: np.ndarray) -> np.ndarray:
    a, d = values[is_aged], values[~is_aged]
    return (a.min(axis=0) > d.max(axis=0)) | (a.max(axis=0) < d.min(axis=0))


def margin(values: np.ndarray, is_aged: np.ndarray) -> np.ndarray:
    """Gap in log2 between the closest aged and closest adult animal.

    Separation is a yes/no, and at n=2 versus 2 a gene can clear it by 0.005
    log2 as easily as by 2.  The margin says which, and it is the difference
    between a result and a coin landing on its edge.
    """
    a, d = values[is_aged], values[~is_aged]
    return np.maximum(a.min(axis=0) - d.max(axis=0), d.min(axis=0) - a.max(axis=0))


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--min-cells", type=int, default=40,
                   help="minimum cells per animal for a population to be tested")
    args = p.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(PROC / "hypothalamus_nuclei.h5ad")
    adata = adata[adata.obs["nucleus_ext"].isin(NUCLEI)].copy()
    counts = counts_matrix(adata)
    genes = adata.var_names.to_numpy()
    animals = adata.obs["animal"].astype(str).to_numpy()
    sections = adata.obs["section"].astype(str).to_numpy()
    pops = (adata.obs["nucleus_ext"].astype(str) + " | "
            + adata.obs["cell_type"].astype(str)).to_numpy()

    animal_is_aged = np.array([a in AGED for a in ANIMALS])
    rows, pop_rows = [], []

    for pop in pd.unique(pops):
        if "unlabelled" in pop:
            continue
        sel = pops == pop
        n_per_animal = {a: int((sel & (animals == a)).sum()) for a in ANIMALS}
        if min(n_per_animal.values()) < args.min_cells:
            continue

        sub, sub_animals, sub_sections = counts[sel], animals[sel], sections[sel]
        cpm_a = pseudobulk(sub, sub_animals, ANIMALS)
        expressed = sub.sum(axis=0) >= MIN_COUNTS
        if expressed.sum() < 30:
            continue

        sep_a = separates(cpm_a, animal_is_aged)
        marg_a = margin(cpm_a, animal_is_aged)
        lfc = cpm_a[animal_is_aged].mean(axis=0) - cpm_a[~animal_is_aged].mean(axis=0)
        # Sign agreement within each batch pair, reported but not required:
        # full separation over four animals already implies it.
        block_ok = np.ones(len(genes), dtype=bool)
        for a, d in BLOCKS.values():
            block_ok &= np.sign(cpm_a[ANIMALS.index(a)] - cpm_a[ANIMALS.index(d)]) \
                == np.sign(lfc)

        # The old criterion, on the same populations, for the comparison.
        keep_s = [s for s in pd.unique(sub_sections)
                  if (sub_sections == s).sum() >= MIN_PER_SECTION]
        sep_s = np.zeros(len(genes), dtype=bool)
        if len(keep_s) == 12:
            cpm_s = pseudobulk(sub, sub_sections, keep_s)
            s_aged = np.array([SECTION_ANIMAL[s] in AGED for s in keep_s])
            s_animals = np.array([SECTION_ANIMAL[s] for s in keep_s])
            sep_s = separates(cpm_s, s_aged)
            for a, d in BLOCKS.values():
                va, vd = cpm_s[s_animals == a], cpm_s[s_animals == d]
                up = (va.min(axis=0) > vd.max(axis=0)) & (cpm_s[s_aged].min(axis=0)
                                                          > cpm_s[~s_aged].max(axis=0))
                dn = (va.max(axis=0) < vd.min(axis=0)) & (cpm_s[s_aged].max(axis=0)
                                                          < cpm_s[~s_aged].min(axis=0))
                sep_s &= up | dn

        rank = pd.Series(np.abs(lfc)[expressed]).rank(pct=True) * 100
        rank_full = np.full(len(genes), np.nan)
        rank_full[expressed] = rank.to_numpy()

        n_exp = int(expressed.sum())
        pop_rows.append({
            "population": pop, "min_cells_per_animal": min(n_per_animal.values()),
            "genes_tested": n_exp,
            "animal_sep_rate_pct": round((sep_a & expressed).sum() / n_exp * 100, 1),
            "section_sep_rate_pct": round((sep_s & expressed).sum() / n_exp * 100, 1),
        })
        for i, g in enumerate(genes):
            if not expressed[i]:
                continue
            rows.append({
                "population": pop, "gene": g,
                "min_cells_per_animal": min(n_per_animal.values()),
                "lfc": round(float(lfc[i]), 3),
                "margin": round(float(marg_a[i]), 3),
                "lfc_pct_rank": round(float(rank_full[i]), 1),
                "animal_separates": bool(sep_a[i]),
                "section_separates": bool(sep_s[i]),
                "blocks_agree": bool(block_ok[i]),
                "genes_tested": n_exp,
                **{f"cpm_{a}": round(float(cpm_a[j, i]), 2)
                   for j, a in enumerate(ANIMALS)},
            })

    res = pd.DataFrame(rows)
    pops_df = pd.DataFrame(pop_rows).sort_values("animal_sep_rate_pct")
    res.to_csv(OUT / "animal_level_sweep.csv", index=False)
    pops_df.to_csv(OUT / "population_calibration.csv", index=False)

    print(f"{pops_df.shape[0]} populations with >={args.min_cells} cells in all "
          f"four animals; {len(res):,} population x gene tests\n")

    print("=== Calibration: how often does each criterion fire? ===")
    print(f"  animal level   {res.animal_separates.mean()*100:5.1f}% of tests   "
          f"(chance for exchangeable genes is 1/6 = 16.7%)")
    print(f"  section level  {res.section_separates.mean()*100:5.1f}% of tests")
    print(f"  the section rule is {res.animal_separates.mean()/max(res.section_separates.mean(),1e-9):.1f}x "
          "stricter than the animal rule")

    strong = res[res.animal_separates & (res.lfc_pct_rank >= TOP_PCT)]
    MIN_MARGIN = 0.25   # ~19%; below this the four animals are effectively tied
    print(f"\n  separating AND in the top {100-TOP_PCT:.0f}% of effects: "
          f"{len(strong)} of {len(res):,} tests "
          f"({len(strong)/len(res)*100:.2f}%; chance ~{16.7*(100-TOP_PCT)/100:.2f}%)")

    print("\n=== What the section rule was discarding ===")
    recovered = res[res.animal_separates & ~res.section_separates
                    & (res.lfc_pct_rank >= TOP_PCT)]
    print(f"  {len(recovered)} strong animal-level hits that the section rule rejected")
    print(f"  {int((res.section_separates & ~res.animal_separates).sum())} section-level "
          "hits that do NOT separate at the animal level (these were false leads)")

    print(f"\n  of those, {int((strong.margin >= MIN_MARGIN).sum())} clear a separation "
          f"margin of {MIN_MARGIN} log2; the rest separate by a hair and are ties")

    with pd.option_context("display.width", 210, "display.max_colwidth", 30):
        print("\n=== Top 25 hits: strong effect AND real daylight between the groups ===")
        top = strong[strong.margin >= MIN_MARGIN]
        top = top.reindex(top["lfc"].abs().sort_values(ascending=False).index)
        print(top.head(25)[["population", "gene", "lfc", "margin", "lfc_pct_rank",
                            "min_cells_per_animal", "section_separates"]]
              .to_string(index=False))

    print("\n=== Galanin system specifically ===")
    foc = res[res.gene.isin(FOCUS) & res.animal_separates]
    foc = foc.reindex(foc["lfc"].abs().sort_values(ascending=False).index)
    if len(foc):
        with pd.option_context("display.width", 210, "display.max_colwidth", 30):
            print(foc[["population", "gene", "lfc", "margin", "lfc_pct_rank",
                       "min_cells_per_animal", "section_separates",
                       "cpm_F536", "cpm_G_073", "cpm_M493", "cpm_M399"]]
                  .to_string(index=False))
    else:
        print("  none separate at the animal level")

    print("\n=== Populations with the cleanest calibration (lowest chance rate) ===")
    print(pops_df.head(10).to_string(index=False))
    print(f"\nWrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
