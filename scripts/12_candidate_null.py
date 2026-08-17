#!/usr/bin/env python3
"""Calibrate the surviving candidates against the whole panel.

11_galanin_search.py applied three filters in sequence: consistent direction in
both blocks, at least 1.5x in the weaker block, and complete separation of aged
from adult sections within each block.  Five galanin-system candidates survived.

The question this script answers is whether surviving all three is rare.  The
same three filters are applied to all 297 panel genes in the same populations,
so the survivor count has a denominator.  A candidate that is one of three
genes out of 297 means something; one of ninety does not.

It also checks each survivor against the rostro-caudal composition axis, which
is what dismantled the Galr1 lead in 10_galr1_nucleus.py.
"""

from __future__ import annotations

import sys
import re
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats as st

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.io import ANIMAL_META, SECTION_ANIMAL, counts_matrix  # noqa: E402

PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "galanin_search"
BLOCKS = {"B1": ("F536", "M493"), "B2": ("G_073", "M399")}
ANIMALS = ["F536", "G_073", "M493", "M399"]
LFC_THRESHOLD = 0.58
MIN_CELLS_PER_SECTION = 15

# The populations whose galanin-system candidates survived all three filters.
POPULATIONS = {
    "ARC | ARC Agrp/Npy": ("cell_type", "ARC", "ARC Agrp/Npy"),
    "DMH | GABA Cacna2d2": ("cell_type", "DMH", "GABA Cacna2d2"),
    "DMH | ml:lat dv:dors": ("spatial", "DMH", ("lat", "dors")),
}
FOCUS = {"ARC | ARC Agrp/Npy": "Gal", "DMH | GABA Cacna2d2": "Gal",
         "DMH | ml:lat dv:dors": "Galr3"}


def _slug(name: str) -> str:
    """Population labels contain '/' and ':', neither safe in a filename."""
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")


def population_mask(adata, spec) -> np.ndarray:
    kind, nucleus, detail = spec
    in_nuc = (adata.obs["nucleus"].astype(str) == nucleus).to_numpy()
    if kind == "cell_type":
        return in_nuc & (adata.obs["cell_type"].astype(str) == detail).to_numpy()
    sub = adata.obs[in_nuc]
    ml = pd.qcut(sub["ml"].abs(), 3, labels=["med", "mid", "lat"], duplicates="drop")
    dv = pd.qcut(sub["dv"], 3, labels=["vent", "mid", "dors"], duplicates="drop")
    keep = (ml.astype(str) == detail[0]) & (dv.astype(str) == detail[1])
    mask = np.zeros(adata.n_obs, dtype=bool)
    mask[np.where(in_nuc)[0][keep.to_numpy()]] = True
    return mask


def gene_table(adata, mask) -> tuple[pd.DataFrame, pd.DataFrame]:
    """CPM per animal and per section for every panel gene in one population."""
    counts = counts_matrix(adata[mask])
    obs = adata.obs[mask]
    genes = adata.var_names.to_numpy()

    def agg(key):
        idx = obs[key].astype(str).to_numpy()
        codes, uniq = pd.factorize(idx)
        summed = np.zeros((len(uniq), counts.shape[1]))
        np.add.at(summed, codes, counts)
        n = np.bincount(codes, minlength=len(uniq))
        cpm = pd.DataFrame(summed / summed.sum(axis=1, keepdims=True) * 1e6,
                           index=uniq, columns=genes)
        return cpm, pd.Series(n, index=uniq)

    per_animal, _ = agg("animal")
    per_section, n_section = agg("section")
    return per_animal, per_section[n_section >= MIN_CELLS_PER_SECTION]


def apply_filters(per_animal: pd.DataFrame, per_section: pd.DataFrame) -> pd.DataFrame:
    lfc = {b: np.log2((per_animal.loc[a] + 1) / (per_animal.loc[d] + 1))
           for b, (a, d) in BLOCKS.items()}
    lfc_b1, lfc_b2 = lfc["B1"].to_numpy(), lfc["B2"].to_numpy()

    animals = pd.Series([SECTION_ANIMAL[s] for s in per_section.index], index=per_section.index)
    groups = animals.map(lambda a: ANIMAL_META[a]["age_group"])
    aged = per_section[groups == "aged"]
    adult = per_section[groups == "adult"]

    sep_overall = (aged.max() < adult.min()) | (aged.min() > adult.max())
    sep_blocks = np.ones(per_section.shape[1], dtype=bool)
    for a, d in BLOCKS.values():
        va, vd = per_section[animals == a], per_section[animals == d]
        sep_blocks &= ((va.max() < vd.min()) | (va.min() > vd.max())).to_numpy()

    consistent = np.sign(lfc_b1) == np.sign(lfc_b2)
    min_abs = np.minimum(np.abs(lfc_b1), np.abs(lfc_b2))
    return pd.DataFrame({
        "lfc_B1": lfc_b1, "lfc_B2": lfc_b2, "mean_lfc": (lfc_b1 + lfc_b2) / 2,
        "consistent": consistent, "min_abs_lfc": min_abs,
        "sep_overall": sep_overall.to_numpy(), "sep_both_blocks": sep_blocks,
        "passes_all": consistent & (min_abs >= LFC_THRESHOLD) & sep_blocks,
    }, index=per_animal.columns)


def ap_axis(adata) -> pd.Series:
    comp = pd.crosstab(adata.obs["section"], adata.obs["cell_type"], normalize="index")
    z = ((comp - comp.mean()) / comp.std().replace(0, 1)).fillna(0)
    u, s, _ = np.linalg.svd(z.to_numpy(), full_matrices=False)
    return pd.Series(u[:, 0] * s[0], index=comp.index)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(PROC / "mbh_atlas_12.h5ad")
    adata = adata[adata.obs["cell_type"] != "Unresolved"].copy()
    pc1 = ap_axis(adata)

    summary = []
    for name, spec in POPULATIONS.items():
        mask = population_mask(adata, spec)
        per_animal, per_section = gene_table(adata, mask)
        if not set(ANIMALS).issubset(per_animal.index):
            print(f"{name}: missing an animal, skipped")
            continue
        res = apply_filters(per_animal, per_section)
        res.to_csv(OUT / f"panel_null_{_slug(name)}.csv")

        gene = FOCUS[name]
        n_pass = int(res["passes_all"].sum())
        rank = int((res["min_abs_lfc"] >= res.loc[gene, "min_abs_lfc"]).sum())
        print(f"\n=== {name}  (focus gene: {gene}) ===")
        print(f"  cells: {int(mask.sum()):,} | sections used: {len(per_section)}")
        print(f"  genes passing ALL three filters: {n_pass} of {len(res)} "
              f"({n_pass/len(res)*100:.1f}% of the panel)")
        print(f"  {gene}: lfc {res.loc[gene,'lfc_B1']:+.2f} / {res.loc[gene,'lfc_B2']:+.2f}, "
              f"rank {rank}/{len(res)} by effect size, passes all = {bool(res.loc[gene,'passes_all'])}")
        others = res[res["passes_all"]].sort_values("min_abs_lfc", ascending=False)
        print(f"  other genes passing: {', '.join(others.index[:14])}"
              f"{' ...' if len(others) > 14 else ''}")

        # Is the effect explained by rostro-caudal position?
        vals = per_section[gene]
        common = vals.index.intersection(pc1.index)
        r, p = st.pearsonr(pc1[common], vals[common])
        groups = pd.Series([ANIMAL_META[SECTION_ANIMAL[s]]["age_group"] for s in common], index=common)
        within = []
        for animal in ANIMALS:
            sel = [s for s in common if SECTION_ANIMAL[s] == animal]
            if len(sel) >= 3:
                within.append(np.polyfit(pc1[sel], vals[sel], 1)[0])
        gap = vals[groups == "aged"].mean() - vals[groups == "adult"].mean()
        swing = max(
            vals[[s for s in common if SECTION_ANIMAL[s] == a]].max()
            - vals[[s for s in common if SECTION_ANIMAL[s] == a]].min()
            for a in ANIMALS
        )
        print(f"  vs rostro-caudal axis: r = {r:+.2f} (p = {p:.3f}); "
              f"aged-adult gap {gap:+.0f} CPM vs largest within-animal swing {swing:.0f} CPM")
        summary.append({
            "population": name, "gene": gene, "n_cells": int(mask.sum()),
            "panel_pass_n": n_pass, "panel_pass_pct": round(n_pass / len(res) * 100, 1),
            "effect_rank": rank, "lfc_B1": round(float(res.loc[gene, "lfc_B1"]), 2),
            "lfc_B2": round(float(res.loc[gene, "lfc_B2"]), 2),
            "r_vs_AP": round(float(r), 2), "p_vs_AP": round(float(p), 3),
            "gap_cpm": round(float(gap)), "max_within_animal_swing_cpm": round(float(swing)),
            "gap_exceeds_swing": bool(abs(gap) > swing),
        })

    out = pd.DataFrame(summary)
    out.to_csv(OUT / "candidate_calibration.csv", index=False)
    print("\n=== Summary ===")
    with pd.option_context("display.width", 250):
        print(out.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
