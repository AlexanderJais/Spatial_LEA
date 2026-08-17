#!/usr/bin/env python3
"""Unbiased search: is ANY Gal or Galr1 population age-regulated?

The targeted analyses looked at pre-named cell types.  This searches wider and,
more importantly, calibrates what it finds.

Search space -- four ways of carving up the MBH, so an effect that is not
aligned to the cell-type labels can still surface:

  A  annotated cell type x nucleus
  B  fine Leiden subclusters (resolution 3.0), which split the named types
  C  the Gal+ and Galr1+ cells themselves, subclustered on their own
  D  spatial bins: nucleus x mediolateral tertile x dorsoventral tertile

Calibration -- the part that makes the answer meaningful.  For every population
the same statistic is computed for **all 297 panel genes**, so ``Galr1`` is
ranked against a null built from the same cells, the same animals and the same
blocked design.  A gene that is "consistent in both blocks" is unremarkable if
half the panel also is; the panel-wide pass rate says how unremarkable.

Robustness -- with three sections per animal, a candidate must also separate at
section level, which is what the rostro-caudal gradient broke last time.  That
is a reproducibility filter, not a p-value: sections within an animal are not
independent, so it cannot rescue n=2 per group.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.io import ANIMAL_META, SECTION_ANIMAL, counts_matrix  # noqa: E402

PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "galanin_search"
TARGETS = ["Galr1", "Galr3", "Gal"]
BLOCKS = {"B1": ("F536", "M493"), "B2": ("G_073", "M399")}
AGED = ["F536", "G_073"]
ADULT = ["M493", "M399"]
ANIMALS = AGED + ADULT

MIN_CELLS_PER_ANIMAL = 40
MIN_TOTAL_COUNTS = 200      # per animal per population, for the gene being tested
LFC_THRESHOLD = 0.58        # 1.5x in the weaker block
SEED = 0


def build_populations(adata) -> dict[str, pd.Series]:
    """Return {scheme: Series of population label per cell}."""
    obs = adata.obs
    schemes: dict[str, pd.Series] = {}

    schemes["A_celltype_nucleus"] = (
        obs["nucleus"].astype(str) + " | " + obs["cell_type"].astype(str)
    )

    sc.tl.leiden(adata, resolution=3.0, key_added="leiden_fine",
                 flavor="igraph", n_iterations=2, random_state=SEED)
    fine = obs["leiden_fine"].astype(str)
    dominant = (
        pd.crosstab(fine, obs["cell_type"]).idxmax(axis=1).astype(str)
    )
    schemes["B_fine_cluster"] = pd.Series(
        [f"c{c} ({dominant[c]})" for c in fine], index=obs.index
    )

    # C: the receptor- and ligand-positive cells in their own right.
    for gene in ("Galr1", "Gal"):
        pos = counts_matrix(adata[:, [gene]]).ravel() > 0
        lab = pd.Series(np.where(pos, f"{gene}+", "other"), index=obs.index)
        sub = obs["nucleus"].astype(str) + " | " + lab
        schemes[f"C_{gene}pos"] = sub.where(pos, other=np.nan)

    # D: geometry, independent of any clustering.
    bins = []
    for nucleus, sub in obs.groupby("nucleus", observed=True):
        ml = pd.qcut(sub["ml"].abs(), 3, labels=["med", "mid", "lat"], duplicates="drop")
        dv = pd.qcut(sub["dv"], 3, labels=["vent", "mid", "dors"], duplicates="drop")
        bins.append(pd.Series(f"{nucleus} | ml:" + ml.astype(str) + " dv:" + dv.astype(str),
                              index=sub.index))
    schemes["D_spatial_bin"] = pd.concat(bins).reindex(obs.index)
    return schemes


def pseudobulk(counts: np.ndarray, groups: pd.Series, animals: pd.Series):
    """Sum counts per (population, animal); returns dict pop -> DataFrame."""
    # Factorising a MultiIndex rather than a joined string: population labels
    # contain arbitrary punctuation, so any separator character is a hazard.
    pairs = pd.MultiIndex.from_arrays(
        [groups.astype(str).to_numpy(), animals.astype(str).to_numpy()]
    )
    codes, uniq = pd.factorize(pairs)
    summed = np.zeros((len(uniq), counts.shape[1]))
    np.add.at(summed, codes, counts)
    sizes = np.bincount(codes, minlength=len(uniq))
    pops = pd.Series([u[0] for u in uniq])
    anim = pd.Series([u[1] for u in uniq])
    return pd.DataFrame(summed), pops, anim, pd.Series(sizes)


def search(adata, schemes) -> pd.DataFrame:
    genes = adata.var_names.to_numpy()
    counts = counts_matrix(adata)
    animals = adata.obs["animal"].astype(str)

    rows = []
    for scheme, labels in schemes.items():
        valid = labels.notna()
        summed, pops, anim, sizes = pseudobulk(counts[valid.to_numpy()],
                                               labels[valid], animals[valid])
        for pop in pops.unique():
            sel = pops == pop
            if set(anim[sel]) != set(ANIMALS):
                continue
            n_by_animal = dict(zip(anim[sel], sizes[sel]))
            if min(n_by_animal.values()) < MIN_CELLS_PER_ANIMAL:
                continue
            mat = summed[sel.to_numpy()]
            mat.index = anim[sel].to_numpy()
            depth = mat.sum(axis=1)
            cpm = mat.div(depth, axis=0) * 1e6

            lfc = {}
            for block, (aged, adult) in BLOCKS.items():
                lfc[block] = np.log2((cpm.loc[aged] + 1) / (cpm.loc[adult] + 1))
            lfc_b1, lfc_b2 = lfc["B1"].to_numpy(), lfc["B2"].to_numpy()
            consistent = np.sign(lfc_b1) == np.sign(lfc_b2)
            min_abs = np.minimum(np.abs(lfc_b1), np.abs(lfc_b2))
            passes = consistent & (min_abs >= LFC_THRESHOLD)
            # The null: how much of the whole panel clears the same bar here.
            pass_rate = float(passes.mean())

            for gene in TARGETS:
                if gene not in genes:
                    continue
                gi = int(np.where(genes == gene)[0][0])
                if float(mat.iloc[:, gi].min()) * 0 + float(mat.iloc[:, gi].sum()) < MIN_TOTAL_COUNTS:
                    continue
                rank = float((min_abs >= min_abs[gi]).mean())  # 0 = biggest effect in panel
                rows.append({
                    "scheme": scheme, "population": pop, "gene": gene,
                    "n_min_animal": int(min(n_by_animal.values())),
                    "cpm_F536": round(float(cpm.loc["F536"].iloc[gi]), 1),
                    "cpm_G_073": round(float(cpm.loc["G_073"].iloc[gi]), 1),
                    "cpm_M493": round(float(cpm.loc["M493"].iloc[gi]), 1),
                    "cpm_M399": round(float(cpm.loc["M399"].iloc[gi]), 1),
                    "lfc_B1": round(float(lfc_b1[gi]), 2),
                    "lfc_B2": round(float(lfc_b2[gi]), 2),
                    "mean_lfc": round(float((lfc_b1[gi] + lfc_b2[gi]) / 2), 2),
                    "consistent": bool(consistent[gi]),
                    "min_abs_lfc": round(float(min_abs[gi]), 2),
                    "panel_pass_rate": round(pass_rate, 3),
                    "panel_rank_pct": round(rank * 100, 1),
                    "passes": bool(passes[gi]),
                })
    return pd.DataFrame(rows)


def section_separation(adata, candidates: pd.DataFrame, schemes) -> pd.DataFrame:
    """Do the aged and adult sections separate completely, in both blocks?

    Reproducibility check, not inference: sections inside one animal are not
    independent, so complete separation cannot substitute for biological n.
    """
    counts = counts_matrix(adata)
    genes = adata.var_names.to_numpy()
    sections = adata.obs["section"].astype(str)
    rows = []
    for _, cand in candidates.iterrows():
        labels = schemes[cand["scheme"]]
        sel = (labels == cand["population"]).fillna(False).to_numpy()
        gi = int(np.where(genes == cand["gene"])[0][0])
        df = pd.DataFrame({
            "section": sections[sel].to_numpy(),
            "g": counts[sel, gi],
            "tot": counts[sel].sum(axis=1),
        })
        per = df.groupby("section").agg(g=("g", "sum"), tot=("tot", "sum"), n=("g", "size"))
        per = per[per["n"] >= 15]
        per["cpm"] = per["g"] / per["tot"] * 1e6
        per["animal"] = [SECTION_ANIMAL[s] for s in per.index]
        per["group"] = per["animal"].map(lambda a: ANIMAL_META[a]["age_group"])
        if per["animal"].nunique() < 4:
            continue
        aged_v, adult_v = per[per.group == "aged"]["cpm"], per[per.group == "adult"]["cpm"]
        sep = aged_v.max() < adult_v.min() or aged_v.min() > adult_v.max()
        block_sep = all(
            (per[per.animal == a]["cpm"].max() < per[per.animal == b]["cpm"].min())
            or (per[per.animal == a]["cpm"].min() > per[per.animal == b]["cpm"].max())
            for a, b in BLOCKS.values()
        )
        rows.append({
            **cand[["scheme", "population", "gene", "mean_lfc", "panel_rank_pct",
                    "panel_pass_rate"]].to_dict(),
            "n_sections": len(per),
            "aged_range": f"{aged_v.min():.0f}-{aged_v.max():.0f}",
            "adult_range": f"{adult_v.min():.0f}-{adult_v.max():.0f}",
            "sections_separate": sep,
            "separate_within_both_blocks": block_sep,
        })
    return pd.DataFrame(rows)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    sc.settings.n_jobs = 4
    adata = sc.read_h5ad(PROC / "mbh_atlas_12.h5ad")
    adata = adata[adata.obs["cell_type"] != "Unresolved"].copy()
    print(f"{adata.n_obs:,} MBH cells, 12 sections\n")

    schemes = build_populations(adata)
    for name, labels in schemes.items():
        print(f"  {name:22s} {labels.dropna().nunique():4d} populations")

    res = search(adata, schemes)
    res.to_csv(OUT / "galanin_search_all.csv", index=False)
    print(f"\nTested {len(res)} (population, gene) combinations across "
          f"{res['population'].nunique()} populations.\n")

    print("=== How often does the WHOLE PANEL clear the same bar? ===")
    print("   (consistent in both blocks and >=1.5x in the weaker one)")
    rate = res.groupby("scheme")["panel_pass_rate"].agg(["mean", "min", "max"]).round(3)
    print(rate.to_string())
    print(f"\n   Panel-wide, {res['panel_pass_rate'].mean()*100:.0f}% of the 297 genes clear it "
          f"in a typical population -- so ~{res['panel_pass_rate'].mean()*297:.0f} genes per "
          f"population pass by chance alone.")

    hits = res[res["passes"]].sort_values("panel_rank_pct")
    print(f"\n=== Galanin-system hits clearing the bar: {len(hits)} ===")
    if len(hits):
        with pd.option_context("display.width", 240, "display.max_colwidth", 40):
            print(hits[["scheme", "population", "gene", "cpm_F536", "cpm_G_073", "cpm_M493",
                        "cpm_M399", "lfc_B1", "lfc_B2", "panel_rank_pct",
                        "panel_pass_rate"]].to_string(index=False))

    print("\n=== Which survive section-level separation? ===")
    if len(hits):
        sep = section_separation(adata, hits, schemes)
        sep.to_csv(OUT / "candidate_section_separation.csv", index=False)
        with pd.option_context("display.width", 240, "display.max_colwidth", 40):
            print(sep.sort_values(["separate_within_both_blocks", "sections_separate"],
                                  ascending=False).to_string(index=False))
        survivors = sep[sep["separate_within_both_blocks"]]
        print(f"\n{len(survivors)} of {len(sep)} candidates separate cleanly within BOTH blocks.")
    else:
        print("  (nothing to test)")
    print(f"\nWrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
