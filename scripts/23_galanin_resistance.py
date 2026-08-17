#!/usr/bin/env python3
"""Test the galanin-resistance hypothesis as a system-level signature.

Resistance predicts a specific shape, the same one that defines leptin
resistance: **ligand rises while receptor stays flat or falls**, so the
ligand-to-receptor balance shifts and the same amount of signalling requires
more ligand.

Tested gene by gene this is hopeless at n=2.  Tested as an ensemble it is not,
because the unit becomes the **population** rather than the animal: across ~50
nucleus x cell-type populations, is the distribution of Gal age-effects shifted
relative to the distribution of Galr1 age-effects, and relative to the panel?
That is a question about a distribution, and 50 populations give it traction
that a single contrast does not.

Four readouts:

  1. ensemble shift    -- Gal vs Galr1 LFC distributions across populations,
                          each calibrated against all 297 panel genes
  2. resistance index  -- per population, LFC(Gal) - LFC(Galr1); its sign
                          consistency, calibrated against random gene pairs
  3. paired test       -- within each population, is Gal above Galr1
  4. spatial coupling  -- local Gal available to each Galr1+ cell, by age

All log fold changes are AP-corrected (21_ap_correction.py method) before use.
The populations share four animals, so they are not independent; the panel
calibration is what makes the comparison interpretable, not a nominal p-value.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc
import statsmodels.formula.api as smf
from scipy import stats as st
from sklearn.neighbors import NearestNeighbors

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.io import ANIMAL_META, SECTION_ANIMAL, counts_matrix  # noqa: E402

PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "resistance"
AP_FILE = REPO / "results" / "ap_matching" / "section_ap_scores.csv"
BLOCKS = {"B1": ("F536", "M493"), "B2": ("G_073", "M399")}
NUCLEI = ("ARC", "ME_3V", "VMH", "DMH", "LHA", "ZI", "DHA_PH")
MIN_CELLS_PER_SECTION = 20
MIN_COUNTS = 150
SEED = 0


def ap_corrected_lfc(mat: pd.DataFrame, totals: pd.Series, ap: pd.Series) -> pd.Series:
    """AP-corrected blocked age effect for every gene in one population."""
    animals = pd.Series([SECTION_ANIMAL[s] for s in mat.index], index=mat.index)
    if animals.nunique() < 4:
        return pd.Series(dtype=float)
    out = {}
    for gene in mat.columns:
        if mat[gene].sum() < MIN_COUNTS:
            continue
        y = np.log2(mat[gene] / totals * 1e6 + 1)
        if y.std() < 1e-9:
            continue
        df = pd.DataFrame({"y": y, "ap": ap[mat.index], "animal": animals})
        try:
            beta = float(smf.ols("y ~ ap + C(animal)", df).fit().params["ap"])
        except Exception:
            continue
        df["yc"] = df["y"] - beta * (df["ap"] - df["ap"].mean())
        per_animal = df.groupby("animal")["yc"].mean()
        out[gene] = float(np.mean([per_animal[a] - per_animal[d] for a, d in BLOCKS.values()]))
    return pd.Series(out)


def build(adata, ap) -> pd.DataFrame:
    """AP-corrected LFC per gene per population."""
    obs = adata.obs
    pops = (obs["nucleus_ext"].astype(str) + " | " + obs["cell_type"].astype(str))
    frames = {}
    for pop in pops.unique():
        sel = (pops == pop).to_numpy()
        if sel.sum() < MIN_CELLS_PER_SECTION * 8:
            continue
        sub = adata[sel]
        counts = counts_matrix(sub)
        sections = sub.obs["section"].astype(str).to_numpy()
        keep = [s for s in pd.unique(sections) if (sections == s).sum() >= MIN_CELLS_PER_SECTION]
        if len({SECTION_ANIMAL[s] for s in keep}) < 4:
            continue
        mat = pd.DataFrame(np.vstack([counts[sections == s].sum(axis=0) for s in keep]),
                           index=keep, columns=adata.var_names)
        lfc = ap_corrected_lfc(mat, mat.sum(axis=1), ap)
        if len(lfc):
            frames[pop] = lfc
    return pd.DataFrame(frames).T


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    adata = sc.read_h5ad(PROC / "hypothalamus_nuclei.h5ad")
    adata = adata[adata.obs["nucleus_ext"].isin(NUCLEI)].copy()
    ap = pd.read_csv(AP_FILE, index_col=0)["ap_score"]

    lfc = build(adata, ap)
    lfc.to_csv(OUT / "ap_corrected_lfc_by_population.csv")
    print(f"AP-corrected age effects for {lfc.shape[1]} genes "
          f"across {lfc.shape[0]} populations\n")

    # ---- 1. ensemble shift -------------------------------------------------
    print("=== 1. Ensemble shift: is Gal displaced upward, and Galr1 not? ===")
    med = lfc.median(axis=0)          # each gene's typical age effect across populations
    frac_up = (lfc > 0).mean(axis=0)  # fraction of populations where it rises
    rows = []
    for gene in ("Gal", "Galr1", "Galr3"):
        if gene not in lfc.columns:
            continue
        rows.append({
            "gene": gene,
            "median_lfc": round(float(med[gene]), 3),
            "pct_populations_up": round(float(frac_up[gene]) * 100, 1),
            "n_populations": int(lfc[gene].notna().sum()),
            "panel_pct_above": round(float((med > med[gene]).mean() * 100), 1),
            "panel_pct_more_consistent": round(
                float((frac_up.sub(0.5).abs() > abs(frac_up[gene] - 0.5)).mean() * 100), 1),
        })
    ens = pd.DataFrame(rows)
    print(ens.to_string(index=False))
    print("  panel_pct_above: % of the 297 genes with a higher median age effect")
    print("  panel_pct_more_consistent: % of genes more one-sided across populations")

    # ---- 2. resistance index ----------------------------------------------
    print("\n=== 2. Resistance index per population: LFC(Gal) - LFC(Galr1) ===")
    both = lfc[["Gal", "Galr1"]].dropna()
    both["resistance_index"] = both["Gal"] - both["Galr1"]
    both = both.sort_values("resistance_index", ascending=False)
    both.round(2).to_csv(OUT / "resistance_index.csv")
    print(both.round(2).head(12).to_string())
    n_pos = int((both["resistance_index"] > 0).sum())
    print(f"\n  index positive in {n_pos}/{len(both)} populations "
          f"({n_pos/len(both)*100:.0f}%), median {both.resistance_index.median():+.2f}")

    # Calibration: random gene pairs from the same populations.
    genes = [g for g in lfc.columns if lfc[g].notna().sum() >= len(both)]
    null_pos, null_med = [], []
    for _ in range(4000):
        a, b = rng.choice(genes, 2, replace=False)
        d = (lfc[a] - lfc[b]).reindex(both.index).dropna()
        if len(d) < 5:
            continue
        null_pos.append((d > 0).mean())
        null_med.append(d.median())
    null_pos, null_med = np.array(null_pos), np.array(null_med)
    p_pos = float((np.abs(null_pos - .5) >= abs(n_pos / len(both) - .5)).mean())
    p_med = float((np.abs(null_med) >= abs(both.resistance_index.median())).mean())
    print(f"  vs random gene pairs: one-sidedness p = {p_pos:.3f}, "
          f"median magnitude p = {p_med:.3f}")

    # ---- 3. paired within population ---------------------------------------
    w = st.wilcoxon(both["Gal"], both["Galr1"])
    print(f"\n=== 3. Paired within population (Gal vs Galr1): "
          f"Wilcoxon W = {w.statistic:.0f}, p = {w.pvalue:.4f} ===")
    print("  Populations are not independent, so read this as a descriptive")
    print("  summary of a consistent within-population ordering.")

    # ---- 4. spatial coupling ----------------------------------------------
    print("\n=== 4. Spatial coupling: Gal available locally to each Galr1+ cell ===")
    spatial = local_ligand(adata)
    spatial.to_csv(OUT / "spatial_coupling.csv", index=False)
    print(spatial.round(3).to_string(index=False))

    _plot(lfc, both, med, spatial)
    ens.to_csv(OUT / "ensemble_shift.csv", index=False)
    print(f"\nWrote {OUT}")
    return 0


def local_ligand(adata, radius=100.0) -> pd.DataFrame:
    """Gal+ neighbours within `radius` of each Galr1+ cell, per section."""
    gi = int(np.where(adata.var_names.to_numpy() == "Galr1")[0][0])
    gl = int(np.where(adata.var_names.to_numpy() == "Gal")[0][0])
    counts = counts_matrix(adata[:, [gi, gl]])
    obs = adata.obs
    rows = []
    for section, idx in obs.groupby("section", observed=True).indices.items():
        xy = obs.iloc[idx][["ml", "dv"]].to_numpy()
        galr1 = counts[idx, 0] > 0
        gal = counts[idx, 1] > 0
        if galr1.sum() < 30 or gal.sum() < 30:
            continue
        nn = NearestNeighbors(radius=radius).fit(xy)
        neigh = nn.radius_neighbors(xy[galr1], return_distance=False)
        frac = np.array([gal[n].mean() if len(n) else np.nan for n in neigh])
        animal = SECTION_ANIMAL[str(section)]
        rows.append({"section": str(section), "animal": animal,
                     "group": ANIMAL_META[animal]["age_group"],
                     "n_galr1": int(galr1.sum()),
                     "mean_pct_gal_neighbours": float(np.nanmean(frac) * 100),
                     "pct_gal_overall": float(gal.mean() * 100)})
    df = pd.DataFrame(rows)
    df["enrichment"] = df["mean_pct_gal_neighbours"] / df["pct_gal_overall"]
    return df


def _plot(lfc, both, med, spatial) -> None:
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.4))

    ax = axes[0]
    ax.hist(med.dropna(), bins=45, color="#C3CCCC", label="panel genes")
    for gene, colour in (("Gal", "#B4531A"), ("Galr1", "#0E5A61"), ("Galr3", "#6A3D9A")):
        if gene in med:
            ax.axvline(med[gene], color=colour, lw=1.8)
            ax.annotate(gene, (med[gene], ax.get_ylim()[1] * .9), color=colour,
                        fontsize=8, style="italic", ha="center")
    ax.axvline(0, color="#555", lw=.8, ls=":")
    ax.set_xlabel("median AP-corrected log2 aged/adult across populations")
    ax.set_ylabel("panel genes")
    ax.set_title("Ligand is displaced up,\nreceptor is not")

    ax = axes[1]
    ax.scatter(both["Galr1"], both["Gal"], s=26, c="#0E5A61", alpha=.8)
    lim = [min(both.min().min(), -1.5), max(both.max().max(), 1.5)]
    ax.plot(lim, lim, color="#555", ls="--", lw=1)
    ax.axhline(0, color="#BBB", lw=.8); ax.axvline(0, color="#BBB", lw=.8)
    ax.set_xlabel("Galr1 log2 aged/adult"); ax.set_ylabel("Gal log2 aged/adult")
    ax.set_title("Per population\n(above the line = ligand outpaces receptor)")

    ax = axes[2]
    vals = both["resistance_index"].sort_values()
    ax.barh(range(len(vals)), vals,
            color=["#B4531A" if v > 0 else "#0E5A61" for v in vals])
    ax.axvline(0, color="#555", lw=.8)
    ax.set_yticks([]); ax.set_xlabel("LFC(Gal) − LFC(Galr1)")
    ax.set_title(f"Resistance index\npositive in {(vals>0).sum()}/{len(vals)} populations")

    ax = axes[3]
    for i, animal in enumerate(["F536", "G_073", "M493", "M399"]):
        sub = spatial[spatial.animal == animal]
        colour = "#B4531A" if ANIMAL_META[animal]["age_group"] == "aged" else "#0E5A61"
        ax.scatter([i] * len(sub), sub["enrichment"], s=40, color=colour)
        if len(sub):
            ax.plot([i - .25, i + .25], [sub["enrichment"].mean()] * 2, color=colour, lw=2)
    ax.set_xticks(range(4))
    ax.set_xticklabels([f"{a}\n{ANIMAL_META[a]['age_group']}" for a in
                        ["F536", "G_073", "M493", "M399"]], fontsize=7)
    ax.set_ylabel("local Gal enrichment around Galr1+ cells")
    ax.set_title("Spatial ligand availability\none point per section")
    ax.axhline(1, color="#999", ls=":", lw=1)

    fig.tight_layout()
    fig.savefig(OUT / "galanin_resistance.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
