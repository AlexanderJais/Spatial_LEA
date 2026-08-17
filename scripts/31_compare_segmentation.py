#!/usr/bin/env python3
"""Baysor against the vendor segmentation, on identical tissue.

The two segmentations see exactly the same molecules -- the Baysor input is the
same hypothalamic window the vendor cells are filtered to -- so every difference
below is segmentation, not tissue, not thresholding and not sampling.

What is being tested.  The vendor cell is a DAPI nucleus dilated by a fixed 5 um.
That is a shape assumption, and measurement showed it fits neither class of cell
in this tissue: neurons come out at 15.8 um across against a real 15-25 um, glia
at 13.9 um against a real 8-12 um.  Under-captured neurons lose the transcripts
in their cytoplasm; over-expanded glia gain their neighbours'.  Both errors blur
cells toward the local average, and the visible symptom is cells carrying
GABAergic and glutamatergic markers at once -- which no neuron in this tissue
should.  Baysor drops the shape assumption and assigns each molecule using local
molecular composition as well as position.

So the report is organised around one question per section:

  A  geometry      -- how big are the cells now, and how many are there
  B  correspondence-- do Baysor cells and DAPI nuclei agree one-to-one
  C  mixed identity-- the marker-level diagnostic, free of any cell-type model
  D  composition   -- does the neuron/glia balance move
  E  Galr1         -- do the receptor claims in results/ survive re-segmentation
  F  age contrast  -- direction only; one section per age is not a test

Cell types are transferred, not re-derived: each Baysor cell takes the label of
the nearest vendor-defined cell-type centroid in log-CPM space.  Using the same
reference for both sides means a composition difference cannot come from having
clustered the two datasets differently.

    python scripts/31_compare_segmentation.py G073_1 M399_3
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial import cKDTree

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.io import ANIMAL_META, SECTION_ANIMAL, counts_matrix  # noqa: E402

SEG = REPO / "data" / "segmentation"
PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "baysor"
GABA = ["Gad1", "Gad2", "Slc32a1"]
GLUT = ["Slc17a6", "Slc17a7"]
MARKER_MIN = 2          # same rule as the published figure: >=2 counts is "positive"
MIN_COUNTS = 10         # vendor QC floor, applied to Baysor cells too


def load_baysor(section: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Per-cell counts, per-cell stats, per-molecule assignments."""
    d = SEG / f"{section}_baysor"
    counts = pd.read_csv(d / "segmentation_counts.tsv", sep="\t", index_col=0).T
    stats = pd.read_csv(d / "segmentation_cell_stats.csv", index_col=0)
    mol = pd.read_csv(d / "segmentation.csv",
                      usecols=["gene", "prior", "ml", "dv", "overlaps_nucleus",
                               "cell", "is_noise", "confidence"])
    return counts, stats, mol


def cell_frame(counts, stats, mol) -> pd.DataFrame:
    """One row per Baysor cell, with position, size and nuclear provenance."""
    assigned = mol[mol["cell"].notna() & (mol["cell"] != "")]
    g = assigned.groupby("cell")
    df = pd.DataFrame({
        "ml": g["ml"].mean(), "dv": g["dv"].mean(),
        "n_mol": g.size(),
        "frac_in_nucleus": g["overlaps_nucleus"].mean(),
    })
    df = df.join(stats[["area", "n_transcripts", "elongation",
                        "avg_assignment_confidence"]], how="left")
    # Which DAPI nucleus, if any, this Baysor cell mostly came from.
    prior = assigned[assigned["prior"] > 0]
    if len(prior):
        top = (prior.groupby(["cell", "prior"]).size().rename("n")
               .reset_index().sort_values("n", ascending=False)
               .drop_duplicates("cell").set_index("cell"))
        df["prior_cell"] = top["prior"]
        df["n_from_prior"] = top["n"]
    df["prior_cell"] = df.get("prior_cell", pd.Series(dtype=float)).fillna(0).astype(int)
    return df.join(counts, how="left")


def reference_profiles(ref) -> pd.DataFrame:
    """Mean log-CPM profile per vendor-annotated cell type."""
    counts = counts_matrix(ref)
    cpm = np.log1p(counts / np.maximum(counts.sum(axis=1, keepdims=True), 1) * 1e4)
    df = pd.DataFrame(cpm, columns=ref.var_names)
    df["cell_type"] = ref.obs["cell_type"].astype(str).to_numpy()
    df = df[~df["cell_type"].str.startswith("unlabelled")]
    return df.groupby("cell_type").mean()


def transfer(counts: pd.DataFrame, profiles: pd.DataFrame) -> pd.Series:
    """Nearest cell-type profile by correlation in log-CPM space."""
    genes = [g for g in profiles.columns if g in counts.columns]
    x = counts[genes].to_numpy(dtype=float)
    x = np.log1p(x / np.maximum(x.sum(axis=1, keepdims=True), 1) * 1e4)
    p = profiles[genes].to_numpy(dtype=float)
    xz = (x - x.mean(1, keepdims=True)) / (x.std(1, keepdims=True) + 1e-9)
    pz = (p - p.mean(1, keepdims=True)) / (p.std(1, keepdims=True) + 1e-9)
    corr = xz @ pz.T / len(genes)
    return pd.Series(profiles.index.to_numpy()[corr.argmax(1)], index=counts.index)


def transmitter(counts: pd.DataFrame) -> pd.Series:
    gaba = counts[[g for g in GABA if g in counts]].sum(axis=1)
    glut = counts[[g for g in GLUT if g in counts]].sum(axis=1)
    out = pd.Series("neither", index=counts.index)
    out[(gaba >= MARKER_MIN) & (glut < MARKER_MIN)] = "GABAergic"
    out[(glut >= MARKER_MIN) & (gaba < MARKER_MIN)] = "glutamatergic"
    out[(gaba >= MARKER_MIN) & (glut >= MARKER_MIN)] = "mixed"
    return out


def summarise(section: str, ref_all) -> dict:
    animal = SECTION_ANIMAL[section]
    print(f"\n{'='*74}\n{section}  ({animal}, {ANIMAL_META[animal]['age_group']})\n{'='*74}")

    counts, stats, mol = load_baysor(section)
    cells = cell_frame(counts, stats, mol)
    gene_cols = [c for c in counts.columns]
    bay = cells[cells[gene_cols].sum(axis=1) >= MIN_COUNTS].copy()

    ven = ref_all[ref_all.obs["section"].astype(str) == section]
    ven_counts = pd.DataFrame(counts_matrix(ven), columns=ven.var_names,
                              index=ven.obs_names)

    # --- A. geometry -------------------------------------------------------
    n_mol_total = len(mol)
    noise = int(mol["is_noise"].sum()) if "is_noise" in mol else 0
    row = {
        "section": section, "animal": animal,
        "age_group": ANIMAL_META[animal]["age_group"],
        "molecules": n_mol_total,
        "noise_pct": round(noise / n_mol_total * 100, 1),
        "baysor_cells": len(bay), "vendor_cells": ven.n_obs,
        "cell_ratio": round(len(bay) / ven.n_obs, 2),
        "baysor_median_counts": int(bay[gene_cols].sum(axis=1).median()),
        "vendor_median_counts": int(np.median(ven_counts.sum(axis=1))),
        "baysor_median_area": round(float(bay["area"].median()), 1),
        "vendor_median_area": round(float(np.median(ven.obs["cell_area"])), 1),
    }
    print("A. geometry")
    print(f"   molecules in window     {n_mol_total:>10,}  ({row['noise_pct']}% called noise)")
    print(f"   cells                   {len(bay):>10,}  vendor {ven.n_obs:,} "
          f"({row['cell_ratio']}x)")
    print(f"   median counts/cell      {row['baysor_median_counts']:>10}  "
          f"vendor {row['vendor_median_counts']}")
    print(f"   median area (um2)       {row['baysor_median_area']:>10}  "
          f"vendor {row['vendor_median_area']}  "
          f"(equiv. diameter {2*np.sqrt(row['baysor_median_area']/np.pi):.1f} vs "
          f"{2*np.sqrt(row['vendor_median_area']/np.pi):.1f} um)")

    # --- B. correspondence with the DAPI nuclei ----------------------------
    with_prior = bay[bay["prior_cell"] > 0]
    dupes = with_prior["prior_cell"].duplicated(keep=False).sum()
    row["pct_cells_with_nucleus"] = round(len(with_prior) / len(bay) * 100, 1)
    row["pct_nuclei_split"] = round(dupes / max(len(with_prior), 1) * 100, 1)
    row["median_frac_in_nucleus"] = round(float(bay["frac_in_nucleus"].median()), 3)
    print("B. correspondence with DAPI")
    print(f"   Baysor cells anchored on a nucleus  {row['pct_cells_with_nucleus']:>5}%")
    print(f"   of those, sharing a nucleus (split) {row['pct_nuclei_split']:>5}%")
    print(f"   median molecules inside the nucleus {row['median_frac_in_nucleus']*100:>5.1f}%"
          "   -- the rest is cytoplasm the 5 um dilation had to guess at")

    # --- C. mixed transmitter identity -------------------------------------
    bt, vt = transmitter(bay[gene_cols]), transmitter(ven_counts)
    def mixed_rate(t):
        neuronal = t[t != "neither"]
        return round((neuronal == "mixed").mean() * 100, 1) if len(neuronal) else np.nan
    row["baysor_mixed_pct"] = mixed_rate(bt)
    row["vendor_mixed_pct"] = mixed_rate(vt)
    print("C. mixed GABAergic/glutamatergic identity  (of cells with any marker)")
    print(f"   vendor {row['vendor_mixed_pct']:>5}%   ->   Baysor {row['baysor_mixed_pct']:>5}%")

    # --- D. composition ----------------------------------------------------
    profiles = reference_profiles(ref_all)
    bay["cell_type"] = transfer(bay[gene_cols], profiles)
    ven_type = transfer(ven_counts, profiles)   # transferred the same way, for fairness
    comp = pd.DataFrame({
        "baysor_pct": bay["cell_type"].value_counts(normalize=True) * 100,
        "vendor_pct": ven_type.value_counts(normalize=True) * 100,
    }).fillna(0).round(2)
    comp["delta"] = (comp["baysor_pct"] - comp["vendor_pct"]).round(2)
    comp.to_csv(OUT / f"{section}_composition.csv")

    glia = ["Astrocyte", "Oligodendrocyte", "OPC", "Microglia", "Tanycyte", "Ependymal"]
    row["baysor_glia_pct"] = round(comp.loc[comp.index.isin(glia), "baysor_pct"].sum(), 1)
    row["vendor_glia_pct"] = round(comp.loc[comp.index.isin(glia), "vendor_pct"].sum(), 1)
    print("D. composition (labels transferred from the same reference on both sides)")
    print(f"   glia + ependymal   vendor {row['vendor_glia_pct']:>5}%  ->  "
          f"Baysor {row['baysor_glia_pct']:>5}%")
    print("   largest shifts:")
    for t, r in comp.reindex(comp["delta"].abs().sort_values(ascending=False).index).head(6).iterrows():
        print(f"     {t:32s} {r.vendor_pct:6.2f}% -> {r.baysor_pct:6.2f}%  ({r.delta:+.2f})")

    # --- E. Galr1 ----------------------------------------------------------
    def galr1_block(counts_df, types, tag):
        pos = counts_df["Galr1"] > 0
        out = {f"{tag}_galr1_pct": round(pos.mean() * 100, 1)}
        t = transmitter(counts_df)
        sig = counts_df["Galr1"]
        gsum, lsum = sig[t == "GABAergic"].sum(), sig[t == "glutamatergic"].sum()
        out[f"{tag}_galr1_gaba_glut"] = round(gsum / max(lsum, 1), 2)
        top = (sig.groupby(types).sum() / max(sig.sum(), 1) * 100).nlargest(4)
        return out, top

    b_stat, b_top = galr1_block(bay[gene_cols], bay["cell_type"], "baysor")
    v_stat, v_top = galr1_block(ven_counts, ven_type, "vendor")
    row.update(b_stat); row.update(v_stat)
    print("E. Galr1")
    print(f"   % of cells Galr1+          vendor {v_stat['vendor_galr1_pct']:>5}%  ->  "
          f"Baysor {b_stat['baysor_galr1_pct']:>5}%")
    print(f"   GABA:Glut share of signal  vendor {v_stat['vendor_galr1_gaba_glut']:>5}x  ->  "
          f"Baysor {b_stat['baysor_galr1_gaba_glut']:>5}x")
    print("   top carrying types (% of all Galr1 signal):")
    for t in dict.fromkeys(list(v_top.index) + list(b_top.index)):
        print(f"     {t:32s} vendor {v_top.get(t, 0):5.1f}%   Baysor {b_top.get(t, 0):5.1f}%")

    # Spatial: does the DMH still carry the receptor?  Nucleus is inherited from
    # the nearest vendor cell, so the parcellation is identical on both sides.
    tree = cKDTree(ven.obs[["ml", "dv"]].to_numpy())
    _, idx = tree.query(bay[["ml", "dv"]].to_numpy())
    bay["nucleus_ext"] = ven.obs["nucleus_ext"].astype(str).to_numpy()[idx]
    nuc_b = (bay.groupby("nucleus_ext")[gene_cols].apply(lambda x: x["Galr1"].sum())
             / max(bay["Galr1"].sum(), 1) * 100)
    ven_nuc = ven.obs["nucleus_ext"].astype(str).to_numpy()
    nuc_v = (pd.Series(ven_counts["Galr1"].to_numpy()).groupby(ven_nuc).sum()
             / max(ven_counts["Galr1"].sum(), 1) * 100)
    print("   Galr1 signal by nucleus (% of section total):")
    for n in nuc_v.sort_values(ascending=False).index[:6]:
        print(f"     {n:32s} vendor {nuc_v.get(n, 0):5.1f}%   Baysor {nuc_b.get(n, 0):5.1f}%")
    for n in nuc_v.index:
        row[f"galr1_pct_{n}_vendor"] = round(float(nuc_v.get(n, 0)), 1)
        row[f"galr1_pct_{n}_baysor"] = round(float(nuc_b.get(n, 0)), 1)

    bay[["ml", "dv", "n_mol", "area", "frac_in_nucleus", "prior_cell",
         "cell_type", "nucleus_ext"]].to_csv(OUT / f"{section}_baysor_cells.csv")
    return row


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("sections", nargs="+")
    args = p.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    ref = sc.read_h5ad(PROC / "hypothalamus_nuclei.h5ad")
    rows = [summarise(s, ref) for s in args.sections]
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "segmentation_comparison.csv", index=False)

    # --- F. age contrast, direction only -----------------------------------
    if len(df) == 2 and df["age_group"].nunique() == 2:
        print(f"\n{'='*74}\nF. age contrast -- one section per age, so direction only\n{'='*74}")
        aged = df[df.age_group == "aged"].iloc[0]
        adult = df[df.age_group == "adult"].iloc[0]
        for col, name in [("baysor_galr1_pct", "Galr1+ cells, Baysor"),
                          ("vendor_galr1_pct", "Galr1+ cells, vendor"),
                          ("baysor_mixed_pct", "mixed-identity rate, Baysor"),
                          ("vendor_mixed_pct", "mixed-identity rate, vendor")]:
            print(f"   {name:32s} aged {aged[col]:6}   adult {adult[col]:6}   "
                  f"diff {aged[col] - adult[col]:+.1f}")
        print("\n   Two sections cannot separate age from animal, section plane or run.")
        print("   These numbers say whether re-segmentation moved the contrast, not")
        print("   whether the contrast is real.")

    print(f"\nWrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
