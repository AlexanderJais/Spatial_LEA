#!/usr/bin/env python3
"""Name the Leiden clusters from their markers, and sanity-check the naming.

The mapping below is written out explicitly rather than derived by a scoring
heuristic, so the evidence for every label is visible and arguable.  Marker
evidence is from results/cluster/cluster_markers.csv (Wilcoxon, res 1.0).

Each entry records the markers the call rests on.  Where the panel cannot
distinguish a population (no Nr5a1 for VMH, no Kiss1 for KNDy), the label says
what the evidence supports and no more -- "VMH-like" rather than "VMH".
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.io import counts_vector  # noqa: E402

PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "annotation"

# cluster -> (label, class, the markers the call rests on)
ANNOTATION = {
    "0":  ("Endothelial",              "Vascular",     "Cldn5, Pecam1, Kdr, Emcn, Adgrl4"),
    "1":  ("VMH-like Glut Rasgrf2",    "Neuron_Glut",  "Rasgrf2, Adcyap1, Slc17a6, Cbln1, Calb1"),
    "2":  ("Mural (pericyte/SMC)",     "Vascular",     "Carmn, Cspg4, Ano1, Car4"),
    "3":  ("VLMC / fibroblast",        "Vascular",     "Dcn, Col1a1, Acta2, Ly6a, Igfbp4"),
    "4":  ("Arachnoid fibroblast",     "Vascular",     "Igfbp5, Dcn, Npnt, Mecom, Sntb1"),
    "5":  ("LHA Hcrt (orexin)",        "Neuron_Glut",  "Hcrt, Pdyn, Cckar, Shisa6"),
    "6":  ("GABA Nts/Gal",             "Neuron_GABA",  "Slc32a1, Gad1/2, Nts, Gal, Plch1"),
    "7":  ("GABA Cartpt/Glp1r",        "Neuron_GABA",  "Cartpt, Glp1r, Gal, Esr1, Prox1, Slc32a1"),
    "8":  ("GABA Gal/Galr1",           "Neuron_GABA",  "Gal, Galr1, Foxp2, Th, Gad1/2, Slc32a1"),
    "9":  ("GABA Cacna2d2",            "Neuron_GABA",  "Gad1/2, Slc32a1, Cacna2d2, Rims3, Cabp7"),
    "10": ("Low-quality / doublet",    "Excluded",     "Prph, Sncg; 12 cells, entropy 0.64"),
    "11": ("Oligodendrocyte",          "Oligo",        "Gjc3, Sox10, Opalin, Dpy19l1"),
    "12": ("Astrocyte",                "Astrocyte",    "Acsbg1, Ntsr2, Slc1a2, Aqp4, Slc39a12"),
    "13": ("Glut Otp/Nrn1",            "Neuron_Glut",  "Slc17a6, Nrn1, Ebf3, Otp, Mdga1"),
    "14": ("ME / meningeal fibroblast", "Vascular",    "Aldh1a2, Fmod, Slc13a4, Igf2, Col1a1, Spp1"),
    "15": ("GABA Otp/Wfs1",            "Neuron_GABA",  "Otp, Slc32a1, Cpne6, Wfs1, Gad1/2"),
    "16": ("Glut Prdm8/Cbln1",         "Neuron_Glut",  "Prdm8, Cbln1, Otp, Adcyap1, Bdnf"),
    "17": ("Microglia",                "Microglia",    "Laptm5, Siglech, Trem2, Cd53, Cd68, Ikzf1"),
    "18": ("VMH-like Glut Calb1",      "Neuron_Glut",  "Slc17a6, Calb1, Adcyap1, Cdh13, Tmem163"),
    "19": ("Glut Bcl11b/Prdm12",       "Neuron_Glut",  "Bcl11b, Prdm12, Calb2, Pthlh, Slc17a6"),
    "20": ("VMH-like Glut Tac1",       "Neuron_Glut",  "Adcyap1, Tac1, Npy2r, Esr1, Slc17a6, Bdnf"),
    "21": ("OPC",                      "OPC",          "Gpr17, Pdgfra, Cspg4, Sox10"),
    "22": ("Glut Synpr/Syt6",          "Neuron_Glut",  "Synpr, Syt6, Cpne4, Slc17a6, Calb2"),
    "23": ("DMH Grp/Ppp1r17",          "Neuron_Glut",  "Grp, Ppp1r17, Pdyn, Foxp2, Prox1, Slc17a6"),
    "24": ("GABA Igf1/Chodl",          "Neuron_GABA",  "Igf1, Chodl, Gad1, Thsd7a, Nwd2"),
    "25": ("GABA Nell1/Pcsk5",         "Neuron_GABA",  "Gad1/2, Slc32a1, Nell1, Grik3, Pcsk5"),
    "26": ("ARC Tac2/Esr1 (KNDy-like)", "Neuron_GABA", "Tac2, Esr1, Gal, Ghsr, Gpr50"),
    "27": ("ARC Pomc",                 "Neuron_Glut",  "Pomc, Cartpt, Esr1, Npy2r"),
    "28": ("ARC Th/Slc6a3 (TIDA)",     "Neuron_GABA",  "Th, Slc6a3, Satb2, Gad1"),
    "29": ("Oxt/Avp magnocellular",    "Neuron_Glut",  "Avp, Oxt, Otp, Fezf2"),
    "30": ("Ependymal",                "Ependymal",    "Cd24a, Spag16, Trp73, Rfx4, Plch1"),
    "31": ("ARC Agrp/Npy",             "Neuron_GABA",  "Npy, Agrp, Otp, Npy2r"),
    "32": ("Tanycyte",                 "Tanycyte",     "Gpr50, Slit2, Cd44, Penk, Cyp1b1"),
}

# Populations whose identity implies an anatomical nucleus.  Used in 05_nuclei.py
# as spatial anchors, and listed here so the link is explicit.
NUCLEUS_ANCHORS = {
    "ARC Agrp/Npy": "ARC", "ARC Pomc": "ARC", "ARC Th/Slc6a3 (TIDA)": "ARC",
    "ARC Tac2/Esr1 (KNDy-like)": "ARC",
    "DMH Grp/Ppp1r17": "DMH",
    "VMH-like Glut Rasgrf2": "VMH", "VMH-like Glut Calb1": "VMH", "VMH-like Glut Tac1": "VMH",
    "LHA Hcrt (orexin)": "LHA",
    "Tanycyte": "ME_3V", "Ependymal": "ME_3V",
}


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(PROC / "mbh_roi_clustered.h5ad")

    key = "leiden_1.0"
    codes = adata.obs[key].astype(str)
    adata.obs["cell_type"] = pd.Categorical(codes.map(lambda c: ANNOTATION[c][0]))
    adata.obs["cell_class"] = pd.Categorical(codes.map(lambda c: ANNOTATION[c][1]))

    n_excluded = int((adata.obs["cell_class"] == "Excluded").sum())
    adata = adata[adata.obs["cell_class"] != "Excluded"].copy()
    print(f"Dropped {n_excluded} cells from the low-quality cluster; {adata.n_obs:,} remain.\n")

    # Sanity check: the class calls must agree with canonical markers, not just
    # with the cluster numbering.  Reported so a wrong call is visible.
    checks = {
        "Neuron_GABA": ["Gad1", "Gad2", "Slc32a1"],
        "Neuron_Glut": ["Slc17a6"],
        "Astrocyte": ["Aqp4", "Slc1a2"],
        "Oligo": ["Sox10", "Opalin"],
        "OPC": ["Pdgfra", "Cspg4"],
        "Microglia": ["Trem2", "Cd53"],
        "Vascular": ["Cldn5", "Dcn"],
        "Tanycyte": ["Gpr50"],
        "Ependymal": ["Spag16"],
    }
    rows = []
    for klass, genes in checks.items():
        sub = adata[adata.obs["cell_class"] == klass]
        rest = adata[adata.obs["cell_class"] != klass]
        for gene in genes:
            if gene not in adata.var_names:
                continue
            rows.append({
                "class": klass, "marker": gene,
                "pct_in_class": round(float((counts_vector(sub, gene) > 0).mean() * 100), 1),
                "pct_elsewhere": round(float((counts_vector(rest, gene) > 0).mean() * 100), 1),
            })
    check = pd.DataFrame(rows)
    check["enriched"] = check.pct_in_class > check.pct_elsewhere
    check.to_csv(OUT / "marker_sanity_check.csv", index=False)
    print("=== Marker sanity check ===")
    print(check.to_string(index=False))
    print(f"\n{int(check.enriched.sum())}/{len(check)} canonical markers enriched in their assigned class.")

    comp = pd.crosstab(adata.obs["cell_type"], adata.obs["animal"])
    frac = pd.crosstab(adata.obs["cell_type"], adata.obs["animal"], normalize="columns") * 100
    table = pd.DataFrame({
        "class": adata.obs.groupby("cell_type", observed=True)["cell_class"].first(),
        "n_cells": comp.sum(axis=1),
        "markers": {v[0]: v[2] for v in ANNOTATION.values()},
    }).join(frac.round(2).add_suffix("_pct"))
    table = table.sort_values(["class", "n_cells"], ascending=[True, False])
    table.to_csv(OUT / "cell_types.csv")

    print("\n=== Annotated cell types (% of each animal's ROI cells) ===")
    with pd.option_context("display.width", 220, "display.max_colwidth", 46):
        print(table.drop(columns="markers").to_string())

    adata.write_h5ad(PROC / "mbh_roi_annotated.h5ad")
    print(f"\nWrote {PROC / 'mbh_roi_annotated.h5ad'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
