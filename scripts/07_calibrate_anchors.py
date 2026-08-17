#!/usr/bin/env python3
"""Calibrate simple marker thresholds for the anatomical anchor populations.

The anatomical frame (08_anatomy.py) has to be built on *whole sections*, before
any ROI exists, so it cannot depend on the ROI clustering.  Landmark detection
therefore uses raw marker thresholds instead of cell types.

This script checks those thresholds against the already-annotated ROI cells,
where the true cell type is known, and reports precision and recall so the
thresholds are chosen on evidence rather than by eye.
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
OUT = REPO / "results" / "anatomy"

# candidate rule -> the annotated cell types it is meant to capture
CANDIDATES = {
    "ventricle_lining": (
        [("Gpr50", 1), ("Gpr50", 2), ("Gpr50", 3), ("Spag16", 2), ("Cd24a", 3)],
        {"Tanycyte", "Ependymal"},
    ),
    "ARC": (
        [("Agrp", 1), ("Agrp", 2), ("Pomc", 2), ("Pomc", 3), ("Pomc", 5), ("Ghrh", 2)],
        {"ARC Agrp/Npy", "ARC Pomc", "ARC Th/Slc6a3 (TIDA)", "ARC Tac2/Esr1 (KNDy-like)"},
    ),
    "VMH": (
        [("Adcyap1", 3), ("Adcyap1", 5), ("Adcyap1", 8), ("Rasgrf2", 2), ("Rasgrf2", 3)],
        {"VMH-like Glut Rasgrf2", "VMH-like Glut Calb1", "VMH-like Glut Tac1"},
    ),
    "DMH": (
        [("Grp", 1), ("Grp", 2), ("Grp", 3), ("Ppp1r17", 2), ("Ppp1r17", 3)],
        {"DMH Grp/Ppp1r17"},
    ),
}


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(PROC / "mbh_roi_annotated.h5ad")
    truth = adata.obs["cell_type"].astype(str).to_numpy()

    rows = []
    for anchor, (rules, target_types) in CANDIDATES.items():
        is_target = np.isin(truth, list(target_types))
        for gene, cut in rules:
            if gene not in adata.var_names:
                continue
            called = counts_vector(adata, gene) >= cut
            tp = int((called & is_target).sum())
            rows.append({
                "anchor": anchor,
                "rule": f"{gene} >= {cut}",
                "n_called": int(called.sum()),
                "n_target": int(is_target.sum()),
                "precision": round(tp / max(called.sum(), 1) * 100, 1),
                "recall": round(tp / max(is_target.sum(), 1) * 100, 1),
            })
    res = pd.DataFrame(rows)
    res["f1"] = (2 * res.precision * res.recall / (res.precision + res.recall).replace(0, np.nan)).round(1)
    res.to_csv(OUT / "anchor_marker_calibration.csv", index=False)

    for anchor in CANDIDATES:
        sub = res[res.anchor == anchor].sort_values("f1", ascending=False)
        print(f"=== {anchor} (target n={sub.n_target.iloc[0]}) ===")
        print(sub[["rule", "n_called", "precision", "recall", "f1"]].to_string(index=False))
        print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
