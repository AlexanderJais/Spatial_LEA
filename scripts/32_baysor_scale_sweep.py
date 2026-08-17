#!/usr/bin/env python3
"""Choose Baysor's scale parameter by measurement rather than by assertion.

``scale`` is the expected cell radius, and it is the one Baysor parameter that
encodes the same kind of assumption we are trying to get rid of.  Picking it by
eye would replace a 5 um dilation with a 6 um guess, so it is swept on a 600 um
patch -- cheap enough to run many settings -- and scored on four things that can
be checked without believing any cell-type model:

  1. cells per DAPI nucleus.  Every cell in a section plane has at most one
     nucleus in that plane, and the nuclei are measured, not inferred.  A setting
     that yields far more or far fewer cells than nuclei is wrong regardless of
     anything else.  **This is the primary criterion.**
  2. molecules assigned.  Molecules dropped as noise are signal thrown away.
  3. mixed GABAergic/glutamatergic rate among cells carrying either marker --
     the artefact being removed.  Read only together with (2), since a setting
     can always look clean by shrinking cells until they hold no markers.
  4. one-to-one agreement with the DAPI nuclei: neither split across two Baysor
     cells nor merged with another nucleus.

    python scripts/32_baysor_scale_sweep.py --patch-csv <molecules.csv>
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

OUT = REPO / "results" / "baysor"
GABA = ["Gad1", "Gad2", "Slc32a1"]
GLUT = ["Slc17a6", "Slc17a7"]
MARKER_MIN = 2
# Radius in um.  4 is about a glial soma, 10 about a large hypothalamic neuron;
# 25%/50% scale_std sets how much of that range one run is allowed to span.
GRID = [(4.0, "25%"), (6.0, "25%"), (6.0, "50%"), (8.0, "25%"), (8.0, "50%"),
        (10.0, "50%")]


def run_one(molecules: Path, scale: float, scale_std: str, workdir: Path) -> Path:
    out = workdir / f"s{scale:g}_{scale_std.rstrip('%')}"
    if (out / "segmentation.csv").exists():
        return out
    out.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    cmd = [env["JULIA_BIN"], f"--project={env['BAYSOR_HOME']}",
           "-e", "using Baysor; Baysor.command_main()", "--",
           "run", "-x", "x", "-y", "y", "-z", "z", "-g", "gene",
           "-m", "30", "-s", str(scale), "--scale-std", scale_std,
           "--prior-segmentation-confidence", "0.5", "--n-clusters", "8",
           "--count-matrix-format", "tsv", "--polygon-format", "none",
           "-o", str(out), str(molecules), ":prior"]
    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if not (out / "segmentation.csv").exists():
        print(proc.stderr[-1500:])
        raise RuntimeError(f"Baysor failed at scale={scale} std={scale_std}")
    print(f"    ({time.time() - t0:.0f}s)", flush=True)
    return out


def score(out: Path, n_prior: int) -> dict:
    counts = pd.read_csv(out / "segmentation_counts.tsv", sep="\t", index_col=0).T
    mol = pd.read_csv(out / "segmentation.csv",
                      usecols=["prior", "cell", "is_noise"])
    assigned = mol[mol["cell"].notna() & (mol["cell"] != "")]

    keep = counts.sum(axis=1) >= 10
    counts = counts[keep]
    gaba = counts[[g for g in GABA if g in counts]].sum(axis=1)
    glut = counts[[g for g in GLUT if g in counts]].sum(axis=1)
    neuronal = (gaba >= MARKER_MIN) | (glut >= MARKER_MIN)
    mixed = (gaba >= MARKER_MIN) & (glut >= MARKER_MIN)

    # How each DAPI nucleus fared: split if its molecules land in >1 Baysor cell
    # that owns a majority of them; merged if a Baysor cell owns two nuclei.
    p = assigned[assigned["prior"] > 0]
    per = p.groupby(["prior", "cell"]).size().rename("n").reset_index()
    dominant = per.sort_values("n", ascending=False).drop_duplicates("prior")
    total = p.groupby("prior").size()
    intact = (dominant.set_index("prior")["n"] / total >= 0.8)
    merged = dominant["cell"].duplicated(keep=False)

    return {
        "cells": int(len(counts)),
        "cells_per_nucleus": round(len(counts) / n_prior, 3),
        "pct_molecules_assigned": round(len(assigned) / len(mol) * 100, 1),
        "median_counts": int(counts.sum(axis=1).median()),
        "pct_mixed": round(mixed.sum() / max(neuronal.sum(), 1) * 100, 1),
        "pct_nuclei_intact": round(float(intact.mean()) * 100, 1),
        "pct_nuclei_merged": round(float(merged.mean()) * 100, 1),
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--patch-csv", required=True, type=Path)
    p.add_argument("--workdir", type=Path,
                   default=Path(os.environ.get("BAYSOR_SWEEP_DIR", "/tmp/baysor_sweep")))
    args = p.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    mols = pd.read_csv(args.patch_csv, usecols=["prior"])
    n_prior = int(mols.loc[mols["prior"] > 0, "prior"].nunique())
    print(f"patch: {len(mols):,} molecules, {n_prior:,} DAPI nuclei\n")

    rows = []
    for scale, std in GRID:
        print(f"  scale={scale:g} um, scale_std={std}", end="", flush=True)
        out = run_one(args.patch_csv, scale, std, args.workdir)
        rows.append({"scale_um": scale, "scale_std": std, **score(out, n_prior)})
        print("   " + "  ".join(f"{k}={v}" for k, v in list(rows[-1].items())[2:]),
              flush=True)
        shutil.rmtree(out / "segmentation.csv", ignore_errors=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "scale_sweep.csv", index=False)
    print(f"\n{df.to_string(index=False)}")

    best = df.iloc[(df["cells_per_nucleus"] - 1).abs().argmin()]
    print(f"\nClosest to one cell per DAPI nucleus: scale={best.scale_um:g} um, "
          f"scale_std={best.scale_std} "
          f"({best.cells_per_nucleus} cells/nucleus, {best.pct_mixed}% mixed, "
          f"{best.pct_molecules_assigned}% of molecules kept)")
    print(f"\nWrote {OUT / 'scale_sweep.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
