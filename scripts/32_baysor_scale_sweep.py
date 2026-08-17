#!/usr/bin/env python3
"""Choose Baysor's scale parameter by measurement rather than by assertion.

``scale`` is the expected cell radius, and it is the one Baysor parameter that
encodes the same kind of assumption we are trying to get rid of.  Picking it by
eye would replace a 5 um dilation with a 6 um guess, so it is swept on a 600 um
patch -- cheap enough to run many settings -- and scored against the DAPI nuclei,
which are measured rather than inferred.

**Primary criteria: the two unambiguous errors.**  Whatever a cell's true shape,
assigning two distinct measured nuclei to one cell is wrong, and scattering one
nucleus's own molecules across several cells is wrong.  Neither needs a
cell-type model, a size assumption or a ground truth beyond the DAPI image.

  merge rate  % of nuclei whose dominant cell also dominates another nucleus
  split rate  % of nuclei with under 80% of their molecules in one cell

**Secondary: is the cell count self-consistent with the cell size?**  The naive
expectation of one cell per nucleus is wrong, because a cell can cross the
section plane while its nucleus does not.  Modelling cells as spheres of
diameter ``d`` with concentric nuclei of diameter ``n``, centres uniform in z, a
cell intersects a section of thickness ``t`` with probability proportional to
``d + t`` and its nucleus with probability proportional to ``n + t``, so

    cells visible per visible nucleus  =  (d + t) / (n + t)

Every term on the right is measured: ``n`` from nucleus_boundaries, ``t`` from
the run's own thickness_of_high_quality_decoded_transcripts, ``d`` from the cell
sizes that same setting produced.  So this is a self-consistency test -- a
setting that invents extra cells must also produce cells large enough to justify
them.  It is secondary because Baysor's area is the hull of the molecules that
happen to fall in the slice, which understates ``d``, biasing the ratio upward
by an unknown amount.  Its use is to rule out the extremes, not to fine-tune.

**Read alongside**: molecules retained (dropping them throws away signal) and
the mixed GABAergic/glutamatergic rate (the artefact being removed).  The mixed
rate must never be read alone -- any setting can look clean by shrinking cells
until they hold no markers at all.

    python scripts/32_baysor_scale_sweep.py --patch-csv <molecules.csv>
    python scripts/32_baysor_scale_sweep.py --patch-csv <...> --rescore-only
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.io import RAW  # noqa: E402

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


def tissue_geometry(section: str) -> tuple[float, float]:
    """Measured nucleus diameter and effective section thickness, in um."""
    ex = json.loads((RAW / section / "experiment.xenium").read_text())
    t = float(ex["thickness_of_high_quality_decoded_transcripts"])
    area = pd.read_parquet(RAW / section / "cells.parquet", columns=["nucleus_area"])
    n = float(np.median(2 * np.sqrt(area["nucleus_area"] / np.pi)))
    return n, t


def score(out: Path, n_prior: int, nucleus_um: float, thickness_um: float) -> dict:
    counts = pd.read_csv(out / "segmentation_counts.tsv", sep="\t", index_col=0).T
    stats = pd.read_csv(out / "segmentation_cell_stats.csv", index_col=0)
    mol = pd.read_csv(out / "segmentation.csv", usecols=["prior", "cell"])
    assigned = mol[mol["cell"].notna() & (mol["cell"] != "")]

    keep = counts.sum(axis=1) >= 10
    counts = counts[keep]
    gaba = counts[[g for g in GABA if g in counts]].sum(axis=1)
    glut = counts[[g for g in GLUT if g in counts]].sum(axis=1)
    neuronal = (gaba >= MARKER_MIN) | (glut >= MARKER_MIN)
    mixed = (gaba >= MARKER_MIN) & (glut >= MARKER_MIN)

    # How each DAPI nucleus fared.  Split: under 80% of its own molecules end up
    # in the one cell that got most of them.  Merged: that cell is also the
    # majority owner of a second nucleus.
    p = assigned[assigned["prior"] > 0]
    per = p.groupby(["prior", "cell"]).size().rename("n").reset_index()
    dominant = per.sort_values("n", ascending=False).drop_duplicates("prior")
    total = p.groupby("prior").size()
    split = (dominant.set_index("prior")["n"] / total < 0.8)
    merged = dominant["cell"].duplicated(keep=False)

    diam = float(np.median(2 * np.sqrt(stats.loc[stats.index.isin(counts.index),
                                                 "area"] / np.pi)))
    expected = (diam + thickness_um) / (nucleus_um + thickness_um)
    observed = len(counts) / n_prior
    return {
        "cells": int(len(counts)),
        "pct_nuclei_merged": round(float(merged.mean()) * 100, 1),
        "pct_nuclei_split": round(float(split.mean()) * 100, 1),
        "median_diam_um": round(diam, 1),
        "cells_per_nucleus": round(observed, 2),
        "expected_per_nucleus": round(expected, 2),
        "excess": round(observed / expected, 2),
        "pct_molecules_assigned": round(len(assigned) / len(mol) * 100, 1),
        "median_counts": int(counts.sum(axis=1).median()),
        "pct_mixed": round(mixed.sum() / max(neuronal.sum(), 1) * 100, 1),
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--patch-csv", required=True, type=Path)
    p.add_argument("--section", default="G073_1", help="section the patch came from")
    p.add_argument("--rescore-only", action="store_true",
                   help="re-score existing runs without re-running Baysor")
    p.add_argument("--workdir", type=Path,
                   default=Path(os.environ.get("BAYSOR_SWEEP_DIR", "/tmp/baysor_sweep")))
    args = p.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    mols = pd.read_csv(args.patch_csv, usecols=["prior"])
    n_prior = int(mols.loc[mols["prior"] > 0, "prior"].nunique())
    nucleus_um, thickness_um = tissue_geometry(args.section)
    print(f"patch: {len(mols):,} molecules, {n_prior:,} DAPI nuclei")
    print(f"measured geometry: nucleus {nucleus_um:.2f} um, "
          f"section {thickness_um:.2f} um\n")

    rows = []
    for scale, std in GRID:
        out = args.workdir / f"s{scale:g}_{std.rstrip('%')}"
        if args.rescore_only:
            if not (out / "segmentation.csv").exists():
                continue
        else:
            print(f"  scale={scale:g} um, scale_std={std}", end="", flush=True)
            out = run_one(args.patch_csv, scale, std, args.workdir)
        rows.append({"scale_um": scale, "scale_std": std,
                     **score(out, n_prior, nucleus_um, thickness_um)})

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "scale_sweep.csv", index=False)
    print(df.to_string(index=False))

    # Rank on the two unambiguous errors only; the rest is context.
    err = df["pct_nuclei_merged"] + df["pct_nuclei_split"]
    best = df.loc[err.idxmin()]
    print(f"\nFewest unambiguous errors: scale={best.scale_um:g} um, "
          f"scale_std={best.scale_std} "
          f"({best.pct_nuclei_merged}% merged + {best.pct_nuclei_split}% split)")
    print(f"  but it produces {best.excess}x the cells its own cell size "
          f"({best.median_diam_um} um) can account for.")
    ok = df[df["excess"].between(0.8, 1.4)]
    if len(ok):
        pick = ok.loc[(ok["pct_nuclei_merged"] + ok["pct_nuclei_split"]).idxmin()]
        print(f"\nFewest errors among the self-consistent settings "
              f"(0.8 <= excess <= 1.4): scale={pick.scale_um:g} um, "
              f"scale_std={pick.scale_std}")
        print(f"  {pick.pct_nuclei_merged}% merged, {pick.pct_nuclei_split}% split, "
              f"{pick.median_diam_um} um cells, {pick.excess}x expected count, "
              f"{pick.pct_mixed}% mixed")
    print(f"\nWrote {OUT / 'scale_sweep.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
