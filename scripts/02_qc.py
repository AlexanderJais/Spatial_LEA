#!/usr/bin/env python3
"""QC, negative-control background calibration, and the ageing sanity check.

Three questions, in the order they should be asked:

1. Is each section technically sound, and what survives filtering?
2. What count level is distinguishable from background?  The panel's negative
   control probes and codewords give a per-section false-detection rate, which
   is what decides whether a 1-count ``Galr1`` cell is real.
3. Does the aged tissue actually look aged?  If the canonical glial ageing
   signature is absent, that reframes everything downstream, so it is checked
   before any effort goes into cell typing.

Writes results/qc/ (tables + figures) and prints a summary.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.io import ANIMAL_META, SECTION_ANIMAL, load_section  # noqa: E402

OUT = REPO / "results" / "qc"
MIN_COUNTS = 10
MIN_GENES = 5

# Canonical glial / inflammatory ageing markers that exist in this panel.
AGEING_UP = ["Gfap", "Cd68", "Trem2", "Spp1", "Cd44", "Ly6a", "Igfbp5", "C4b" ]
GAL_SYSTEM = ["Galr1", "Galr3", "Gal"]


def qc_table(sections) -> tuple[pd.DataFrame, dict]:
    rows, objects = [], {}
    for section in sections:
        adata = load_section(section)
        roi = adata[adata.obs["in_roi"]].copy()
        objects[section] = roi

        n_ctrl_probe = int((adata.uns["control_types"] == "Negative Control Probe").sum())
        n_ctrl_cw = int((adata.uns["control_types"] == "Negative Control Codeword").sum())
        obs = roi.obs

        # Background per probe per cell: what a gene measuring nothing would read.
        neg_probe = obs["negative_control_probe_counts"].to_numpy()
        neg_cw = obs["negative_control_codeword_counts"].to_numpy()
        bg_per_probe = neg_probe.sum() / (n_ctrl_probe * roi.n_obs)
        bg_per_cw = neg_cw.sum() / (n_ctrl_cw * roi.n_obs)

        keep = (obs["gene_counts"] >= MIN_COUNTS) & (obs["n_genes"] >= MIN_GENES) & (obs["nucleus_count"] == 1)
        rows.append(
            {
                "section": section,
                "animal": obs["animal"].iloc[0],
                "group": obs["age_group"].iloc[0],
                "batch": obs["batch"].iloc[0],
                "roi_cells": roi.n_obs,
                "median_counts": float(np.median(obs["gene_counts"])),
                "median_genes": float(np.median(obs["n_genes"])),
                "median_cell_area": round(float(np.median(obs["cell_area"])), 1),
                "median_nuc_area": round(float(np.median(obs["nucleus_area"])), 1),
                "neg_probe_rate": round(float(neg_probe.sum() / obs["gene_counts"].sum()), 5),
                "bg_counts_per_probe_per_cell": round(float(bg_per_probe), 4),
                "bg_counts_per_codeword_per_cell": round(float(bg_per_cw), 4),
                "multinucleate_pct": round(float((obs["nucleus_count"] > 1).mean() * 100), 2),
                "kept": int(keep.sum()),
                "kept_pct": round(float(keep.mean() * 100), 1),
            }
        )
        roi.obs["qc_pass"] = keep.to_numpy()
    return pd.DataFrame(rows), objects


def pseudobulk_cpm(objects, qc_only=True) -> pd.DataFrame:
    """Per-animal CPM over ROI cells -- the unit of replication is the animal."""
    frames = {}
    for section, roi in objects.items():
        sub = roi[roi.obs["qc_pass"]] if qc_only else roi
        counts = np.asarray(sub.X.sum(axis=0)).ravel()
        frames[SECTION_ANIMAL[section]] = pd.Series(counts * 1e6 / counts.sum(), index=sub.var_names)
    return pd.DataFrame(frames)


def detection_rate(objects, genes, qc_only=True) -> pd.DataFrame:
    out = {}
    for section, roi in objects.items():
        sub = roi[roi.obs["qc_pass"]] if qc_only else roi
        vals = {}
        for gene in genes:
            if gene in sub.var_names:
                col = np.asarray(sub[:, gene].X.todense()).ravel()
                vals[gene] = float((col > 0).mean() * 100)
        out[SECTION_ANIMAL[section]] = pd.Series(vals)
    return pd.DataFrame(out)


def blocked_lfc(cpm: pd.DataFrame) -> pd.DataFrame:
    """Within-batch log2 fold changes, aged - adult, one per block."""
    aged = {m: a for a, m in ANIMAL_META.items() for m in [m["batch"]] if a}  # placeholder
    blocks = {}
    for animal, meta in ANIMAL_META.items():
        blocks.setdefault(meta["batch"], {})[meta["age_group"]] = animal
    result = {}
    for batch, pair in sorted(blocks.items()):
        result[f"lfc_{batch}"] = np.log2((cpm[pair["aged"]] + 1) / (cpm[pair["adult"]] + 1))
    out = pd.DataFrame(result)
    out["mean_lfc"] = out.mean(axis=1)
    out["consistent"] = np.sign(out.iloc[:, 0]) == np.sign(out.iloc[:, 1])
    return out


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    sections = list(SECTION_ANIMAL)

    qc, objects = qc_table(sections)
    qc.to_csv(OUT / "section_qc.csv", index=False)
    print("\n=== Section QC (MBH ROI cells) ===")
    print(qc.to_string(index=False))

    cpm = pseudobulk_cpm(objects)
    cpm.to_csv(OUT / "pseudobulk_cpm.csv")
    lfc = blocked_lfc(cpm)
    lfc.to_csv(OUT / "blocked_lfc.csv")

    print("\n=== Ageing signature (pseudobulk CPM, aged - adult per block) ===")
    present = [g for g in AGEING_UP if g in cpm.index]
    view = pd.concat([cpm.loc[present].round(1), lfc.loc[present].round(2)], axis=1)
    print(view.to_string())
    n_consistent = int(lfc.loc[present, "consistent"].sum())
    n_up = int(((lfc.loc[present, "mean_lfc"] > 0) & lfc.loc[present, "consistent"]).sum())
    print(f"\n{n_consistent}/{len(present)} ageing markers agree in direction across both blocks; "
          f"{n_up} consistently UP in aged.")

    print("\n=== Galanin system (pseudobulk CPM) ===")
    gal_present = [g for g in GAL_SYSTEM if g in cpm.index]
    print(pd.concat([cpm.loc[gal_present].round(1), lfc.loc[gal_present].round(2)], axis=1).to_string())

    print("\n=== Detection rate, % of QC-passing ROI cells with >=1 count ===")
    det = detection_rate(objects, gal_present + present)
    det.to_csv(OUT / "detection_rate.csv")
    print(det.round(2).to_string())

    # Background reference: the strongest single negative-control probe.
    print("\n=== Background reference ===")
    for section, roi in objects.items():
        sub = roi[roi.obs["qc_pass"]]
        ctrl = np.asarray(sub.obsm["control_counts_matrix"].todense())
        is_probe = roi.uns["control_types"] == "Negative Control Probe"
        probe_det = (ctrl[:, is_probe] > 0).mean(axis=0) * 100
        galr1 = np.asarray(sub[:, "Galr1"].X.todense()).ravel()
        print(
            f"  {section:8s} Galr1 detected in {float((galr1>0).mean()*100):5.2f}% of cells | "
            f"neg-control probes: median {np.median(probe_det):.3f}%, max {probe_det.max():.3f}% | "
            f"ratio vs max probe = {float((galr1>0).mean()*100)/max(probe_det.max(),1e-9):.0f}x"
        )

    _plot_qc(qc, objects, cpm, present)
    print(f"\nWrote tables and figures to {OUT}")
    return 0


def _plot_qc(qc, objects, cpm, ageing_genes) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.5))
    order = qc.sort_values(["batch", "group"])["section"].tolist()
    colors = {"aged": "#B4531A", "adult": "#0E5A61"}

    ax = axes[0, 0]
    for section in order:
        roi = objects[section]
        grp = roi.obs["age_group"].iloc[0]
        ax.hist(np.log10(roi.obs["gene_counts"] + 1), bins=60, histtype="step",
                color=colors[grp], alpha=.85, label=f"{section} ({grp})")
    ax.axvline(np.log10(MIN_COUNTS + 1), color="k", ls="--", lw=1)
    ax.set_xlabel("log10(gene counts + 1)"); ax.set_ylabel("cells"); ax.set_title("Counts per ROI cell")
    ax.legend(fontsize=7)

    ax = axes[0, 1]
    ax.bar(range(len(qc)), qc["neg_probe_rate"] * 100,
           color=[colors[g] for g in qc["group"]])
    ax.set_xticks(range(len(qc))); ax.set_xticklabels(qc["section"], rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("% of counts"); ax.set_title("Negative control probe rate")

    ax = axes[0, 2]
    ax.bar(range(len(qc)), qc["kept_pct"], color=[colors[g] for g in qc["group"]])
    ax.set_xticks(range(len(qc))); ax.set_xticklabels(qc["section"], rotation=45, ha="right", fontsize=8)
    ax.set_ylim(0, 100); ax.set_ylabel("% cells"); ax.set_title(f"Passing QC (>={MIN_COUNTS} counts, 1 nucleus)")

    ax = axes[1, 0]
    for section in order:
        roi = objects[section]
        ax.scatter(roi.obs["roi_x"], roi.obs["roi_y"], s=.4, alpha=.25,
                   color=colors[roi.obs["age_group"].iloc[0]], rasterized=True)
    ax.set_aspect("equal"); ax.invert_yaxis()
    ax.set_title("MBH ROI cells, common frame"); ax.set_xlabel("µm"); ax.set_ylabel("µm")

    ax = axes[1, 1]
    idx = np.arange(len(ageing_genes))
    width = .2
    for i, animal in enumerate(cpm.columns):
        grp = ANIMAL_META[animal]["age_group"]
        ax.bar(idx + i * width, cpm.loc[ageing_genes, animal], width,
               color=colors[grp], alpha=.55 + .45 * (i % 2), label=f"{animal} ({grp})")
    ax.set_xticks(idx + 1.5 * width); ax.set_xticklabels(ageing_genes, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("CPM"); ax.set_title("Ageing markers, per animal"); ax.legend(fontsize=7)

    ax = axes[1, 2]
    gal = [g for g in GAL_SYSTEM if g in cpm.index]
    idx = np.arange(len(gal))
    for i, animal in enumerate(cpm.columns):
        grp = ANIMAL_META[animal]["age_group"]
        ax.bar(idx + i * width, cpm.loc[gal, animal], width,
               color=colors[grp], alpha=.55 + .45 * (i % 2), label=f"{animal} ({grp})")
    ax.set_xticks(idx + 1.5 * width); ax.set_xticklabels(gal, fontsize=9)
    ax.set_ylabel("CPM"); ax.set_title("Galanin system, per animal"); ax.legend(fontsize=7)

    fig.tight_layout()
    fig.savefig(OUT / "qc_overview.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
