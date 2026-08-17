#!/usr/bin/env python3
"""Variance components and detectable effect size for this design.

Three sections per animal make the variance decomposable, which is what turns
"n=2 is underpowered" into a number the paper can quote:

  sigma_animal   between animals within a group -- the variance the age test
                 must beat, and the only one adding biological replicates fixes
  sigma_section  between sections within an animal -- anatomical sampling, fixed
                 by AP matching, not by more animals
  sigma_cell     between cells within a section -- already negligible at these
                 cell counts

From those, the minimum detectable effect is computed for a blocked design at
n = 2, 3, 5, 8, 10 animals per group, both as-is and with AP matching removing
the section term.  This is what justifies the size of the next cohort.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats as st

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spatial_lea.io import ANIMAL_META, SECTION_ANIMAL, counts_matrix  # noqa: E402

PROC = REPO / "data" / "processed"
OUT = REPO / "results" / "power"
READOUTS = [
    ("Galr1", "GABA Gal/Galr1", "Galr1 in DMH Gal/Galr1 neurons"),
    ("Gal", "ARC Agrp/Npy", "Gal in ARC Agrp/Npy neurons"),
    ("Gal", "GABA Cacna2d2", "Gal in DMH GABA Cacna2d2 neurons"),
]
N_PER_GROUP = (2, 3, 5, 8, 10)
N_SIM = 4000
ALPHA = 0.05
SEED = 0


def section_values(adata, cell_type: str, gene: str) -> pd.DataFrame:
    sel = (adata.obs["cell_type"].astype(str) == cell_type).to_numpy()
    sub = adata[sel]
    gi = int(np.where(adata.var_names.to_numpy() == gene)[0][0])
    counts = counts_matrix(sub)
    df = pd.DataFrame({
        "section": sub.obs["section"].astype(str).to_numpy(),
        "g": counts[:, gi],
        "tot": counts.sum(axis=1),
    })
    per = df.groupby("section").agg(g=("g", "sum"), tot=("tot", "sum"), n=("g", "size"))
    per = per[per["n"] >= 15]
    per["log_cpm"] = np.log2(per["g"] / per["tot"] * 1e6 + 1)
    per["animal"] = [SECTION_ANIMAL[s] for s in per.index]
    per["group"] = per["animal"].map(lambda a: ANIMAL_META[a]["age_group"])
    return per


def variance_components(per: pd.DataFrame) -> dict:
    """Nested one-way decomposition on log2 CPM: animal within group, section within animal."""
    # Section-level residual around each animal's mean.
    resid = per.groupby("animal")["log_cpm"].transform(lambda x: x - x.mean())
    df_section = len(per) - per["animal"].nunique()
    var_section = float((resid**2).sum() / max(df_section, 1))

    animal_means = per.groupby(["group", "animal"])["log_cpm"].mean().reset_index()
    dev = animal_means.groupby("group")["log_cpm"].transform(lambda x: x - x.mean())
    df_animal = len(animal_means) - animal_means["group"].nunique()
    n_sec = per.groupby("animal").size().mean()
    # Animal means already carry section noise; subtract it to isolate the
    # biological component rather than double-counting.
    var_animal_obs = float((dev**2).sum() / max(df_animal, 1))
    var_animal = max(var_animal_obs - var_section / n_sec, 0.0)

    return {
        "sd_section": np.sqrt(var_section),
        "sd_animal": np.sqrt(var_animal),
        "sd_animal_observed": np.sqrt(var_animal_obs),
        "sections_per_animal": n_sec,
        "df_animal": df_animal,
    }


def min_detectable(sd_animal: float, sd_section: float, n_sections: float,
                   n_per_group: int, rng) -> float:
    """Smallest |log2 FC| a paired/blocked t-test would call at alpha, 80% power."""
    se_unit = np.sqrt(sd_animal**2 + sd_section**2 / max(n_sections, 1))
    lo, hi = 0.0, 8.0
    for _ in range(40):
        mid = (lo + hi) / 2
        aged = rng.normal(mid, se_unit, size=(N_SIM, n_per_group))
        adult = rng.normal(0.0, se_unit, size=(N_SIM, n_per_group))
        diff = aged - adult                      # blocked: paired within block
        if n_per_group < 2:
            return np.nan
        sd = diff.std(axis=1, ddof=1)
        # A simulated group with zero spread would give an infinite t; guard it
        # so the search converges on the variance, not on a numerical artefact.
        t = np.divide(diff.mean(axis=1), sd / np.sqrt(n_per_group),
                      out=np.zeros(N_SIM), where=sd > 0)
        crit = st.t.ppf(1 - ALPHA / 2, n_per_group - 1)
        power = float((np.abs(t) > crit).mean())
        lo, hi = (lo, mid) if power >= 0.8 else (mid, hi)
    return (lo + hi) / 2


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    adata = sc.read_h5ad(PROC / "hypothalamus_domains.h5ad")

    rows, comps = [], []
    for gene, cell_type, label in READOUTS:
        per = section_values(adata, cell_type, gene)
        if per["animal"].nunique() < 4:
            print(f"{label}: fewer than 4 animals available, skipped")
            continue
        vc = variance_components(per)
        obs_lfc = (per[per.group == "aged"]["log_cpm"].mean()
                   - per[per.group == "adult"]["log_cpm"].mean())
        comps.append({"readout": label, **{k: round(v, 3) for k, v in vc.items()},
                      "observed_lfc": round(float(obs_lfc), 2)})

        print(f"\n=== {label} ===")
        print(f"  sd between sections within an animal : {vc['sd_section']:.3f} log2")
        print(f"  sd between animals within a group    : {vc['sd_animal']:.3f} log2 "
              f"(observed {vc['sd_animal_observed']:.3f}, {vc['df_animal']} df)")
        print(f"  observed aged-adult difference       : {obs_lfc:+.2f} log2")
        for n in N_PER_GROUP:
            mde_now = min_detectable(vc["sd_animal"], vc["sd_section"],
                                     vc["sections_per_animal"], n, rng)
            mde_matched = min_detectable(vc["sd_animal"], 0.0, 1, n, rng)
            rows.append({"readout": label, "n_per_group": n,
                         "mde_log2": round(mde_now, 2),
                         "mde_fold": round(2**mde_now, 2),
                         "mde_log2_ap_matched": round(mde_matched, 2),
                         "mde_fold_ap_matched": round(2**mde_matched, 2),
                         "detectable_now": bool(abs(obs_lfc) >= mde_now)})
            flag = "  <- observed effect detectable" if abs(obs_lfc) >= mde_now else ""
            print(f"    n={n:2d}/group  min detectable {mde_now:5.2f} log2 "
                  f"({2**mde_now:5.2f}x) | AP-matched {mde_matched:5.2f} "
                  f"({2**mde_matched:5.2f}x){flag}")

    pd.DataFrame(comps).to_csv(OUT / "variance_components.csv", index=False)
    power = pd.DataFrame(rows)
    power.to_csv(OUT / "min_detectable_effect.csv", index=False)

    # A between-animal SD estimated from 2 animals per group has 2 df and can
    # land on zero by chance, as it does for Gal in ARC Agrp/Npy.  A zero there
    # is not evidence of zero variance, so the cohort recommendation uses the
    # pooled estimate across readouts instead of the per-readout one.
    pooled_sd = float(np.sqrt(np.mean([c["sd_animal"] ** 2 for c in comps])))
    pooled_sec = float(np.sqrt(np.mean([c["sd_section"] ** 2 for c in comps])))
    per_readout = ", ".join(f"{c['sd_animal']:.3f}" for c in comps)
    print("\n=== Pooled variance estimate (used for the recommendation) ===")
    print(f"  sd between animals  : {pooled_sd:.3f} log2   (per-readout: {per_readout})")
    print(f"  sd between sections : {pooled_sec:.3f} log2")
    boundary = [c["readout"] for c in comps if c["sd_animal"] == 0]
    if boundary:
        print(f"  boundary (zero) estimates, 2 df: {'; '.join(boundary)} "
              f"-- not evidence of zero variance")
    pooled_rows = []
    for n in N_PER_GROUP:
        pooled_rows.append({
            "n_per_group": n,
            "mde_log2": round(min_detectable(pooled_sd, pooled_sec, 3, n, rng), 2),
            "mde_log2_ap_matched": round(min_detectable(pooled_sd, 0.0, 1, n, rng), 2),
        })
    pooled_df = pd.DataFrame(pooled_rows)
    pooled_df["mde_fold"] = (2 ** pooled_df.mde_log2).round(2)
    pooled_df["mde_fold_ap_matched"] = (2 ** pooled_df.mde_log2_ap_matched).round(2)
    pooled_df.to_csv(OUT / "min_detectable_effect_pooled.csv", index=False)
    print(pooled_df.to_string(index=False))

    print("\n=== Cohort size needed for the surviving candidate ===")
    gal = power[power.readout.str.startswith("Gal in ARC")]
    if len(gal):
        obs = next(c["observed_lfc"] for c in comps if c["readout"].startswith("Gal in ARC"))
        ok = pooled_df[pooled_df.mde_log2 <= abs(obs)]
        ok_m = pooled_df[pooled_df.mde_log2_ap_matched <= abs(obs)]
        print(f"  observed effect {obs:+.2f} log2 ({2**abs(obs):.2f}x)")
        print(f"  using the pooled between-animal SD ({pooled_sd:.3f} log2):")
        print(f"    as sampled       : n >= {ok.n_per_group.min() if len(ok) else '>10'} per group")
        print(f"    with AP matching : n >= {ok_m.n_per_group.min() if len(ok_m) else '>10'} per group")
    print(f"\nWrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
