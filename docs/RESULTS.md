# Results — Galr1 in the ageing MBH (cohort A, 2 aged vs 2 adult)

Pipeline run on the four ROI sections: `scripts/02_qc.py` → `03_cluster.py` →
`04_annotate.py` → `05_galr1.py` → `06_section_matching.py`.
Tables and figures under `results/`.

**Headline: Galr1 shows no systematic age difference across MBH cell types at
this n.** Cross-block agreement is at chance level. The most consistent signal
in the galanin system is in the ligand `Gal`, not the receptor — and one
striking-looking composition result is probably an artefact of where the
sections were cut. Details below.

---

## 1. The data is technically excellent

| | F536 (aged) | G_073 (aged) | M493 (adult) | M399 (adult) |
|---|---|---|---|---|
| MBH ROI cells | 11 614 | 8 592 | 10 291 | 10 366 |
| passing QC | 99.9% | 100.0% | 99.8% | 99.7% |
| median counts/cell | 201 | 269 | 243 | 245 |
| median genes/cell | 72 | 82 | 78 | 79 |
| negative-control probe rate | 0.012% | 0.012% | 0.013% | 0.012% |
| multinucleate segmentations | 0% | 0% | 0% | 0% |

40 783 cells survive QC (≥10 counts, ≥5 genes, single nucleus).

**Galr1 is measured well above background.** It is detected in 15.3–19.0% of
ROI cells, versus 0.07–0.09% (median) and 0.55–0.80% (max) for the 27
negative-control probes — a **23–30× margin**. Galr3 is detected in 5.2–6.7%,
`Gal` in 30–42%. The receptor question is answerable with this data.

## 2. A clean 32-type MBH atlas

Leiden at resolution 1.0 gives 33 clusters, animal-mixing entropy median 0.99
(1.0 = perfectly mixed), and all canonical markers land where they should —
**16/16 class markers enriched in their assigned class**. One 12-cell cluster
was dropped as low quality.

Every expected MBH population is recovered, including `ARC Agrp/Npy`,
`ARC Pomc`, `ARC Th/Slc6a3` (TIDA), `ARC Tac2/Esr1` (KNDy-like),
three `VMH-like` glutamatergic types, `DMH Grp/Ppp1r17`, `LHA Hcrt`,
`Oxt/Avp` magnocellular, tanycytes, ependyma, and full glial and vascular sets.

**One cluster has `Galr1` itself as a defining marker** (LFC 2.9): a GABAergic
`Gal`/`Galr1`/`Foxp2`/`Th` population, 1 705 cells, **53–67% Galr1⁺**. This is
the cell type where galanin autoreceptor-like signalling is concentrated, and it
is the natural focus of the biology.

## 3. Galr1: no systematic age effect

Per-animal detection rate and depth-normalised expression, within-block log2
fold changes (B1 = F536/M493, B2 = G_073/M399), gated at ≥30 cells and ≥15
expressing cells per animal so fold changes are not driven by near-zero
denominators.

| gene | metric | cell types tested | agree across blocks | binomial p |
|---|---|---|---|---|
| Galr1 | detection rate | 16 | 8 | 0.60 |
| Galr1 | counts per 10k | 16 | 10 | 0.23 |
| Galr3 | detection rate | 8 | 3 | 0.86 |
| Galr3 | counts per 10k | 8 | 3 | 0.86 |
| Gal | counts per 10k | 29 | 11 | 0.93 |

Agreement is what two coin flips would give. **Nothing in Galr1 or Galr3 passes
a shortlist of "consistent in both blocks and ≥1.5× in the weaker block."**

The closest thing to a Galr1 signal is a modest, consistent *decrease* in the
`GABA Gal/Galr1` population itself — the one type where the measurement is most
precise:

| | F536 aged | G_073 aged | M493 adult | M399 adult | LFC B1 | LFC B2 |
|---|---|---|---|---|---|---|
| detection % | 53.4 | 64.2 | 63.9 | 67.2 | −0.26 | −0.07 |
| counts / 10k | 54.5 | 54.7 | 67.3 | 71.9 | −0.31 | −0.39 |

Consistent in direction in both blocks and both metrics, ~22–24% lower in aged.
Worth carrying forward as a hypothesis — but the aged animals do not separate
from the adults beyond within-group spread, so it is a lead, not a finding.

## 4. The ligand moves more than the receptor

Three effects pass the shortlist, all in `Gal`:

| cell type | metric | LFC B1 | LFC B2 | direction |
|---|---|---|---|---|
| Tanycyte | detection | −0.87 | −0.89 | **down in aged (~1.8×)** |
| Tanycyte | counts/10k | −0.67 | −1.35 | **down in aged (~2×)** |
| ME / meningeal fibroblast | counts/10k | −0.68 | −0.76 | down in aged |
| ARC Agrp/Npy | counts/10k | +1.02 | +0.62 | **up in aged (~1.8×)** |

Tanycyte `Gal` down and ARC Agrp/Npy `Gal` up, each consistent across two
independent blocks, is a more interesting shape than a flat receptor change: it
suggests the ageing shift is in galanin *availability* and *source*, not in
receptor abundance. It is still n=2 vs 2 and needs replication.

## 5. The ageing signature is present but modest

`Gfap` (+0.24), `Cd68` (+0.25) and `Trem2` (+0.30 log2) rise consistently in
both blocks. `Spp1`, `Cd44` and `Igfbp5` disagree between blocks, and `Ly6a`
falls. A partial glial/microglial signature is what 15–16-month animals should
look like — consistent with the middle-aged framing, and an argument for adding
an 18–24-month group.

## 6. Important caveat — one striking result is probably geometry

`DMH Grp/Ppp1r17` looks like a dramatic age effect: 1.57% and 1.22% of ROI cells
in the aged animals versus 0.46% and 0.10% in the adults — consistent in both
blocks, up to 16× range. **It should not be reported as an ageing phenotype
without more work.**

Composition correlation between sections is **r = 0.934 within group vs 0.905
between group** — group differences are barely above section-to-section noise.
And the AP-informative populations do not shift coherently:

| population | AP hint | F536 | G_073 | M493 | M399 |
|---|---|---|---|---|---|
| Oxt/Avp magnocellular | rostral | 0.33 | 0.44 | 0.54 | 0.62 |
| ARC Pomc | mid | 1.91 | 2.37 | 1.88 | 1.96 |
| VMH-like Rasgrf2 | mid | 4.22 | 5.74 | 5.23 | 5.41 |
| DMH Grp/Ppp1r17 | caudal | 1.57 | 1.22 | 0.46 | 0.10 |
| LHA Hcrt | caudal-lateral | 1.41 | 1.38 | 2.63 | 2.61 |

If the aged sections were simply cut more caudally, DMH **and** Hcrt should both
rise. DMH rises 16× while Hcrt *falls* — so this is not a clean AP shift but a
combination of cut level and ROI placement. The ROI rectangles differ in area by
40% (3.68–5.15 mm²) and are hand-drawn, so how much DMH versus LHA each box
captures is partly a drawing decision.

**The fix is the nucleus assignment step**: express each population as a
fraction of cells *within its own nucleus* rather than of the whole ROI, which
removes box geometry from the denominator entirely.

---

## What this means for the paper

The atlas, the QC and the Galr1 measurement are solid and publication-grade.
The age contrast is not yet: for Galr1 specifically the honest statement is
*"no detectable systematic change across MBH cell types in 2 vs 2 middle-aged
versus adult male mice,"* with the `GABA Gal/Galr1` decrease and the tanycyte
`Gal` decrease as the leads worth pursuing.

Next, in order of value:

1. **Nucleus assignment** (plan §5) — removes the ROI-geometry confound and is
   the only way to answer "is it the DMH?" properly.
2. **The other 6 sections** of these same animals — each animal has 2–3
   sections; using them stabilises per-animal estimates and directly tests
   whether the DMH result survives a different cut.
3. **Cohort B** as independent replication.
4. **More animals.** At n=2 vs 2 the exact blocked permutation floor is p=0.5
   two-sided; a real Galr1 effect of the size seen here would need n≈5–8 per
   group to be demonstrable.
