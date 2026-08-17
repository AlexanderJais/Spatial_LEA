# Results — Galr1 in the ageing MBH

Cohort A: 4 male mice (2 aged ~15–16 mo, 2 adult ~7 mo), **12 sections, 3 per
animal**, anatomically registered MBH.

Pipeline: `02_qc` → `03_cluster` → `04_annotate` → `07_calibrate_anchors` →
`08_anatomy` → `09_mbh_atlas` → `10_galr1_nucleus`. Tables under `results/`.

---

## Headline

**Galr1 shows no age difference that survives the additional sections.** The one
lead from the 4-section analysis — a Galr1-high GABAergic population — turned
out to sit in the **DMH**, which is the anatomy the study was asking about, and
its apparent ~20% decrease in aged mice is **inseparable from a steep
rostro-caudal gradient**: a single animal's own three sections span more Galr1
than the two age groups differ by.

That is a negative result for the receptor, and a positive result about study
design: **AP-matched sections are required**, and the data now says by how much.

### Two corrections to the earlier 4-section report

| earlier claim | with 12 sections |
|---|---|
| "Tanycyte `Gal` down ~2× in aged — the tightest, most reproducible effect" | **Does not replicate.** Detection LFC B1 −0.40 / B2 +0.11; expression B1 +0.56 / B2 −0.02 — inconsistent in both metrics. |
| "The DMH Grp/Ppp1r17 abundance result is probably geometry" | **Not geometry.** Per-nucleus normalisation made it *larger* (LFC +2.30 / +5.05). But it is still not safe — see §4. |

---

## 1. Anatomical registration replaces the drawn boxes

Each section now carries a coordinate frame built from its own anatomy: `ml` =
signed distance from the third-ventricle midline, `dv` = distance dorsal from
the ventral surface. The ventricle is traced from tanycytes (`Gpr50`) and
ependyma (`Spag16`), seeded at the arcuate centroid; the arcuate also fixes
which way is ventral, so mounting angle stops mattering.

ARC, VMH and DMH are fitted once in that shared frame and applied to all 12
sections. The geometry was recovered, not imposed:

| nucleus | centre \|ml\| | centre dv | anatomy |
|---|---|---|---|
| ARC | 166 µm | 188 µm | ventral midline |
| VMH | 366 µm | 529 µm | dorsolateral to ARC |
| DMH | 336 µm | 1172 µm | dorsal to VMH |

Validation against annotated cells: 99.5% of ARC Agrp/Npy, 92.7% of ARC Pomc,
92.2% of VMH-like Rasgrf2 and 81.1% of DMH Grp/Ppp1r17 land in their expected
nucleus. KNDy/Tac2 neurons do not (7.5%) — a known limitation, not tuned away.

**MBH cells are now balanced across animals** (ARC 3479–3925, DMH 3015–3662,
VMH 5205–5734 per animal), which the 40%-varying rectangles were not.

M399_2, the section the lab flagged as distorted, was independently flagged by
the frame — ventricle elongation 1.25 versus 2.9–5.8 elsewhere, and an inverted
ventral direction before arcuate seeding.

## 2. The Galr1 population is a DMH population

The cluster with `Galr1` as a defining marker — GABAergic, `Gal`/`Galr1`/
`Foxp2`/`Th`, 2 338 cells across 12 sections — is **96.6% DMH** (2 259 DMH, 61
VMH, 18 ARC). So the answer to "is it DMH neurons?" is yes: that is where MBH
Galr1 is concentrated, in a GABAergic galanin-co-expressing population.

Galr1 detection in it is high and well measured: 46–83% of cells depending on
section, against 0.07–0.09% for negative-control probes.

## 3. But the age effect is confounded with rostro-caudal position

Per-animal means look convincing — aged 51.3% and 70.5%, adult 73.5% and 80.9%,
a 16.3-point gap, consistent in both blocks. The per-section values are not:

| animal | group | section 1 | section 2 | section 3 | swing |
|---|---|---|---|---|---|
| M399 | adult | 81.5 | 80.3 | 80.9 | 1.1 |
| M493 | adult | 79.3 | 73.7 | 67.5 | 11.8 |
| F536 | aged | 57.1 | 50.8 | 45.9 | 11.1 |
| G_073 | aged | **80.4** | 71.1 | 59.9 | **20.5** |

`G073_1` at 80.4% is an *aged* section above every adult section except M399's.
Galr1 detection correlates **r = +0.85 (p < 0.001)** with a composition axis
(PC1 of cell-type proportions across sections), and that correlation holds
**within** animals, where age is constant and any slope must be anatomical.

**One animal's own sections span 20.5 points; the age groups differ by 16.3.**

The axis also separates the groups perfectly (aged PC1 −4.74…−0.69, adult
−0.10…+5.44, no overlap), so at n=2 per group age and anatomical position are
fully collinear and cannot be separated. Adjusting for PC1 flips the effect from
−16.3 to +4.1 points, which is over-adjustment rather than evidence of absence.
Either way the claim is not supportable from this data.

## 4. DMH Grp/Ppp1r17 abundance — bigger, still not safe

As a share of DMH cells: F536 17.7%, G_073 7.9% (aged) versus M493 3.6%,
M399 0.2% (adult). Consistent and large in both blocks. But M399 has almost no
Grp⁺ DMH neurons at all (0.23%), and the within-group spread is 2–15×. `Grp`
marks a rostro-caudally restricted DMH subpopulation, so this is the same AP
problem as §3 in a different guise. Per-nucleus normalisation fixed the *box*
confound; it cannot fix the *plane* confound.

## 5. What did hold up: the ligand, not the receptor

Four `Gal` effects are consistent in both blocks and ≥1.5× in the weaker one:

| nucleus | cell type | LFC B1 | LFC B2 |
|---|---|---|---|
| DMH | GABA Cacna2d2 | +1.14 | +1.92 |
| VMH | Astrocyte | +0.96 | +1.15 |
| ARC | Microglia | +0.96 | +0.70 |
| ARC | ARC Agrp/Npy | +0.79 | +0.66 |

`Gal` **up** in aged across neurons, astrocytes and microglia in three different
nuclei, while the receptor stays flat, is a more coherent shape than anything
seen in Galr1 — but it has not been checked against the AP gradient the way §3
was, and must be before it is trusted.

## 6. Ageing signature

`Gfap` (+0.24), `Cd68` (+0.25), `Trem2` (+0.30 log2) up consistently in both
blocks; `Spp1`, `Cd44`, `Igfbp5` inconsistent, `Ly6a` down. A partial glial
signature, as expected at 15–16 months.

---

## What this means for the study

The atlas, the anatomical registration and the Galr1 measurement are solid. The
biological claim is not there yet, and the reason is now specific and fixable
rather than vague:

1. **Match sections on rostro-caudal position.** This is the single highest-value
   change. The composition PC1 built here is already a usable AP proxy; with the
   DAPI images now available, Allen CCF registration would give a calibrated one.
   Without matching, a 20-point Galr1 swing is available for free.
2. **Re-test `Gal`** against that axis, the way Galr1 was tested in §3.
3. **More animals.** With n=2 per group, any axis that happens to separate the
   four animals is perfectly collinear with age. n≈5–8 per group breaks that.
4. **Consider the DMH Gal/Galr1 population as the target** for orthogonal
   validation. It is a real, well-defined, Galr1-rich DMH population and it is
   the right cell type to count by RNAscope in AP-matched sections.
