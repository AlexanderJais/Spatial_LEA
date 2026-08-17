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

## 6. Unbiased search — is ANY galanin population regulated?

`11_galanin_search.py` carves the MBH four ways so an effect that is not aligned
to the cell-type labels can still surface: annotated type × nucleus (82
populations), fine Leiden subclusters (54), the `Gal`⁺ and `Galr1`⁺ cells in
their own right (6), and pure spatial bins — nucleus × mediolateral tertile ×
dorsoventral tertile (27). 179 (population, gene) combinations tested.

**The calibration is what makes the answer interpretable.** For every population
the same statistic is computed for all 297 panel genes, so the galanin genes are
ranked against a null built from the same cells, animals and blocked design.

> **~9% of the panel clears "consistent in both blocks and ≥1.5× in the weaker
> one" in a typical population** — about 27 genes per population by chance.
> 13 of 179 galanin tests cleared it, *fewer* than the ~16 expected. The galanin
> system as a whole is not enriched for age regulation.

Adding the third filter — complete separation of aged from adult sections within
both blocks — leaves 5 candidates, and `12_candidate_null.py` calibrates those
against the panel too:

| population | gene | genes passing all 3 filters | effect rank | gap vs largest within-animal swing | r vs AP axis |
|---|---|---|---|---|---|
| ARC \| ARC Agrp/Npy | **`Gal`** | **9 / 297 (3.0%)** | 27/297 | **+1727 vs 1480 CPM — gap wins** | −0.77 (p=0.003) |
| DMH \| GABA Cacna2d2 | `Gal` | 17 / 297 (5.7%) | 19/297 | +6426 vs 4105 CPM — gap wins | −0.82 (p=0.001) |
| DMH \| lateral-dorsal bin | `Galr3` | 16 / 297 (5.4%) | 34/297 | +159 vs 224 CPM — **swing wins** | −0.64 (p=0.024) |

### The answer

**No `Galr1` population survives.** Not one, in any of the four carvings.

**The best galanin candidate in the dataset is `Gal` in ARC Agrp/Npy neurons** —
roughly 1.7× higher in aged, found independently by the annotated-cell-type and
fine-cluster carvings, one of only 9 genes out of 297 to pass all three filters
in that population, and the group gap exceeds the largest within-animal section
swing. That last point is what the DMH Galr1 lead failed.

It is still not clean. `Gal` in these neurons correlates r = −0.77 with the
rostro-caudal axis, so age and position remain collinear at n=2 per group. The
honest status is: **the strongest galanin-system candidate this data can
produce, and the right one to test in a properly powered, AP-matched cohort** —
not an established ageing effect.

The other 8 genes passing in ARC Agrp/Npy are `Cldn5`, `Pglyrp1`, `Igf1`,
`Penk`, `Prss35`, `Pomc`, `Fign`, `Kcnmb2` — a plausible ageing set, which is
mild independent support that this population is genuinely changing rather than
the filter being noise.

## 8. How big is the DMH, really?

Short answer: **the delineation was too small, and even the improved one is at
the low end.**

The Gaussian model fits each nucleus to marker-defined anchors. ARC had 1 856
anchor cells across 4 cell types and VMH 5 525 across 3 — but **DMH had 343
cells from a single type**, `DMH Grp/Ppp1r17`. `Grp` marks a restricted DMH
subpopulation, so the ellipse is the extent of the Grp⁺ core, not the nucleus.

`13_spatial_domains.py` re-derives the domains from tissue structure instead:
every cell is described by the cell-type composition of its 30 nearest spatial
neighbours, and those descriptions are clustered into 14 niches over the whole
hypothalamic window (164 943 cells). Domain 13 is the DMH — `GABA Gal/Galr1`
23%, `DMH Grp/Ppp1r17` 9%, sitting dorsal to the VMH domains.

| definition | cells/section | \|ml\| p5–95 | dv p5–95 | size |
|---|---|---|---|---|
| Gaussian ellipse (Grp-anchored) | 1 103 | 136–542 | 985–1364 | 1084 × **380** µm |
| Domain 13 (tissue-derived) | 1 037 | 71–513 | 868–1337 | 1026 × **469** µm |
| published mouse DMH | — | ~400–600/side | — | ~800–1200 × **500–700** µm |

The tissue-derived domain is 23% taller and captures more of the relevant
populations (`GABA Gal/Galr1` 68.6% vs 51.1%; `DMH Grp` 85.7% vs 78.5%), but
**both fall short dorsoventrally**, and the two definitions agree only at
**Jaccard 0.44** — 37% of the domain lies outside the ellipse and 41% of the
ellipse outside the domain. That disagreement brackets the real uncertainty.

The domain is also far less stable section to section (CV 31% vs 10%, ranging
259–1 322 cells), which is partly genuine: the DMH's size really does change
with rostro-caudal level. That is a worse denominator but a truer boundary.

**Neither method can settle this, because both infer the DMH from its contents.**
With the DAPI images now available, registration to the Allen CCF would give a
boundary defined by anatomy rather than by cell type, and is the way to close
this properly.

## 9. AP matching — and what it proves

`14_ap_matching.py` scores each section's rostro-caudal level from **morphometry**
— third-ventricle length, tissue half-width at four dorsoventral levels, midline
extent, VMH position and width, window cell count. Morphometry rather than
cell composition on purpose: matching on composition would risk regressing out
the very age differences being tested.

The score orders sections by cut order in **4 of 4 animals**, and unlike the
composition axis, the aged and adult ranges now **overlap**.

**But a four-animal matched set is not available.** F536's sections score +2.2
to +4.2 while the other three animals span −2.1 to −0.1 — that block was cut at
a different level. The tightest one-per-animal set still spans 3.50 of the full
6.23. Within blocks, however:

| block | best pair | AP gap |
|---|---|---|
| B2 (G_073 aged vs M399 adult) | **G073_1 vs M399_3** | **0.04** |
| B1 (F536 aged vs M493 adult) | F536_2 vs M493_3 | 2.36 |

For scale, adjacent sections within one animal differ by ~0.51. So **B2 pairs
almost exactly; B1 cannot be matched at all.**

### The discrimination

Computing each candidate across *all 18 possible section pairings* and
correlating the effect size with the AP gap separates artifact from signal:

| candidate | effect vs AP mismatch | best-matched pair | poorly-matched mean | verdict |
|---|---|---|---|---|
| `Galr1`, DMH GABA Gal/Galr1 | **r = −0.69 (p = 0.002)** | −0.17 | −0.62 | **artifact** — the effect grows with mismatch |
| `Gal`, ARC Agrp/Npy | r = +0.12 (p = 0.63) | **+0.74** | +0.68 | **survives** — indifferent to matching |
| `Gal`, DMH GABA Cacna2d2 | r = −0.10 (p = 0.70) | +0.66 | +0.36 | survives, noisier |

This is the cleanest result in the study. The DMH Galr1 "effect" is a direct
function of how badly the compared sections are matched — at near-perfect
matching it collapses from −0.62 to −0.17, and its detection-rate version
reverses sign. **`Gal` in ARC Agrp/Npy neurons does the opposite**: it is
unchanged by matching, positive in every one of the 18 pairings, and largest at
the best-matched pair.

### Practical answer

Yes, a matched core set exists — **but only within block B2**. The usable
comparison is `G073_1` vs `M399_3`, matched to 0.04 AP units. Block B1 should be
reported separately, or F536 excluded and replaced. Any future cohort should
collect sections spanning a wider AP range per animal so matching is possible
across all animals, not just two.

## 11. Extended parcellation — and how big the DMH really is

The three-nucleus Gaussian model was replaced by a 14-domain parcellation
derived from tissue structure, then mapped onto named nuclei by
nucleus-diagnostic marker z-scores (`18_extended_nuclei.py`). Four proposed
domain groupings were tested rather than assumed:

| proposal | verdict | evidence |
|---|---|---|
| VMH = domains 7 + 4 | **confirmed** | both Slc17a6-dominant (z +2.13, +1.76); 7 is dorsomedial (\|ml\| 246, Rasgrf2 58%), 4 ventrolateral (\|ml\| 443, Calb1 41%) |
| ZI = domains 10 + 0 | **0 yes, 10 no** | 0 is GABAergic (z +1.77, Cacna2d2 48%) as ZI must be; 10 is glutamatergic (Slc17a6 z +1.23, ZI z −0.10) → dorsal hypothalamic area |
| ARC = domains 12 + 9 | **9 yes, 12 is ME/3V** | 9 scores ARC z +2.96 (Agrp, Pomc, Ghrh, Tac2); 12 is the 3V floor, 64% tanycyte. Legitimate as a combined "ARC+ME" unit if labelled so |
| DMH = domains 13 + 8 | **13 yes, 8 mostly LHA** | 8 scores LHA z +2.11, sits at \|ml\| 742 (p90 1011), and is 64.5% Hcrt⁺ laterally vs 40.8% medially |

Domain 8 is one niche by composition but two structures by anatomy, so it is
split at \|ml\| = 500 µm: the medial fringe (Grp⁺ 18.6% vs 11.9% laterally) joins
the DMH, the lateral portion becomes the LHA.

### The DMH size question, answered

| definition | cells/section | size |
|---|---|---|
| Gaussian ellipse (Grp-anchored) | 1 103 | 1084 × **380** µm |
| domain 13 alone | 1 037 | 1026 × **469** µm |
| domain 13 + medial fringe of 8 | 1 254 | 1000 × **477** µm |
| published mouse DMH | — | ~800–1200 × **500–700** µm |

Adding the medial fringe of domain 8 grows the DMH by 21% in cells but barely in
height (469 → 477 µm), because domain 8 occupies the same dorsoventral band. So
the DMH really does come out at ~480 µm in these sections — at the low end of,
but not far outside, the published range for a single coronal plane. The
mediolateral width (1000 µm bilateral) matches published values well.

### Full parcellation

| nucleus | cells | per section | width | height | between-animal CV |
|---|---|---|---|---|---|
| ARC | 6 639 | 553 | 688 µm | 310 µm | **10.2%** |
| ME/3V | 8 050 | 671 | 214 µm | 1278 µm | **14.8%** |
| VMH | 21 686 | 1 807 | 1194 µm | 521 µm | **11.2%** |
| DMH | 15 053 | 1 254 | 1000 µm | 477 µm | 35.4% |
| LHA | 16 080 | 1 340 | 2204 µm | 714 µm | 30.3% |
| ZI | 17 703 | 1 475 | 2227 µm | 548 µm | 49.2% |
| DHA/PH | 17 340 | 1 445 | 2876 µm | 990 µm | 33.5% |

**The CV column decides which nuclei can be compared across animals.** ARC, ME
and VMH are sampled evenly (10–15%); DMH, LHA, ZI and DHA vary 30–49% because
their cross-sectional area changes steeply with rostro-caudal level. M399 in
particular contributes 1 807 DMH cells against 4 100–4 800 for the other three,
and 7 162 ZI cells against 2 044 for F536. **Age comparisons in the AP-limited
nuclei are confounded with section level before any gene is examined.**

### Galanin system in the new territory

Extending to LHA, ZI and DHA/PH added no robust effect. Nine candidates cleared
the blocked filters; applying the AP test and the sampling-balance criterion
leaves **one**:

| nucleus | cell type | gene | LFC B1 | LFC B2 | best-matched | same-sign pairings | r vs AP gap |
|---|---|---|---|---|---|---|---|
| ME/3V | ARC Agrp/Npy | `Gal` | +1.18 | +1.23 | **+1.46** | 18/18 | −0.41 (p = 0.09) |

This is the same finding as §6, localised more precisely: the `Gal` increase in
AgRP neurons is **strongest in the ones bordering the median eminence** (~2.3×
there versus ~1.7× in the arcuate as a whole), and it is *larger* at the
best-matched pairing than on average — the opposite of an artifact.

Two candidates were caught by the AP test and are anatomy, not age: `Galr3` in
DHA/PH glutamatergic neurons (r = +0.56, p = 0.016) and in DHA/PH overall
(r = +0.64, p = 0.004). Several others pass the AP correlation but collapse at
the best-matched pairing (DHA/PH `Gal`: +2.02/+2.79 overall, +0.37 matched) and
sit in nuclei with 30–49% sampling CV, so they are not pursued.

## 12. Ageing signature

`Gfap` (+0.24), `Cd68` (+0.25), `Trem2` (+0.30 log2) up consistently in both
blocks; `Spp1`, `Cd44`, `Igfbp5` inconsistent, `Ly6a` down. A partial glial
signature, as expected at 15–16 months.

---

## What this means for the study

The atlas, the anatomical registration and the Galr1 measurement are solid. The
biological claim is not there yet, and the reason is now specific and fixable
rather than vague:

0. **Test `Gal` in ARC Agrp/Npy neurons in a properly powered, AP-matched
   cohort.** It is the one galanin-system candidate that survived every filter
   the data supports. Everything else below serves this.
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
