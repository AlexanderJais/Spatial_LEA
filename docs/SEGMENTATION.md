# Re-segmentation with Baysor

Every result in this repository was computed on the vendor segmentation: a DAPI
nucleus dilated by a fixed 5 µm. This document tests whether that choice was
load-bearing.

**Answer: no.** Baysor changes the objects a little and the conclusions not at
all. The one thing it does improve — cells carrying GABAergic and glutamatergic
markers at once — it improves by about a tenth, not by the amount the raw
comparison first suggested.

Scripts `29`–`34`, figure `results/figures/baysor_segmentation.png`.

---

## Why it was worth checking

Measuring the vendor objects showed the 5 µm dilation fits neither cell class in
this tissue:

| | nucleus | vendor cell | real soma |
|---|---|---|---|
| neurons | 8.4 µm | 15.8 µm | 15–25 µm — **under-captured** |
| glia | 6.5 µm | 13.9 µm | 8–12 µm — **over-expanded** |

Under-capturing neurons loses the transcripts in their cytoplasm; over-expanding
glia imports the neighbours'. Both errors push every cell toward the local
average, and the visible symptom is the ~25–36% of marker-positive cells that
read as GABAergic *and* glutamatergic — which no neuron in this tissue should.

The dilation also governs most of the data. **Only 39% of molecules in the
hypothalamic window fall inside a nucleus.** The other 61% was assigned by
geometry alone.

## What was run

`transcripts.parquet` for the closest AP-matched cross-age pair — G073_1 (aged,
AP −1.38) and M399_3 (adult, AP −1.33) — cropped to the same anatomical window
every existing analysis uses. The crop validates against the vendor's own cell
count: 13,269 DAPI nuclei against 13,191 vendor cells for G073_1, 14,640 against
14,497 for M399_3.

Prior = the DAPI nuclei at confidence 0.5, zero elsewhere. That keeps the
measured nuclei and leaves the cytoplasmic extent — the part we do not trust —
to be inferred.

Baysor v0.7.1 built from source (`scripts/setup_baysor.sh`); 4.9M and 5.4M
molecules, ~65 min per section on 4 threads.

### Choosing `scale`, and a criterion I got wrong first

`scale` is the expected cell radius — the one Baysor parameter carrying the same
kind of assumption as the 5 µm dilation, so it was swept rather than asserted.

My first criterion was **wrong**: I ranked settings by cells-per-DAPI-nucleus
approaching 1.0. That picks 10 µm, which merges 29% of measured nuclei into
shared cells and leaves the mixed rate at 37%, no better than the vendor. The
error is that a cell can cross the section plane while its nucleus does not, so
1.0 was never the target.

Replaced with two errors that need no model — **two measured nuclei in one
cell**, and **one nucleus scattered across several** — plus a self-consistency
check from measured geometry (nucleus 7.12 µm, effective section thickness
5.54 µm from the run's own metadata): expected cells per visible nucleus is
`(d + t) / (n + t)`, with `d` taken from the cell sizes that setting produced.

| scale | merged | split | cell Ø | cells/nucleus | vs expected | mixed |
|---|---|---|---|---|---|---|
| 4 µm | 2.9% | 9.7% | 8.1 µm | 2.69 | **2.50×** | 19.3% |
| **6 µm, std 50%** | **10.5%** | **13.3%** | **11.3 µm** | **1.68** | **1.26×** | **25.5%** |
| 8 µm | 21.2% | 17.2% | 13.9 µm | 1.29 | 0.84× | 31.5% |
| 10 µm | 29.4% | 17.7% | 14.8 µm | 1.10 | 0.69× | 37.2% |

4 µm wins on raw error rate only because it barely leaves the nucleus — 8.1 µm
cells against a 7.12 µm nucleus, 2.5× more cells than its own size can account
for, half the counts per cell. It returns nuclei relabelled as cells. 6 µm is
the only self-consistent setting.

---

## Results

### 1. The improvement is real, and smaller than it first looks

Baysor assigns **99% of molecules against the vendor's 88%**, and produces 1.85
cells per DAPI nucleus at 87 median transcripts against the vendor's 218.

Comparing all Baysor cells to all vendor cells gives 25.6% → 17.8% mixed
identity. **That comparison is not fair.** The mixed rate rises steeply with
sequencing depth — from 1.4% at 30 transcripts per cell to 41% above 500 — and
Baysor's cell set includes ~11,000 nucleus-free objects at a median of 37
transcripts, which dilute the rate for arithmetic reasons.

Two ways of removing the confound:

| | G073_1 | M399_3 |
|---|---|---|
| vendor | 25.6% | 26.3% |
| **Baysor, cells on a measured nucleus** (189/196 vs 218/219 transcripts) | **23.3%** | **24.0%** |
| Baysor, all cells (87 transcripts) | 17.8% | 18.6% |

Within bins of equal transcript count, Baysor is lower in 6 of 7 bins for
G073_1 and 5 of 7 for M399_3, losing only in the two shallowest bins.

**So: a ~9% relative reduction in mixed identity, not the ~30% the raw numbers
imply.** Re-segmentation improves the problem; it does not solve it. Roughly a
quarter of marker-positive neurons still read as both transmitters, which points
at the 5 µm dilation being a smaller contributor than assumed — more likely
genuine transcript diffusion, optical bleed between adjacent cells, or real
co-expression the panel cannot resolve.

### 2. The extra cells are mostly fragments

44% of Baysor cells sit on no DAPI nucleus. Geometry predicts some of this — a
10.7 µm cell in a 5.54 µm section with a 7.12 µm nucleus gives 1.28 cells per
nucleus — but the observed 1.85 is more than that accounts for.

| | with a nucleus | without |
|---|---|---|
| n cells (M399_3) | 14,860 | 11,577 |
| median transcripts | 196 | **37** |
| median diameter | 13.7 µm | 7.8 µm |

**This matters because it created a false result.** Across all Baysor cells,
astrocytes appear to gain +7.5 to +9.8 percentage points of the population.
Restricted to cells with a nucleus — the only objects matched to what a vendor
cell is — the shift is **+0.7 to +1.2 pp**, and every other cell type moves less
than 1.3 pp. The astrocyte gain was 34% of the nucleus-free fragments being
labelled astrocyte.

**Composition is unchanged.**

### 3. Every Galr1 conclusion survives

This is the part that matters for the project.

| | vendor | Baysor |
|---|---|---|
| DMH share of Galr1 signal (G073_1) | 35.0% | 34.2% |
| LHA | 17.5% | 17.9% |
| VMH | 10.2% | 10.1% |
| ARC | 0.5% | 0.5% |
| GABA:Glut ratio of Galr1 signal | 3.74× | 3.85× |
| top carrier: GABA Gal/Galr1 | 23.3% | 23.7% |
| GABA Nts/Gal | 17.6% | 17.9% |

Nothing moves by more than a percentage point. The receptor's DMH/LHA
localisation, its near-absence from the ARC, and its 2.5–3.9× bias toward
GABAergic cells were **not artefacts of the dilation**.

The per-cell detection rate falls from 14.0% to 8.8% — expected, since Baysor
cells are smaller and carry fewer transcripts each. That is a change in the
denominator, not in the biology.

### 4. No age signal either way

Galr1⁺ cells: aged 8.8% vs adult 8.8% under Baysor; 14.0% vs 13.7% under the
vendor. Mixed-identity rate: 17.8% vs 18.6% (Baysor), 25.6% vs 26.3% (vendor).

One section per age cannot separate age from animal, section plane or run. These
numbers say re-segmentation did not move the contrast — not that the contrast is
real.

---

## What this changes

**Nothing needs redoing.** The vendor segmentation was good enough for every
claim in `RESULTS.md`. The Galr1 localisation results are robust to the
segmentation method, which is worth one sentence in a methods section and a
supplementary figure.

The honest caveat to carry forward: ~24% of marker-positive cells still carry
both transmitter markers after re-segmentation, so **single-cell transmitter
identity in this dataset remains unreliable at the individual-cell level**.
Population-level statements — which nucleus, which cell type, what fraction of
signal — are unaffected, because the mixing is roughly symmetric and does not
move the aggregates.

Two sections were run, not twelve. Extending to all 12 would cost ~13 hours of
compute and would not change these conclusions, since the two sections agree
closely on every measure and neither differs from the vendor in any way that
matters.
