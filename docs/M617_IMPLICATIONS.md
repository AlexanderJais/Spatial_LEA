# What this dataset says for the M617 experiment

M617 is a Galr1-selective agonist. For a treatment study the useful question is
not "what changed with age" — it is **which cells carry the receptor, what they
are, and where they sit**. That is descriptive: 102 551 cells across seven
nuclei, no age contrast, and none of the n=2 limitation that constrains
everything else in this project.

Source: `22_m617_target_atlas.py`, tables in `results/m617_atlas/`.

---

## 1. Where M617 can act

**17.7% of MBH-region cells are Galr1-positive**, and the receptor is highly
concentrated: the **top 5 populations carry 53% of all Galr1 signal, the top 10
carry 69%.**

| population | transmitter | cells/section | % Galr1⁺ | share of all Galr1 |
|---|---|---|---|---|
| **DMH GABA Gal/Galr1** | mixed/GABA | 245 | **68.2%** | **17.6%** |
| VMH-like Glut Rasgrf2 | glutamatergic | 540 | 43.4% | 12.3% |
| LHA GABA Nts/Gal | GABAergic | 295 | 45.5% | 10.5% |
| ZI GABA Cacna2d2 | GABAergic | 707 | 14.3% | 6.5% |
| DHA/PH Glut Otp/Nrn1 | mixed | 651 | 14.3% | 6.2% |
| DMH GABA Nell1/Pcsk5 | mixed | 68 | 62.8% | 3.9% |
| DMH Glut Prdm8/Cbln1 | mixed | 34 | **86.0%** | 3.9% |
| DMH GABA Nts/Gal | mixed | 55 | 61.6% | 3.2% |

By nucleus:

| nucleus | % cells Galr1⁺ | Galr1⁺ cells/section | share of all Galr1 |
|---|---|---|---|
| **DMH** | 35.8% | 449 | **37.2%** |
| **LHA** | 25.0% | 336 | **23.0%** |
| VMH | 20.4% | 368 | 17.6% |
| DHA/PH | 11.1% | 160 | 10.5% |
| ZI | 11.0% | 162 | 10.2% |
| ME/3V | 3.0% | 20 | 0.9% |
| **ARC** | **3.3%** | **18** | **0.7%** |

**The ARC is effectively Galr1-free.** It carries 0.7% of the region's receptor
signal — so M617 will barely act on the arcuate directly, including on AgRP
neurons. Any arcuate phenotype after M617 is more likely to be a network
consequence of DMH/LHA inhibition than a direct effect. (This also means my
earlier focus on `Gal` in AgRP neurons, while statistically defensible, was the
wrong thing to emphasise for your experiment — that finding is about the ligand,
in a nucleus the drug scarcely touches.)

## 2. What the drug will do — Galr1 is Gi-coupled

An agonist **inhibits** the cell carrying the receptor, so the transmitter
identity of Galr1⁺ cells decides whether M617 silences a circuit or releases it.

Using a strict call (≥2 counts of one class, zero of the other; the permissive
call puts 40% of cells in "mixed", which is transcript bleed-through from the
5 µm nucleus-expansion segmentation rather than real co-expression):

| class | cells | carries |
|---|---|---|
| GABAergic | 28 359 (27.7%) | **29.4% of Galr1** |
| glutamatergic | 17 133 (16.7%) | 12.0% of Galr1 |
| ambiguous | 22 569 (22.0%) | 35.7% of Galr1 |

Among confidently-called neurons, **Galr1 sits on GABAergic cells about 2.5× more
than on glutamatergic ones**. The dominant predicted action of M617 in this
region is therefore **disinhibition of downstream targets**, not direct
silencing. Worth keeping in mind when interpreting a cFos or behavioural
readout: the responding cells may be ones that never carried the receptor.

## 3. Three mechanistic features worth designing around

**The top target is an autoreceptor configuration.** DMH GABA Gal/Galr1 neurons
are **68.2% Galr1⁺ and 94.9% `Gal`⁺** — they make galanin and carry its
receptor. M617 there mimics an endogenous autoinhibitory loop, and these cells
already sit in high local ligand, so the agonist has less headroom than at a
naive receptor. If M617 has a weaker effect than expected in the DMH, saturation
is a candidate explanation, not a failure of delivery.

**Galr3 co-expression flags where selectivity matters.** Galr3 is also
Gi-coupled and also on the panel:

| population | % Galr1⁺ | % Galr3⁺ | % both |
|---|---|---|---|
| DMH Glut Prdm8/Cbln1 | 86.0 | 59.9 | **55.7** |
| LHA Glut Prdm8/Cbln1 | 49.2 | 29.5 | 23.6 |
| VMH-like Glut Rasgrf2 | 43.4 | 22.5 | 11.1 |

DMH Glut Prdm8/Cbln1 is the population where any residual Galr3 activity of
M617 would be least distinguishable from its Galr1 action. It is small
(34 cells/section) but the most receptor-dense population in the region.

**Endogenous tone varies enormously between targets.** Among the top targets,
`Gal` positivity ranges from 82–95% (DMH Gal/Galr1, LHA Nts/Gal, DMH Nts/Gal) to
**9.9% in VMH-like Glut Rasgrf2**. The VMH target is receptor-rich but
ligand-poor — the closest thing here to an unoccupied receptor, and the place an
exogenous agonist should produce the largest change from baseline.

## 4. The age question, framed usefully

The Galr1⁺ fraction of each nucleus, per animal:

| nucleus | aged mean | adult mean | ratio |
|---|---|---|---|
| DMH | 33.0% | 38.4% | 0.86 |
| LHA | 27.4% | 22.6% | 1.21 |
| VMH | 20.8% | 20.1% | 1.03 |
| ZI | 14.6% | 9.8% | 1.49 |
| DHA/PH | 13.7% | 8.5% | 1.61 |
| ARC | 3.1% | 3.6% | 0.86 |
| ME/3V | 3.2% | 2.7% | 1.19 |

Ratios sit near 1 in every nucleus, and the two that deviate most (ZI, DHA/PH)
are the two with the worst sampling balance between animals (CV 49% and 34%), so
those are the least trustworthy numbers in the table.

**The practically useful statement:** *the size of the drug-target population is
broadly comparable between adult and middle-aged mice in every nucleus M617 can
act on.* So if M617 produces a **different response** in aged animals, this
dataset argues the explanation is unlikely to be "there are fewer receptors" —
it points downstream, to receptor coupling, effector state, or circuit context.
That is a genuinely useful constraint to have before the experiment reads out,
and it is a claim this design can support, unlike a per-gene age difference.

---

## Practical suggestions for the M617 study

1. **Read out in the DMH and LHA first.** They carry 60% of the region's Galr1
   between them. The ARC is not where the drug acts.
2. **Expect disinhibition.** Most receptor sits on GABAergic neurons, so plan
   readouts that can detect increased activity in cells that are *not*
   Galr1-positive.
3. **Use `Gal`-negative, receptor-rich populations as the sensitive test case.**
   VMH-like Glut Rasgrf2 (43% Galr1⁺, 10% `Gal`⁺) has the most headroom.
4. **Include a Galr3 control** if any readout centres on DMH Glut Prdm8/Cbln1.
5. **If aged and adult respond differently, look downstream of the receptor** —
   the target population itself is comparable.
6. **AP-match the sections** you compare, in the treatment study too. The
   variance components here show section-to-section anatomical variation
   (SD 0.265 log2) exceeds between-animal biological variation (0.152) — that
   applies to a treatment contrast just as much as to an age contrast.

---

# The galanin-resistance hypothesis

Resistance predicts a specific shape — the one that defines leptin resistance:
**ligand rises while receptor stays flat or falls**, so the ligand-to-receptor
balance shifts and the same signalling needs more ligand. Tested gene by gene
that is hopeless at n=2; tested as an **ensemble** it is not, because the unit
becomes the population rather than the animal. `23_galanin_resistance.py` runs
it over 64 populations with AP-corrected effects, calibrated against all 293
panel genes.

## Six tests, all negative

| test | result |
|---|---|
| **Ensemble shift** — is `Gal` displaced up and `Galr1` not? | `Gal` median LFC **+0.007** (48th percentile of the panel), `Galr1` **−0.018** (55th). Both sit in the middle of the panel distribution. |
| **Resistance index** — LFC(`Gal`) − LFC(`Galr1`) per population | Positive in **12/29 populations (41%)**, median **−0.15**. Against random gene pairs from the same populations: one-sidedness p = 0.58, magnitude p = 0.50 — indistinguishable from noise. |
| **Paired within population** | Wilcoxon p = 0.73. |
| **Spatial coupling** — local `Gal` available to each `Galr1`⁺ cell | aged 1.09–1.35, adult 1.09–1.42. Completely overlapping. |
| **Activity set-point** — failing inhibition should raise activity in `Galr1`⁺ cells | `Fos` aged/adult ratio **0.73**, `Arc` 0.84, `Bdnf` 0.93, `Nrn1` 0.97. If anything the drift is *downward* — the opposite of disinhibition. |
| **Autoreceptor coverage** — do `Gal`⁺ cells keep their own receptor? | 27.4% aged vs 26.4% adult, ratio **1.04**. No loss. |

**There is no transcriptional signature of galanin resistance in this dataset.**

The ensemble tests carry more weight than the per-gene ones. Their unit is the
population (n = 29–52), not the animal, so they are far better powered than the
~11-fold minimum detectable effect that limits a single age contrast — and they
are calibrated against random gene pairs rather than a nominal p-value. This is
a more informative negative than the per-gene nulls elsewhere in the project.

## Why that does not refute the hypothesis

**Pharmacological resistance is usually post-transcriptional, and Xenium sees
transcripts.** The mechanisms that actually produce a blunted agonist response —
receptor internalisation, GRK/β-arrestin desensitisation, G-protein uncoupling,
RGS upregulation, GIRK channel changes — are largely invisible at the mRNA
level, and the few with transcriptional components (`Rgs*`, `Grk*`, `Arrb*`,
`Kcnj*`) are not on this 297-gene panel at all.

So a null here means something narrower but still useful: **if galanin
resistance exists in these animals, it is not because receptor transcript is
lost, nor because the ligand-to-receptor transcript balance has shifted.**

## What this gives the M617 experiment

This is the constructive form of the negative, and it is a genuine constraint to
have in hand *before* the treatment reads out:

> If M617 produces a blunted response in aged mice, this dataset argues the
> cause is **not** reduced receptor availability, **not** loss of receptor from
> galanin-producing cells, and **not** a ligand–receptor imbalance in
> transcript abundance. The mechanism would have to sit downstream of receptor
> expression — coupling, desensitisation, effector state or circuit context.

That is exactly where the field usually finds resistance, and it narrows the
follow-up considerably.

### What would actually test resistance

1. **Dose–response to M617** — a rightward EC50 shift in aged animals is the
   definition of resistance, and it is the experiment you are already running.
   A Bmax/Kd separation by receptor binding would then say whether it is
   receptor number or receptor affinity/coupling.
2. **Proximal signalling readout, not transcript** — pERK, pCREB, or GIRK
   currents in `Galr1`⁺ DMH neurons after M617. This is the layer where
   resistance lives.
3. **A panel that can see the machinery** — if a future spatial run is planned,
   adding `Rgs2`/`Rgs4`, `Grk2`/`Grk5`, `Arrb1`/`Arrb2` and GIRK subunits
   (`Kcnj3`, `Kcnj6`, `Kcnj9`) would make the desensitisation axis visible.
   The current panel simply cannot address it.
4. **Read out where the receptor is** — DMH and LHA carry 60% of it, ARC 0.7%.

## One positive finding worth keeping

`Galr1`⁺ cells sit in **`Gal`-enriched neighbourhoods in all 12 sections**
(enrichment 1.09–1.42, every section above 1). Ligand and receptor are
spatially organised rather than independently distributed — the receptor is
positioned where the ligand is made. That is a Tier 1 descriptive result,
replicated in every section and independent of age, and it is a reasonable
figure for the paper regardless of how the resistance question resolves.
