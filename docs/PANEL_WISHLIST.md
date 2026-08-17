# Supplementary panel wishlist

Targets to add in a future Xenium run. Compiled from the gaps this project ran
into: the current `mBrain_50g` design (247 predesigned + 50 custom = 297 genes)
has **zero coverage of the Galr1 signal-transduction and desensitisation
machinery**, which is the layer where pharmacological resistance actually lives
and the one the resistance analysis could not reach.

Machine-readable version: [`panel_wishlist.csv`](../config/panel_wishlist.csv).

**Budget note.** 10x custom add-on panels are typically 50 or 100 targets. Tier 1
alone is 14 genes and answers the resistance question; Tiers 1–3 are 27 and cover
the galanin system completely; Tiers 1–5 are 46 and fit a 50-plex add-on with
room to spare.

---

## Tier 1 — the resistance axis (14 genes)

Galr1 is Gi/Go-coupled. Its response is set by four things downstream of the
receptor, none of which is currently measurable. **If M617 shows a blunted
response in aged mice and receptor transcript is unchanged — which is what this
project found — the explanation is in this tier.**

| gene | role | why it matters for resistance |
|---|---|---|
| **`Kcnj3`** (GIRK1) | effector channel | Galr1 hyperpolarises via GIRK. Reduced GIRK = weaker response at unchanged receptor |
| **`Kcnj6`** (GIRK2) | effector channel | dominant neuronal GIRK subunit |
| **`Kcnj9`** (GIRK3) | effector channel | modulates GIRK trafficking and desensitisation |
| **`Grk2`** (`Adrbk1`) | receptor kinase | phosphorylates agonist-occupied GPCRs; upregulation is a classic resistance mechanism |
| **`Grk3`** (`Adrbk2`) | receptor kinase | main partner to GRK2 for Gi-coupled receptors |
| **`Grk5`** | receptor kinase | slower, membrane-anchored arm |
| **`Arrb1`** | β-arrestin 1 | desensitisation and internalisation |
| **`Arrb2`** | β-arrestin 2 | as above; the pair determines receptor availability at the membrane |
| **`Rgs4`** | Gα GAP | broadly hypothalamic; accelerates Gi signal termination |
| **`Rgs7`** | Gα GAP (R7) | with Gβ5, sets GIRK deactivation kinetics |
| **`Rgs9`** | Gα GAP (R7) | as above |
| **`Rgs2`** | Gα GAP | activity-regulated; links neuronal state to GPCR gain |
| **`Rgs20`** | Gαo-selective GAP | the most Galr1-relevant RGS by coupling |
| **`Gnb5`** | Gβ5 | obligate partner of R7-family RGS; without it the RGS reads are hard to interpret |

## Tier 2 — transduction context (6 genes)

Distinguishes "fewer transducers" from "faster termination".

`Gnao1` (dominant neuronal Gαo, Galr1's principal transducer) · `Gnai1` ·
`Gnai2` · `Adcy5` (Gi target) · `Pde10a` · `Prkaca`

## Tier 3 — completing the galanin system (7 genes)

**`Galr2` is the single highest-value addition on this list.** It is
**Gq-coupled and excitatory**, the opposite sign to Galr1 and Galr3, and it is
entirely absent from the current panel. Any galanin study without it is reading
half the system — and an apparent loss of galanin effect could be a shift in the
Galr1:Galr2 balance rather than resistance at Galr1.

| gene | role |
|---|---|
| **`Galr2`** | Gq-coupled galanin receptor — **missing, and it inverts the sign of the response** |
| `Galp` | galanin-like peptide, arcuate-specific, a second endogenous ligand |
| `Spx` | spexin, endogenous Galr2/Galr3 ligand (low expressor — flag with 10x) |
| `Kcnj5` | GIRK4, completes the channel family |
| `Rgs17`, `Rgs10`, `Rgs19` | further Gi-family GAPs if budget allows |

## Tier 4 — nucleus delineation (10 genes)

This project could not delineate the DMH properly: the panel has no
DMH-restricted transcription factor, so the nucleus had to be inferred from
`Grp` (a restricted subpopulation) and came out ~480 µm dorsoventrally against a
published 500–700. These would fix that, and would also give the AP-registration
a marker-based anchor instead of a purely morphometric one.

| gene | nucleus it resolves |
|---|---|
| **`Nr5a1`** (SF-1) | **VMH — the definitive marker, currently absent** |
| **`Tbx3`** | **ARC — definitive, currently absent** |
| `Sim1` | PVN / DMH |
| `Foxb1` | mammillary / premammillary, and an AP landmark |
| `Rax` | tanycytes, ventricular zone |
| `Nkx2-1` | ventral hypothalamic territory |
| `Lhx6` | zona incerta — would have settled the ZI call directly |
| `Trh` | PVN, DMH |
| `Pmch` | LHA — currently the LHA rests on `Hcrt` alone |
| `Kiss1` | ARC KNDy — the population the current model assigns worst (7.5%) |

## Tier 5 — M617 response readout (5 genes)

The panel has `Fos`, `Arc` and `Bdnf`. For a drug-response experiment a wider
immediate-early set gives both a faster and a slower window, and `Npas4` is
specific to activity rather than to stress.

`Npas4` · `Egr1` · `Junb` · `Fosb` · `Nr4a1`

## Tier 6 — ageing and neuroinflammation (5 genes)

The ageing signature here rested on `Gfap`, `Cd68` and `Trem2` and came out
partial. Complement and MHC-I are the better-established ageing readouts.

`C1qa` · `C1qb` · `B2m` · `Cd74` · `Ctss`

---

## Suggested build

| build | tiers | genes | answers |
|---|---|---|---|
| minimum | 1 | 14 | is resistance downstream of the receptor? |
| **recommended** | 1 + 3 + 4 | **31** | resistance axis, full galanin system, proper nuclei |
| 50-plex add-on | 1–5 | 46 | the above plus drug-response readout |
| 100-plex | 1–6 + spares | 51 + | all of it, with room for replicates of low expressors |

## Practical notes for panel design

- **Flag low expressors with 10x at design time.** `Galr2`, `Spx`, `Galp` and
  `Rgs20` are expressed at low levels; 10x will model detection efficiency and
  may advise against or suggest alternative probe sets. Better to know before
  the run than to discover a target is undetectable afterwards.
- **Symbol aliases.** `Grk2` = `Adrbk1` and `Grk3` = `Adrbk2` in older
  annotations; supply both so the design tool resolves them.
- **Keep the existing 297 unchanged** if you want the new cohort comparable to
  this one — a custom add-on preserves the current targets, whereas a redesign
  would break the cross-cohort comparison this project already set up.
- **Segmentation.** If the next run is on a newer Xenium chemistry, prefer
  multimodal segmentation over 5 µm nucleus expansion. The 40% "mixed
  GABA/glutamate" rate in this dataset is transcript bleed-through, and it
  degrades exactly the co-expression calls a receptor-machinery panel depends on.
