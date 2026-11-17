# Is the lead design's scoring allele-neutral, and does it matter?

An offline audit of lead [Design #3](../design-3/README.md). It proposes no
antigen, assembles no construct, runs no inference and changes no shared file.

Design #3 scores a candidate peptide from its predicted affinity **in nM** —
`score-config.json` records vaxrank `logistic_normalized(350, 150)` with a
5000 nM cutoff — and gates eligibility at 500 nM. MHCflurry also reports the
percentile of each affinity among random peptides *for that allele*. The two
are not interchangeable, so the audit asks two questions in order: is the nM
scale allele-neutral, and if not, does that explain anything about the
selection?

## 1. The nM scale is not allele-neutral

Recovered from 233,562 recorded affinity/percentile pairs across 13 committed
prediction tables (~44,000 distinct affinity points per allele, 0 monotonicity
inversions over all 47,016 raw rows per allele):

| Allele | Percentile at 500 nM | nM at 2% rank |
|---|---:|---:|
| HLA-A\*01:01 | 0.70 | 2609 |
| HLA-B\*08:01 | 0.73 | 2482 |
| HLA-B\*15:08 | 0.81 | 2310 |
| HLA-C\*07:02 | 1.85 | 569 |
| HLA-C\*01:02 | 1.85 | 560 |
| HLA-A\*33:03 | 2.33 | 389 |

The 500 nM gate is 3.3-fold stricter in rank terms for HLA-A\*01:01 than for
HLA-A\*33:03. Read the other way, the conventional 2% weak-binder threshold is
2,609 nM for HLA-A\*01:01 and 389 nM for HLA-A\*33:03, a 6.7-fold spread. The
six alleles fall into two clusters, not six distinct levels: A\*01:01, B\*08:01
and B\*15:08 near 0.7–0.8%, and C\*07:02, C\*01:02 and A\*33:03 near 1.8–2.3%.

The effect is easiest to see on peptides that have nothing to do with this
design. Across the **38,927 arbitrary peptides** in these tables that were
scored against all six alleles, the allele holding the best *nM* prediction is
HLA-A\*33:03 for **54.8%** of them, and the locus split is A 57.5 / B 5.1 /
C 37.4. Restated as percentiles the same 38,927 peptides give **A 34.5 /
B 33.5 / C 32.0** — near-even, as a calibrated rank should be. That is the
sharpest form of this section's claim, and it is what `locus-background.tsv`
is for.

## 2. It does not explain the design's concentration on HLA-C

Of the 30 selected regions, 23 have their best predicted binder on an HLA-C
allele. It would be tempting to read that off the table above as a scale
artefact. **It is not.** Three checks, all in `locus-coverage.tsv`:

- **The best allele per region is unchanged by the rescaling.** Counting which
  allele holds each region's best target peptide gives the identical
  distribution on the nM scale and on the percentile scale — 1, 3, 2, 1, 14, 9
  for A\*01:01, A\*33:03, B\*08:01, B\*15:08, C\*01:02, C\*07:02 respectively,
  so A 4 / B 3 / C 23 either way, for all 30 regions.
- **HLA-C leads at every threshold tested.** 50 nM on the allele-neutral scale
  is the 0.107th percentile; regions reaching it are A\*01:01 2, A\*33:03 2,
  B\*08:01 1, B\*15:08 3, **C\*01:02 6, C\*07:02 6**. The lead is *widest* at
  the design's own 500 nM gate (24 regions against 16 for the best HLA-A or
  HLA-B allele) and *narrowest* at the allele-neutral and 2% rank gates (20
  against 19, and 24 against 21). So correcting the scale narrows HLA-C's lead
  without removing it.
- **The offset does not even track locus.** The most lenient allele at 500 nM
  is **HLA-A\*33:03**, more lenient than either HLA-C allele, yet it covers 16
  regions against HLA-C's 23 and 24. The largest part of the 6.7-fold nM spread
  is *within* HLA-A, which is why correcting it redistributes between HLA-A and
  HLA-B rather than away from HLA-C.

What the correction does change is coverage *counts* at the gate, and mostly by
redistributing within HLA-A and HLA-B: moving from the 500 nM gate to its
allele-neutral equivalent (1.209% rank) takes A\*01:01 from 11 regions to 17
and B\*08:01 from 13 to 19, while A\*33:03 falls from 16 to 12. The
max/min spread across alleles narrows from 2.2-fold to 1.7-fold. HLA-C stays on
top throughout.

### How much of the 23/30 is this selection, and how much is background?

A share means nothing without a baseline, so `locus-background.tsv` computes
three, all from tables already in the manifest. Share of groups whose best
target peptide sits on HLA-C, by percentile:

| Population | Groups | HLA-C share |
|---|---:|---:|
| Single arbitrary peptides | 38,927 | 32.0% |
| Random groups of 29 arbitrary peptides (size-matched to a region's target set) | 4,000 | 43.2% |
| The candidate regions this selection did **not** take | 162 | 51.9% |
| **The 30 selected regions** | 30 | **76.7%** |

23 of 30 is well above every baseline — against the size-matched 43.2%,
P(≥23/30) ≈ 2.5e-04 — so the selection really is HLA-C-skewed beyond chance.
But roughly two thirds of the 76.7% is already present before the selection is
applied: grouping 29 peptides together lifts HLA-C from 32% to 43%, and the
eligible candidate pool lifts it again to 52%. Only the remaining ~25
percentage points are attributable to this selection.

**So the HLA-C concentration is not a unit artefact — no change of scale
removes it — but it is mostly inherited from the candidate pool and from how
this predictor behaves on peptide groups of this size, rather than created by
the selection. Either way the design's class-I presentation leans on the locus
with the lowest cell-surface density, which is a substantive risk the
arithmetic cannot explain away.**

## 3. What the selected set depends on

**Order is broadly stable, but only on one chosen neutral scale.** Re-scoring
the **23 regions that carry a recorded `original_weighted_score`** with the
per-allele offset removed — same peptides, windows, flanks, logistic, gate and
expression weights — leaves F5, KLF5 and PRAME first, second and third, moves
9 of 23 regions not at all, and moves nothing more than six places (SF1 4th →
10th). But no scale is uniquely "neutral": `reference-scale-sensitivity.tsv`
repeats the exercise against each single allele's nM scale instead of the
geometric mean. Against the geometric mean those scales move 8 to 12 of the 23
regions by up to **8** places
(`reference_scale_max_rank_change_vs_geometric_mean`); against the **recorded**
order they reach **10** (`reference_scale_max_rank_change_vs_recorded`). The
stability above is a property of the geometric-mean scale, not of the design.

The other 7 of the 30 regions carry no recorded weighted score — the five
reference-derived Guardant additions, H2BC12 (expression deliberately left
missing) and the CDH1::NRG1 in-frame product. They are reported with
contribution sums only, as the design left them, rather than being silently
ranked alongside the rest.

**Structural dependence is gate-specific**, so all three gates are reported
rather than one preferred:

| Gate | Regions covered | Worst single-allele loss | Haplotype-loss range | Regions resting on one allele |
|---|---:|---:|---:|---|
| 500 nM (the design's) | 29/30 | 27 | 24–28 | HSPG2, PLEC, SYT12, TGFBR2 |
| 2% rank | 30/30 | 29 | 27–28 | CDH1::NRG1, PLEC, SYT12, TGFBR2 |
| 0.5% rank | 23/30 | 19 | 15–22 | DHX36, HSF1, HSPG2, KMT2D, NBPF14, TGFBR2 |

Only **TGFBR2 p.N515fs** rests on a single allele under every gate, on
HLA-A\*33:03. That does not make the other regions safer, and the comparison
has to be read with the zero-coverage column in mind: a region can leave the
single-allele set either by gaining a second allele **or by failing the gate
entirely**. PLEC and SYT12 leave it between the 2% and 0.5% gates the second
way. Regions inside no allele's gate are CDH1::NRG1 at 500 nM; none at 2%; and
ARHGAP40, CDH1::NRG1, H2BC12, PLEC, RAF1, SYT12 and VPS11 at 0.5% — all listed
in `summary.json` under `regions_inside_no_allele_gate`.

Crossing a gate removes a region from a count, not its peptides.
`single-allele-dependence.tsv` carries the next-best affinity and the next-best
percentile for each single-allele region, each labelled with the allele it
belongs to, because in 5 of the 9 distinct regions those are *different*
alleles. HSPG2's next-best by affinity is HLA-A\*33:03 at 547.0 nM, the 2.46th
percentile; its next-best by percentile is HLA-B\*08:01 at 0.944%, which is
754.6 nM. Pairing the two columns as one allele would misstate that percentile
2.6-fold — the same confusion between the scales that this audit is about.

Haplotype loss (HLA LOH removes a linked A/B/C triple) was enumerated over all
four possible phasings of the recorded alleles, eight candidate haplotypes.
Phase was never measured and no LOH assay was run, so these are hypotheses
about what a loss would cost.

**Two regions cross the 500 nM gate when the offset is removed.** PLEC
p.R2171W goes 413 → 789 nM, and CDH1::NRG1 — recorded as
`admission_exception = DRIVER_JUNCTION_ABOVE_500NM` at ~509 nM — goes
509 → 246 nM and would clear the gate outright. Both flips depend on holding
the threshold at 500 nM on the neutral scale, which is the 1.209th percentile
and neither conventional rank threshold: **under a 2% rank gate nothing flips
and all 30 regions are covered.** PLEC remains the weakest region on every
measure tested (last or second-last in all three rankings, no target peptide
at ≤0.5% rank on any allele), but the gate crossing itself is not robust.

**HLA-C surface density is the open risk.** HLA-C protein runs 4–5-fold below
HLA-A/B in total cellular protein and 13–18-fold at the surface in normal
lymphocytes ([Apps et al. 2015](https://doi.org/10.4049/jimmunol.1403234)).
Since §2 shows the HLA-C concentration is real, `hla-c-sensitivity.tsv` scans
weights of 1.00, 0.25, 0.20 and 0.06 on top of the allele-neutral scale.
Ranks move up to 9 places; NBPF14 falls furthest and RXRA p.A241S rises most.
This is a sensitivity parameter, not a measurement: allele-specific surface
density was not measured for Yash, and the bundled vendor expression is
gene-level RNA, neither allele-resolved nor a surface-protein measurement.

**The splice gene set is unchanged** — F5, HLA-B, HSF1, NBPF14 and KMT2D top
the 167 eligible windows under all three scales, with only the order differing.
This is weaker evidence than it looks: there are only 10 candidate genes, the
fifth-to-sixth score ratio is 2.24, and the three scales are monotone
reparametrisations of the same affinities, so agreement is close to expected.

## Run it

```sh
python designs/allele-calibration/run.py            # recompute and check audit/
python designs/allele-calibration/run.py --write    # refresh audit/
```

Python 3.11 or 3.12 (the repository's supported versions; also byte-identical
here on 3.13 and 3.15), standard library plus this repository's `src/`. No
install, no download, no model, no network, about two seconds. With no flag it
recomputes every table and fails if any committed file differs, including line
endings.

`calibration.py` is a standalone helper any other design can import to get an
allele-comparable scale without installing a predictor:

```python
from calibration import Calibration
calibration = Calibration('.')                       # reads recorded tables only
calibration.percentile('HLA-C*01:02', 500.0)         # -> 1.849
calibration.affinity_at('HLA-A*01:01', 2.0)          # -> 2609.2 nM
calibration.neutral_affinity(0.5)                    # -> 160.9 nM, allele-neutral
```

### What the script checks before it re-scores anything

- **The objective constants** are read from `score-config.json` and asserted to
  be midpoint 350, width 150, floor 1e-5, cutoff 5000 in affinity mode.
- **The per-peptide score** is checked against every one of the 4,692 rows of
  `affinity-contributions.tsv`: max absolute difference **1.11e-16**. The 1e-5
  floor also accounts for the break between the largest scored affinity
  (2,075.52 nM) and the smallest unscored one (2,124.00 nM).
- **The per-region weights** implied by the recorded scores are checked against
  the recorded expression formula for 179 of the 185 regions carrying a
  recorded score: max relative difference **4.95e-07**. That residual is *not*
  explained by rounding of the recorded proxies — rows sharing one proxy carry
  different residuals, and the proxy's own precision bounds the effect roughly
  eightfold smaller — so it is reported as consistent with single-precision
  arithmetic in the recorded pipeline, with the formula itself confirmed at the
  per-peptide level above.
- The remaining 6 regions have a recorded score but a blank
  `supplied_gene_value`; their weight is recovered from the recorded score and
  contribution sum, and the implied proxies are listed in `summary.json`. They
  are *not* joined to another design's candidate pool, where one of them
  (SPATA31A6) is filed under a different gene symbol.
- It aborts unless the plan holds exactly 30 regions, and unless the six
  audited alleles match `data/hla/class-I.tsv`. It also calls this
  repository's `require_assertions`, so it refuses to run under
  `-O`/`PYTHONOPTIMIZE` — though that is a repository convention rather than
  protection for anything here, since all seven checks above raise rather than
  assert and `-O` would not remove them.

## Results

| File | Contents |
|---|---|
| [summary.json](audit/summary.json) | Headline numbers, the checks above, scope and limitations. |
| [allele-calibration.tsv](audit/allele-calibration.tsv) | Per allele and peptide length: percentile at 50/500 nM and nM at the 0.5%/2% ranks, with allele-level columns alongside. |
| [locus-coverage.tsv](audit/locus-coverage.tsv) | Per allele: coverage at the nM gates, at their allele-neutral rank equivalents and at both conventional rank gates, plus how many regions have their best affinity and best percentile there. |
| [locus-background.tsv](audit/locus-background.tsv) | The locus split for single arbitrary peptides, size-matched random groups of them, the non-selected candidate regions and the selected 30 — the baseline §2 is measured against. |
| [coverage-matrix.tsv](audit/coverage-matrix.tsv) | Each selected region × each allele: best target affinity, percentile, allele-neutral equivalent, gate and rank verdicts. |
| [escape-robustness.tsv](audit/escape-robustness.tsv) | Six single-allele and eight candidate haplotype losses, at each of the three gates, with the regions each would drop. |
| [single-allele-dependence.tsv](audit/single-allele-dependence.tsv) | Regions inside exactly one allele's gate, per gate, with next-best affinity and percentile. |
| [region-rescoring.tsv](audit/region-rescoring.tsv) | All 30 regions under the recorded, allele-neutral and HLA-C-down-weighted scores, flagged by whether they are in the recorded-score population. |
| [reference-scale-sensitivity.tsv](audit/reference-scale-sensitivity.tsv) | How far the order moves when each single allele's scale is used as the neutral one instead of the geometric mean. |
| [splice-window-rescoring.tsv](audit/splice-window-rescoring.tsv) | All 167 eligible splice windows re-ranked under the three scales, with per-gene order. |
| [hla-c-sensitivity.tsv](audit/hla-c-sensitivity.tsv) | Region ranks across HLA-C weights 1.00/0.25/0.20/0.06. |
| [inputs.manifest.json](audit/inputs.manifest.json) | SHA-256 and size of all 19 recorded inputs read. |

[methods.md](methods.md) gives the derivation, the sign conventions and the
reasoning, including why the allele-neutral scale was defined as it was.

## Scope and limits

Predicted binding only. No processing, presentation, immunogenicity or T-cell
measurement; no experiment; no peptide validated. MHCflurry percentiles are
internal to the pinned model bundle — they normalise its nM scale across
alleles, they are not measured ligand yields, and the 0.5% and 2% thresholds
are NetMHCpan conventions applied to a different predictor whose reference
distribution is its own. Rank normalisation assumes each allele presents a
comparable share of the peptide universe, which is itself an approximation, so
part of the measured offset may be real repertoire breadth rather than a scale
artefact. Target peptides, 25-aa windows, flanks and expression proxies are the
recorded ones; only the per-allele nM scale and the HLA-C weight vary.

**What this audit cannot do.** It cannot propose replacement antigens. The
annotated candidate pool it re-scores — the 192 rows of
`region-annotations.tsv` — has already passed the 500 nM gate, so the *variants*
an allele-neutral gate would have admitted are not among them.

One bundled table is genuinely pre-gate:
`designs/design-3/inputs/predictions/all-windows.tsv` holds 15 of 125 rows above
500 nM, up to 21,470 nM. It does not help here, because those 125 rows are 25
alternative 25-aa window placements for each of the same 5 Guardant variants,
all of which were selected anyway — so it offers different placements, not new
antigens. It does contain a worked example of the effect, though: two RAF1
placements recorded at 523.4 nM on HLA-A\*01:01 sit at the 0.726th percentile,
which this audit's neutral scale restates as 252.5 nM — inside both the
allele-neutral 500 nM gate and the 2% rank gate, while outside the recorded one.

Finding genuinely new candidates means re-running
`tools/candidate_annotations/run.py` in its documented fresh mode over the
pre-gate pool, which needs the recorded Modal inputs and `--execute`.

## Suggested next steps

1. **Treat the HLA-C concentration as real.** §2 is the load-bearing result:
   23 of 30 regions depend on the locus with the lowest surface density, no
   change of scale alters that, and it is above every baseline measured here —
   while being mostly inherited from the candidate pool rather than created by
   the selection, which is where a fix would have to act. Measuring
   HLA-A/B/C allele-specific expression
   and HLA LOH from the published tumour RNA and matched-normal DNA BAMs would
   replace the weight scan here with a measurement, and is the single most
   useful follow-up.
2. **TGFBR2 p.N515fs** is the only region resting on one allele under every
   gate tested. If HLA-A\*33:03 were lost or down-regulated it would contribute
   no predicted class-I binder at any of the three gates here; its next-best
   allele is HLA-B\*08:01 at 3,937 nM / 2.67% rank, which the recorded
   objective scores as exactly zero (the 1e-5 floor bites at 2,091 nM).
3. **PLEC p.R2171W** is the weakest region on every measure tested and worth
   reconsidering in a future selection — but on its ranking and its lack of any
   ≤0.5%-rank binder, not on the 500 nM gate crossing, which is not robust.
4. For future designs, report percentile alongside nM. The %Rank convention
   exists because it normalises across MHC molecules ([Reynisson et al.
   2020](https://doi.org/10.1093/nar/gkaa379); [Nielsen & Andreatta
   2016](https://doi.org/10.1186/s13073-016-0288-x), both for NetMHCpan rather
   than for MHCflurry), and §3 shows how much the single-allele-dependency set
   changes with the threshold.
