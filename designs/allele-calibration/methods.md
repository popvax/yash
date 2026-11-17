# Methods

Everything below runs from files already committed to this repository. No
model, no download, no inference, no cloud job. `audit/inputs.manifest.json`
records the SHA-256 and size of all 19 inputs read.

## 1. The recorded objective, re-implemented and checked

`designs/initial/inputs/mutations/rank/score-config.json` records the scoring
rule as vaxrank's

```
EpitopeConfig(logistic_epitope_score_midpoint=350.0,
              logistic_epitope_score_width=150.0,
              min_epitope_score=1e-05,
              binding_affinity_cutoff=5000.0,
              scoring_mode='affinity')
resolved_score_expr = (affinity.value < 5000.0)
                      * affinity.value.logistic_normalized(350.0, 150.0)
```

`run.py` parses those four constants out of that record and aborts if they
differ from the ones it uses, so the formula below is pinned to the recorded
configuration rather than hard-coded on trust. `logistic_normalized` is
normalised at zero, so the per-peptide contribution is

```
L(a) = 1 / (1 + exp((a - 350) / 150))
score(a) = L(a) / L(0)      zeroed when score < 1e-5, or when a >= 5000 nM
```

Checked against all 4,692 rows of
`designs/initial/inputs/mutations/rank/affinity-contributions.tsv`, this
reproduces every recorded `contribution` to a maximum absolute difference of
**1.11e-16**. The `min_epitope_score` floor accounts for the observed break
between the largest scored affinity (2,075.52 nM) and the smallest unscored one
(2,124.00 nM). This check runs inside `run.py`; the file is in the manifest.

### Region scores and the expression weight

A region's score is the contribution sum over its target peptide/allele pairs
times an expression weight. `score-config.json` records the weight as
`sqrt(supplied gene value)` for splice, CTA/self and fusion regions, and
`sqrt(supplied gene value * rna_alternate_fraction)` for mutation regions, the
fraction coming from `designs/design-3/selection/selected.tsv`.

For ranking, the weight is taken as

```
weight = original_weighted_score / recorded_contribution_sum
```

so the "recorded" column of `audit/region-rescoring.tsv` is the design's own
number, not a reconstruction of it. `run.py` then checks that implied weight
against the formula above. Of the 192 annotated region rows, 185 carry a
recorded `original_weighted_score`; 179 of those also carry a
`supplied_gene_value` and are checked against the formula, with a maximum
relative difference of **4.95e-07**.

That residual is reported as consistent with single-precision arithmetic in the
recorded pipeline, not as rounding of the recorded proxies. Rounding does not
fit it: the worst row is a non-plan KMT2D splice window with no
`rna_alternate_fraction` at all, the five-decimal precision of its
`supplied_gene_value` bounds the effect on the score about eightfold smaller
than the observed residual, and rows sharing one proxy carry different
residuals. The formula itself is confirmed at the per-peptide level above, where
it is exact to 1.11e-16.

The remaining 6 regions — LRRC36, PRAME, SPATA31A6, TEKT5, YBX2 and the
CDH1--NRG1 fusion — have a recorded score but a deliberately blank
`supplied_gene_value`. Their weight comes from the division above and the
implied proxies are listed in `audit/summary.json` under
`weight_implied_from_recorded_score_only`.

They are **not** joined to the Design #1 or Design #2 candidate pools, even
though those pools do carry an `expression_proxy` matching all six implied
values to at least six significant figures. The reason is that the join is
unsafe on the key that looks natural. In both pools the SPATA31A6 region id
`REG_a391cc51fbd08e062be3` appears on **four** rows, under the gene symbols
SPATA31A1, SPATA31A3, SPATA31A6 and SPATA31A7, with four different proxies
(3.44875, 5.60661, 13.60587 and 4.86063). Only the SPATA31A6 row carries the
value this design used. A join by region id would therefore pick one of four
proxies for that region depending on iteration order — a gene-symbol join
happens to be correct for all six here, but the ambiguity is real and
`AGENTS.md` flags gene assignment as something to treat carefully. Dividing
the recorded score by the recorded contribution sum avoids the question
entirely and is exact.

Seven of the 30 plan regions carry no recorded weighted score — the five
reference-derived Guardant additions (AR, DHX36, ERBB4, KMT2D p.R2789Q, RAF1),
H2BC12 whose expression is deliberately left missing, and the CDH1::NRG1
in-frame product. They are reported with contribution sums and coverage only,
and are flagged `in_recorded_score_population = False` rather than ranked
alongside the rest.

## 2. Per-allele calibration from recorded predictions

MHCflurry returns an affinity in nM and the percentile of that affinity among
random peptides for the same allele. Thirteen committed tables carry both for
the six case alleles, all produced by the same pinned MHCflurry 2.3 bundle:

```
data/shared-rna/splicing/scoring/affinities.tsv.gz
designs/initial/construct/terminal-inference/predictions.tsv
designs/initial/construct/protein/boundary-predictions.tsv
designs/initial/construct/selected-inference/predictions.tsv
designs/initial/selection/boundary-predictions.tsv
designs/initial/inputs/mutations/rank/peptide-HLA-affinity.tsv
designs/design-2/construct/protein/boundary-predictions.tsv
designs/design-2/construct/selected-inference/predictions.tsv
designs/design-3/construct/protein/boundary-predictions.tsv
designs/design-3/construct/selected-inference/predictions.tsv
designs/design-3/inputs/predictions/peptide-HLA-affinity.tsv
designs/design-3/review/additional-mutation-affinities.tsv
designs/design-3/review/candidate-annotations/predictions.tsv.gz
```

Together they give 233,562 distinct peptide/allele percentiles, 47,016 rows and
about 44,000 distinct affinity points per allele.

Two checks establish that the pooled points define a single curve rather than a
scatter, because one alone would not. `Calibration.inversions()` counts
adjacent pairs where a higher affinity carries a lower percentile, over the raw
rows rather than the affinity-keyed curve, and reports 0 for all six alleles.
That test cannot by itself detect two rows that share an affinity but disagree
on percentile, since the rows are sorted by `(affinity, percentile)` and such a
pair is already in ascending order. `Calibration.tied_affinity_disagreements()`
is the test for that case: it groups the raw rows by affinity and counts the
affinities whose recorded percentiles are not all equal. It reports 0 for all
six alleles, and the largest relative percentile spread inside any single
affinity is 2.1e-16, which is floating-point representation noise.

The relation is close to, but not exactly, length-independent. The per-length
columns of `audit/allele-calibration.tsv` differ — most visibly HLA-B\*08:01 at
50 nM, where 8, 9, 10 and 11-mers give 0.0905, 0.0905, 0.2415 and 0.1527, and
HLA-A\*33:03 at 50 nM, 0.3917 for 8-mers against 0.2917 for the rest. This does
not affect any result here: all 8,994 target peptide/allele rows of
`peptide-annotations.tsv.gz` have a **recorded** percentile
(`target_rows_needing_interpolated_percentile` is 0), so the per-length curves
are reported for inspection and never relied on.

Interpolation, where used, is linear in log affinity and log percentile between
neighbouring recorded points. Each curve is inverted the same way to read an
affinity at a given percentile. `audit/allele-calibration.tsv` carries both
per-length and allele-level columns so the two are never confused.

## 3. The allele-neutral scale, and its arbitrariness

To separate the allele preference of the nM scale from the antigens, each
prediction is restated on one common scale while keeping its percentile:

```
neutral_nM(peptide, allele) = geometric mean over the six case alleles
                              of nM_allele(percentile(peptide, allele))
```

The geometric mean is used because affinity is log-distributed. Nothing else
changes: the peptides, their percentiles, the 25-aa windows, the flanks, the
logistic, the 500 nM gate and the expression weights are all the recorded ones.
Summed over the 30 regions the contribution mass moves from 221.57 to 194.31, a
12.3% change, so the comparison is not dominated by a global shift.

**This scale is a choice, and the resulting order is not invariant to it.**
Because the scoring logistic is non-linear, a monotone reparametrisation of nM
re-weights peptides and changes the sums.
`audit/reference-scale-sensitivity.tsv` repeats the whole re-scoring with each
single allele's nM scale used as the neutral one. Against the geometric mean,
8 to 12 of the 23 ranked regions move and the largest change is 4 to 8 places;
against the **recorded** order the largest change reaches **10** places. So the
"nothing moves more than six places" result in the README is a property of the
geometric-mean scale, and is reported as such.

## 4. Coverage, locus concentration and loss scenarios

A region counts as covered by an allele when at least one of its recorded
target peptides reaches the threshold for that allele. Three gates are carried
through every structural result:

- `recorded_500nm` — the design's own 500 nM affinity gate.
- `rank_2p0_percent` and `rank_0p5_percent` — the conventional weak and strong
  rank thresholds. These are NetMHCpan's defaults, applied here to MHCflurry's
  own `affinity_percentile`, whose reference distribution is its own.

Reporting all three matters because the answer depends on the gate: the set of
regions inside exactly one allele's gate is {HSPG2, PLEC, SYT12, TGFBR2} at
500 nM, {CDH1::NRG1, PLEC, SYT12, TGFBR2} at 2% and {DHX36, HSF1, HSPG2,
KMT2D, NBPF14, TGFBR2} at 0.5%. Only TGFBR2 appears under all three.

To ask whether the nM offset explains the concentration on HLA-C,
`audit/locus-coverage.tsv` compares the loci at matched stringency. The
allele-neutral equivalents of 50 nM and 500 nM are found by bisecting the
neutral scale, giving the 0.107th and 1.209th percentiles, and coverage is
counted at those ranks as well as at the two conventional ones. The table also
counts, per allele, how many regions have their best affinity there and how
many have their best percentile there; these two counts are identical for all
six alleles, which is the strongest single form of the negative result.

A share also needs a baseline, so `audit/locus-background.tsv` computes the
locus split for four populations: single arbitrary peptides (every peptide in
the calibration tables scored against all six alleles), random groups of those
peptides sized to the mean number of distinct target peptides per selected
region, the candidate regions this selection did not take, and the selected 30.
Sampling is seeded (`BACKGROUND_SEED = 20261117`, 4,000 draws) so the table is
reproducible byte for byte. The same table carries the per-allele share of
best-affinity calls on arbitrary peptides, which is the clearest demonstration
of §1: HLA-A\*33:03 holds the best nM prediction for 54.8% of 38,927 arbitrary
peptides, and the locus split A 57.5 / B 5.1 / C 37.4 on the nM scale becomes
A 34.5 / B 33.5 / C 32.0 on the percentile scale.

`audit/escape-robustness.tsv` enumerates, at each gate:

- **Single-allele loss**: each of the six alleles removed in turn.
- **Haplotype loss**: HLA loss of heterozygosity removes a linked HLA-A, HLA-B
  and HLA-C allele together. Phase was never measured here, so all four A/B/C
  pairings of the recorded alleles are enumerated, giving eight candidate
  haplotypes, and every one is reported. These are hypotheses about what a loss
  would cost, not evidence that any loss occurred; no LOH assay was run.

A region can also leave the single-allele set by failing a gate entirely rather
than by gaining a second allele, which a count of covered regions does not
distinguish. `summary.json` therefore lists
`regions_inside_no_allele_gate` per gate: CDH1::NRG1 at 500 nM, none at 2%, and
ARHGAP40, CDH1::NRG1, H2BC12, PLEC, RAF1, SYT12 and VPS11 at 0.5%. PLEC and
SYT12 drop out of the single-allele set between the 2% and 0.5% gates that way.

Crossing a gate removes a region from a count, not its peptides.
`audit/single-allele-dependence.tsv` therefore carries the next-best affinity
and the next-best percentile for every single-allele region. These two need not
belong to the same allele — in 5 of the 9 distinct regions they do not, for
instance HSPG2, whose next-best allele by affinity is HLA-A\*33:03 at 547.0 nM
(2.46%) while its next-best by percentile is HLA-B\*08:01 at 0.944%
(754.6 nM). Each figure is therefore paired with the allele it came from, and
`next_best_alleles_agree` flags the rows where they coincide. Reading the two
columns as one allele's pair would misstate HSPG2's next-best percentile by
2.6-fold, which is the same confusion between the two scales that this audit
is about.

## 5. HLA-C weight scan

HLA-C protein is expressed below HLA-A and HLA-B. Apps et al. measured, in
normal peripheral blood lymphocytes of one haplotype, HLA-A at 5-fold and
HLA-B at 4-fold the level of HLA-C by mass spectrometry of total cellular
protein, and HLA-A at 13–16-fold and HLA-B at 17–18-fold by flow cytometry of
surface protein. `audit/hla-c-sensitivity.tsv` multiplies every HLA-C
contribution by 1.00 (the recorded assumption), 0.25 and 0.20 (the reciprocals
of the total-cellular ratios) and 0.06 (inside the 1/18–1/13 = 0.056–0.077
surface range), on top of the allele-neutral scale, and reports each region's
rank at each weight.

This is a sensitivity parameter, not a measurement. Allele-specific surface
density was not measured for this tumour, and the published ratios come from
normal cells of a different haplotype. The bundled vendor expression is
gene-level RNA, so it is neither allele-resolved nor a measurement of surface
protein, and the audit does not claim a mechanism for the protein-level
difference — Apps et al. report protein measurements and attribute the gap
between their two methods to surface-versus-total measurement and antibody
specificity, not to transcription.

## 6. Sign conventions

Rank changes are `recorded_rank - alternative_rank`, so a positive value means
the region moves up the order. In `audit/hla-c-sensitivity.tsv` the change is
`rank at weight 1.00 - rank at weight 0.06`, both computed on the
allele-neutral scale, so those columns isolate the weight effect rather than
the calibration effect.

## 7. Verification

`python designs/allele-calibration/run.py` with no flag rebuilds every table
and JSON file from the recorded inputs and fails if any committed file differs,
including line endings. It aborts unless the plan contains exactly 30 regions,
unless the audited alleles match the recorded class-I calls in
`data/hla/class-I.tsv`, unless the objective constants match
`score-config.json`, unless the per-peptide score reproduces
`affinity-contributions.tsv`, and unless at least 185 regions carry a recorded
score whose implied weight matches the recorded formula. Runtime is about two
seconds, and output is byte-identical on Python 3.11, 3.12, 3.13 and 3.15.

It also calls `yash_vaccine.io.require_assertions`, so it refuses to run under
`-O` or `PYTHONOPTIMIZE`. That follows this repository's convention rather than
protecting anything here: every check listed above raises `ValueError`, which
`-O` does not strip, and neither module in this directory contains an `assert`.

## References

- Apps R, Meng Z, Del Prete GQ, Lifson JD, Zhou M, Carrington M. Relative
  expression levels of the HLA class-I proteins in normal and HIV-infected
  cells. *J Immunol.* 2015;194(8):3594–3600.
  [doi:10.4049/jimmunol.1403234](https://doi.org/10.4049/jimmunol.1403234)
- O'Donnell TJ, Rubinsteyn A, Laserson U. MHCflurry 2.0: improved pan-allele
  prediction of MHC class I-presented peptides by incorporating antigen
  processing. *Cell Syst.* 2020;11(1):42–48.e7.
  [doi:10.1016/j.cels.2020.06.010](https://doi.org/10.1016/j.cels.2020.06.010)
  — the predictor whose recorded outputs this audit reuses, and the source of
  the `affinity_percentile` values the calibration is built from.
- Reynisson B, Alvarez B, Paul S, Peters B, Nielsen M. NetMHCpan-4.1 and
  NetMHCIIpan-4.0: improved predictions of MHC antigen presentation by
  concurrent motif deconvolution and integration of MS MHC eluted ligand data.
  *Nucleic Acids Res.* 2020;48(W1):W449–W454.
  [doi:10.1093/nar/gkaa379](https://doi.org/10.1093/nar/gkaa379) — %Rank
  normalises prediction scores across MHC molecules; 0.5% and 2% are the
  default class-I strong and weak binder thresholds. For NetMHCpan, not for
  MHCflurry.
- Nielsen M, Andreatta M. NetMHCpan-3.0; improved prediction of binding to MHC
  class I molecules integrating information from multiple receptor and peptide
  length datasets. *Genome Med.* 2016;8:33.
  [doi:10.1186/s13073-016-0288-x](https://doi.org/10.1186/s13073-016-0288-x) —
  percentile ranks rather than affinity thresholds for ligand identification.
  Also for NetMHCpan.
