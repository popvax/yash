# Presentation re-ranking of the Design #2 candidate peptides

An analysis of the regions Design #2 already selected and the pool it selected them from.
It adds two reported columns, processing and presentation, and measures how far the
ranking moves. It changes no selection, applies no threshold, and proposes no construct.

## Design and rationale

### The question

Design #2 ranks class I candidates on predicted binding affinity and admits them at a
500 nM gate. Affinity answers one question: given that this peptide reaches the groove,
how tightly does it sit there. It does not answer whether the peptide reaches the groove
at all. Getting there means surviving the proteasome, being transported, and being
trimmed to the right ends, and none of that is in an affinity number.

MHCflurry already models the second half. `Class1PresentationPredictor` combines the
affinity prediction with an antigen processing score through a two-feature logistic
regression, and the coefficients for that regression ship inside the
`models_class1_presentation` 2.3.0 bundle this repository already pins. The processing
models are 22 of the 44 files in
`designs/initial/inputs/checkpoints/provenance/construct-model-files.json`, and the
coefficients are in `weights.csv` of that same bundle. So the ranking Design #2 could
have used is sitting unused inside the download it already made.

The maintainers' own README asks for immunogenicity prediction. This is the smallest
honest step in that direction, because it needs no new model, no new pin and no new hash.

### What a processing score is worth, in numbers

The shipped coefficients, read out of `weights.csv` and recorded in
`presentation-run.json`:

```
with_flanks:     presentation = sigmoid(-7.2029 + 11.4267 * affinity_score + 4.0259 * processing_score)
without_flanks:  presentation = sigmoid(-6.9760 + 11.6033 * affinity_score + 3.3885 * processing_score)
```

`affinity_score` is `from_ic50(nM) = clip(1 - log(nM) / log(50000), 0, 1)`, so 20 nM is
0.723 and 500 nM is 0.426.

Going from the worst processing score to the best is worth 4.026 logits. Dividing by the
affinity coefficient, that is 0.3523 of affinity score, and inverting `from_ic50` that is
a **45.2-fold** change in nM. A peptide the processing model likes is treated as though it
bound 45 times tighter than it does.

Worked, with the numbers this run actually produced:

| Peptide | Gene | Affinity nM | Processing | Presentation |
|---|---|---:|---:|---:|
| `LSPPAPPAPL` | SPATA31E1 | 16.2 | 0.045 | 0.812 |
| `GKIHSIGDAF` | LRRC36 | 482.1 | 0.952 | 0.822 |

The first binds 30 times tighter than the second and is ranked slightly below it, because
the second looks far more likely to be produced. By affinity alone, `LSPPAPPAPL` is the
24th best admitted peptide in the pool and `GKIHSIGDAF` is the 2,326th of 2,363. By
presentation they are 706th and 674th. A 2,302-place gap becomes a 32-place one.

### Flanks

The processing model reads the residues on either side of the peptide, because the ends
are where cleavage happens. MHCflurry takes 5 residues on each side, not 15, and the two
coefficient rows exist because a prediction made without flanks is a different prediction,
not a worse one. The N and C flanks are both-or-neither: the predictor takes
`n_flanks` and `c_flanks` together or not at all.

Two sources of flanking context are in the repository, and the script prefers the longer:

- The isovar reconstructed protein contexts in
  `designs/design-2/inputs/rna/primary-proteins.tsv`, which extend past the 25mer window.
  432 peptides got their flanks from there.
- The 25mer window itself, for the 3,534 peptides with no longer context recorded.

3,344 of 3,966 peptides got both flanks and were scored on the `with_flanks` row, 622 got
one or neither and were scored on `without_flanks`. Only 1,513 peptides, 38.15%, have the
full 5 on both sides; the rest sit near a window edge and get what the window holds. That
is a real limitation of scoring from recorded windows rather than from full transcripts,
and it is why the run also scores every flanked peptide a second time with no flanks at
all and reports the difference. Median absolute change in presentation score 0.0077,
maximum 0.5086; median absolute change in processing score 0.1110, maximum 0.6148. The
flanks move individual peptides and barely move the distribution.

### Method

1. Read every recorded candidate from `designs/design-2/inputs/integrated/candidate-pool.tsv`:
   211 windows across CTA/SELF, SPLICE, MUTATION and FUSION, keyed on `region_id`.
2. Take the peptides each row records as scored, so the set matches what Design #2 ranked
   rather than a fresh enumeration.
3. Attach 5+5 flanks from the longest available context.
4. Predict affinity, processing and presentation for all 6 of the patient's class I
   alleles with `Class1PresentationPredictor`, in two groups by coefficient row.
5. Compare the affinity ranking with the presentation ranking, inside the 500 nM gate.

3,966 peptides, 24,522 peptide and allele rows, 36.9 seconds on CPU.

One window is out of scope. `REG_1f7babf47bbafeebeedc`, the CDH1::NRG1 fusion labelled
AG24, records no target span and no target peptides in its pool row, so the peptides that
cover its junction cannot be identified from that row. Enumerating the whole window
instead returns a 60.4 nM self peptide against the recorded 508.95 nM, which is a
different quantity, so the window is skipped and recorded as skipped in
`presentation-run.json` rather than scored wrongly.

### Results

**The two rankings agree on the big picture and disagree on individual peptides.**
Spearman correlation between affinity and presentation across the 2,363 admitted peptides
is **0.699**. Of the top 50 by affinity, **16** are in the top 50 by presentation.

**Processing is a real filter on a third of what passed.** **317** of the 2,363 admitted
peptides have a processing score below 0.1. These are peptides the gate let through on
affinity that the processing model thinks are unlikely to be produced at all.

**The gate itself holds up.** Median presentation score among admitted peptides is
**0.666**. Of the 1,603 peptides rejected at 500 nM, only **16** score above that median.
This is a null result and it is the most useful one here: presentation does not rescue the
rejected tail, so the 500 nM gate is not discarding a population of well-processed
peptides. The disagreement is all inside the gate, about ordering, not at the boundary.

**Region ordering moves a lot.** Across the 203 regions with scored peptides, region-level
Spearman is **0.629**. **174** regions move more than 5 places and **132** move more than
20.

The 23 selected regions in scope, ordered by the affinity Design #2 ranked them on:

| Label | Gene | Source | Best affinity nM | Affinity rank | Best presentation | Presentation rank | Move |
|---|---|---|---:|---:|---:|---:|---:|
| AG18 | PRAME | CTA/SELF | 15.3 | 16 | 0.963 | 73 | −57 |
| AG20 | LRRC36 | CTA/SELF | 21.0 | 42 | 0.985 | 16 | +26 |
| AG22 | TEKT5 | CTA/SELF | 23.7 | 60 | 0.972 | 52 | +8 |
| AG15 | MUC6 | SPLICE | 24.2 | 68 | 0.807 | 165 | −97 |
| AG13 | F5 | SPLICE | 24.2 | 69 | 0.968 | 63 | +6 |
| AG03 | AQP3 | MUTATION | 29.4 | 94 | 0.972 | 55 | +39 |
| AG21 | YBX2 | CTA/SELF | 29.8 | 96 | 0.980 | 33 | +63 |
| AG19 | SPATA31A6 | CTA/SELF | 33.6 | 114 | 0.897 | 134 | −20 |
| AG07 | PPIC | MUTATION | 36.5 | 121 | 0.859 | 151 | −30 |
| AG08 | KCTD19 | MUTATION | 36.9 | 122 | 0.904 | 132 | −10 |
| AG01 | KLF5 | MUTATION | 37.0 | 124 | 0.958 | 84 | +40 |
| AG06 | TGFBR2 | MUTATION | 38.0 | 127 | 0.769 | 177 | −50 |
| AG02 | SF1 | MUTATION | 41.5 | 136 | 0.872 | 148 | −12 |
| AG14 | MUC5B | SPLICE | 44.9 | 150 | 0.780 | 173 | −23 |
| AG04 | RXRA | MUTATION | 58.1 | 173 | 0.876 | 145 | +28 |
| AG17 | HLA-B | SPLICE | 66.2 | 177 | 0.969 | 60 | +117 |
| AG23 | CDH1--NRG1 | FUSION | 70.2 | 178 | 0.915 | 124 | +54 |
| AG05 | HSPG2 | MUTATION | 81.1 | 182 | 0.468 | 193 | −11 |
| AG09 | VPS11 | MUTATION | 115.6 | 187 | 0.842 | 154 | +33 |
| AG10 | ARHGAP40 | MUTATION | 138.2 | 191 | 0.433 | 195 | −4 |
| AG12 | SYT12 | MUTATION | 243.7 | 196 | 0.752 | 178 | +18 |
| AG16 | MUC3A | SPLICE | 260.8 | 197 | 0.330 | 197 | 0 |
| AG11 | PLEC | MUTATION | 413.3 | 201 | 0.274 | 198 | +3 |

A positive move means the region ranks better by presentation than by affinity. AG24 is
absent for the reason given above.

Reading that table from both ends. The three weakest selections by affinity, AG11 PLEC at
413.3 nM, AG16 MUC3A at 260.8 nM and AG10 ARHGAP40 at 138.2 nM, are also the three
lowest by presentation, 0.274, 0.330 and 0.433. Presentation does not argue for keeping
them. On the other side, AG17 is the largest single move in the table at 117 places, and
it is the HLA-B splice region, which is already the one selection in the cassette worth
a separate question for other reasons. AG18 PRAME drops 57 places while still scoring
0.963, which is the honest reading of a drop near the top: the regions above it are
crowded, and a drop there means less than a drop in the middle.

The two directions have a single cause, and it is visible in the extremes. The five
peptides that fall furthest all have processing scores between 0.011 and 0.020. The five
that rise furthest all have processing scores between 0.844 and 0.952. Nothing subtler
than that is going on.

## Sequences and results

No sequence is proposed, assembled or changed. No cassette, no construct, no new antigen.

Three output files, all in this directory:

- [`presentation.tsv`](presentation.tsv), 24,522 rows, one per peptide and allele:
  `region_id`, `source`, `gene`, `peptide`, `length`, `allele`, `offset_in_window`,
  `flank_context`, `n_flank`, `c_flank`, `full_5_5_flanks`, `coefficient_row`,
  `affinity_nm`, `affinity_percentile`, `processing_score`, `presentation_score`,
  `presentation_percentile`, and the two no-flank sensitivity columns.
- [`region-rank-comparison.tsv`](region-rank-comparison.tsv), 203 rows, one per region:
  best affinity and best presentation with the peptide and allele that produced each,
  the two ranks and the change.
- [`presentation-run.json`](presentation-run.json), the run record: counts, the shipped
  coefficients as read from the bundle, the flank breakdown, the sensitivity pass, the
  ranking statistics, the skipped window with its reason, the parity comparison, and
  library versions.

Keys are `region_id`. AG labels appear in this README only, to read against the recorded
cassette.

## Reproduction instructions

### Inputs, all already in the repository

- `designs/design-2/inputs/integrated/candidate-pool.tsv`, the candidate windows and their
  recorded peptides.
- `designs/design-2/inputs/rna/primary-proteins.tsv`, the isovar protein contexts used for
  flanks.
- `config/case.json`, the patient's six class I alleles: HLA-A\*01:01, HLA-A\*33:03,
  HLA-B\*08:01, HLA-B\*15:08, HLA-C\*01:02, HLA-C\*07:02.
- `designs/initial/inputs/checkpoints/provenance/construct-model-files.json`, the pinned
  model manifest, used to check the bundle.
- `designs/design-2/inputs/predictions/mutations/best-peptide-per-HLA.tsv` and
  `designs/design-2/inputs/predictions/fusions/affinities.tsv`, as parity checks only.

Nothing in `designs/design-2/` is read for anything other than input and comparison, and
nothing there is written.

### The one external download, which is the bundle already pinned

```sh
pip install 'mhcflurry==2.3.14'
mhcflurry downloads fetch models_class1_presentation --release 2.3.0
MODELS="$(mhcflurry downloads path models_class1_presentation --release 2.3.0)/models"
python designs/presentation-reranking/run.py --models "$MODELS"
```

Note the path ends at `models`, not `models/affinity_predictor`: the presentation
predictor needs the processing models and `weights.csv` beside the affinity models.

Direct URL from `mhcflurry/downloads.yml` at tag 2.3.14:
`https://github.com/openvax/mhcflurry/releases/download/2.3.0/models_class1_presentation.20260928.tar.bz2`,
160,356,140 bytes compressed.

The script hashes all 44 files of the extracted bundle against the pinned manifest. They
are identical, in order, byte for byte, and the run record states that as
`models_match_pinned_bundle: true`. The same weights that produced the recorded
predictions produce these, and nothing new is pinned.

### Dependencies

What I ran: mhcflurry 2.3.14, numpy 2.5.4, torch 2.15.0 on CPU, Python 3.12.3. Those are
not the recorded versions, which `construct-inference.json` gives as torch 2.14.1+cu130
and numpy 2.4.6 on a GPU. The parity check below exists because of that gap. Spearman is
computed with numpy only, no scipy.

### Parity against the recorded predictions

Every number here is new inference. Nothing is a relabelled cached prediction, and
`presentation-run.json` records both as `new_inference: true` and
`cached_predictions_relabelled: false`.

To show the affinity half is not drifting, the run compares its own affinities against the
repository's recorded tables wherever the peptide, allele and region coincide:

| Recorded table | Rows compared | Rows outside this scan | Max relative difference |
|---|---:|---:|---:|
| `mutations/best-peptide-per-HLA.tsv` | 72 | 24 | 1.2e-06 |
| `fusions/affinities.tsv` | 72 | 132 | 1.1e-06 |

Six significant figures of agreement on a different accelerator and a different numpy.
The processing and presentation columns have no recorded counterpart to compare against,
because this repository never computed them.

## Interpretation

### What the evidence covers

A local replay of nothing, and new inference on everything. The affinity column is a
reproduction of a quantity this repository already has, agreeing to 1.2e-06. The
processing and presentation columns are new quantities computed from already-pinned
weights. No experimental result of any kind is involved.

Within that, three claims are supported. First, presentation reorders the admitted pool
substantially but not randomly, rho 0.699 at peptide level and 0.629 at region level.
Second, 317 admitted peptides carry a processing score below 0.1, so affinity alone
admitted 13% of the pool on a criterion the processing model disagrees with. Third,
presentation does not rescue the rejected tail: 16 of 1,603 rejected peptides beat the
admitted median, which supports keeping the 500 nM gate exactly where it is.

### What it does not cover

No threshold is applied and no selection is changed. Nothing here says the cassette should
be different. A region's presentation rank is not a statement that it is a better or worse
vaccine target, because presentation predicts presentation, and immunogenicity needs more
than that, including the agretopicity-style question of whether the repertoire was
tolerized against something nearly identical.

Class I only. MHCflurry has no class II model, so the class II gap in this repository is
untouched by this work.

The flank limitation is the main one to weigh. 38.15% of peptides have full 5+5 context;
the rest are scored on what the recorded window holds. The sensitivity pass bounds the
cost at a median 0.0077 of presentation score, but the maximum is 0.5086, so individual
peptides near a window edge can be meaningfully mis-scored in either direction. Full
transcript context would fix this and is not reconstructable from the recorded windows
alone.

MHCflurry ships no percentile calibration for the processing score, so processing is
reported as a raw score and is not comparable across alleles in the way an affinity
percentile is. The presentation percentile that MHCflurry does provide is carried in the
table unchanged.

One window, AG24, is out of scope, so the fusion source is represented by AG23 alone.

### Open questions I cannot answer from this evidence

Whether the 317 low-processing admitted peptides would actually fail to be presented. The
processing model is trained on mass spectrometry ligands, so a low score means
under-represented in that data, which is evidence and not proof of absence.

Whether re-ranking would have changed the cassette. The optimizer in
`designs/initial/construct/optimize.py` takes affinity contributions as part of a larger
objective, and substituting presentation there is a different piece of work with a
different review burden. This analysis deliberately stops before touching it.

### Unfinished work

Nothing in this directory is incomplete, but two extensions are obvious and were left
out on purpose: a presentation-weighted objective for the optimizer, and the full
Łuksza and Balachandran quality model, of which neither presentation nor the
agretopicity analysis in `designs/agretopicity/` is more than a component.
