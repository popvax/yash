# ND-only normal-reference screening and limited strand checks

This analysis implements additional check 1 and limited check 4 on the existing
patient-RNA splice pool. It introduces no distinct-pattern requirement or extra
substitution filter. The initial design uses this ND-only selection. These files record
sequence eligibility; assembly and RNA encoding live under designs/initial/.

## Method

Search exact GRCh38 junction endpoints across every resolvable UID in the two
supplied normal-control matrices, regardless of gene assignment and strand.
Every eligible splice event must now be **NOT_DETECTED_IN_SEARCHED_REFERENCE
(ND)**: no exact coordinate entry in either searched catalog. Any catalog match,
including an all-zero entry or an opposite-strand assignment, is REVIEW.
Unresolved coordinates cannot qualify as ND. Missing means stay null rather than
being imputed as measured zeros.

The [selection policy](selection-policy.json) records this offline admission rule
separately from the original [measurement policy](checks-policy.json). Measurement
files, executed worker snapshots and receipts retain their original hashes. The
original per-tissue mean threshold of 3 remains in those measurement records;
**a low mean no longer admits a detected junction**. Counts, sample denominators
and tissue means remain available as diagnostics. Same-coordinate same-strand
aliases use the per-sample maximum, avoiding duplicate counting.

Measure library orientation from paired reads at annotated junctions, independently
of motif/XS. For each full-window witness set, require a strict majority agreeing
with the nominated gene's direction. Tied, opposite, conflicting or unevaluable
evidence is REVIEW. Motif disagreement is report-only. The
[policy](checks-policy.json) and [stage guide](../../../workflow/stages/04_region_integration/patient_splice/CHECKS.md)
state every operational threshold and limitation.

All existing sequences, frame anchors, support requirements, native review,
reference exclusions, affinity predictions and expression scores are reused
unchanged. PASS here means passing these specific added requirements; it does
not establish tumour specificity, translation or presentation.

## Results

The search scanned 2,476,734 GTEx catalog rows and 5,403,219 TCGA-adjacent-normal
rows. Coordinates resolved for 2,463,061 and 5,379,440 rows, respectively; 13,673
and 23,779 rows could not be resolved. Counts cover 2,629 GTEx samples in 51
source tissue categories and 705 TCGA-adjacent-normal samples in 21 categories.
The search used these processed references, not raw normal RNA alignments.

All 4,929 input events received normal-reference decisions. All 39,269
reconstructed windows received limited strand decisions. Among the 3,804 windows
previously eligible for primary scoring:

| Added decision | Windows held for review |
|---|---:|
| Normal check alone | 2,729 |
| Strand check alone | 159 |
| Both checks | 164 |
| Remain eligible | 752 |

The 752 remaining windows cover **42 events across 22 genes**. All 22 genes have
a positive supplied expression value. Ranking uses the original scores and best
passing window per event, then per gene. Across the complete 4,929-event input
pool, 1,527 are ND; other existing reconstruction, native QC, strand, reference
and affinity requirements still apply before an event can enter the 42-event pool.

| Rank | Gene | Weighted score | Best target affinity, nM | Normal-reference finding |
|---|---|---:|---:|---|
| 1 | F5 | 272.94 | 24.24 | Not detected in searched reference |
| 2 | MUC5B | 86.20 | 44.87 | Not detected in searched reference |
| 3 | MUC6 | 68.73 | 24.16 | Not detected in searched reference |
| 4 | MUC3A | 45.34 | 260.85 | Not detected in searched reference |
| 5 | HLA-B | 44.19 | 66.16 | Not detected in searched reference |

These are the five splice selections used by the initial design.
SETD2 is held for review because its coordinates are detected in normal
references: 745/2,629 GTEx samples and 1/705 TCGA-adjacent-normal samples. Its
highest tissue mean of 0.85 does not permit admission under the ND rule. HLA-B
is fifth in the five-gene selection; its complete window has two agreeing
fragment names and no opposing names. Those names do not establish independent
original molecules.

ND describes junction-catalog membership, not demonstrated absence of a peptide
from healthy tissues. A peptide could be produced by another normal junction,
an unrepresented transcript or an inherited alternative sequence. This policy
adds no normal-tissue immunopeptidome comparison or experimental presentation test.

The excluded CFH event is REVIEW because of normal-reference detection, including
137/137 GTEx liver samples. Other MUC3A and MUC5B events have opposite-strand
normal-coordinate matches and remain REVIEW. The final selected events for those genes are ND.

Library read-1 orientation is antisense in 1,252/1,275 sampled fragments (98.20%)
across 200 genes. The limited strand rule is evaluated separately from normal-reference admission. MUC6
has two agreeing and one opposing full-window witness; that minority is reported.
All five selected windows have canonical motifs pointing opposite the
nominated gene. All five pass the limited strand check: motif disagreement is
not an additional exclusion rule. SETD2's noncanonical motif remains recorded
but is not the reason it fails the new ND criterion.

See [selected rows](selected.tsv), [gene ranking](gene-ranking.tsv),
[previous-window decisions](previous-selection-checks.tsv),
[every event](event-checks.tsv), [summary](summary.json) and
[verification](verification.json). The normal positive-sample observation totals
are source-sample observations, not independent-person or molecule counts.

## Reproduce the completed analysis offline

From the repository root, in the installed project environment:

```sh
python workflow/stages/04_region_integration/patient_splice/checks.py apply \
  --source data/splice-evidence --output data/splice-evidence/normal-screen
python workflow/stages/04_region_integration/patient_splice/checks.py verify \
  --source data/splice-evidence --output data/splice-evidence/normal-screen
```

Apply first installs the current selection-policy.json (or checks an identical
existing copy), then verifies the complete original splice reconstruction,
reference and prediction evidence. It recomputes normal/strand gates and the
selection using included measurements and unchanged affinity/expression scores.
Verify uses the recorded selection policy, rebuilding and comparing all ten
output artifacts byte for byte. No Modal login or new inference is needed.

The [independent check](independent-check.json) verifies the ND decision directly
from all recorded coordinate-match lists, the unchanged tissue-count artifact,
each window's strand eligibility, unchanged source sequence/score fields and
the event/gene rankings. Raw measurement files and receipts are unchanged by admission replay.

## Repeat the BAM/reference measurements

Use the authenticated existing Modal account and research-lens-workspace volume.
Read the [stage guide](../../../workflow/stages/04_region_integration/patient_splice/CHECKS.md).
The completed run used remote source
/data/results/TN25-112103/patient-rna-splices-20261103-final and result root
/data/results/TN25-112103/splice-normal-strand-20261104. Use a new result root
and local output directory for a new attempt.

Both jobs were CPU-only. The [normal receipt](normal-receipt.json) and
[strand receipt](strand-receipt.json) record input hashes and software versions.
normal-state.json and strand-state.json retain the actual sandbox handles.
execution-sources/ contains the controller/worker text whose hashes match the
executed jobs. The checked source BAM, annotation, genome and coordinate adapter
match those used for the original reconstructed evidence.

## Scope retained

The starting pool is the 4,929 events already passing the original pooled normal
screens. This analysis tightens admission within that pool; it does not recover
potential candidates rejected upstream. It adds no checks 2/3, DNA filter,
splice-fraction weighting, sequence substitution, new HLA model or construct
optimization. Detected and uncertain cases remain available for review rather
than being deleted from the evidence.
