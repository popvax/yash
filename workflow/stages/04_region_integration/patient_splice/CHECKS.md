# Additional normal-reference and limited strand checks

This is an added admission stage for the existing patient-RNA splice workflow.
It reuses the reconstructed sequences, native review, exact-reference exclusions,
MHCflurry affinities and Vaxrank/expression scores. It applies only the two added
checks below before selecting the best passing window per event and gene.

This stage records normal-reference measurements and admission alongside the
source reconstruction. The initial design uses its passing selections. The
measurement and reconstruction files retain their own hashes. The complete evidence is in
[data/splice-evidence/normal-screen](../../../../data/splice-evidence/normal-screen/README.md).

## Check 1: coordinate-aware normal reference

Resolve every UID in the supplied GTEx and TCGA-adjacent-normal matrices using
the full Alt91 exon table and the original MSNAF coordinate adapter. Match the
exact GRCh38 exon-boundary endpoints across all gene assignments and both
strands. The same physical coordinates on the opposite strand are recorded as
an interpretation requiring review, not asserted to encode the same peptide.

For aliases on the same strand, use the per-sample maximum count rather than
adding alternate annotations of the same gap. Decode the source tissue labels,
including the legacy H5AD categorical representation. Retain sample denominators,
positive sample counts, mean counts and maximum counts separately by source
tissue category. TCGA categories are preserved as supplied; they are not all
translated into a harmonized organ ontology.

The current admission rule requires NOT_DETECTED_IN_SEARCHED_REFERENCE (ND)
for every selected splice event. Any exact coordinate entry in either catalog,
on either strand and under any resolved gene assignment, is REVIEW. This also
holds catalog entries with all-zero sample counts: a known entry is distinct
from the ND label. Unresolved coordinates cannot qualify as ND. Detected events
remain in the evidence and are excluded from automatic primary selection.

The original measurement policy retained MSNAF's normal-mean threshold of 3,
applied separately by tissue, plus opposite-strand review. Those measurements,
policy and receipts remain unchanged. The separately hashed
[selection policy](selection-policy.json) tightens admission using the existing
counts; a low mean no longer allows a detected junction into selection. Tissue
means remain diagnostic outputs. No new cloud measurement or HLA inference is
claimed. There is no added tumour/normal count ratio or comparison of supplied
transformed expression with GTEx counts.

Missing coordinate matches are labelled NOT_DETECTED_IN_SEARCHED_REFERENCE.
Their normal mean stays null. Catalog rows containing only zeros remain distinct
from absent rows. Search scope is limited to resolvable rows of the two supplied
processed matrices; unresolved-row totals are recorded. Raw normal reads, newer
GTEx releases and the entire upstream rejected junction pool are not reprocessed.
Normal junction detection does not establish normal peptide translation or HLA
presentation. ND is limited to the references and resolvable rows actually searched.

## Limited check 4: transcription orientation

Infer library read-1 orientation from unambiguous annotated CDS-to-CDS junctions
using paired BAM flags. Sample at most 20 fragments per selected gene, at most
200 genes and 5,000 intervals. The configurable diagnostic requires at least
100 fragments across 20 genes and a dominant orientation fraction >=90%; these
are operational confidence requirements, not a biological universal.

For every complete reconstructed window, recover the named junction/boundary
witnesses and classify their orientation under that measured library convention.
Require a strict majority supporting the nominated gene direction. Unrecovered,
conflicting or unevaluable witnesses and tied/opposite majorities are REVIEW.
A minority in the other direction is reported rather than an automatic failure.

Splice motif and XS tags are report-only. In particular, a motif pointing opposite
the gene cannot by itself remove a window. STAR's motif-derived XS tags are not
counted as an independent orientation measurement.

No distinct-alignment-pattern minimum, UMI/molecule-count claim, additional
substitution check, DNA-variant gate, new frame filter or sequence replacement is
introduced. Existing fragment-name and native-QC requirements remain intact.

## Run on a completed evidence directory

Use the repository environment and a new local and remote result directory.
The two CPU workers run independently on the existing Modal volume. Cloud calls
and downloads require --execute. Without it, measure/status/collect only print
plans. The source root must be a complete recorded patient-RNA splice run.

```sh
python workflow/stages/04_region_integration/patient_splice/checks.py measure \
  --source data/splice-evidence --output designs/my-splice-checks \
  --remote-source /data/results/TN25-112103/patient-rna-splices-20261103-final \
  --remote-root /data/results/TN25-112103/splice-checks-NEW --execute

python workflow/stages/04_region_integration/patient_splice/checks.py status \
  --source data/splice-evidence --output designs/my-splice-checks --execute

python workflow/stages/04_region_integration/patient_splice/checks.py collect \
  --source data/splice-evidence --output designs/my-splice-checks --execute

python workflow/stages/04_region_integration/patient_splice/checks.py apply \
  --source data/splice-evidence --output designs/my-splice-checks

python workflow/stages/04_region_integration/patient_splice/checks.py verify \
  --source data/splice-evidence --output designs/my-splice-checks
```

Collect reads the actual recorded sandbox handles and remote paths; it downloads
only completed jobs and verifies worker/input/output hashes. Repeated launches
are blocked. --retry-failed allows retries only after a confirmed nonzero exit.
A changed measurement input or measurement policy requires a new result directory.
Offline admission is recorded separately as selection-policy.json; apply installs
the current ND policy without relabelling or altering measurement receipts. An
existing, different admission policy requires a new result directory. Status and
collect do not start substitute jobs when a process is merely slow.

To apply the included completed measurements offline, use the last two commands
with --output data/splice-evidence/normal-screen. No Modal login, patient BAM download
or new model inference is needed for those commands. The apply stage first runs
the independent verifier for the entire source evidence, then adds the two gates.

## Outputs and verification

- event-checks.tsv: all input events, detection labels, reasons and window counts.
- normal-tissues.tsv.gz: counts, denominators and matched UIDs by source tissue/strand.
- window-checks.tsv.gz: every reconstructed window's added decisions and orientation.
- all-windows.tsv.gz: original sequences/scores plus added checks and eligibility.
- event-ranking.tsv, gene-ranking.tsv, selected.tsv and selected.fasta: passing selections.
- previous-selection-checks.tsv: distinguish an old window surviving from its gene surviving with another event.
- normal-counts.json.gz and strand-counts.json.gz: measurements for offline replay.
- Receipts, source-verification.json and verification.json: input/model/worker hashes and complete-scope checks.
- selection-policy.json: the current ND-only admission rule, hashed independently
  of the policy used when the measurements were collected.

The implementation is in [check_measurements.py](check_measurements.py),
[checks.py](checks.py), [checks-policy.json](checks-policy.json),
[selection-policy.json](selection-policy.json) and
[splice_checks.py](../../../../src/yash_vaccine/splice_checks.py).
