# Design #2 — lead design with matched-normal DNA

This lead design uses Yash’s matched-normal white-blood-cell DNA alongside tumour
DNA and RNA. It contains **24 distinct 25-aa regions: 12 mutation, 5 splice,
5 CTA/self and 2 fusion regions**, assembled into a **724-aa protein**.

**All 12 eligible mutation regions are included.** Expression-weighted scores
order comparison tables; they do not impose a top-N mutation cutoff. Positive
available expression remains an eligibility requirement. Splice and CTA
selection still use expression-weighted ranking. [Methods](methods.md) document
all gates, source-specific scoring and limitations.

Four splice candidates rely on additional RNA changes; see the
[sequence details](methods.md#7-integrate-recorded-rna-source-candidates).

[Design #1](../initial/README.md), generated without matched-normal DNA, remains
available for comparison: 29 regions / 852 aa. Seven of its mutation targets were
detected in normal DNA and removed; AQP3 and TGFBR2 were added. Twenty-two region
sequences are retained, including all twelve RNA-source sequences.

## Verify or reproduce

Install from the repository root with Python 3.11 or 3.12:

```sh
python -m pip install -e '.[mrna]'
yash-vaccine verify
yash-vaccine reproduce --output work/design-2-replay
yash-vaccine compare designs/initial designs/design-2
```

`verify` checks the recorded evidence, source snapshot, selected regions, protein,
known nucleotide segments and GenBank annotations. It is read-only.

`reproduce` reselects regions from the included candidate pool, solves the joint
order/linker problem again using recorded edge costs, independently recounts the
chosen junctions against included predictions, and reruns synonymous encoding.
The selected TSV/FASTA, protein/components/order and nucleotide outputs must match
byte for byte. This is an **offline replay of recorded measurements**; it does not
repeat raw-read discovery or model inference, download data or start cloud jobs.

Explicit entry points are also available:

```sh
yash-vaccine reproduce designs/design-2 --output work/lead-replay
python designs/design-2/verify_results.py
python designs/design-2/run.py --output work/another-lead-replay
```

Choose a new output directory for each replay. [Raw-read processing and Modal
controllers](pipeline/README.md) have separate inputs, runtimes and execution
requirements.

## Current results

- [Selection](selection/selected.tsv), [sequences](selection/selected.fasta), [DNA policy](selection/paired-policy.json).
- [Every-antigen review](review/antigen-review.md), [sortable evidence table](review/antigen-review.tsv), [completed matched-normal reassessment](review/final/report.md).
- [Protein](construct/protein/cassette.fasta), [order](construct/protein/antigen-order.tsv), [components/linkers](construct/protein/components.tsv).
- [Coding DNA](construct/mrna/cds.fasta), [known RNA segment](construct/mrna/mrna-segments.fasta), [annotated DNA segment](construct/mrna/template-segments.gb).
- [Junction check](construct/final-check.json), [hydropathy](construct/protein/hydrophobicity.json), [nucleotide/GenBank checks](construct/mrna/verification.json).
- [Execution audit](execution/completion-audit.json), [native receipts and local path map](execution/README.md).

PolyA and any downstream extension are **unspecified**. The nucleotide files are
known upstream segments, not complete molecular sequences. The coding sequence
has 2,172 sense nucleotides. Hard translation/sequence checks and all 64 GenBank
features pass. The encoder’s final paired-bracket heuristic is 21 nt versus its
14-nt target; that soft target was not retained. Predictions and heuristics do
not establish processing, presentation or a T-cell response.

## Evidence and data access

Independent paired Mutect2 discovery produced 5,210 native-PASS alleles, then
342 protein-changing variants and 99 primary DNA candidates. RNA review covered
101 inputs and produced 17 leading protein contexts; 12 regions met the final
mutation criteria. The fingerprint matched, and 97.23% of the assessed normal
coding-plus-flank bases met GQ≥20/depth≥20. That is not whole-genome callability.
The normal instrument is unconfirmed. See [methods](methods.md) for exact scope.

[Matched-normal data](../../data/matched-normal/README.md) provides the final
normal BAM/index, scoped GVCF, paired calls, QC and exact references through the
[documented downloader](../../tools/r2/README.md), without credentials.
Original normal FASTQs are excluded.

| Location | Contents |
|---|---|
| selection/ | Canonical chosen antigens, sequences and paired DNA policy |
| construct/ | Protein, partial nucleotide segments, annotations and assembly evidence |
| exports/ | Input/evidence tables for replay and independent analyses |
| review/ | Current individual assessments and final matched-normal reassessment |
| review/historical/ | Superseded intermediate normal-only checkpoint |
| pipeline/ | Raw-read/Modal controllers, workers and recorded input configuration |
| execution/ | Immutable native receipts, source-adaptation notes and canonical-path mapping |

Methods, sequences and evidence are design-specific. Native receipts retain their
original run paths and dates; [recorded source snapshots](../../provenance/snapshots/README.md)
preserve execution-time code separately from the maintained implementation.
