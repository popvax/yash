# Patient-RNA splice evidence

This shared evidence reconstructs splice-source candidates from Yash's actual
RNA bases before HLA prediction or selection. Reference annotation supplies a
coding-frame anchor, rather than replacement bases for his sequence.

## Find the evidence

| Files | Contents |
|---|---|
| [policy.json](policy.json) | Read, frame, reference, binding and expression rules |
| [Input events](reconstruction/input-events.tsv), [outcomes](reconstruction/event-audit.tsv) | Every normal-filtered event, including rejected and unresolved events |
| [Patient windows](reconstruction/patient-windows.jsonl.gz) | Coding sequences, amino acids, alignment paths, fragment names and stop-free upstream frame witnesses |
| [Reference transcripts](reconstruction/reference-transcripts.tsv) | Complete reference CDS-to-protein translation checks |
| [Native review](review/receipt.json) | LENS review of actual patient coding sequences |
| [Exact-self check](reference-check/receipt.json) | Full annotated-reference inventory and two complete exact-match scans |
| [Affinities](scoring/affinities.tsv.gz), [all windows](scoring/all-windows.tsv.gz) | Fresh six-allele class-I predictions and scored candidates |
| [Event ranking](scoring/event-ranking.tsv), [gene ranking](scoring/gene-ranking.tsv) | Best qualifying window per event and gene |
| [Top five](scoring/selected.tsv), [sequences](scoring/selected.fasta) | Splice regions used by Design #1 |
| [Independent BAM check](read-verification/verification.json) | Every window and named fragment checked through an independent alignment-array implementation |
| [Inference receipt](scoring/receipt.json) | Inputs, model hashes, package versions and executed worker |

Start with all 4,929 junction UIDs passing both native MSNAF normal-control
screens, before reference translation or HLA filtering. Require actual 75-nt
haplotypes in at least two qualifying fragment names. Overlapping mates must
agree; different fragments are never assembled into a synthetic consensus.
A verified GENCODE 37 upstream CDS anchors the frame, and observed codons from
that anchor through the join must be high-quality and stop-free.

Native LENS event and peptide review precedes admission. Exact reference-matching
join peptides do not contribute to scoring. The complete native Ensembl 103
inventory includes annotated protein isoforms and their source identities.
Fresh MHCflurry affinity predictions are checked against the OpenVax adapter.
An admitted region needs at least one nonself native-PASS peptide at <=500 nM.
Native Vaxrank affinity sums are multiplied by the square root of the supplied
gene value; no splice-fraction or raw-read abundance multiplier is used.

Rank windows, then the best per event and gene. Select five positive-expression
genes, one region per gene, with no sequence-identity filter. A window is a
25-aa region, an event is an RNA join, and a gene may contribute several events.

## Verify offline

After installing the repository, run from its root:

```sh
python workflow/stages/04_region_integration/patient_splice/replay.py \
  --results data/splice-evidence --output work/splice-verification.json
```

This checks the input partition and hashes, translations, fragment/path/frame
witnesses, independent BAM-check receipt, exact-self decisions, complete allele
coverage, affinity contributions, expression multiplication and gene ranking.
It starts no cloud job, rereads no BAM and makes no new model predictions.
The normal initial-design reproduction includes the same verification.

For raw-read processing and fresh inference, follow the
[upstream instructions](../../workflow/stages/04_region_integration/patient_splice/README.md).
Use matching staged inputs and authenticated Modal access. Explicit launchers
print plans unless --execute is provided.

A supported local sequence and stop-free frame do not establish a complete
translated isoform, tumour specificity, peptide presentation or a T-cell response.
Fragment names are not UMI molecule counts. Linear CIGAR gaps cannot establish
circular/back-splice joins; unresolved and nonlinear coordinates remain in the audit.
