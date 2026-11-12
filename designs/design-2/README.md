# Design #2 — lead design

This design uses Yash’s **matched-normal white-blood-cell DNA**, tumour DNA,
tumour RNA and supplied gene-expression values. It contains **24 regions of
25 amino acids: 12 mutation, 5 splice, 5 CTA/self and 2 CDH1::NRG1 fusion
regions**, assembled into a **724-aa protein**.

[Design #1](../initial/README.md) was made without matched-normal DNA. Seven of
its mutation targets were detected in normal DNA and excluded here; AQP3 and
TGFBR2 were added. All twelve splice, CTA and fusion regions were retained.

## How we made it

1. **Call mutations using paired DNA.** We aligned and quality-checked the normal
   reads, confirmed the tumour–normal sample match, and ran paired Mutect2
   discovery. Native filtering, population-frequency checks and direct normal
   evidence left 99 primary DNA candidates.
2. **Reconstruct RNA-supported proteins.** Isovar examined the tumour RNA and
   reconstructed altered protein contexts. Candidate peptides were checked
   against annotated reference proteins and sequence possibilities derived from
   normal DNA. MHCflurry predicted binding to Yash’s six HLA-I alleles; Vaxrank
   selected 25-aa windows by combined binding score among RNA-supported windows.
   **All 12 qualifying mutation regions were included**, with no top-N cutoff.
   Positive supplied expression was required; expression-weighted ranks were
   descriptive.
3. **Integrate the RNA sources.** We reused the recorded LENS splice, CTA/self
   and fusion evidence, rechecked altered targets against the normal-derived
   sequences, and freshly predicted fusion binding. Expression-weighted ranking
   selected five splice genes and five CTA genes; two fusion regions were
   included. Splice joins had no detections in either searched normal catalog.
   The in-frame fusion retained its explicit exception to the 500-nM cutoff.
4. **Rebuild and check the protein.** Patient-supported HLA-derived SEC/MITD
   trafficking sequences were fixed. Fresh predictions supported joint antigen
   order/linker optimization to reduce unintended junction binding, followed by
   independent junction and hydropathy checks.
5. **Encode the RNA.** [PopRNA](construct/mrna/optimizer.py) chose synonymous
   codons for the selected protein using human codon usage, motif and repeat
   avoidance, and RNA-folding heuristics. We verified that the encoded protein
   remained unchanged and retained the existing UTRs and stop codons. Sequence
   constraints and GenBank annotations were checked; **PolyA and downstream
   extensions remain unspecified**.

[Detailed methods](methods.md) document the thresholds and limitations, including
the additional RNA changes in four splice candidates. Binding predictions do
not establish tumour-cell presentation or a T-cell response.

## Results and evidence

| What | Files |
|---|---|
| Selected antigens | [Table](selection/selected.tsv) · [Sequences](selection/selected.fasta) |
| Protein | [FASTA](construct/protein/cassette.fasta) · [Order](construct/protein/antigen-order.tsv) · [Components](construct/protein/components.tsv) |
| Nucleotide segments | [RNA](construct/mrna/mrna-segments.fasta) · [Coding DNA](construct/mrna/cds.fasta) · [GenBank](construct/mrna/template-segments.gb) |
| Antigen review | [Individual assessments](review/antigen-review.md) · [Evidence table](review/antigen-review.tsv) |
| Supporting inputs | [Design #2 inputs](inputs/README.md) · [Shared RNA evidence](../../data/shared-rna/README.md) · [Sequencing data](../../data/README.md) |

## Verify or reproduce

From the repository root, with Python 3.11 or 3.12:

```sh
python -m pip install -e '.[mrna]'
yash-vaccine verify
yash-vaccine reproduce --output work/design-2-replay
```

`verify` checks recorded evidence and outputs. `reproduce` reselects regions,
re-solves order/linker optimization from recorded costs, checks junctions and
reruns RNA encoding. Selected sequences and nucleotide outputs must match byte
for byte. Neither command downloads reads or starts cloud jobs.

See [usage](../../docs/usage.md) for alternatives and comparisons, or the
[pipeline guide](pipeline/README.md) for new raw-read processing and inference.
