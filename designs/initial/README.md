# Design #1 — without matched-normal DNA

This design used Yash’s tumour DNA, tumour RNA and supplied gene-expression values.
It contains **29 regions of 25 amino acids: 17 mutation, 5 splice, 5 CTA/self and
2 CDH1::NRG1 fusion regions**, assembled into an **852-aa protein**.

It is retained for comparison. **[Design #2 is the lead](../design-2/README.md)**
and adds matched-normal DNA to the mutation analysis.

## How we made it

1. **Find mutations and reconstruct their proteins.** Mutect2 called variants
   independently from tumour DNA. Quality and population-frequency filters
   narrowed the calls; Isovar reconstructed altered protein sequences supported
   by tumour RNA. Without matched-normal DNA, inherited variants could remain.
2. **Add RNA-source candidates.** RNA-only LENS supplied splice, fusion and
   cancer–testis antigen (CTA) candidates. Splice windows required observed RNA
   sequences in an annotated coding frame, strand support and no detected join
   in either searched normal-reference catalog. CTA candidates used reference
   proteins from genes with positive supplied expression.
3. **Score and select regions.** MHCflurry predicted 8–11-aa peptide binding to
   Yash’s six HLA-I alleles; Vaxrank combined these contributions within 25-aa
   windows. Expression-weighted ranking selected five splice genes and five CTA
   genes. All 17 qualifying mutation regions were included after hydropathy
   screening, alongside two fusion regions. The in-frame fusion was an explicit
   exception to the usual 500-nM binding cutoff.
4. **Assemble the protein.** We added patient-supported HLA-derived SEC/MITD
   trafficking sequences, jointly optimized antigen order and linker choices to
   reduce unintended junction binding, and checked the complete protein’s
   junctions and hydropathy.
5. **Encode the RNA.** [PopRNA](construct/mrna/optimizer.py) chose synonymous
   codons for the selected protein using human codon usage, motif and repeat
   avoidance, and RNA-folding heuristics. We verified that the encoded protein
   remained unchanged and retained the existing UTRs and stop codons. Sequence
   constraints and GenBank annotations were checked; **PolyA and downstream
   extensions remain unspecified**.

[Detailed methods](methods.md) explain the filters, expression proxies,
exceptions and limitations. These binding predictions do not establish
tumour-cell presentation or a T-cell response.

## Results and evidence

| What | Files |
|---|---|
| Selected antigens | [Table](selection/selected.tsv) · [Sequences](selection/selected.fasta) |
| Protein | [FASTA](construct/protein/cassette.fasta) · [Order](construct/protein/antigen-order.tsv) · [Components](construct/protein/components.tsv) |
| Nucleotide segments | [RNA](construct/mrna/mrna-segments.fasta) · [Coding DNA](construct/mrna/cds.fasta) · [GenBank](construct/mrna/template-segments.gb) |
| Antigen review | [Individual assessments](antigen-review.md) |
| Supporting inputs | [Design #1 inputs](inputs/README.md) · [Shared RNA evidence](../../data/shared-rna/README.md) · [Sequencing data](../../data/README.md) |

## Reproduce

From the repository root, with Python 3.11 or 3.12:

```sh
python -m pip install -e '.[mrna]'
yash-vaccine reproduce designs/initial --output work/design-1-replay
```

This replays selection, protein assembly, junction checks and RNA encoding from
bundled evidence and recorded predictions, and checks the outputs against the
recorded design. It downloads no reads and starts no cloud jobs.

See [usage](../../docs/usage.md) for alternatives and comparisons, or
[upstream processing](../../docs/upstream.md) to rerun from sequencing data.
