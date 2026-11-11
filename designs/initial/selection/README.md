# Antigen selection

Select 29 distinct 25-aa regions: all 17 segmented mutation regions passing the hydropathy screen,
one region from each of five splice genes, five CTA/self genes and two fusion regions.
Splice windows require actual patient RNA bases, complete-window support in at
least two qualifying fragment names, a verified annotated upstream frame and
no observed intervening stop codon. They also require ND in both searched normal
catalogs and a strict majority of full-window witnesses agreeing with the nominated
gene under the measured library orientation.

[Selected regions](selected.tsv) · [Sequences](selected.fasta) · [Recipe](design.json) ·
[Selection audit](selection-audit.tsv) · [Splice verification](splice-rna-verification.json) ·
[Gene ranking](splice-gene-ranking.tsv) · [Methods](../methods.md)

The compact pool contains 223 entries: 18 mutation, 42 splice, 161 CTA/self and 2 fusion.
Mutation DNA eligibility is calculated from native calls and population observations
using the [research policy](../../../data/mutation-evidence/policy.json). It retains
native GATK quality filters, with no extra median MQ/BQ thresholds or ancestry-group
sample-size minimum. [All DNA decisions](../../../data/mutation-evidence/eligibility/decisions.tsv)
are included; reproduction verifies them before mutation-window scoring.
Native affinity contributions are independently recounted from the recorded
model predictions and exact-reference matches. Region scores use the square
root of the documented expression proxy. Splice weighting uses gene expression
alone, with no splice fraction or raw-read abundance multiplier. The in-frame
fusion is an explicit 508.95 nM driver-junction inclusion, with no fabricated
expression score. See [fusion evidence](../../../data/fusion-evidence/README.md).

## Construct-compatibility gate

All 223 nominated regions receive a [hydropathy screen](hydrophobicity-screen.tsv)
before source quotas. A complete 19-aa Kyte–Doolittle mean strictly above 1.6
excludes a region. PTCH1 is the sole exclusion; upstream DNA/RNA evidence remains
intact. This is a conservative heuristic, not a topology validation. See the
[policy, independent calculation and full-chain checks](../review/hydrophobicity/README.md).
The solver also excludes incompatible internal joins and independently rescans the
assembled antigen/linker cassette.

## Optional internal-cassette helper

The recipe runner also constructs a cassette without SEC/MITD for comparisons
using that helper. The final encoded **852-aa protein** is jointly assembled
with the patient-supported routing elements in [construct/](../construct/README.md).

```sh
yash-vaccine reproduce designs/initial/selection
yash-vaccine init my-candidate --from designs/initial/selection
```

Different selections need matching predictions and checks. The optional helper
does not automatically apply SEC/MITD or nucleotide encoding to another design.
