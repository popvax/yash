# Class II typing and presentation-aware scoring

**Status: an analysis, not a finished construct.** It adds a missing input
(class II HLA typing) and re-scores Design #2's own candidates with a
presentation model. No cassette is assembled here, and nothing below has been
tested in a laboratory.

## The question

Design #2's methods say what it leaves out:

> This design uses class-I affinity predictions; class II, presentation and
> immunogenicity predictions do not contribute.

The README's "other directions to explore" names MHC class II and
immunogenicity first. Both are blocked on the same missing fact: `config/case.json`
records six class I alleles and **no class II typing at all**, so no
class II-aware selection is currently possible.

This directory does two things:

1. **Types Yash's class II alleles** from the matched-normal WGS already
   published with the project, and validates the method by recovering his known
   class I type. This is the part that unblocks future work.
2. **Re-scores the 17 mutant contexts** with MHCflurry's presentation model
   (affinity + antigen processing) instead of affinity alone, and reports where
   that disagrees with Design #2.

## 1. Class II typing

`hla/type_hla.py` pulls reads over the MHC from the 143 GB matched-normal BAM by
HTTP range request (~17 s, no full download), reduces them to 31-mers, and
scores every IMGT/HLA genomic allele by k-mer coverage, then picks the pair that
best explains the locus. Depth over the MHC is 40–95x.

### Validation first

Run on class I, the method must reproduce the alleles Caris
(`config/case.json`) and Guardant360 independently report. It does, **6 of 6**:

| Locus | Called | Known | |
|---|---|---|---|
| HLA-A | A\*01:01 / A\*33:03 | A\*01:01, A\*33:03 | match |
| HLA-B | B\*08:01 / B\*15:08 | B\*08:01, B\*15:08 | match |
| HLA-C | C\*01:02 / C\*07:02 | C\*01:02, C\*07:02 | match |

### Class II result

| Locus | Called | Coverage | Pair explains |
|---|---|---|---|
| HLA-DRB1 | DRB1\*03:01 / DRB1\*08:02 | 0.854 / 0.842 | 0.854 |
| HLA-DRB3/4/5 | DRB3\*02:02 present; DRB4, DRB5 absent | 0.850 / 0.216 / 0.116 | |
| HLA-DQA1 | DQA1\*05:01 / DQA1\*04:01 | 0.975 / 0.967 | 0.910 |
| HLA-DQB1 | DQB1\*02:01 / DQB1\*04:02 | 0.947 / — | 1.000 |
| HLA-DPA1 | DPA1\*01:03 / DPA1\*02:01 | 1.000 / 1.000 | 1.000 |
| HLA-DPB1 | DPB1\*02:01 / DPB1\*13:01 | 1.000 / 1.000 | 0.998 |

The calls resolve into two textbook haplotypes, which is a second check the
method never optimises for — each class II locus was typed independently:

```
Haplotype 1 : DRB1*03:01 - DRB3*02:02 - DQA1*05:01 - DQB1*02:01   (DR3-DQ2, DR52)
Haplotype 2 : DRB1*08:02 -      none  - DQA1*04:01 - DQB1*04:02   (DR8-DQ4)
```

Every DQ allele called is the one in known linkage with the DRB1 allele on its
haplotype, and the DRB3-present / DRB4-DRB5-absent pattern is exactly what a
DR3 + DR8 genotype predicts. Haplotype 1 together with A\*01:01-B\*08:01 is the
common 8.1 ancestral haplotype.

```sh
pip install pysam
cd hla
python type_hla.py                          # all loci -> hla_calls.json
python type_hla.py --loci A,B,C --out ctrl  # class I validation only
```

Caveats: 2-field resolution; the pair search uses a commonness tie-break, so a
genuinely rare allele could be reported as its common near-neighbour. Within
DQA1\*05 the subtypes encode near-identical alpha chains, so the group matters
more than the digits. The BAM is aligned to the GRCh38 primary assembly with no
alt contigs, so DRB3/4/5 reads land on DRB1; DRB3\*02:02 is therefore reported as
present rather than as a fully resolved independent call. None of this is a
substitute for a clinical typing laboratory, and it should be confirmed before
any clinical use.

## 2. Presentation-aware scoring

`presentation/score_presentation.py` reuses Design #2's own mutant protein
contexts and re-scores every 8–11mer overlapping a mutated residue with
MHCflurry 2.3.14's presentation predictor, which adds an antigen-processing
model trained on eluted ligands to the affinity model Design #2 used alone.

The script enumerates **604 peptides**, the same number Design #2's methods
report, which is a check that the peptide set was reproduced rather than
re-invented.

### What changes

**FES p.R804W is the clearest gap.** Design #2 excluded it for one stated
reason: "Only 15 aa available", too short for its fixed 25-aa window. By
presentation it ranks **4th of 17**, above **9 of the 12 regions that were
included**:

```
AYEPGQWPSF   HLA-C*07:02   94.5 nM   presentation 0.902
```

The exclusion was a window-format constraint, not a judgement about the
antigen. The limit is the length of the RNA-reconstructed context, not of the
FES protein, so a shorter region, or one padded with reference FES sequence
around the mutation, would recover it. That is worth deciding deliberately
rather than by default.

**The two affinity-gate rejects stay rejected.** SERPINA1 p.S22Y (0.111) and
CTDSP2 p.R46C (0.038) rank last by presentation too. Presentation and affinity
agree where Design #2 excluded on affinity, which is reassuring about both.

**Ranking within the cassette moves a lot.** Spearman rho between Design #2's
affinity rank and the presentation rank, over the 12 selected regions, is
**0.53**. HSPG2 falls from rank 5 to 12, PLEC from 11 to 15, while AQP3 rises
from 3 to 1. Design #2 included all 12 regardless of rank, so membership does
not change — but this matters for any design with less cassette space, and for
antigen order.

**The 25-aa window itself differs in 11 of 16 contexts.** For the same
mutation, the window maximising summed presentation is usually not the window
maximising summed affinity.

```sh
pip install 'mhcflurry==2.3.14' && mhcflurry-downloads fetch models_class1_presentation
cd presentation && python score_presentation.py --repo-root ../../..
```

Results: [`region-presentation.tsv`](presentation/region-presentation.tsv)
(17 contexts), [`peptide-presentation.tsv`](presentation/peptide-presentation.tsv)
(604 peptides).

## 3. What is not done here

**No class II binding predictions were run.** Every class II predictor I could
find is licence-encumbered: NetMHCIIpan and the IEDB class II tools need an
accepted academic licence, and MixMHC2pred's licence takes effect on download.
Accepting one on the project's behalf, or adding an academic-only dependency to
a public repository that maintains `THIRD_PARTY_NOTICES.md`, is the maintainers'
call, not mine.

So `classII/prepare_netmhciipan.py` does everything up to the binary. It writes
the 11 class II molecules in NetMHCIIpan allele syntax (three DR, plus cis and
trans DQ and DP heterodimers), the 24 antigen regions, the full 724-aa cassette
— scoring the whole cassette also catches epitopes that exist only across
antigen junctions — and a `run.sh`:

```sh
cd classII && python prepare_netmhciipan.py --repo-root ../../..
cd job && sh run.sh        # requires netMHCIIpan on PATH
```

The question that job answers is the one worth answering: how many of the 24
regions contain a CD4 epitope on any of Yash's class II molecules, and would a
different window or antigen set raise that number. A cassette selected purely
on class I may still contain class II epitopes by chance; whether this one does
is currently unknown, and now measurable.

## Incidental cross-check: RXRA

The Guardant360 report (21 Oct 2025) lists the RXRA variant as **A214S**, while
this repository and Caris call it **A241S**. The repository is right: the 24
residues preceding the mutation in `primary-proteins.tsv` align to positions
217–240 of RXRA (UniProt P19793), which puts the mutated residue at 241, and the
wild-type residue there is alanine. Position 214 is aspartate. No change needed;
recorded because the discrepancy is visible to anyone comparing the two reports.

Of the twelve mutation antigens, RXRA p.A241S and TGFBR2 p.N515Kfs\*26 are the
two that sit on the Guardant360 panel, and both were still detected in the
October 2025 plasma, rising with tumour fraction (12.1% and 11.7%). The other
ten genes are not on that panel, so their absence from the report says nothing
either way.

## Scope of the evidence

Everything here is prediction and genotyping from existing sequencing data. No
new sequencing, no laboratory work, no immunogenicity measurement. Binding and
presentation predictions do not establish that a peptide is presented on Yash's
tumour cells or that a T cell responds to it. The class II typing should be
confirmed by a clinical laboratory before it is relied on.

## Reproducing

| Input | Source |
|---|---|
| Matched-normal WGS BAM | `data/catalog/access.json` public base URL, `raw/TN25-112103/normal/normal.recal.bam` (range requests only) |
| IMGT/HLA allele sequences | `github.com/ANHIG/IMGTHLA` `Latest/fasta/*_gen.fasta`, fetched by the script |
| Mutant protein contexts | `designs/design-2/inputs/rna/primary-proteins.tsv` |
| Affinity ranking compared against | `designs/design-2/inputs/predictions/mutations/all-ranked-regions.tsv` |
| Antigen regions and cassette | `designs/design-2/construct/protein/` |

Tools: Python 3.12, pysam 0.24.1, mhcflurry 2.3.14 with the
`models_class1_presentation` bundle (2.3.0). MHCflurry needs numpy 1.26.4 with
scipy 1.13.1 and scikit-learn 1.5.2; newer scipy builds fail against that numpy.
`yash-vaccine verify` still passes unchanged.
