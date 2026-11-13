# Class II typing and presentation-aware scoring

**Status: an analysis, not a finished construct.** It adds a missing input
(class II HLA typing) and re-scores Design #2's own candidates with a
presentation model. No cassette is assembled here, and nothing below has been
tested in a laboratory.

## The question

Design #2's methods say what it leaves out:

> This design uses class-I affinity predictions; class II, presentation and
> immunogenicity predictions do not contribute.

The repository root README's "other directions to explore" names MHC class II
and immunogenicity first. Both are blocked on the same missing fact:
`config/case.json` records six class I alleles and **no class II typing at
all**, so no class II-aware selection is currently possible.

This directory does two things:

1. **Types Yash's class II alleles** from the matched-normal WGS already
   published with the project, validated against two unrelated public donors.
   This is the part that unblocks future work.
2. **Re-scores the 17 mutant contexts** with MHCflurry's presentation model
   (affinity + antigen processing) instead of affinity alone, and reports where
   that disagrees with Design #2.

## 1. Class II typing

`hla/type_hla.py` pulls reads over the MHC from the 143 GB matched-normal BAM by
HTTP range request (~17 s, no full download), reduces them to 31-mers, and
scores every IMGT/HLA genomic allele by k-mer coverage, then picks the pair that
best explains the locus. Depth over the MHC is 40–95x.

### Validation

Run on Yash's class I, the method reproduces all six alleles that Caris
(`config/case.json`) and Guardant360 report:

| Locus | Called | Reported | |
|---|---|---|---|
| HLA-A | A\*01:01 / A\*33:03 | A\*01:01, A\*33:03 | match |
| HLA-B | B\*08:01 / B\*15:08 | B\*08:01, B\*15:08 | match |
| HLA-C | C\*01:02 / C\*07:02 | C\*01:02, C\*07:02 | match |

**That check is weaker than it looks, in two ways.** The tie-breaks in the pair
search were settled while looking at this answer, so it is in-sample. And
neither comparator is clinical histocompatibility typing: Caris labels its
result "Predicted … HLA Genotype Results", and Guardant's is likewise an
inferred aid. They agree with each other, which is worth something, but neither
is a gold standard.

So the method was also run against two unrelated public donors with published
types, in `validation/`:

| Donor | Data | Truth source | Result |
|---|---|---|---|
| NA12878 / HG001 | GIAB 30x, GRCh37 | Cellosaurus CVCL_7526 (PubMed 27792722) | **16/16** over 8 loci, class I and II |
| HG002 / NA24385 | GIAB 60x, GRCh38 | clinical typing in arXiv 2209.09155 | **5/6** over 3 class I loci, at 2 fields |

NA12878 **found a real defect**, which is the main reason to keep it here. The
original candidate-admission margin of 0.05 excluded any second allele sitting
more than 0.05 below a perfectly covered first one — exactly what happens when
a divergent allele loses reads to mismapping. It lost DRB1\*03:01 and
DQB1\*02:01 and returned false homozygotes: 14/16. Both failures showed a
**low `pair explains`** (0.47 and 0.55) where every correct locus scored
0.875–1.000, so the metric flags its own failure mode. Widening the margin to
0.15 or more gives 16/16, and **Yash's eight loci are identical at every margin
from 0.05 to 0.30**, so the defect never affected his calls. The margin is now
0.20.

Having been used to diagnose that, NA12878 is a development set. HG002 is the
clean test: no parameter was chosen using it, and it is scored at full 2-field
resolution. Its one miss is C\*12:24 against a true C\*12:03 — the right allele
group, the wrong subtype — and it came with a *high* explains (0.908). So
`pair explains` catches a missed second allele but does **not** catch a subtype
error.

Taken together across both donors: **22/22 at allele-group (1-field)
resolution, 10/11 where 2-field truth is available.**

```sh
cd validation && python validate_held_out.py      # both donors, ~3 min
```

### Class II result

| Locus | Called | Coverage | Pair explains |
|---|---|---|---|
| HLA-DRB1 | DRB1\*03:01 / DRB1\*08:02 | 0.854 / 0.842 | 0.854 |
| HLA-DRB3/4/5 | DRB3\*02:02 present; DRB4, DRB5 absent | 0.850 / 0.216 / 0.116 | |
| HLA-DQA1 | DQA1\*05:01 / DQA1\*04:01 | 0.975 / 0.967 | 0.910 |
| HLA-DQB1 | DQB1\*02:01 / DQB1\*04:02 | 0.899 / 0.947 | 1.000 |
| HLA-DPA1 | DPA1\*01:03 / DPA1\*02:01 | 1.000 / 1.000 | 1.000 |
| HLA-DPB1 | DPB1\*02:01 / DPB1\*13:01 | 1.000 / 1.000 | 0.998 |

Depth per locus is recorded in [`hla/hla_calls.json`](hla/hla_calls.json) and
ranges from 40x (DRB1) to 95x (HLA-A).

The calls resolve into two haplotypes that are standard in immunogenetics, and
each class II locus was typed independently:

```
Haplotype 1 : DRB1*03:01 - DRB3*02:02 - DQA1*05:01 - DQB1*02:01   (DR3-DQ2, DR52)
Haplotype 2 : DRB1*08:02 -      none  - DQA1*04:01 - DQB1*04:02   (DR8-DQ4)
```

Every DQ allele called is the one in known linkage with the DRB1 allele on its
haplotype, and the DRB3-present / DRB4-DRB5-absent pattern is what a DR3 + DR8
genotype predicts. Haplotype 1 together with A\*01:01-B\*08:01 is the common 8.1
ancestral haplotype.

This coherence is corroborating but not independent evidence. The pair search
breaks near-ties toward better-documented alleles, and common alleles are
exactly the ones that form recognised haplotypes, so the method is biased
toward this outcome. What it does rule out is a call pattern that no real
haplotype would produce.

```sh
pip install pysam
cd hla
python type_hla.py                          # all loci -> hla_calls.json
python type_hla.py --loci A,B,C --out ctrl  # class I validation only
```

Caveats, in the order they would bite:

- **Subtype digits are the weak point.** The held-out donors were right at
  allele-group level every time (22/22) but missed one 2-field subtype in
  eleven (C\*12:24 for C\*12:03), and NA12878's B and DQB1 calls landed on rare
  subtypes of the correct groups. Treat the first field as solid and the second
  as probable. For binding prediction this mostly does not matter — within
  DQA1\*05, for instance, the subtypes encode near-identical alpha chains — but
  it is not a clinical-grade result.
- The pair search breaks near-ties toward better-documented alleles, so a
  genuinely rare allele can be reported as its common near-neighbour.
- `pair explains` detects a missed second allele but not a wrong subtype, so a
  high value is not a guarantee of correctness.
- The BAM is aligned to the GRCh38 primary assembly with no alt contigs, so
  DRB3/4/5 reads land on DRB1. DRB3\*02:02 is reported as present rather than
  as a fully resolved independent call.
- None of this substitutes for a clinical typing laboratory, and it should be
  confirmed there before any clinical use.

## 2. Presentation-aware scoring

`presentation/score_presentation.py` reuses Design #2's own mutant protein
contexts and re-scores every 8–11mer overlapping a mutated residue with
MHCflurry 2.3.14's presentation predictor, which adds an antigen-processing
model trained on eluted ligands to the affinity model Design #2 used alone.

The script enumerates **604 peptides**, matching the count Design #2's methods
report. On its own that count proves little — it is forced, being every 8–11mer
overlapping a mutation across the 17 published contexts. The real check is that
all **96** of Design #2's `best-peptide-per-HLA` rows appear among the 604, and
where the best allele agrees the predicted affinity matches to under 0.01 nM.
So this is the same peptide set, re-scored, not a different one.

### What changes

**FES p.R804W is the clearest gap.** Design #2 excluded it for one stated
reason: "Only 15 aa available", too short for its fixed 25-aa window. By
presentation it ranks **4th of 17**, above **9 of the 12 regions that were
included**:

```
AYEPGQWPSF   HLA-C*07:02   94.5 nM   presentation 0.902
```

The exclusion was a window-format constraint, not a judgement about the
antigen. FES is not missing expression either (19.5 in the supplied expression
table), and Design #2 scored its peptides before dropping it — the 604 only
reaches 604 because FES's 20 peptides are included; the other 16 contexts give
584.

The 15 aa is the length of the RNA-reconstructed context, not of the FES
protein. A shorter region, or one padded with reference FES sequence around the
mutation, *might* therefore recover this antigen — but both would break rules
Design #2 holds to (a fixed 25-aa window, and RNA-supported sequence only), so
this is a suggestion to weigh, not a fix to apply.

**The two affinity-gate rejects stay rejected.** SERPINA1 p.S22Y (0.111) and
CTDSP2 p.R46C (0.038) rank last by presentation too. Presentation and affinity
agree where Design #2 excluded on affinity, which is reassuring about both.

**Ranking within the cassette shifts, but agrees more than it disagrees.**
Spearman rho between Design #2's `primary_rank` and the presentation rank, over
the 12 selected regions, is **0.71**. Ranking within those 12, HSPG2 falls from
5 to 10 and PLEC from 11 to 12, while AQP3 rises from 3 to 1. Design #2
included all 12 regardless of rank, so membership does not change; this would
matter for a design with less cassette space.

Two caveats on that comparison. `primary_rank` is not an affinity ranking: it
orders by affinity-contribution-sum × sqrt(expression proxy), so part of the
disagreement is expression weighting rather than presentation-versus-affinity.
And Design #2 states that this score "does not determine a top-N mutation
cutoff or construct order", so the shift does not imply its cassette order is
wrong.

**The 25-aa window itself differs in 11 of the 16 contexts where a window
exists.** (FES is the 17th; its 15-aa context admits no 25-aa window at all.)
For the same mutation, the window maximising summed presentation is usually not
the one Design #2 chose.

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
October 2025 plasma at cfDNA variant allele fractions of 12.1% and 11.7%, up
from 0.5% and 0.7% in March, alongside a rise in tumour fraction from 1.4% to
31.9%. (Those are variant allele fractions, not tumour fraction.) The other ten
genes are not on that panel, so their absence from the report says nothing
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
| Held-out validation BAMs | GIAB: NA12878 `NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam` (GRCh37); HG002 `HG002.GRCh38.60x.1.bam` |
| Held-out truth | Cellosaurus CVCL_7526 (PubMed 27792722); arXiv 2209.09155 |
| Affinity ranking compared against | `designs/design-2/inputs/predictions/mutations/all-ranked-regions.tsv` |
| Antigen regions and cassette | `designs/design-2/construct/protein/` |

Every number in this README is recomputable from the committed TSV/JSON files;
`validation/held-out.json` records the two public-donor runs.

Tools: Python 3.12, pysam 0.24.1, mhcflurry 2.3.14 with the
`models_class1_presentation` bundle (2.3.0). MHCflurry needs numpy 1.26.4 with
scipy 1.13.1 and scikit-learn 1.5.2; newer scipy builds fail against that numpy.
`yash-vaccine verify` still passes unchanged.
