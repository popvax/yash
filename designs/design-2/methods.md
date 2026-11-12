# Design #2 — detailed methods

Design #2 is the lead: **24 distinct 25-amino-acid regions in a 724-aa protein**,
comprising 12 mutation, 5 splice, 5 cancer–testis antigen (CTA)/self and 2 fusion
regions. It uses Yash’s tumour DNA and RNA, matched-normal white-blood-cell DNA,
and supplied gene-expression values. [Design #1](../initial/methods.md) used no
matched-normal DNA.

## Workflow at a glance

| Stage | Recorded result |
|---|---|
| Prepare normal DNA | 603,687,217 read pairs; tumour–normal fingerprint match |
| Discover variants from both DNA samples | 146,864 emitted ALT alleles; 5,210 native PASS |
| Apply coding annotation and DNA policy | 342 protein-changing variants; 99 primary DNA candidates |
| Reconstruct tumour-RNA proteins | 101 inputs including 2 review-tier variants; 17 usable leading contexts |
| Screen mutation windows | 12 qualifying regions; all 12 retained |
| Add recorded RNA-source candidates | 5 splice, 5 CTA/self and 2 fusion regions |
| Assemble and encode | 724-aa protein; checked RNA/DNA segments with PolyA unspecified |

Follow [normal preparation](#1-normal-alignment-and-sample-identity),
[normal evidence](#2-normal-genotypes-and-explicit-normal-evidence),
[paired discovery](#3-independent-paired-discovery),
[DNA eligibility](#4-additional-dna-research-policy),
[RNA reconstruction](#5-fresh-tumour-rna-reconstruction),
[binding and selection](#6-personal-normal-self-sequences-and-class-i-ranking),
[RNA-source candidates](#7-integrate-recorded-rna-source-candidates),
[assembly](#8-joint-assembly-and-fresh-independent-checks) and
[RNA encoding](#9-rna-encoding-genbank-and-replay).

A **region** is a 25-aa vaccine block; a **peptide** here is an 8–11-aa sequence
scored for HLA-I binding. One region can contain several peptides, each scored
against six alleles. These counts are not interchangeable. **Fragment** means
an alignment read-name group, including its paired mates, rather than a
UMI-confirmed original molecule.

The [data guide](../../data/README.md) provides the sequencing files and small
tables; [inputs/](inputs/README.md) organizes this design’s processed evidence.

## 1. Normal alignment and sample identity

We aligned the normal reads, marked duplicates and recalibrated base qualities
before paired calling. The sequencing instrument is unknown; the original
header was retained and GATK used `--default-platform UNKNOWN`.

| Measurement | Result |
|---|---:|
| Compressed input size | 130,524,712,208 bytes |
| Read pairs / read length | 603,687,217 / 150 nt |
| Reported Q30 | 96.13% |
| Primary mapping / proper pairing | 99.99% / 99.31% |
| Marked primary duplicates | 2.27% |
| Tumour–normal fingerprint LOD | 48.887978; required ≥5 |

Source files passed full-object CRC64NVME and SHA-256 checks on the persistent
Modal volume. fastp was used for QC only; it removed no reads or bases.
BWA-MEM2, fixmate, coordinate sorting and duplicate marking produced the indexed
alignment.

The compatible GRCh38 analysis FASTA has SHA-256
`c3e61ad6d770db22a6f27fdbcd3dba93d0d80f649018c8887ec803f12ea90450`.
Caris’s exact custom-masked alignment reference is unavailable.

Native GATK base-quality score recalibration (BQSR) used dbSNP/Mills known sites.
Training was scattered over primary nuclear intervals. Application used disjoint
whole-contig groups and a separate no-coordinate tail, followed by native
gathering. A 16,111-record comparison established exact recalibrated-record
equivalence with whole-file application. Gathering preserved all 1,212,214,833
BAM records, including supplementary and unmapped records. Final BAM SHA-256:
`1d864150551f6d47260050c52105c88b7b9f4120d005be33656094767326ce18`.

Fingerprinting retained all 933 markers in 240 blocks. The Broad map header was
adapted to the validated shared BAM dictionary after coordinate/reference-allele
checks; marker evidence was unchanged. CalculateContamination estimated 0 for
the normal and 0.000203642 for the tumour. These are cross-sample contamination
estimates, not tumour-purity measurements.

Evidence: [alignment and QC records](execution/README.md) ·
[pipeline and runtime instructions](pipeline/README.md).

## 2. Normal genotypes and explicit normal evidence

### Coding genotypes

HaplotypeCaller GVCFs cover primary-nuclear annotated CDS, start and stop intervals
plus 100-bp flanks. **73,917,616 of 76,021,670 bases (97.23%)** had GQ≥20 and
depth≥20, using `MIN_DP` for reference blocks. This is coding-plus-flank
callability, not whole-genome callability. Gathering produced no scatter conflicts.

### Evidence at each candidate allele

We combined native normal genotypes with direct fragment evidence from the
recalibrated BAM. Eligible normal reads use MQ≥20 and focal BQ≥20.

| Classification | Required evidence |
|---|---|
| `NORMAL_REFERENCE_SUPPORTED` | GQ≥20 genotype excludes the target allele, at least 20 usable fragment names and fewer than 3 high-quality ALT fragments |
| `PRESENT_IN_NORMAL` | Confident target-containing genotype with observed alternate support |
| Unresolved | Insufficient coverage, incomplete/complex genotype or conflicting evidence |

Missing normal calls are never treated as measured absence. Presence in blood
cells does not by itself distinguish an inherited variant from every other origin.

### Reassessment of Design #1

All 17 mutation targets from Design #1 were checked, including targets not emitted
by paired discovery. Seven were present in normal DNA: **KMT2D G3189R, PARP3 K37Q,
CAST I309V, CCNI P285S, PPP2R1A K280N, NRXN2 D203H and DOCK8 G1251R**. The remaining
ten had adequate normal reference evidence and passed the paired DNA policy.
None remained unresolved because of missing normal coverage.

Evidence: [17-target audit](review/final/current-17-final-audit.tsv). It separates
failure to rediscover an exact allele from rejection by a native filter.

## 3. Independent paired discovery

GATK **Mutect2 4.6.2.0** called variants from both BAMs with the normal sample
explicitly identified. Discovery covered 24 primary nuclear chromosomes in
**143 shards**. Neither Caris’s VCF nor the 17 earlier targets nominated alleles.
The vendor capture BED was unavailable, so the plan used chromosome-wide
intervals. No panel of normals was supplied.

The calling germline resource was corrected QC-PASS gnomAD v4.1 genomes.
Every emitted allele also received gnomAD v4.1 genome/exome review. Before final
filtering, its population prior was set to the maximum of the recorded native
prior and available QC-PASS global cohort frequencies. This corrected priors
in emitted calls; it did not rerun local assembly. Missing or QC-filtered
frequencies remained unknown, and lookup failures stopped processing.
Conflicting duplicate/scatter measurements were held for review.

Native **FilterMutectCalls** combined sample-wide calling statistics, the learned
F1R2 orientation model, contamination estimates and allele-fraction segments.
We selected `FALSE_DISCOVERY_RATE` with a **0.05 model target**, rather than
`OPTIMAL_F_SCORE`. The target is not an experimentally measured error rate.

The result contained **143,552 raw records, 146,864 normalized ALT alleles and
5,210 native-PASS alleles**, with no scatter conflicts.

Evidence: [calling/filtering audits](review/README.md) ·
[executed commands and completion records](execution/README.md).

## 4. Additional DNA research policy

Of the native-PASS alleles, 471 overlapped coding/splice-boundary annotations and
342 had predicted protein-changing consequences. We applied the following
additional research eligibility rules:

| Rule | Requirement |
|---|---|
| Native site and effective allele filters | PASS |
| Tumour depth / alternate AD | ≥20 / ≥5 |
| Model tumour fraction | ≥5% |
| Germline quality (`GERMQ`) | ≥20 |
| Population frequency | Exclude if either QC-PASS gnomAD cohort reports global or recorded maximum-frequency ancestry-group AF≥0.1%; no group-size minimum |
| Direct normal evidence | `NORMAL_REFERENCE_SUPPORTED`, as defined above |
| Context review | No unresolved sex/ploidy, low-complexity/duplication, population-prior, normal-disagreement or exact Caris FP/low-quality flag |

Native GATK handles mapping/base-quality filtering; we added no tumour alternate
median MQ/BQ cutoff. Caris annotations were comparison/review evidence, not
nomination or proof of somatic status. Unknown population observations were not
converted to AF=0. The numerical gates are research choices, not validated
probabilities of tumour specificity.

This left **99 primary variants in 98 genes**. Two variants with unresolved normal
evidence entered RNA review as a separate tier and could not enter primary
selection.

Evidence: [exact policy](selection/paired-policy.json) ·
[all variant assessments and rejections](inputs/mutations/variant-evidence.tsv).

## 5. Fresh tumour RNA reconstruction

**Isovar 1.57.1**, Varcode **10.16.0** and Pyensembl **2.27.3** examined the original
tumour RNA BAM at all 101 review variants. Isovar assembled mutant-bearing RNA,
matched compatible transcript frames and translated local protein contexts.

Read eligibility excluded duplicate, secondary and supplementary alignments,
required MQ≥20 and required NH=1 when present. Mate names were unioned as
fragments. No additional three-fragment minimum was imposed. Repeat-shifted
indel representations were explicitly checked for equivalence.

| RNA outcome | Inputs |
|---|---:|
| Local coverage, no alternate allele | 58 |
| No eligible local coverage | 7 |
| Alternate allele observed, no protein reconstructed | 17 |
| Requires Isovar filter review | 2 |
| Usable leading protein context | 17 |
| **Total** | **101** |

The 17 leading contexts came from 797 total hypotheses. All 17 were DNA-primary;
16 also had unambiguous gene attribution. The receipt field
`primary_DNA_and_RNA_eligible=16` includes that attribution condition. Lack of
RNA reconstruction or coverage was not called biological absence.

Evidence: [RNA outcomes](inputs/rna/rna-evidence.tsv) ·
[reconstructed proteins](inputs/rna/primary-proteins.tsv) ·
[RNA inputs and representation audit](inputs/rna/).

## 6. Personal-normal self sequences and class-I ranking

### Reference and personal-normal comparison

Normal coding models combined nearby variants, retained measured local phase
blocks and enumerated unphased combinations as **possible self sequences**, with
a 4,096-combination cap. Frameshift models retained downstream coding variants.
Compatible full HLA reference proteins came from the earlier measured two-field
typing; this was not new full normal HLA typing.

There were **112,353 local/HLA sequence possibilities** and **8,987 unresolved
model records**. Exact altered-target hits were conservatively excluded. A
negative search covers these bounded models and annotated references, not every
possible personal protein. None of the 604 mutation peptides below matched the
searched reference or normal models.

### Binding predictions and window scores

MHCflurry **2.3.14**, affinity-model bundle **2.3.0**, freshly scored 604 unique
mutation-overlapping 8–11-aa peptides against six alleles:
**A*01:01, A*33:03, B*08:01, B*15:08, C*01:02 and C*07:02** — **3,624 peptide–HLA
pairs**. This design uses class-I affinity predictions; class II, presentation
and immunogenicity predictions do not contribute.

Vaxrank **3.46.0** scored **262 fully supported 25-aa windows**. We summed native
per-allele affinity contributions and checked arithmetic parity. The highest
sum selected the window within each context; native centering order broke ties.
**Expression did not choose the mutation window.**
The [native contribution curve](../initial/methods.md#4-score-and-choose-25-aa-regions)
converts affinity to window contributions; the separate 500-nM rule below
determines admission.

For descriptive ranking, we also calculated:

```text
RNA ALT fraction = ALT fragment names / union(REF, ALT, OTHER fragment names)
Expression proxy = supplied gene value × RNA ALT fraction
Weighted score   = affinity contribution sum × sqrt(expression proxy)
```

Gene abundance came from supplied transformed expression, not RNA read counts.
Its calibration as ordinary TPM is unresolved. This score orders comparison
tables; it does not determine a top-N mutation cutoff or construct order.

### Final mutation eligibility

A region needed primary DNA eligibility, unambiguous gene assignment, positive
available expression, a supported 25-aa window and at least one non-reference
target peptide with affinity **≤500 nM**. **All 12 passing regions were kept.**
Expression availability was an eligibility requirement; its magnitude and the
ALT fraction did not remove lower-ranked qualifying regions.

| Leading context not selected | Reason |
|---|---|
| FES | Only 15 aa available |
| SERPINA1; CTDSP2 | Best eligible affinity >500 nM |
| H2BC12 | No supplied expression after approved/previous/alias-symbol checks |
| AL358113.1 | Missing expression and ambiguous gene attribution |

Missing expression remained missing rather than zero.

Evidence: [mutation scores and HLA summaries](inputs/predictions/mutations/) ·
[candidate pool](inputs/integrated/candidate-pool.tsv) ·
[individual antigen review](review/antigen-review.md).

## 7. Integrate recorded RNA-source candidates

We reused the recorded RNA-only LENS discovery and coordinate-aware normal splice
screen. Shared measurements retain their checksums; this was not a new raw-RNA
discovery run. [Design #1’s RNA-source methods](../initial/methods.md#3-add-the-rna-only-lens-sources)
explain their discovery, reconstruction and native review in full.

| Source | Admission and selection in Design #2 |
|---|---|
| Splice | Complete observed 75-nt window; ≥2 full-window fragment names; stop-free verified upstream frame; strand consistency; no exact junction entry in either searched normal catalog. Rank with supplied gene expression; one region per selected gene. |
| CTA/self | Retain the reference-protein hypotheses, native review and expression-weighted selection of five genes. These intentionally include self sequences. |
| Fusion | Recheck altered targets against personal-normal models and freshly infer affinities with original-model parity checks; retain two distinct breakpoint/translation contexts. |

Splice/fusion target scoring was rechecked against the personal-normal models.
Splice/CTA affinity measurements were reused. The original event/region tie-break
was retained. Five splice and five CTA regions remained unchanged. **ND** means
no exact resolvable junction entry in the searched catalogs; it does not show
absence from all healthy tissues or establish translation. Low support and
motif disagreement remain visible.

### Additional RNA changes in the splice sequences

For **F5, MUC5B, MUC3A and MUC6**, reads show both an unusual splice join and an
extra base change. The extra change makes the selected peptides different from
the main reference proteins. [Sequence checks](review/splice-sequence-check.json)
record those comparisons and the targeted MUC6 normal-DNA check, independently
reproducing the reference comparison in [issue #1](https://github.com/popvax/yash/issues/1).

To repeat the checks, obtain Ensembl 103 proteins from the
[provider guide](../../docs/external-resources.md), install `pysam==0.23.3`, then run:

```sh
python designs/design-2/pipeline/check_splice_sequences.py \
  --reference-proteins /path/to/Homo_sapiens.GRCh38.pep.all.fa.gz \
  --check-normal --output work/splice-sequence-check.json
```

`--check-normal` reads the indexed MUC6 normal-BAM region over public HTTPS.
Normal DNA does not verify tumour-specific splicing. For CTA/self targets, it
also does not establish normal-tissue expression or the exact translated patient
isoform; **SPATA31A6 has particularly weak local RNA evidence**.

### The two fusion contexts

The alternative frameshift-tail context uses its own exact-breakpoint RNA
fraction. The in-frame driver context has a **508.95-nM inclusion exception** and
no borrowed expression score. Its NAT sequon is a potential glycosylation site;
occupancy was not modelled. [Every-region review](review/antigen-review.md)
records the sequence and normal-tissue annotation issues for all sources.

## 8. Joint assembly and fresh independent checks

All 24 regions passed the 19-aa Kyte–Doolittle screen, which rejects a region
if any window mean is >1.6. Fixed routing elements were the **HLA-C*01:02 leader
plus CS (26 aa)** and **HLA-B*15:08 terminal domain (55 aa)**. Their original
tumour-DNA/RNA support is in [allele evidence](../initial/construct/allele-evidence/README.md);
it was not relabelled as new normal typing or measured cleavage.

We evaluated all **11 linker choices** for every directed internal and terminal
option: **6,600 options, 96,352 unique peptides and 578,112 fresh peptide–HLA
predictions**. Internal high-hydropathy options were excluded before optimization.

OR-Tools CP-SAT used one worker and seed 31. It solved these objectives in order,
preserving each earlier optimum:

1. Minimize junction peptide–HLA rows at affinity percentile rank ≤0.5%.
2. Minimize rows at rank ≤2%.
3. Maximize the lowest junction percentile rank.
4. Minimize inserted linker residues.

All four phases were proven optimal within this fixed search space. This reduces
predicted unwanted-junction binding; it does not optimize cleavage or intended
antigen presentation.

The protein contains **600 antigen + 26 SEC + 55 MITD + 43 linker = 724 aa**.
Independent whole-protein enumeration matched the pairwise objective. A separate
fresh prediction run checked **6,132 selected-junction peptide–HLA pairs** with
unchanged threshold classifications.

| Selected-junction result | Peptide–HLA rows |
|---|---:|
| Affinity rank ≤0.5% | 1 |
| Affinity rank ≤2% | 26 |
| Affinity ≤500 nM | 24 |
| Affinity ≤50 nM | 0 |

The thresholds overlap; these are not counts of distinct antigens. No unexpected
high-hydropathy antigen/linker windows remained. Routing and boundary overlaps
were reported separately, without experimental topology validation.

Evidence: [assembly and solver record](construct/protein/construct.json) ·
[all option costs](construct/protein/edge-costs.tsv) ·
[independent junction check](construct/final-check.json) ·
[full-chain hydropathy](construct/protein/hydrophobicity.json).

## 9. RNA encoding, GenBank and replay

### PopRNA settings and checks

[PopRNA](construct/mrna/optimizer.py) chose synonymous codons for the fixed protein.
It used human codon usage, cutoff **0.10**, seed **0** and ViennaRNA **2.7.2**.
Its sequential heuristics address configured motifs, repeats, alternative-frame
stops, codon distribution and RNA folding. The existing UTRs, **TAA–TAG–TGA** stops
and known template context were retained. Exact settings and implementation
changes are in [backbone.json](construct/mrna/backbone.json).

All hard translation, motif, codon and component checks passed. Strict GenBank
round-tripping preserved **64 features and 24 antigen identities**. Coding GC was
**56.08%**, and no repeated 12-nt pairs remained. The final paired-bracket run
was **21 nt**, exceeding the heuristic’s **14-nt target** after sequential repairs.
That soft target was not retained; these heuristics do not establish an
expression or stability benefit.

### What the nucleotide files specify

| Sequence | Length |
|---|---:|
| Sense coding sequence | 2,172 nt |
| Coding sequence including first stop | 2,175 nt |
| Known RNA segment | 2,315 nt |
| Known DNA-context segment | 4,256 bp |

**PolyA and downstream extensions remain unspecified.** Complete molecular
lengths, cap chemistry and nucleotide modifications are not assigned. The files
represent known segments rather than complete RNA/plasmid molecules.

Evidence: [RNA builder](construct/mrna/run.py) ·
[sequence manifest](construct/mrna/construct.json) ·
[nucleotide and GenBank checks](construct/mrna/verification.json).

### Verify or replay locally

From the repository root, with Python 3.11 or 3.12:

```sh
python -m pip install -e '.[mrna]'
yash-vaccine verify
yash-vaccine reproduce --output work/design-2-replay
yash-vaccine compare designs/initial designs/design-2
```

`verify` checks recorded evidence and outputs. `reproduce` reselects regions,
re-solves assembly using recorded edge costs, recounts chosen junctions and
reruns RNA encoding. The recorded sequence outputs must match byte for byte.
These commands start no discovery, model inference, downloads or cloud jobs.

For new raw-read processing or inference, use the [pipeline guide](pipeline/README.md)
and [upstream guide](../../docs/upstream.md). [Native receipts](execution/stage-receipts.json),
[worker exits](execution/native-execution.json) and [artifact identities](execution/artifacts.json)
distinguish executed inference from reused evidence and offline replay. The
[data catalog](../../data/README.md#sequencing-data) provides the final normal
alignment, scoped normal genotypes, paired calls and QC. Original normal FASTQs
remain on the persistent Modal volume and are excluded from public downloads.

This design has computational sequence, binding and compatibility checks;
it includes no experimental tumour-cell presentation or T-cell-response result.
