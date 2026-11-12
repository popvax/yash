# Design #1 — detailed methods

Design #1 (case TN25-112103) used Yash’s tumour DNA, tumour RNA and supplied gene-expression values,
without matched-normal DNA. It contains **29 distinct 25-amino-acid regions in an
852-aa protein**: 17 mutation, 5 splice, 5 cancer–testis antigen (CTA)/self and
2 CDH1::NRG1 fusion regions. [Design #2](../design-2/methods.md) is the matched-normal
lead; this design is retained for comparison.

## Workflow at a glance

| Stage | Recorded result |
|---|---|
| Discover tumour DNA variants | 888,490 normalized ALT alleles; 114,386 native PASS |
| Apply protein-change annotation and DNA policy | 1,060 protein-changing variants; 105 primary DNA candidates |
| Reconstruct tumour-RNA proteins | 23 usable primary contexts; 18 qualifying mutation regions |
| Integrate RNA-source regions | 223 entries: 18 mutation, 42 splice events, 161 CTA/self genes and 2 fusion regions |
| Screen hydropathy and select | All 17 compatible mutation regions, 5 splice genes, 5 CTA genes and 2 fusion regions |
| Assemble and encode | 852-aa protein; checked RNA/DNA segments with PolyA unspecified |

Follow [DNA discovery](#1-discover-and-filter-tumour-dna-variants),
[RNA reconstruction](#2-reconstruct-rna-supported-mutant-protein-contexts),
[LENS sources](#3-add-the-rna-only-lens-sources),
[binding scores](#4-score-and-choose-25-aa-regions),
[expression](#5-apply-the-expression-multiplier),
[selection](#6-screen-construct-compatibility-and-select-29-regions),
[assembly](#7-add-patient-supported-routing-sequences-and-optimize-all-joins),
[RNA encoding](#8-encode-and-annotate-the-rnadna-segments) and
[reproduction](#verification-and-upstream-execution).

A **region** is a 25-aa vaccine block; a **peptide** here is an 8–11-aa sequence
scored against an HLA-I allele. One region can supply several peptides. A
**fragment** groups alignments by read name, including paired mates; it is not
a UMI-confirmed original molecule. **ND** denotes no exact junction entry in the
searched normal catalogs, as defined below.

The [data guide](../../data/README.md) provides sequencing access; [inputs/](inputs/README.md)
organizes this design’s evidence. [Usage](../../docs/usage.md) covers alternative analyses.

## 1. Discover and filter tumour DNA variants

### Independent discovery and native filtering

GATK **Mutect2 4.6.2.0** discovered SNVs and indels independently from the tumour
BAM, without Caris VCF nomination or an `--alleles` list. Calling used compatible
GRCh38 coordinates across 24 primary nuclear chromosomes in **143 intervals**.
The vendor capture BED and exact custom-masked alignment reference were
unavailable; calling therefore used chromosome-wide intervals. No matched-normal
DNA or assay-matched panel of normals was available.

The run emitted **879,912 raw records / 888,490 normalized ALT alleles**. Its
calling germline resource used gnomAD v4.1 genomes. Before final filtering,
every emitted allele was checked against QC-PASS gnomAD v4.1 genomes and exomes.
The population prior became the maximum of the native prior and available global
cohort frequencies. Cohort AC/AN values were not summed. Unknown or QC-filtered
observations remained unknown. This updated priors in emitted calls; it did not
repeat local assembly with a different resource.

Native **FilterMutectCalls** combined sample-wide calling statistics, the learned
read-orientation model, contamination estimates and tumour allele-fraction
segments from GetPileupSummaries/CalculateContamination. Segmentation adjusts the
expected inherited heterozygous fraction for local imbalance; it is not a
matched-normal comparison or definitive copy-number call.

We chose `FALSE_DISCOVERY_RATE` with a **0.05 model target**, rather than the
`OPTIMAL_F_SCORE` default. The target is not an experimentally measured error
rate. The segmented result contained **114,386 native-PASS alleles**, including
**1,060 protein-changing variants**.

### Additional DNA eligibility

| Rule | Requirement |
|---|---|
| Native site and allele filters | PASS |
| DNA depth / alternate AD | ≥20 / ≥5 |
| Model variant fraction | ≥5% |
| Germline quality (`GERMQ`) | ≥20 |
| Population frequency | Exclude if either QC-PASS gnomAD cohort reports global or recorded maximum-frequency ancestry-group AF≥0.1%; no group-size minimum |

Mapping/base quality was handled by native GATK, including its long-indel
mapping-quality exception. No extra median MQ/BQ gate was added. Ancestry screening
used recorded `AF_grpmax` regardless of group size.

Sex/ploidy ambiguity, low complexity or segmental duplication, high allele
fraction and insufficient population evidence triggered review rather than
primary admission. Overlapping Caris `FP`, `lowQ` or `lowDP` flags also triggered
review. Caris supplied annotations, not nomination. An absent population record
was distinguished from a lookup failure; a technical lookup failure stopped
processing.

These research rules left **105 primary DNA candidates**. Native PASS and these
extra gates do not prove somatic status; rare inherited variants can remain
without matched-normal DNA.

Evidence: [policy](inputs/mutations/policy.json) ·
[eligibility decisions](inputs/mutations/eligibility/decisions.tsv) ·
[policy verification](inputs/mutations/eligibility/verification.json) ·
[source measurements](inputs/mutations/annotate/variant-evidence.tsv) ·
[segmented commands](inputs/mutations/filter/segmented/commands.json) ·
[filtering statistics](inputs/mutations/filter/segmented/filtering-stats.tsv) ·
[pipeline](../../workflow/stages/01_mutation_calling/discovery/README.md).

## 2. Reconstruct RNA-supported mutant protein contexts

**105 primary DNA candidates → 23 usable RNA-supported contexts → 18 qualifying
mutation regions → 17 retained after hydropathy screening.**

**Isovar 1.57.1** examined the tumour RNA BAM at each nominated allele, assembled
mutant-bearing RNA, found compatible transcript frames and translated local
protein contexts. We used the primary usable context after native filters.
Alternative reconstructions and unsuccessful outcomes remain in the evidence.
RNA review also included DNA review-tier variants; only primary segmented-DNA
candidates could enter mutation selection.

Eligible alignments had mapping quality ≥20, excluded duplicate, secondary and
supplementary records, and required NH=1 when present. Mate names were unioned
as fragments. No additional three-fragment minimum was imposed; supporting reads
still had to allow native Isovar reconstruction containing the alteration.

The ALT fraction was `ALT fragment names / union(REF, ALT, OTHER fragment names)`.
Using a union avoids double-counting overlapping categories. Counts supplied
sequence support and local fractions; gene abundance came from the separate
expression file. Missing reconstruction could reflect coverage or assay limits,
so it was not called biological absence.

Evidence: [RNA outcomes](inputs/mutations/rna/rna-evidence.tsv) ·
[protein contexts](inputs/mutations/rna/primary-proteins.tsv) ·
[RNA worker](../../workflow/stages/01_mutation_calling/discovery/review/rna.py) ·
[evidence guide](inputs/mutations/README.md).

## 3. Add the RNA-only LENS sources

LENS **1.9.2** and MSNAF **1.2.1** used original-FASTQ discovery checkpoints:
splice counts, STAR-Fusion coding effects and Salmon transcript identifiers.
FusionInspector realigned the original FASTQs. Native read review used the
**vendor RNA BAM**, because the original LENS STAR BAM was not preserved.

### CTA/self proteins

The curated list contained **279 genes**: 206 with positive supplied expression,
37 with zero and 36 missing. Positive genes entered without an additional
abundance-percentile cutoff. Their 588 coding transcripts produced 568
methionine-starting reference proteins. Native review, reference filtering and
the HLA seed gate retained **404 isoform contexts from 161 genes**.

Exact matches to non-CTA reference proteins were excluded; matches within the
curated CTA family were intentionally allowed. These are reference self-protein
hypotheses. Gene-level expression does not establish the translated isoform.

Normal-tissue names and medullary thymic epithelial-cell (mTEC) mean/median
annotations were retained but did not affect scores. Normal-tissue names were
missing for 47 of 161 genes. An unavailable mTEC p95 was not substituted with a
mean or represented as a completed p95 filter. Among selected genes, SPATA31A6
had no tissue name and TEKT5 had the literal label `Actin`, which is not an organ.
These annotations did not establish healthy-tissue absence.

### Splice discovery and patient-RNA reconstruction

MSNAF’s initial `maxmin` screen required normal mean <3 and tumour-count excess
>20 against each supplied GTEx/ENCODE and TCGA control. We started with **all
4,929 UIDs** passing both decisions, before translation or HLA filtering. The
source workflow filled missing control entries with zero. That discovery rule
was followed by the stricter coordinate-aware admission screen below.

Patient-window reconstruction required all of the following:

1. **Resolved coordinates and observed bases.** Recover coordinates with the
   native MSNAF functions and cross-check 3,656 recorded translated events.
   Use primary, nonduplicate QC-pass alignments with MAPQ≥20 and NH=1 when
   present. Retain query insertions and separate supported haplotypes; exclude
   soft clips and reference-filled gaps. Do not build a consensus across names.
2. **Consistent mates.** Overlapping mates could extend one fragment only when
   at least 10 aligned overlap bases established a single placement. High-quality
   overlapping bases had to agree.
3. **An observed coding path.** Verify each GENCODE 37 reference CDS against its
   reference protein. Anchor in at least 10 upstream annotated CDS bases, then
   propagate the frame through actual query insertions/deletions. All complete
   codons from anchor to join had to be high-quality and stop-free. Matching a
   modulo-three grouping after an observed stop did not qualify.
4. **A complete window.** Enumerate every 25-aa window crossing the join. All
   75 bases needed quality ≥20 in one fragment or consistently overlapping mates.
   At least **two distinct fragment names** had to support the identical full
   nucleotide haplotype. Each genomic path retained its supporting names and
   stop-free upstream-anchor witnesses.
5. **Native review and reference exclusion.** Apply native event metrics, flags,
   penalties, PASS threshold, coding-kmer and cysteine rules to patient sequences.
   Direct full-window evidence replaced reference-only sequence consistency.
   Nonlinear/back-splice and unresolved-coordinate outcomes stayed in the audit;
   ordinary linear CIGAR gaps were not evidence of circular products.

Join-overlapping 8–11-aa peptides were compared with the complete native Vaxrank
Ensembl 103 inventory: **89,878 distinct protein sequences / 104,806 source
records**, including alternative isoforms. GTF/FASTA associations were revalidated;
the full multi-pattern scan was checked against exhaustive reference k-mers.
Exact reference matches could not contribute to target scores.

Fresh six-allele affinities were checked against OpenVax’s `MHCflurry_Affinity`
adapter. Peptides had to cross the join or contain its split codon; ordinary
flanks could not score. An independent BAM pass verified every window, fragment
and upstream path through pysam query/reference-position arrays.

### Exact-coordinate normal screening

We searched exact GRCh38 exon-boundary endpoints across all resolvable GTEx and
TCGA-adjacent-normal UIDs, regardless of gene label or strand.

| Catalog | Samples / source categories | Resolved rows / total |
|---|---|---|
| GTEx | 2,629 / 51 | 2,463,061 / 2,476,734 |
| TCGA-adjacent normal | 705 / 21 | 5,379,440 / 5,403,219 |

Same-coordinate, same-strand aliases used the per-sample maximum, not a sum.
Positive counts, tissue means and maxima remained diagnostic outputs.

Final admission required **`NOT_DETECTED_IN_SEARCHED_REFERENCE` (ND)**: no exact
coordinate entry in either catalog. Any matching entry was held for review,
including an opposite-strand assignment or an all-zero row. Unresolved coordinates
could not qualify as ND. Missing entries retained null means, not measured zeros;
a low tissue mean could not rescue a detected event. The separate admission
policy was hashed without changing measurement receipts.

Unresolved catalog rows remain outside this searchable universe. ND does not
establish absence from every healthy tissue, exclude another transcript encoding
the same peptide or constitute normal-tissue immunopeptidomics.

### Library orientation and native review

Orientation was inferred from paired BAM flags at annotated CDS junctions,
independently of motifs and XS. Read 1 was antisense in **1,252/1,275 fragments
across 200 genes (98.20%)**. The orientation baseline required ≥100 fragments,
≥20 genes and a dominant fraction ≥90%.

A strict majority of each full-window witness set had to agree with the nominated
gene under this convention. Ties, conflicts and unresolved evidence were held
for review; minority opposition was reported. Motif and XS disagreement stayed
report-only. No distinct-alignment-pattern minimum or extra substitution
exclusion was added.

Native review covered support, mapping/duplicate flags, clipping, bias, anchors,
motifs and coding consistency where implemented. Missing metrics stayed unknown;
repeat/mappability tracks and matched-normal evidence were unavailable. No splice
rows entered the stricter native DNA-linked high-confidence tier. RNA-only
admission used the policy described here.

Evidence: [patient-RNA reconstruction and inference](../../data/shared-rna/splicing/README.md) ·
[normal-screen evidence](../../data/shared-rna/splicing/normal-screen/README.md) ·
[ND policy](../../data/shared-rna/splicing/normal-screen/selection-policy.json) ·
[screen implementation](../../src/yash_vaccine/splice_checks.py) ·
[execution guide](../../workflow/stages/04_region_integration/patient_splice/CHECKS.md).

### Fusion

Eight coding contexts were recovered from STAR-Fusion. FusionInspector support
was matched by **exact breakpoint pair**, not gene names alone. Native RNA PASS,
reference checks and normal-fusion annotation filtering left one primary-eligible
event, **CDH1 → NRG1**, with a predicted frameshift context. Altered-tail peptides
could score without crossing the physical join; a full-length fusion protein was
not experimentally demonstrated.

We also explicitly included the in-frame driver junction at
**chr16:68819425:+ → chr8:32727949:+**. Its call reported 200 junction reads and
one spanning fragment. The 25-aa window was centred on the junction codon.
Model-file-verified inference evaluated 38 distinct peptides × 6 alleles
(**228 pairs**), excluding eight reference-matching peptides. The best retained
peptide, **ATDNATSTST**, had predicted affinity **508.95 nM**, rank **0.7111%**, for
HLA-A*01:01.

This was an explicit driver-junction exception to the general 500-nM cutoff.
It received no invented expression score or fraction borrowed from the other
breakpoint. Both contexts remain translation/presentation hypotheses.

Evidence: [fusion evidence and admission](../../data/shared-rna/fusions/README.md).

### Executed LENS workflow and adaptations

Supplied expression replaced capture-derived CTA abundance gates, not read QC.
Current ranking did not use native LENS linear-affinity or splice-percentage
weights. Coordinate contracts, gene-ID normalization and exact-breakpoint
corrections remain explicit patches/tests.

Code: [stage runner](../../workflow/stages/02_lens/canonical_run.sh) ·
[expression adapter](../../workflow/stages/02_lens/prepare_expression.py) ·
[CTA review](../../workflow/stages/02_lens/run_cta_review.sh) ·
[coordinate contract](../../workflow/stages/02_lens/test_splice_contract.py) ·
[exact fusion review](../../workflow/stages/02_lens/run_fusion_review_exact.sh) ·
[eligibility audit](../../workflow/stages/02_lens/finalize_candidates.py) ·
[native verification](../../workflow/stages/02_lens/verify_final.py).

## 4. Score and choose 25-aa regions

All sources used **MHCflurry 2.3.14**, affinity-model bundle **2.3.0**, for
**A*01:01, A*33:03, B*08:01, B*15:08, C*01:02 and C*07:02**. Class II, processing,
presentation and immunogenicity models did not contribute.

Vaxrank **3.46.0** converted target-overlapping 8–11-aa peptide/HLA affinities into
contributions and summed them within each window. Its normalized logistic curve
used midpoint **350 nM**, width **150 nM**, and zero contributions at affinity
≥5,000 nM or contribution <0.00001. No mutant/reference affinity-ratio gate was
applied. Weaker peptides could contribute inside a window admitted by a stronger
peptide.

Ordinarily a window needed a source-eligible target with predicted affinity
**≤500 nM**, applying the relevant reference policy. The in-frame fusion’s
508.95-nM exception did not change that gate for other candidates.

| Source | Window nomination and recorded pool |
|---|---|
| Mutation | Score every available alteration-overlapping 25-aa window after exact reference-peptide exclusion. Choose the maximum binding sum, breaking ties by native centering order. Require primary DNA eligibility, unambiguous gene, positive expression and a target ≤500 nM: **18 regions in 18 genes**. |
| CTA/self | Highest affinity-sum eligible window per gene: **161 regions**. |
| Splice | **39,269 supported windows / 1,564 events** before HLA filtering; native review, exact-self exclusion and binding admission leave **3,804 / 182 events / 80 genes**. ND and strand admission leave **752 / 42 events / 22 genes**. |
| Fusion | One qualifying frameshift-context window plus the explicitly included in-frame driver window: **2 regions**. |

The compact pool therefore held **223 entries**. Each splice entry represented
an event; final selection still chose one per gene. [Mutation replay](../../src/yash_vaccine/mutations.py),
[splice replay](../../src/yash_vaccine/splice_evidence.py) and the LENS contribution
recount checked native scores against saved affinities and reference exclusions.

## 5. Apply the expression multiplier

```text
Weighted score = native affinity contribution sum × sqrt(expression proxy)
```

| Source | Expression proxy |
|---|---|
| Mutation | Supplied gene value × ALT-supporting RNA fragment fraction |
| Splice | Supplied assigned primary-gene value; no splice-fraction multiplier |
| CTA/self | Supplied gene value |
| Frameshift fusion | Left-partner gene value × bounded exact-breakpoint fusion fraction |
| In-frame fusion | Explicit inclusion; no expression-weighted rank |

The source labels its values TPM and records
`WTSPercentileTransformationVersion=1.0.0`. We preserved them without asserting
calibration as conventional TPM. Raw RNA counts supplied sequence/QC evidence
and local fractions, not the gene-abundance multiplier. Both fusion partners’
measurements were retained; the documented proxy used the left partner.

All 22 splice genes had positive supplied expression. Missing expression remained
missing and could not enter a weighted quota. Expression-weighted scores ranked
splice and CTA candidates and described mutation comparisons; **there was no
top-N mutation cutoff**. Mutation windows were nominated by binding sum before
this multiplier was applied.

## 6. Screen construct compatibility and select 29 regions

### Hydropathy admission

All **223 source-nominated regions** were screened in complete **19-aa
Kyte–Doolittle windows**. Any mean **strictly >1.6** excluded the region. No
padding or rounding preceded the decision; nonstandard residues or short inputs
stopped the helper. This membrane-like-segment heuristic uses the
[Kyte–Doolittle scale](https://doi.org/10.1016/0022-2836(82)90515-0)
([ExPASy values](https://web.expasy.org/protscale/pscale/Hydropath.Doolittle.html));
it is not a calibrated topology predictor or validated vaccine threshold.

**PTCH1 V1083M was the sole exclusion**, maximum KD19 2.121; VPS11 was next at
0.689. All 25 recorded PTCH1 windows also failed, with the least hydrophobic
maximum 1.747. The screen did not substitute a different nominated window.
The underlying DNA/RNA evidence and 223-row pool were retained.

Evidence: [complete region screen](selection/hydrophobicity-screen.tsv) ·
[independent calculation](review/hydrophobicity/verification.json) ·
[screen scope](review/hydrophobicity/README.md).

### Source counts and selected sequences

We included **all 17 compatible mutation regions**, the top **five CTA genes**,
the best event/window from each of the top **five splice genes**, and **both fusion
regions**. Splice ranking considered every qualifying window, then the best per
event and gene. No sequence-identity filter was applied. The source counts were
design choices, not experimentally established optimal counts.

| Splice gene | Selected 25-aa sequence | Full-window fragments | Weighted score | Agreeing/opposing strand names |
|---|---|---:|---:|---|
| F5 | `ALGQMPIFPDLSHTTLSLDFSQTNL` | 2 | 272.94 | 2/0 |
| MUC5B | `LTTTATTPAATISTVTPSSALGTTH` | 2 | 86.20 | 2/0 |
| MUC6 | `PSSFSMAKTSTSLPYHTSSTHHPEV` | 3 | 68.73 | 2/1 |
| MUC3A | `TSQPTTITPSSVGFSGSLPMMTDLT` | 3 | 45.34 | 3/0 |
| HLA-B | `ITQRKWEAAREAEQRSLPGGHVRGV` | 2 | 44.19 | 2/0 |

Every selected splice event was ND in both catalogs. These fragment counts do
not demonstrate independent original molecules or rule out another normal
transcript producing the same peptide. Subsequent sequence interpretation is
explained in the lead design’s [splice sequence checks](../design-2/methods.md#additional-rna-changes-in-the-splice-sequences).

The CTA genes were **PRAME, SPATA31A6, LRRC36, YBX2 and TEKT5**. Mutation genes,
in descriptive score order, were **KLF5, SF1, KMT2D, RXRA, HSPG2, PARP3, CAST,
CCNI, PPIC, PPP2R1A, KCTD19, VPS11, ARHGAP40, NRXN2, PLEC, DOCK8 and SYT12**.
All selected regions were distinct. `region_id` is stable across selections;
AG01-style labels are local to a selection.

The independent BAM pass covered **39,269 windows and 617,268 named
fragment/window observations**, including stop-free upstream paths. Offline
replay recounted native PASS, reference exclusion, admission and rank.

Evidence: [selected table](selection/selected.tsv) · [FASTA](selection/selected.fasta) ·
[splice gene ranking](selection/splice-gene-ranking.tsv) ·
[offline splice audit](selection/splice-rna-verification.json) ·
[full patient-RNA evidence](../../data/shared-rna/splicing/README.md).

## 7. Add patient-supported routing sequences and optimize all joins

### Routing elements

| Element | Patient sequence | Role |
|---|---|---|
| SEC | HLA-C*01:02 residues 1–26 | Annotated 24-aa signal peptide plus two native residues |
| MITD | HLA-B*15:08 residues 308–362 | Membrane anchor and trafficking tail |

Tumour DNA/RNA typing supported both alleles. The read audit checked every
native base and RNA-supported codon. Cleavage and trafficking were not measured
in this construct. [Architecture](construct/architecture.json) and
[allele evidence](construct/allele-evidence/README.md) pin the sequences.

### Search space and predictions

We fixed the 29 antigen windows and jointly selected their order and one of
[11 linker choices](../../config/linkers.json) at every join, including
SEC→first antigen and last antigen→MITD: **8,932 internal options + 638 terminal
options = 9,570 options across 30 selected joins**.

Internal options were rejected if any complete 19-aa KD mean in the
antigen–linker–antigen sequence exceeded 1.6. All options among the retained
29 regions passed; the rule remained enforced. Every region also had to pass
alone. The final full-chain screen was unmasked.

Peptides touching a linker were enumerated at lengths 8–11 aa, including
linker-only peptides; direct joins had to cross the boundary. SEC supplied the
initiating methionine. Model files, packages, native Vaxrank enumeration and
OpenVax/mhctools adapter agreement were checked. Only affinity was used.

| Prediction set | Unique peptides | Six-allele pairs |
|---|---:|---:|
| Full verified internal-cassette cache | 135,344 | 812,064 |
| Queries required by the internal-cassette helper | 127,115 | 762,690 |
| Fresh terminal predictions | 14,040 | 84,240 |
| Full routed-assembly query union | 134,034 | 804,204 |

The internal-helper queries were a checked subset of the cache. Terminal
inference was separate; the routed union includes internal and terminal options.
The sets overlap, so their counts cannot simply be added. These are all-option
queries, not just the junctions selected in the final protein.

### Joint objective and final check

OR-Tools CP-SAT **9.14.6206**, one worker and seed 31, solved these objectives in
order while preserving each earlier optimum:

1. Minimize position-specific junction peptide–HLA rows at rank ≤0.5%.
2. Minimize rows at rank ≤2%.
3. Maximize the lowest junction percentile rank.
4. Minimize inserted linker residues.

All four phases were proven optimal within the fixed antigen/linker/model space.
The protein contained **26 SEC + 725 antigen + 46 linker + 55 MITD = 852 aa**.
A separate fresh run checked all **30 joins: 1,204 unique peptides × 6 alleles =
7,224 pairs**. A whole-protein scan confirmed positions and threshold agreement.

| Selected-junction result | Peptide–HLA rows |
|---|---:|
| Affinity rank ≤0.5% | 0 |
| Affinity rank ≤2% | 33 |
| Affinity ≤500 nM | 38 |

The thresholds overlap. These are binding rows, not distinct antigens or measured
presentation. The objective did not optimize cleavage, intended-target processing
or immunogenicity.

The independent screen found no high windows outside the fixed routing elements:
**753 windows, maximum KD19 1.000**. Ten high windows wholly inside MITD and six
crossing its boundary were separately reported, without topology validation.

Evidence: [assembly code](construct/run.py) · [solver](construct/optimize.py) ·
[protein](construct/protein/cassette.fasta) · [components](construct/protein/components.tsv) ·
[junctions](construct/protein/junctions.tsv) ·
[solver/prediction record](construct/protein/construct.json) ·
[independent junction check](construct/final-check.json).

## 8. Encode and annotate the RNA/DNA segments

### PopRNA settings

[PopRNA](construct/mrna/optimizer.py) chose synonymous codons for the fixed
852-aa protein using the corrected supplied heuristic.

| Setting | Value |
|---|---|
| Human codon table | `h_sapiens_9606` |
| Codon-frequency cutoff / seed | 0.10 / 0 |
| Folding engine | ViennaRNA 2.7.2, 37 °C |
| Folding passes | Up to 100 iterations per pass |

Sequential steps addressed motif repair, first-codon selection,
alternative-frame-stop placement, paired-bracket-run reduction, codon-distribution
balancing and repeat removal. [backbone.json](construct/mrna/backbone.json)
records exact settings, source hash and implementation changes.

The builder checked exact translation and every protein component. It retained
the 5′ UTR/Kozak, **TAA–TAG–TGA** stops, 3′ UTR and upstream template context.
Fixed UTRs were outside the coding motif constraint. PolyA length, sequence and
architecture, downstream extensions, cap chemistry and modifications remained
unspecified.

### Exported segments and sequence checks

| Result | Value |
|---|---:|
| Coding DNA including first stop | 2,559 nt |
| Known RNA segment | 2,699 nt |
| Known upstream DNA-context segment | 4,640 bp |
| Coding GC / standard CAI | 56.2989% / 0.77396802 |
| Forbidden motifs / repeated 12-nt pairs | 0 / 0 |
| Longest homopolymer | 7 nt |
| Final paired-bracket run / heuristic target | 28 nt / 14 nt |

Exports end after the 3′ UTR without an invented tail or complete molecule
length. Coding motifs and below-cutoff codons were excluded; the known template
had no XbaI site. Later repairs did not retain the folding heuristic’s 14-nt
target. Bracket-run length is not a continuous-stem or ensemble-accessibility
objective, and these checks do not establish improved expression or stability.

[GenBank code](construct/mrna/genbank.py) derives features from selected regions,
routing-element sources and pinned backbone references. Antigens retain gene,
category, event and region identities. Upstream context includes a
reference-matching Tn903 aphA1 CDS and pUC/pMB1-like origin with its documented
insertion. Unknown intervals and incomplete topology are labelled. Round-tripping
checks sequence, coordinates, features and translations; PolyA has no assigned
sequence or length.

Evidence: [sequence manifest](construct/mrna/construct.json) ·
[nucleotide checks](construct/mrna/verification.json) ·
[nucleotide builder](construct/mrna/run.py).

## Additional mutation checks

### Filter-model sensitivity

The selected global priors agreed with retained gnomAD records; some rare
population alleles remained below the 0.1% exclusion threshold. Two native
FilterMutectCalls diagnostics used either segmentation with higher observed
population-group priors, or the original priors without optional segmentation.
The segmented model remained the selection policy.

KLF5, SF1, RXRA, PPIC, PPP2R1A, KCTD19, VPS11, ARHGAP40, PLEC, DOCK8 and SYT12
passed research gates in both diagnostics. KMT2D, HSPG2, PARP3, CCNI, PTCH1 and
NRXN2 were more segmentation-dependent. CAST fell below the extra GERMQ gate
with higher priors while retaining segmentation. These were sensitivity checks,
not inherited-status assignments. [Commands and measurements](inputs/mutations/sensitivity/README.md).

### Selected-mutation read review

All 18 pre-hydropathy mutation candidates received a fresh read-only DNA/RNA base
check; 17 entered the construct. Each alternate was seen in both datasets and
both read orientations under the supplementary quality policy. Paired orientations
do not establish independent molecules. KCTD19, NRXN2 and SYT12 had two alternate
RNA fragment names each, also supporting their longer Isovar contexts. No extra
minimum was added.

Every region was compared with compatible annotated proteins; PPIC retained its
documented background difference. This review added evidence, not a new selection
filter or somatic-status assignment. See [code, counts and receipts](review/mutation-read-check/README.md)
and [individual assessments](antigen-review.md).

## Verification and upstream execution

From the repository root, with Python 3.11 or 3.12:

```sh
python -m pip install -e '.[mrna]'
yash-vaccine reproduce designs/initial --output work/design-1-replay
```

This rebuilds selection, the internal-cassette checkpoint, full routed assembly,
recorded junction verification and nucleotide encoding. It compares results with
[candidate.json](candidate.json). Selection uses the maintained implementation;
protein assembly runs verified recorded source in a temporary input tree so
historical code hashes remain accurate. Outputs must match the recorded design.
No patient reads are downloaded and no new model inference or cloud job starts.

Offline replay recalculates eligibility, scores and selection from recorded
measurements; it does not call variants, query new population data or reread
BAMs. New DNA/RNA analysis and prediction require separately staged inputs and
the documented runtime. [Assembly execution](construct/README.md),
[RNA encoding](construct/mrna/README.md) and [upstream processing](../../docs/upstream.md)
provide commands. No experimental tumour-cell presentation or T-cell-response
result is included.

This is one design approach. The shared data also support MHC class II,
immunogenicity and presentation-aware linker studies.
