# Methods for the initial vaccine design

This initial design for case TN25-112103 combines **29 distinct 25-aa regions**:
17 mutation regions, one region from each of five splice genes, five CTA/self
genes and two CDH1::NRG1 fusion regions. Patient-supported SEC/MITD routing sequences
and jointly selected linkers form one **852-aa protein**. PolyA and any downstream
extension are unspecified; the nucleotide files contain the known upstream segment.

This is one approach to vaccine design for Yash. The shared data support other
antigen selections, MHC class II analysis, immunogenicity prediction and linker
optimization for processing and presentation. Commands are in [usage](../../docs/usage.md),
and sequencing access is described in [data](../../data/README.md).

## 1. Discover and filter tumour DNA variants

We used GATK **Mutect2 4.6.2.0** to discover SNVs and indels independently from
the tumour DNA BAM, without Caris VCF nomination or an `--alleles` list. The
reference is compatible GRCh38. Calling spans 24 primary nuclear chromosomes
in **143 parallel intervals**. The vendor capture-region BED and exact masked
alignment reference were unavailable, so calling used chromosome-wide intervals.
No matched-normal DNA or assay-matched panel of normals was available.

The discovery set contains **879,912 raw records / 888,490 normalized ALT alleles**.
The calling resource uses gnomAD v4.1 genomes. Every emitted allele is also checked
against QC-passing gnomAD v4.1 genomes and exomes before final filtering. The
population prior is the maximum of the native prior and those cohorts' global
allele frequencies. Cohort AC/AN values are not summed. Unknown or QC-filtered
observations are not treated as measured zero frequencies. This updates priors
in the emitted calls before filtering; it does not repeat local assembly with a
new calling resource.

Native **FilterMutectCalls** uses merged sample-wide calling statistics, the learned
read-orientation model, contamination estimates and the **tumour-segmentation**
file from GetPileupSummaries/CalculateContamination. Segmentation adjusts the
expected inherited heterozygous allele fraction for local imbalance; it is an
estimate from tumour data, not matched-normal confirmation or a definitive copy-number call.
The selected threshold strategy is `FALSE_DISCOVERY_RATE` with a **0.05 model target**.
This is an explicit choice rather than the `OPTIMAL_F_SCORE` default, and it
is not an experimentally established 5% error rate for the candidate list.

The segmented callset has **114,386 native-PASS alleles**, including **1,060
protein-changing variants**. Primary research eligibility additionally requires:

| Rule | Requirement |
|---|---|
| Native site and allele filters | PASS |
| DNA depth / alternate AD | At least 20 / at least 5 |
| Model variant fraction | At least 5% |
| Germline quality GERMQ | At least 20 |
| QC-PASS population AF | Below 0.1% globally and for the recorded maximum-frequency ancestry group, in either gnomAD cohort; no minimum group allele count |

Mapping/base quality is handled by native GATK filtering, including its long-indel
mapping-quality exception; no additional median MQ/BQ thresholds are applied.
Ancestry-group frequency screening uses the recorded `AF_grpmax` regardless of
that group's sample size. Fixed depth, support,
AF, GERMQ and population exclusions are additional research choices, without a
validated sensitivity estimate. Sex-chromosome, low-complexity/segmental-duplication,
high-AF/ploidy-ambiguous and insufficient-population cases require review rather
than entering the primary pool. An overlapping Caris `FP`, `lowQ` or `lowDP` flag
also triggers review; Caris is used for annotations, not nomination. Missing
population lookups are distinguished from technical lookup failure, which stops
processing. These rules leave **105 primary DNA candidates**.

Neither native PASS nor these extra gates confirms somatic status. Rare inherited
variants remain possible without matched-normal DNA. The full rules are in the
[research policy](../../data/mutation-evidence/policy.json), with
[DNA eligibility](../../data/mutation-evidence/eligibility/decisions.tsv),
[policy verification](../../data/mutation-evidence/eligibility/verification.json),
[source allele measurements](../../data/mutation-evidence/annotate/variant-evidence.tsv),
[segmented commands](../../data/mutation-evidence/filter/segmented/commands.json),
[filtering statistics](../../data/mutation-evidence/filter/segmented/filtering-stats.tsv)
and [pipeline source](../../workflow/stages/01_mutation_calling/discovery/README.md).

Offline reproduction calculates DNA eligibility directly from the recorded native
calls and population observations using this policy. It verifies the resulting
decisions before mutation-window scoring, without new variant calling or population
lookups. Source execution records are documented in the
[evidence guide](../../data/mutation-evidence/README.md).

## 2. Reconstruct RNA-supported mutant protein contexts

**105 primary DNA candidates → 23 with usable RNA-supported protein contexts →
18 qualifying mutation regions before construct-compatibility screening → 17 retained.**

Isovar **1.57.1** examines the tumour RNA BAM at each nominated genomic allele,
assembles mutant-bearing sequence, finds compatible transcript reading frames
and translates local protein contexts. We use the primary usable context after
native filters; alternate reconstructions and unsuccessful outcomes remain in
the recorded RNA evidence. RNA review also includes DNA review-tier variants,
but only segmented primary-DNA candidates can enter the mutation selection.

Eligible alignment evidence uses mapping quality ≥20, excludes duplicate,
secondary and supplementary alignments, and requires NH=1 where NH is present.
Mates share one fragment identifier. The alternate fraction divides alternate
fragment IDs by the union of reference, alternate and other fragment IDs; it does
not sum potentially overlapping categories or claim UMI molecule counting.
No additional three-fragment minimum is imposed. Supporting fragments must
allow native Isovar protein reconstruction containing the alteration.

A lack of usable reconstruction can reflect local coverage or assay limitations,
not biological absence. Read counts provide sequence support and local allele
fractions; gene abundance comes from the separately supplied expression file.

[RNA evidence](../../data/mutation-evidence/rna/rna-evidence.tsv),
[protein contexts](../../data/mutation-evidence/rna/primary-proteins.tsv) and
[RNA worker](../../workflow/stages/01_mutation_calling/discovery/review/rna.py)
record the implementation and outcomes.

## Mutation-model sensitivity

The selected global population priors agree with the retained gnomAD records.
Several rare population alleles remain below the policy’s 0.1% exclusion threshold.
Two native GATK FilterMutectCalls diagnostics retain either segmentation with
higher observed population-group priors, or the original priors without optional
segmentation. The segmented model remains the selection policy. KLF5, SF1, RXRA
PPIC, PPP2R1A, KCTD19, VPS11, ARHGAP40, PLEC, DOCK8 and SYT12 pass the research gates in both diagnostics. KMT2D, HSPG2, PARP3, CCNI,
PTCH1 and NRXN2 are more dependent on segmentation; CAST also falls below the extra
GERMQ gate with higher priors while segmentation is retained. These diagnostics
are model sensitivity checks, not inherited-status assignments. [Commands and
measurements](../../data/mutation-evidence/sensitivity/README.md) are retained.

## Supplementary selected-mutation review

All 18 mutation candidates before hydropathy screening received a fresh read-only DNA/RNA base check; 17 enter the construct.
Every alternate is observed in both datasets and both read orientations under the
supplementary quality policy; paired orientations are not independent molecules.
KCTD19, NRXN2 and SYT12 have only two alternate RNA fragment names each, also
supporting their longer Isovar context. No extra minimum was added. Every region
is compared with compatible annotated protein contexts; PPIC retains its documented
background difference. This review adds evidence, not a somatic-status assignment
or a new eligibility filter. [Code, counts and receipts](review/mutation-read-check/README.md)
and [individual assessments](antigen-review.md) retain those limits.

## 3. Add the RNA-only LENS sources

The RNA workflow used LENS **1.9.2**, MSNAF **1.2.1** and original-FASTQ discovery
checkpoints: the splice-count matrix, STAR-Fusion
coding-effect calls and Salmon transcript IDs. FusionInspector realigned the
original RNA FASTQs. Native read review used the **vendor RNA BAM**, because the
original LENS STAR BAM had not been preserved. This distinction is recorded.

### CTA/self proteins

LENS's curated list contains 279 genes. Supplied expression was positive for
206, zero for 37 and missing for 36. We admitted positive genes without a
further abundance-percentile cutoff: 588 coding transcripts produced 568
methionine-starting reference proteins. Native RNA review, reference filtering
and the HLA seed gate left **404 isoform contexts from 161 genes**. Gene-level
expression does not establish which reference isoform is translated.

Exact non-CTA reference matches are excluded. Matches within the curated CTA
family are intentionally allowed: these are self-protein hypotheses. Normal
tissue names and mTEC mean/median annotations remain in the evidence, but do
not enter the current score. Normal-tissue names are missing for 47 of the 161
genes; absent annotations do not establish absence of healthy-tissue expression.
An unavailable mTEC p95 was not replaced by a mean or represented as a completed
p95-based filter.
Among the selected genes, SPATA31A6 has no normal-tissue name, and TEKT5 carries
the literal source label `Actin`, which cannot be interpreted as an organ.
These source annotations require clarification before use as a normal-organ
screen; neither affected this design's score.

### Splice events

MSNAF's `maxmin` screen requires normal mean <3 and tumour-count excess >20
against each supplied GTEx/ENCODE and TCGA control. Start with **all 4,929 UIDs**
passing both decisions, before reference translation or HLA filtering. Missing
control entries are represented as zero by the source workflow; this does not
establish biological absence in healthy tissue.

**Final splice admission requires ND in both searched normal catalogs.** The
initial MSNAF screen defines the discovery pool, not the final normal-reference
decision. The additional coordinate-aware stage described below is mandatory.

Recover coordinates with the original MSNAF functions and cross-check all 3,656
recorded translated events. Reconstruct from actual tumour RNA query bases,
using primary nonduplicate QC-pass alignments, MAPQ >=20 and NH=1 when present.
Do not use soft clips, fill gaps from the reference or combine different
fragments into a consensus. Retain query insertions and separately supported
nucleotide haplotypes. Consistently overlapping mates may extend a fragment;
at least ten aligned overlap bases must establish a single placement.

GENCODE 37 supplies a local upstream frame. First verify each reference CDS
against its matching reference protein. Anchor actual RNA translation in at
least ten upstream annotated CDS bases. Propagate that frame along the query,
including local insertions/deletions, and require all observed complete codons
between the anchor and the join to be high-quality and **stop-free**. A region
after an observed termination codon does not qualify merely by matching the
same modulo-three grouping.

Enumerate all 25-aa windows crossing the join. All 75 nucleotide bases must have
quality >=20 in a junction fragment or its consistently overlapping mate. Require
at least **two distinct fragment names** supporting the same complete nucleotide
haplotype; fragment names are not UMI molecule counts. High-quality overlapping
mates must agree. Record every genomic path, its own supporting fragment names
and its stop-free coding-anchor witnesses.

Apply native LENS event metrics, flags, penalties and PASS threshold, plus its
peptide coding-kmer and cysteine rules on the actual patient coding sequence.
Direct full-window RNA evidence replaces the reference-only sequence-consistency
check. Keep native motif handling, without an added direction exclusion.
Nonlinear/back-splice and unresolved coordinates remain in the outcome audit;
ordinary linear CIGAR gaps do not establish circular products.

Compare join-overlapping 8–11-aa peptides with the complete native Vaxrank
Ensembl 103 reference inventory: **89,878 distinct protein sequences and 104,806
protein-source records**, including annotated alternative isoforms. Revalidate
every association against the current GTF/protein FASTA. Independently compare
the full multi-pattern scan with exhaustive k-mer enumeration of every reference
protein. Exact self matches cannot contribute to target scoring.

Predict fresh six-allele affinities and check every result against the OpenVax
`MHCflurry_Affinity` adapter. Contributing peptides must cross the join or contain
its split codon; ordinary flanks do not score as splice targets. The complete
pool and read/inference receipts are included in
[patient-RNA splice evidence](../../data/splice-evidence/README.md). An independent
BAM pass verifies every reconstructed window and named fragment, including the
observed stop-free upstream path, through pysam's query/reference-position arrays.

#### Exact-coordinate normal screening and read orientation

Search exact GRCh38 exon-boundary endpoints across every resolvable UID in the
supplied GTEx and TCGA-adjacent-normal matrices, regardless of gene assignment
and strand. This covers 2,629 GTEx samples in 51 source tissue categories and
705 TCGA-adjacent-normal samples in 21 categories. Same-coordinate same-strand
aliases use the per-sample maximum rather than being added. Retain positive
sample counts, tissue means and maxima as diagnostics.

Every admitted splice event must be **NOT_DETECTED_IN_SEARCHED_REFERENCE (ND)**:
no exact coordinate entry in either searched catalog. Any matching entry,
including an opposite-strand assignment or a catalog row with all-zero counts,
is held for review. Unresolved coordinates cannot qualify as ND. Missing entries
retain null means rather than measured zeros. A low normal-tissue mean does not
permit a detected event into selection. The separate admission policy is hashed
without altering the original measurement receipts.

Library orientation is inferred from paired BAM flags at annotated CDS junctions,
independently of splice motifs and XS. Read 1 is antisense in 1,252/1,275 sampled
fragments across 200 genes (98.20%). The operational orientation check requires
at least 100 fragments across 20 genes and a dominant fraction of at least 90%.
For each complete-window witness set, a strict majority must agree with the
nominated gene under that library convention; ties, conflicting or unresolved
evidence are held for review. Minority opposition is reported. Motif and XS
disagreement remain report-only. No distinct-alignment-pattern minimum or extra
substitution exclusion is added.

The normal search resolves 2,463,061/2,476,734 GTEx rows and
5,379,440/5,403,219 TCGA-normal rows; unresolved catalog rows remain outside the
searchable universe. ND therefore describes this search, not demonstrated absence
from all healthy tissues or measured absence of an HLA-presented peptide. A
different normal transcript could also produce the same peptide. No experimental
normal-tissue immunopeptidome comparison is included.

[Normal-screen evidence](../../data/splice-evidence/normal-screen/README.md),
[admission policy](../../data/splice-evidence/normal-screen/selection-policy.json),
[implementation](../../src/yash_vaccine/splice_checks.py) and
[execution guide](../../workflow/stages/04_region_integration/patient_splice/CHECKS.md)
record the complete screen and offline replay.

### Fusion

Eight coding contexts were recovered from the original STAR-Fusion calls.
FusionInspector support was matched by **exact breakpoint pair**, not merely
gene names. Native RNA PASS, reference checks and the normal-fusion annotation
filter left one primary-eligible event: **CDH1 → NRG1**. The admitted context is
a predicted frameshift product, not a demonstrated full-length fusion protein.
Relevant altered-tail peptides can score without spanning the physical join.

The design also explicitly includes the in-frame CDH1::NRG1 junction at
chr16:68819425:+ → chr8:32727949:+. Its source call reports 200 junction reads
and one spanning fragment. The 25-aa window is centred on the recorded junction
codon. A model-file-verified screen predicts 38 distinct 8–11-aa peptides across
six alleles (228 pairs); eight reference-matching peptides are excluded. The best
remaining peptide, ATDNATSTST, is predicted at **508.95 nM**, rank **0.7111%**, for
HLA-A*01:01. Its inclusion is an explicit driver-junction exception to the general
500 nM admission gate. It receives no invented expression score and no RNA fraction
borrowed from the other breakpoint. Both fusion products remain hypotheses of
translation/presentation. [Fusion evidence](../../data/fusion-evidence/README.md).

### Native review versus current ranking

Read checks cover support, mapping/duplicate flags, clipping, bias, anchors,
motifs and coding-sequence consistency where implemented. Missing metrics remain
unknown; repeat/mappability tracks and matched-normal evidence were unavailable.
RNA-only splice rows are retained under the stated policy; the stricter native
DNA-linked high-confidence tier contains no splice rows.

The supplied expression adapter replaces capture-derived CTA abundance gates;
it does not replace read-level QC. The exon-boundary coordinate contract,
gene-ID normalization and exact-breakpoint fusion corrections are preserved as
explicit patches/tests. Current ranking does **not** use native LENS's
linear-affinity or splice-percentage weights.

Executed code:
[canonical stage runner](../../workflow/stages/02_lens/canonical_run.sh),
[expression adapter](../../workflow/stages/02_lens/prepare_expression.py),
[CTA review](../../workflow/stages/02_lens/run_cta_review.sh),
[splice coordinate contract](../../workflow/stages/02_lens/test_splice_contract.py),
[exact fusion review](../../workflow/stages/02_lens/run_fusion_review_exact.sh),
[eligibility and audit](../../workflow/stages/02_lens/finalize_candidates.py),
[native verification](../../workflow/stages/02_lens/verify_final.py).

## 4. Score and choose 25-aa regions

All sources use **MHCflurry 2.3.14**, model bundle **2.3.0**, for Yash's six
HLA-I alleles: **A*01:01, A*33:03, B*08:01, B*15:08, C*01:02, C*07:02**.
HLA-II, processing/presentation and immunogenicity models do not contribute.

A region normally must contain a non-reference target peptide with predicted affinity
**≤500 nM**, under the relevant source-specific reference policy. Native
Vaxrank **3.46.0** converts target-overlapping 8–11-aa peptide/HLA affinities to
contributions and sums them within each window. Its default normalized logistic
curve has midpoint 350 nM and width 150 nM; contributions are zero at ≥5,000 nM
or below 0.00001. Thus weaker binders inside an admitted window can contribute.
No mutant-versus-wild-type affinity-ratio gate is applied. The explicitly included
in-frame driver junction is the documented 508.95 nM exception, not a change to
this gate for the rest of the pool.

- **Mutation:** score every available 25-aa window overlapping the altered interval,
  excluding peptides with an exact annotated reference-protein match. Select the
  exact maximum binding-contribution window; native centring order breaks ties.
  Require positive supplied expression, unambiguous gene assignment, primary
  segmented DNA eligibility and at least one retained target ≤500 nM. This leaves
  **18 regions in 18 genes**.
- **CTA:** keep the highest affinity-sum eligible window per CTA gene. There are
  **161 CTA/self regions**.
- **Fusion:** retain the qualifying frameshift-context window and explicitly
  include the in-frame driver-junction window described above: **two regions**.
- **Splice:** reconstruct 39,269 supported RNA windows from 1,564 events before
  HLA filtering. Native review, exact-self exclusion and affinity admission leave
  3,804 windows across 182 events and 80 genes. Required ND-only normal screening
  and the limited strand check leave **752 windows, 42 events and 22 genes**.
  Rank passing windows, then retain the best per event and gene.

The compact pool contains **223 entries**: 18 mutation, 161 CTA/self, 42 splice
and 2 fusion. [Mutation replay](../../src/yash_vaccine/mutations.py) and
[splice replay](../../src/yash_vaccine/splice_evidence.py) independently recount
native window scores from the saved affinities and reference matches.
The LENS contribution recount separately checks its recorded CTA/fusion inputs.

## 5. Apply the expression multiplier

**Weighted score = native affinity contribution sum × √(expression proxy).**

| Source | Expression proxy |
|---|---|
| Mutation | Supplied gene value × alternate-supporting RNA fragment fraction |
| CTA/self | Supplied gene value |
| Splice | Supplied assigned primary-gene value; no splice fraction |
| Frameshift fusion | Left-partner gene value × bounded fusion fraction; both raw partner measurements retained |
| In-frame fusion | Explicit inclusion; no invented expression-weighted score |

The source table labels its values TPM but carries
`WTSPercentileTransformationVersion=1.0.0`; we preserve the numbers without
asserting they are calibrated conventional TPM. Capture-derived read counts
supply sequence/QC evidence rather than the gene-expression multiplier.

At the final splice-gene level, **all 22 genes** have positive supplied expression.
Missing expression remains missing and cannot enter the weighted quota. Splices
use their assigned gene value without a read-fraction multiplier. Both fusion
partner measurements remain in the evidence; scoring uses the documented
left-partner proxy. No pseudocount is added.

## 6. Screen construct compatibility and select 29 regions

Screen every one of the **223 source-nominated regions** using uniform, complete
**19-aa Kyte–Doolittle windows**. Reject a region if any mean is **strictly greater
than 1.6**. No padding or rounding precedes the decision; nonstandard residues or
short inputs stop this helper. This is a conservative membrane-like-segment
heuristic from [Kyte and Doolittle (1982)](https://doi.org/10.1016/0022-2836(82)90515-0),
using the [ExPASy scale](https://web.expasy.org/protscale/pscale/Hydropath.Doolittle.html).
It is not a calibrated topology predictor or a validated vaccine threshold.

**PTCH1 V1083M is the sole exclusion**: maximum KD19 2.121, compared with 0.689
for the next-highest nominated region (VPS11). All 25 recorded PTCH1 mutation-spanning
windows also fail, with the least hydrophobic maximum at 1.747. The source adapters'
window nomination remains unchanged; this gate does not silently choose substitute
windows. DNA/RNA measurements, somatic eligibility and the 223-row pool are retained.
See the [complete screen](selection/hydrophobicity-screen.tsv),
[independent calculation](review/hydrophobicity/verification.json) and
[scope and limitations](review/hydrophobicity/README.md).

Starting from the complete patient-RNA splice pool, rank every qualifying window
using the native affinity sum times the square root of the supplied gene value.
Take the best window per event and the best event per gene. No sequence-identity
filter is applied. The independent BAM check covers all **39,269 windows**
and **617,268 named fragment/window observations**,
including the observed stop-free upstream paths. The offline replay independently
recounts every native-PASS window, reference exclusion, admission decision and rank.

Include **all 17 mutation regions passing the hydropathy screen**, top **five CTA genes**, the **two fusion regions**,
and the highest-ranked region from each of the top **five distinct splice genes**.
The mutation count includes every region passing both upstream eligibility and construct compatibility, rather than a score-rank cutoff. The other source counts remain design choices, not experimentally established optimal counts.

| Splice gene | Selected 25-aa sequence | Full-window fragments | Weighted score |
|---|---|---:|---:|
| F5 | `ALGQMPIFPDLSHTTLSLDFSQTNL` | 2 | 272.94 |
| MUC5B | `LTTTATTPAATISTVTPSSALGTTH` | 2 | 86.20 |
| MUC6 | `PSSFSMAKTSTSLPYHTSSTHHPEV` | 3 | 68.73 |
| MUC3A | `TSQPTTITPSSVGFSGSLPMMTDLT` | 3 | 45.34 |
| HLA-B | `ITQRKWEAAREAEQRSLPGGHVRGV` | 2 | 44.19 |

Every selected splice event is ND in both searched catalogs. Full-window strand
support is respectively 2/0, 2/0, 2/1, 3/0 and 2/0 agreeing/opposing fragment
names. These are not demonstrated independent original molecules. ND does not
establish that the peptide cannot arise through another healthy-tissue transcript.

The CTA genes are **PRAME, SPATA31A6, LRRC36, YBX2 and TEKT5**.
The mutation genes, in score order, are **KLF5, SF1, KMT2D, RXRA, HSPG2,
PARP3, CAST, CCNI, PPIC, PPP2R1A, KCTD19, VPS11, ARHGAP40, NRXN2,
PLEC, DOCK8 and SYT12**. Every selected sequence is a distinct 25-aa region.

`region_id` identifies a sequence across selections. AG01-style labels are local
to a selection. Compare independent designs by stable region IDs.

[Selected regions](selection/selected.tsv), [FASTA](selection/selected.fasta),
[splice gene ranking](selection/splice-gene-ranking.tsv),
[offline splice audit](selection/splice-rna-verification.json) and
[full patient-RNA evidence](../../data/splice-evidence/README.md) are included.

## 7. Add patient-supported routing sequences and optimize all joins

SEC uses Yash's **HLA-C*01:02 residues 1–26**, including its annotated 24-aa
signal peptide plus the next two native residues. MITD uses **HLA-B*15:08
residues 308–362**, including membrane anchor and trafficking tail. Tumour DNA
and RNA typing support both alleles; the read-level audit checks every native
base and RNA-supported codon. Cleavage and trafficking have not been measured
in this construct. [Architecture](construct/architecture.json) and
[allele evidence](construct/allele-evidence/README.md) pin these sequences.

Keep all 29 antigen windows fixed. Jointly choose antigen order and one of
[11 linker choices](../../config/linkers.json) at all **30 joins**, including
SEC→first antigen and last antigen→MITD. The search contains **8,932 directed
internal pair/linker options** and **638 terminal options**, totaling **9,570**.
Before optimization, reject an internal option if any complete 19-aa KD mean in
its antigen–linker–antigen sequence exceeds 1.6. All options among the retained
29 regions pass; the feasibility rule is still enforced. The full assembled chain
receives an independent unmasked screen after optimization. It has zero high
windows outside the fixed routing elements (753 windows, maximum KD19 1.000).
Ten high windows wholly inside MITD and six crossing its boundary are separately
reported for review, not asserted to have validated topology. No antigen can be
exempted by touching a routing element because each must also pass alone.

Enumerate all 8–11-aa peptides touching a linker; direct joins must cross the
boundary. Linker-only peptides are included. SEC supplies the initiating M.
The verified cache contains **135,344 internal-cassette query peptides × six
alleles = 812,064 pairs**. Current helper requirements are **127,115 peptides ×
six = 762,690 pairs**, a checked subset. Fresh terminal inference covers **14,040
peptides × six = 84,240 pairs**. Checks pin model files and package versions,
verify native Vaxrank junction enumeration and confirm OpenVax/mhctools adapter
agreement. Only the affinity submodel is used.

OR-Tools CP-SAT **9.14.6206**, one worker and seed 31, solves four lexicographic
objectives in order:

1. Minimize position-specific junction peptide–HLA hits at rank ≤0.5%.
2. Holding that optimum, minimize hits at rank ≤2%.
3. Holding both counts, maximize the lowest junction percentile rank.
4. Holding all three, minimize inserted linker residues.

All four phases are proven optimal within the fixed antigen/linker/model space.
The resulting protein has **0 hits at rank ≤0.5%, 33 at rank ≤2% and 38 at predicted
affinity ≤500 nM**. These overlapping counts describe peptide–HLA binding rows.
The chain contains **26 SEC + 725 antigen + 46 linker + 55 MITD = 852 aa**.
This objective does not optimize cleavage, intended-target processing,
presentation probability or immunogenicity.

[Assembly code](construct/run.py), [solver](construct/optimize.py),
[protein](construct/protein/cassette.fasta), [coordinates](construct/protein/components.tsv),
[junctions](construct/protein/junctions.tsv) and
[solver/prediction record](construct/protein/construct.json) accompany the design.

## 8. Encode and annotate the RNA/DNA segments

PopRNA encodes the fixed **852-aa protein** with human codon usage
(`h_sapiens_9606`), a **0.10** codon-frequency cutoff and seed **0**. The corrected
supplied heuristic performs synonymous motif repair, first-codon selection,
alternative-frame-stop placement, paired-bracket-run reduction, codon-distribution
balancing and repeat removal. It uses ViennaRNA **2.7.2** minimum-free-energy
folding at **37 °C**, with up to 100 iterations in each folding pass. The exact
settings, source SHA-256 and implementation changes are in
[backbone.json](construct/mrna/backbone.json).

The builder independently verifies exact translation and every protein component.
It preserves the 5′ UTR/Kozak, TAA–TAG–TGA stop block, 3′ UTR and upstream DNA
context used in the RNA Foundry backbone. PolyA length, sequence and architecture,
any downstream extension, cap chemistry and nucleotide modifications remain
unspecified.

The coding DNA including its first stop is **2,559 nt**. The known RNA segment is
**2,699 nt** and its upstream DNA-template record is **4,640 bp**. Exports end after
the 3′ UTR, without an invented tail or complete molecule length. The coding
sequence excludes the configured motifs and codons below the cutoff; the known
upstream template contains no XbaI site. Fixed UTRs are outside the coding motif
constraint and are not altered.

The coding region has **56.2989% GC**, standard **CAI 0.77396802**, zero forbidden
motif occurrences and zero repeated 12-nt pairs. Its longest homopolymer is 7 nt.
The final paired-bracket run is **28 nt**: later stages do not retain the folding
heuristic’s 14-nt target. This run-length heuristic is not a continuous-stem or
ensemble-accessibility objective. Sequence checks do not imply experimentally
improved expression, stability or immunogenicity. Metrics and verification are in
[the sequence manifest](construct/mrna/construct.json) and
[nucleotide checks](construct/mrna/verification.json).

[GenBank annotation code](construct/mrna/genbank.py) derives features from the
selected regions, native routing-element sources and pinned backbone references.
Antigens carry gene/category, event and region identities. The upstream context
has a reference-matching Tn903 aphA1 CDS and a pUC/pMB1-like origin with its
documented insertion. Unknown intervals and incomplete topology are labelled.
The partial GenBank record is parsed back to check sequence, features, coordinates
and translations. PolyA is labelled without assigning a sequence or length.

## Verification and upstream execution

`yash-vaccine reproduce` rebuilds selection, the internal-cassette checkpoint,
full SEC–antigens–MITD assembly, recorded junction verification and nucleotide
encoding from included evidence. It compares outputs with [candidate.json](candidate.json).
It downloads no patient reads and starts no model inference or cloud job.

A separate fresh run checks all **30 selected junctions**: **1,204 unique peptides ×
six alleles = 7,224 peptide–HLA pairs**. A complete protein scan verifies
positions and threshold agreement with the optimization inputs. All selected
regions occur exactly once, including the exact patient-supported SEC/MITD
sequences. Nucleotide checks verify translation, fixed context, constraints and
GenBank annotations. Fresh inference is distinct from offline cached replay.

[Assembly execution](construct/README.md), [final check](construct/final-check.json),
[nucleotide builder](construct/mrna/run.py) and [upstream guide](../../docs/upstream.md)
provide commands and provenance. Raw DNA/RNA reprocessing needs separately
staged data, references and the recorded runtime; it is not performed by local
reproduction. No tumour-cell presentation or T-cell response experiment is included.
