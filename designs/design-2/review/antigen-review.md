# Review of all 24 selected regions

Reviewed against the completed matched-normal outputs. This is a computational evidence review, not an experimental or clinical validation. No new selection rule is introduced here.

The ten retained mutation regions now have direct normal-reference evidence and paired native PASS calls. Seven old targets are present in white-blood-cell DNA and are removed. AQP3 and TGFBR2 are added. All twelve RNA-source sequences remain unchanged; their earlier regional reviews remain relevant.

Counts use fragment names, not UMI-confirmed independent molecules. Mutation RNA ALT/total describes local allele support; the longer protein-context count is different. Supplied expression retains the vendor's transformation and is not asserted to be calibrated TPM.

All 12 mutation regions meeting the evidence, expression-availability, binding and construct-compatibility criteria are included. Their expression-weighted scores order comparison tables, not a top-N inclusion list or the assembled construct. Expression magnitude and RNA ALT fraction do not remove lower-ranked qualifying mutation regions. Splice and CTA selections still use expression-weighted ranking.

| ID | Source | Gene / change | Best affinity (nM) | Normal REF / ALT | RNA ALT / total | Longer-context fragments | Supplied gene value |
|---|---|---|---:|---:|---:|---:|---:|
| AG01 | MUTATION | KLF5 p.D418N | 37.04 | 50 / 0 | 30 / 116 | 21 | 1297.240 |
| AG02 | MUTATION | SF1 p.R218Q | 41.47 | 57 / 0 | 21 / 66 | 10 | 157.412 |
| AG03 | MUTATION | AQP3 p.C11Y | 29.42 | 46 / 0 | 2 / 4 | 2 | 19.636 |
| AG04 | MUTATION | RXRA p.A241S | 58.08 | 42 / 0 | 23 / 42 | 18 | 23.015 |
| AG05 | MUTATION | HSPG2 p.R3952L | 81.11 | 46 / 0 | 5 / 9 | 4 | 167.615 |
| AG06 | MUTATION | TGFBR2 p.N515fs | 38.03 | 40 / 0 | 85 / 206 | 67 | 55.743 |
| AG07 | MUTATION | PPIC p.P187S | 36.51 | 47 / 0 | 9 / 41 | 6 | 72.916 |
| AG08 | MUTATION | KCTD19 p.G920R | 36.94 | 49 / 0 | 2 / 7 | 2 | 9.912 |
| AG09 | MUTATION | VPS11 p.D151A | 115.60 | 32 / 0 | 7 / 23 | 4 | 15.707 |
| AG10 | MUTATION | ARHGAP40 p.R438W | 138.22 | 49 / 0 | 9 / 20 | 8 | 6.272 |
| AG11 | MUTATION | PLEC p.R2171W | 413.26 | 26 / 0 | 8 / 16 | 8 | 92.286 |
| AG12 | MUTATION | SYT12 p.P391L | 243.74 | 56 / 0 | 2 / 8 | 2 | 7.922 |
| AG13 | SPLICE | F5  | 24.24 | — | — | — | 470.598 |
| AG14 | SPLICE | MUC5B  | 44.87 | — | — | — | 284.038 |
| AG15 | SPLICE | MUC6  | 24.16 | — | — | — | 65.815 |
| AG16 | SPLICE | MUC3A  | 260.85 | — | — | — | 656.582 |
| AG17 | SPLICE | HLA-B  | 66.16 | — | — | — | 185.958 |
| AG18 | CTA/SELF | PRAME  | 15.30 | — | — | — | 17.014 |
| AG19 | CTA/SELF | SPATA31A6  | 33.56 | — | — | — | 13.606 |
| AG20 | CTA/SELF | LRRC36  | 21.01 | — | — | — | 9.152 |
| AG21 | CTA/SELF | YBX2  | 29.75 | — | — | — | 15.478 |
| AG22 | CTA/SELF | TEKT5  | 23.74 | — | — | — | 4.845 |
| AG23 | FUSION | CDH1--NRG1  | 70.20 | — | — | — | — |
| AG24 | FUSION | CDH1::NRG1  | 508.95 | — | — | — | — |

## Individual assessments

### AG01 · KLF5 p.D418N

`TWEGCDWRFARSNELTRHYRKHTGA` · `REG_e16d2a87273645d8d7b4` · max KD19 -1.342.

Substantial DNA and RNA support, with 50 reference-only normal fragments. Its large supplied expression value raises the exported comparison score, but all 12 eligible mutation regions are included regardless of that rank. This is not a measured peptide-abundance ranking.

### AG02 · SF1 p.R218Q

`PSPEPIYNSEGKQLNTREFRTRKKL` · `REG_1dfce9e906fa21c20196` · max KD19 -1.426.

Substantial RNA allele and longer-context support. DNA and RNA alternate fractions are similar. The normal sample now supports the somatic interpretation that was unresolved without it.

### AG03 · AQP3 p.C11Y

`MGRQKELVSRYGEMLHIRYRLLRQA` · `REG_c310fc71d7b0fb823535` · max KD19 -0.232.

Newly selected C11Y. DNA evidence is substantial, but the RNA estimate is only 2 alternate of 4 fragments; both support the longer context. The apparent 50% RNA fraction is poorly sampled. No extra three-fragment rule was introduced.

### AG04 · RXRA p.A241S

`DMPVERILEAELSVEPKTETYVEAN` · `REG_8687690f602539a7060e` · max KD19 -0.237.

Good local DNA/RNA support. RNA alternate fraction is higher than DNA model fraction, which can reflect allele-specific expression or sampling; the discrepancy is not itself evidence of an invalid sequence.

### AG05 · HSPG2 p.R3952L

`LPALTNTHHELLLDVEFKPLAPDGV` · `REG_9954d136ed6d4ba9f759` · max KD19 0.121.

Strong DNA evidence and adequate normal reference support. RNA is limited to 5 alternate of 9 fragments, with 4 supporting the longer context. Rare population observation is retained, rather than labelled absent.

### AG06 · TGFBR2 p.N515fs

`EPPFGSKVREHPCVESMKDKRVERS` · `REG_f4ca2bba236d14d54d3f` · max KD19 -0.958.

Newly selected frameshift region. Eighty-five alternate RNA fragments, including 67 for the longer context, give substantially stronger sequence support than the other newly selected region AQP3. Transcript-dependent protein numbering is not a separate genomic mutation. This is independent paired rediscovery, not a claim that the clinical frameshift was previously unknown.

### AG07 · PPIC p.P187S

`SIELQATDGHDRSLTSCSIINSGKI` · `REG_e6aeb364b2cc7e909162` · max KD19 -0.026.

The RNA-derived window also retains the previously reviewed N190S background difference. The focal P187S is supported somatic here. The best peptide ATDGHDRSL includes the focal change and ends before that background difference. Do not describe the complete window as having exactly one difference from reference or claim that every background residue has been separately classified.

### AG08 · KCTD19 p.G920R

`MDLLIQRGLSRSVSYSILRKYLQED` · `REG_9dde223e9bf58444173b` · max KD19 0.258.

Only two alternate RNA fragment names out of seven; both support the longer context. The matched normal now resolves focal somatic support, but does not improve that small RNA sample.

### AG09 · VPS11 p.D151A

`HENLNFMAIGFTAGSVTLNKGDITR` · `REG_7bb9945d495d7fb8ce3f` · max KD19 0.689.

Seven alternate RNA fragments and four for the longer context provide modest sequence support. No exact match in the searched normal/reference sequence models was found for the scored mutant peptides; that is a scoped computational negative.

### AG10 · ARHGAP40 p.R438W

`TVMAPNLFLHQGWPPKLPKGKEKQL` · `REG_7176562878cc11bcbfaa` · max KD19 -0.211.

Nine alternate RNA fragments and eight for the longer context. Rare population observations remain below the policy cutoff, while direct normal reference evidence supports the focal somatic interpretation.

### AG11 · PLEC p.R2171W

`RRLRERAEQESAWQLQLAQEAAQKR` · `REG_f0a90a31c4c6e8899224` · max KD19 -1.058.

Eight alternate RNA fragments support the longer context. The best target prediction is approximately 413 nM, relatively near the 500-nM gate. It qualifies computationally without establishing tumour-cell display or a T-cell response.

### AG12 · SYT12 p.P391L

`GRGDNVGHVIIGLSASGMGTTHWNQ` · `REG_f93fe5a69c414329a31d` · max KD19 0.511.

Only two alternate RNA fragments out of eight; both support the longer context. Adequate normal evidence supports the somatic call but does not eliminate uncertainty from sparse RNA sampling.

### AG13 · F5

`ALGQMPIFPDLSHTTLSLDFSQTNL` · `REG_ae4e341fdea4a4a9e453` · max KD19 0.363.

RNA reads support the join and extra S-to-F change that makes the scored peptides differ from the main reference. Recorded support is 5 fragments at the changed base and 2 across the full 25-aa window. The earlier targeted check reported 5 / 1,441 RNA fragments at the base and 0 alternate / 1,300 reference tumour-DNA fragments.

Recorded complete-window fragments: 2; strand agreement/opposition 2/0. The splice motif disagrees with the nominated gene direction. The selected window passed the recorded read-orientation rule and annotated upstream-frame reconstruction. ND means not detected at these coordinates in the two searched normal catalogs, not absent from every healthy tissue. New normal DNA does not settle RNA direction, splicing or translation.

### AG14 · MUC5B

`LTTTATTPAATISTVTPSSALGTTH` · `REG_ea6a38a75f1c7d884bee` · max KD19 0.500.

Two full-window RNA fragments support the join and extra S-to-I change that makes the scored peptides differ from the main reference. Low-complexity mapping and molecule independence remain unresolved. The recorded event/region tie-break is retained.

Recorded complete-window fragments: 2; strand agreement/opposition 2/0. The splice motif disagrees with the nominated gene direction. The selected window passed the recorded read-orientation rule and annotated upstream-frame reconstruction. ND means not detected at these coordinates in the two searched normal catalogs, not absent from every healthy tissue. New normal DNA does not settle RNA direction, splicing or translation.

### AG15 · MUC6

`PSSFSMAKTSTSLPYHTSSTHHPEV` · `REG_35621fc5dd92cf1c215d` · max KD19 -0.426.

RNA reads support the join and extra T-to-M change: 9 fragments at the changed base and 3 across the full window. Two full-window fragments agree with the nominated strand; one opposes it. Earlier targeted counts were 9 / 269 RNA fragments and 8 alternate / 5,583 reference tumour-DNA fragments. The selected sequence also matches alternate-reference translation ENSP00000487059.1. At the tested position (chr11:1017576), Yash's normal DNA showed 150 fragments with the usual base and 0 with the methionine-producing change.

Recorded complete-window fragments: 3; strand agreement/opposition 2/1. The splice motif disagrees with the nominated gene direction. The selected window passed the recorded read-orientation rule and annotated upstream-frame reconstruction. ND means not detected at these coordinates in the two searched normal catalogs, not absent from every healthy tissue. New normal DNA does not settle RNA direction, splicing or translation.

### AG16 · MUC3A

`TSQPTTITPSSVGFSGSLPMMTDLT` · `REG_d58dfb6f2f0b754ec2b9` · max KD19 0.532.

Three full-window RNA fragments support the join and extra I-to-F change that makes the scored peptides differ from the main reference. All three agree with the nominated strand. The repetitive mucin context still complicates alignment attribution.

Recorded complete-window fragments: 3; strand agreement/opposition 3/0. The splice motif disagrees with the nominated gene direction. The selected window passed the recorded read-orientation rule and annotated upstream-frame reconstruction. ND means not detected at these coordinates in the two searched normal catalogs, not absent from every healthy tissue. New normal DNA does not settle RNA direction, splicing or translation.

### AG17 · HLA-B

`ITQRKWEAAREAEQRSLPGGHVRGV` · `REG_fadc531c5476f4384f0d` · max KD19 -0.879.

The first 15 residues ITQRKWEAAREAEQR occur in compatible HLA-C reference proteins. They are not the scored tumour-specific target: scoring uses peptides overlapping the nominated junction at window residues 18 (one based). Best target EQRSLPGGHVR is approximately 66 nM on A*33:03 and has no match in the searched personal-normal models. Two full-window fragments support a hypothesis in a highly homologous HLA region, not a validated novel HLA protein.

Recorded complete-window fragments: 2; strand agreement/opposition 2/0. The splice motif disagrees with the nominated gene direction. The selected window passed the recorded read-orientation rule and annotated upstream-frame reconstruction. ND means not detected at these coordinates in the two searched normal catalogs, not absent from every healthy tissue. New normal DNA does not settle RNA direction, splicing or translation.

### AG18 · PRAME

`SYEDIHGTLHLERLAYLHARLRELL` · `REG_8389833c12d42dfd778f` · max KD19 0.047.

Intentional reference-protein target. The recorded regional audit has at least two matching fragments per codon and two covering YLHARLREL, but no single merged pair spans all 25 residues. No PRAME IHC or patient tumour-cell peptide measurement was found in that review. Do not use the transcript-level 85-fragment count as peptide support.

### AG19 · SPATA31A6

`SLSASQPPEPSLPLEHPSPEPPALF` · `REG_a391cc51fbd08e062be3` · max KD19 -0.637.

Weakest local CTA sequence evidence: 24 of 25 codons have no qualifying reference-matching RNA fragment in the recorded regional audit; neither the full region nor the best peptide is covered. The best peptide is shared with other SPATA31A-family reference proteins. Positive supplied gene expression does not resolve this regional/isoform attribution. Kept under the existing CTA policy, not newly confirmed by normal DNA.

### AG20 · LRRC36

`GSYAGKIHSIGDAFRNFKNLRSLDL` · `REG_73d812fbc5c9353da35b` · max KD19 -0.226.

Recorded regional evidence is sparse: one matching fragment at 24 codons, none at one; one supports SYAGKIHSI and none establishes the whole 25-aa region. The source testis label is not evidence of testis-exclusive expression.

### AG21 · YBX2

`PNRRKSRRFIPRPPSVAPPPMVAEI` · `REG_436ff2a1af597e184676` · max KD19 -0.089.

One fragment pair matches the full selected window, two match VAPPPMVAEI, and each codon has at least two matching fragment names in the recorded regional check. This is local RNA support, not protein presentation. Normal-tissue restriction remains a separate question.

### AG22 · TEKT5

`FTIDDTLQTLKLRLRETQDTLQLLV` · `REG_faac00cfb1186184893a` · max KD19 -0.089.

Two fragment pairs match the full selected window and TIDDTLQTL in the recorded regional check. The literal source normal-tissue label Actin is not an organ and cannot support a tissue-restriction claim; it is preserved as a provenance error. Normal protein/tissue validation remains unresolved.

### AG23 · CDH1--NRG1

`EPDTFMEQKITYIYIHHWDKPSCKM` · `REG_98875bfa1d09b7fab791` · max KD19 -0.558.

Alternative predicted frameshift-tail context at the recorded breakpoint, distinct from the explicitly included in-frame driver junction. Exact-breakpoint review reports 62 junction-plus-spanning units and 473 unfused left-partner units; weighting uses 62/(62+473) times supplied CDH1 expression. Other generic baseline columns are not this denominator. Translation of this alternative context and detailed saved artifact metrics remain unresolved.

### AG24 · CDH1::NRG1

`STYTALIIATDNATSTSTTGTSHLV` · `REG_1f7babf47bbafeebeedc` · max KD19 0.184.

Explicit in-frame driver-junction region. The recorded exact-breakpoint call reports 200 junction reads and one spanning fragment. Its best non-reference prediction is approximately 508.95 nM, an explicit inclusion exception above 500 nM. The NAT motif is a potential N-linked glycosylation sequon; occupancy and its effect on processing/presentation were not modelled. This region and the alternative tail are two hypotheses from the same gene pair, not two independent drivers.

## Construct-level findings

The 724-aa protein contains each selected region exactly once, 43 linker residues, a 26-aa SEC element and a 55-aa MITD. Fresh inference covered all 6,600 order/linker options (578,112 peptide/HLA queries), and a separate final-junction run covered 6,132 queries. All four recorded optimizer objectives were proven optimal for this finite search and its unwanted-junction binding objective. This does not optimize cleavage or intended antigen presentation.

The final junction scan still has one row at affinity percentile ≤0.5%, 26 at ≤2%, and 24 at affinity ≤500 nM; none is ≤50 nM. Passing verification means agreement with the objective and independent scan, not absence of predicted junction binders. There are no high-hydropathy windows confined to antigens/linkers; SEC/MITD overlaps remain separately reported and are not topology-validated.

The fixed protein was synonymously encoded with PopRNA (human usage, 0.10 cutoff, seed 0). Translation, forbidden motifs, codon cutoff and all 64 GenBank feature annotations pass their hard checks. Coding GC is 56.08%. The final folding heuristic is 21 consecutive paired brackets versus a target of 14, so that heuristic target is **not retained** after the sequential repairs. Repeated 12-nt pairs are zero. These structure/codon heuristics are not expression or stability measurements. PolyA and the downstream extension remain unspecified.

## Evidence

- [All 17 original mutation reassessments](final/current-17-final-audit.tsv), including seven not retained.
- [Selected regions](../selection/selected.tsv), [DNA evidence](../inputs/mutations/variant-evidence.tsv), [RNA evidence](../inputs/rna/rna-evidence.tsv), and [per-allele best mutation peptides](../inputs/predictions/mutations/best-peptide-per-HLA.tsv).
- [Splice sequence comparisons and focal MUC6 normal-DNA evidence](splice-sequence-check.json), reproduced with [the check script](../pipeline/check_splice_sequences.py).
- [Earlier regional observations for unchanged RNA sources](../../initial/antigen-review.md): retained as historical measurements, without importing its old tumour-only somatic assignments.
- [Fresh final junction verification](../construct/final-check.json) and [nucleotide/GenBank verification](../construct/mrna/verification.json).
