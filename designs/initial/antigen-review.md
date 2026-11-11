# Evidence for the selected antigens

The design contains 29 distinct 25-aa regions, including 17 mutation regions. Mutation eligibility uses segmented tumour-only GATK filtering; all selected splice events require ND in both searched normal catalogs. Every region also passes the separate hydropathy compatibility screen. CTA regions intentionally use reference proteins. Evidence categories describe observations, not confirmed tumour specificity or presentation.

| Region | Source | Gene / alteration | Best affinity, nM | Max KD19 | Evidence |
|---|---|---|---:|---:|---|
| AG01 | MUTATION | KLF5 D418N | 37.04 | -1.342 | Stronger sequence evidence |
| AG02 | MUTATION | SF1 R218Q | 41.47 | -1.426 | Stronger sequence evidence |
| AG03 | MUTATION | KMT2D G3189R | 18.03 | 0.089 | High priority germline review |
| AG04 | MUTATION | RXRA A241S | 58.08 | -0.237 | Stronger sequence evidence |
| AG05 | MUTATION | HSPG2 R3952L | 81.11 | 0.121 | Limited RNA and germline uncertainty |
| AG06 | MUTATION | PARP3 K37Q | 42.23 | -0.042 | Germline model dependence |
| AG07 | MUTATION | CAST I309V | 88.60 | -0.395 | Near threshold germline confidence |
| AG08 | MUTATION | CCNI P285S | 198.05 | -0.605 | High priority germline review |
| AG09 | MUTATION | PPIC P187S | 36.51 | -0.026 | Patient haplotype annotation needed |
| AG10 | MUTATION | PPP2R1A K280N | 97.82 | 0.332 | RNA-supported mutation hypothesis |
| AG11 | MUTATION | KCTD19 G920R | 36.94 | 0.258 | Limited RNA sampling |
| AG12 | MUTATION | VPS11 D151A | 115.60 | 0.689 | RNA-supported mutation with modest sampling |
| AG13 | MUTATION | ARHGAP40 R438W | 138.22 | -0.211 | RNA-supported mutation hypothesis |
| AG14 | MUTATION | NRXN2 D203H | 103.30 | 0.411 | Limited RNA and near-threshold germline confidence |
| AG15 | MUTATION | PLEC R2171W | 413.26 | -1.058 | Eligible but weaker predicted binding |
| AG16 | MUTATION | DOCK8 G1251R | 253.07 | 0.568 | Near-threshold germline confidence |
| AG17 | MUTATION | SYT12 P391L | 243.74 | 0.511 | Limited RNA sampling |
| AG18 | SPLICE | F5  | 24.24 | 0.363 | ND junction; RNA-supported local sequence |
| AG19 | SPLICE | MUC5B  | 44.87 | 0.500 | ND junction; RNA-supported local sequence |
| AG20 | SPLICE | MUC6  | 24.16 | -0.426 | ND junction; RNA-supported local sequence |
| AG21 | SPLICE | MUC3A  | 260.85 | 0.532 | ND junction; RNA-supported local sequence |
| AG22 | SPLICE | HLA-B  | 66.16 | -0.879 | ND junction; RNA-supported local sequence |
| AG23 | CTA/SELF | PRAME  | 15.30 | 0.047 | Locally supported CTA hypothesis |
| AG24 | CTA/SELF | SPATA31A6  | 33.56 | -0.637 | Selected CTA region lacks RNA coverage |
| AG25 | CTA/SELF | LRRC36  | 21.01 | -0.226 | Sparse local CTA support |
| AG26 | CTA/SELF | YBX2  | 29.75 | -0.089 | Locally supported CTA hypothesis |
| AG27 | CTA/SELF | TEKT5  | 23.74 | -0.089 | Locally supported CTA with annotation gap |
| AG28 | FUSION | CDH1--NRG1  | 70.20 | -0.558 | Predicted alternative fusion tail |
| AG29 | FUSION | CDH1::NRG1  | 508.95 | 0.184 | Explicit in-frame driver-junction inclusion |

## AG01 · KLF5 D418N

Sequence: `TWEGCDWRFARSNELTRHYRKHTGA`. Maximum KD19: -1.342.

The selected N is exactly the D418N substitution in a reference-matching KLF5 context. The RNA allele fraction is 25.9%, close to the DNA fraction of 29.4%. The alternate is present in both DNA/RNA read orientations, with median alternate base quality 37 in the independent check. It is unobserved in the two queried gnomAD cohorts, which is not a measured population frequency of zero. It passes segmented DNA eligibility.

This is one of the stronger mutation candidates on the available evidence. Its high rank is substantially helped by the large supplied gene value. Normal DNA is still needed to establish somatic origin; no tumour-cell peptide or T-cell measurement establishes the predicted target.

## AG02 · SF1 R218Q

Sequence: `PSPEPIYNSEGKQLNTREFRTRKKL`. Maximum KD19: -1.426.

The R218Q sequence matches the nominated transcript context, with 21 alternate RNA fragments and 10 fragments supporting the chosen longer Isovar context. DNA and RNA fractions are similar, approximately 34% and 32%. The alternate has balanced forward/reverse evidence in the fresh base check. It passes segmented DNA eligibility; the queried gnomAD cohorts do not contain this allele.

This is another comparatively well-supported mutation hypothesis. The different transcript protein coordinates do not represent different mutations. Tumour specificity remains a matched-normal question.

## AG03 · KMT2D G3189R

Sequence: `GHTAEKASFGATRGPPAHLLTPSPL`. Maximum KD19: 0.089.

The sequence change is real in the tumour material: DNA is 2,869/5,750 alternate and RNA is 86/177 alternate, with 74 fragments supporting the reconstructed protein context. However, these approximately 50% fractions are also compatible with an inherited heterozygote. The segmented model passes it with GERMQ 93.

The QC-passing gnomAD exome record contains 46 alternate alleles, all in the recorded South Asian group: 46/83,412, or 0.0551%. That is below our 0.1% exclusion threshold, so the rules behaved as written. It is nevertheless meaningful evidence for an inherited alternative. I would confirm this allele in normal DNA before treating it as a somatic target. Strong predicted binding does not resolve that question.

## AG04 · RXRA A241S

Sequence: `DMPVERILEAELSVEPKTETYVEAN`. Maximum KD19: -0.237.

The A241S change and its selected flanks match compatible RXRA reference transcripts. There are 230/658 alternate DNA reads, 23/42 alternate RNA fragments and 18 fragments supporting the reconstructed protein context. The fresh read check finds high-quality alternate bases on both orientations. It passes segmented DNA eligibility, and the two queried gnomAD cohorts do not contain the allele.

RNA has a higher alternate fraction than DNA, approximately 55% versus 35%. That can reflect allelic expression and sampling; it does not by itself invalidate the reconstruction or establish tumour specificity. This remains one of the better-supported research mutation candidates.

## AG05 · HSPG2 R3952L

Sequence: `LPALTNTHHELLLDVEFKPLAPDGV`. Maximum KD19: 0.121.

DNA support is substantial, but the expression fraction rests on only 5 alternate and 4 reference RNA fragments; the longer protein context has 4 supporting fragments. The mutation is at position 12 of this 25-aa window, rather than the usual position 13, and the sequence/highlight correctly reflect that offset.

The allele occurs once in QC-passing gnomAD genomes. Segmented GERMQ is 93. The high supplied gene value should not obscure the limited RNA sampling. I would retain an explicit limited-RNA and germline-uncertainty flag and resolve normal DNA before prioritizing it highly.

## AG06 · PARP3 K37Q

Sequence: `EEDPFRSTAEALQAIPAEKRIIRVD`. Maximum KD19: -0.042.

The K37Q sequence has strong local support: 81 alternate RNA fragments, 75 fragments for the longer protein context and 1,094 alternate DNA reads. The alternate is well represented in both read orientations. The main uncertainty is inherited origin, not whether the altered RNA exists.

The QC-passing gnomAD exome record contains 23 alternate alleles, including 22/85,540 in the recorded South Asian group (0.0257%). Segmented GERMQ is 44. Normal confirmation should precede a somatic claim.

## AG07 · CAST I309V

Sequence: `GITGPPADSSKPVGPDDAIDALSSD`. Maximum KD19: -0.395.

The I309V sequence and its flanks are consistent with the RNA reconstruction, supported by 30 alternate RNA fragments and 22 longer-context fragments. The native DNA and fresh read evidence do not show an obvious one-orientation artifact.

Segmented GERMQ is 22, close to the research gate of 20. The allele occurs in both gnomAD cohorts, with a largest recorded group frequency of 0.0414%. This is a conservative I-to-V change. Exact-reference exclusion does not establish that T cells distinguish it from the wild-type peptide. I would resolve inherited status and assess mutant-versus-self recognition rather than prioritize it on gene expression alone.

## AG08 · CCNI P285S

Sequence: `DKGVFRLHPSSVSGPDFSKDNSKPE`. Maximum KD19: -0.605.

The P285S sequence is supported by 41 alternate RNA fragments, including 31 supporting the longer reconstructed context. DNA and RNA fractions are approximately 46% and 44%. There is good two-orientation base evidence.

Segmented GERMQ is 77. QC-passing exomes contain 56 alternate alleles globally; the recorded South Asian count is 53/86,256 (0.0614%). It passes our frequency threshold but remains a serious inherited-alternative concern. Its best predicted affinity, about 198 nM, passes the recorded gate.

## AG09 · PPIC P187S

Sequence: `SIELQATDGHDRSLTSCSIINSGKI`. Maximum KD19: -0.026.

P187S is supported by 9 alternate RNA fragments and 6 longer-context fragments. It passes segmented DNA eligibility, although GERMQ is more modest than the strongest candidates. The allele is unobserved in the queried gnomAD cohorts.

The 25-aa window contains a second difference from the annotated reference: N190S. This is not an invented residue: 8 RNA fragment names spanning both genomic sites support the P187S-plus-N190S haplotype, and 30 support reference P187 with the same N190S background. The current DNA record also contains a phased-site link to the nearby site. The best predicted peptide ATDGHDRSL includes P187S and stops before N190S. Keep the RNA-derived sequence, but annotate its background difference and do not silently label it a single-difference reference window. The status of the background allele has not been established by this review.

## AG10 · PPP2R1A K280N

Sequence: `TELQKAVGPEITNTDLVPAFQNLMK`. Maximum KD19: 0.332.

K280N has 2,426 alternate DNA reads at depth 5,050, 91 alternate RNA fragments out of 203, and 60 fragments supporting the longer Isovar protein context. Its selected window differs from the compatible annotated reference only at K→N. Both DNA and RNA independently show high-quality alternate bases in both read orientations. The best eligible binding prediction is 97.82 nM. The allele is unobserved in the queried gnomAD cohorts, which is not a measured zero frequency. Segmented GERMQ is 43 and remains 43 with higher group priors; the diagnostic without segmentation also retains the research gate. Near-heterozygous DNA balance still does not establish somatic origin.

## AG11 · KCTD19 G920R

Sequence: `MDLLIQRGLSRSVSYSILRKYLQED`. Maximum KD19: 0.258.

G920R has 134 alternate DNA reads at depth 419. Two alternate RNA fragment names out of seven support the allele, and both support the longer Isovar context. The substitution is at window position 19, not the usual position 13; the selected frame and context agree with compatible annotation. The fresh check finds the alternate on both read orientations with median base quality 37. Three QC-passing gnomAD exome alleles are recorded, all in the South Asian group (3/86,234); the group frequency remains below 0.1%. Segmented GERMQ is 80, or 68 with higher group priors. The 36.94 nM best binding prediction is strong within this computational screen, but does not compensate for the small RNA sample or establish immunogenicity.

## AG12 · VPS11 D151A

Sequence: `HENLNFMAIGFTAGSVTLNKGDITR`. Maximum KD19: 0.689.

D151A has 78 alternate DNA reads at depth 205, seven alternate RNA fragments out of 23, and four fragments supporting the longer reconstructed context. The selected 25-aa window matches the compatible reference apart from D→A. The fresh check finds high-quality alternate bases on both orientations. Segmented GERMQ is 34 and remains 34 with the higher-prior diagnostic; the no-segmentation diagnostic also passes. The allele is unobserved in the queried gnomAD cohorts. Best predicted affinity is 115.60 nM. Limited context sampling and unconfirmed inherited status remain attached to this candidate.

## AG13 · ARHGAP40 R438W

Sequence: `TVMAPNLFLHQGWPPKLPKGKEKQL`. Maximum KD19: -0.211.

R438W has 196 alternate DNA reads at depth 596, nine alternate RNA fragments out of 20, and eight longer-context fragments. Its selected reference comparison contains only the intended R→W substitution. The fresh check finds median alternate base quality 37 and both read orientations in DNA and RNA. Rare QC-passing observations occur in both gnomAD cohorts; the highest recorded group frequency is approximately 0.0227%, below the policy cutoff. Segmented GERMQ is 87, or 72 with higher group priors, and the no-segmentation diagnostic passes. Its best binding prediction is 138.22 nM. These data establish an expressed sequence hypothesis, not a confirmed somatic or immunogenic target.

## AG14 · NRXN2 D203H

Sequence: `LLGSQGLRGATAHPLCAPARNPCAN`. Maximum KD19: 0.411.

D203H has 142 alternate DNA reads at depth 293 but only two alternate RNA fragment names out of eight; both support the longer Isovar context. The context differs from the compatible annotated reference only at D→H. The fresh base check observes the alternate on both read orientations, with median base quality 37. The allele is unobserved in the queried gnomAD cohorts. Segmented GERMQ is 22 and remains 22 with higher group priors. The diagnostic without segmentation is native PASS with GERMQ 19, below the additional research gate; the selected workflow nevertheless uses segmentation. Best predicted affinity is 103.30 nM. Retain the sparse-RNA and germline-model uncertainty.

## AG15 · PLEC R2171W

Sequence: `RRLRERAEQESAWQLQLAQEAAQKR`. Maximum KD19: -1.058.

R2171W has 65 alternate DNA reads at depth 232, eight alternate RNA fragments out of 16, and eight longer-context fragments. The selected context differs from compatible reference only at R→W, with high-quality alternate bases on both orientations in the fresh check. QC-passing gnomAD exomes contain 75 alternate alleles and genomes contain two; the largest recorded group frequency is approximately 0.0207%, below 0.1%. Segmented GERMQ is 57, or 51 with higher group priors, and the no-segmentation diagnostic passes. Its best eligible affinity is 413.26 nM: it passes the gate but is the weakest of the 18 best-affinity values. Binding and inherited status remain predictions/uncertainties rather than measured tumour presentation.

## AG16 · DOCK8 G1251R

Sequence: `DEEQERAGAINQNVALAIAGNNFNL`. Maximum KD19: 0.568.

G1251R has 124 alternate DNA reads at depth 290, eight alternate RNA fragments out of 12, and six longer-context fragments. The substitution is at window position 6; the selected context matches the annotated reference except for G→R. Fresh DNA and RNA bases have median alternate quality 37 and both read orientations. Segmented GERMQ is 23 and remains 23 with higher group priors; the no-segmentation diagnostic passes. The allele is unobserved in the queried gnomAD cohorts. Best predicted affinity is 253.07 nM. Its high RNA alternate fraction is based on 12 fragments and does not establish a somatic mutation.

## AG17 · SYT12 P391L

Sequence: `GRGDNVGHVIIGLSASGMGTTHWNQ`. Maximum KD19: 0.511.

P391L has 290 alternate DNA reads at depth 845, but only two alternate RNA fragment names out of eight; both support the longer reconstructed context. The window matches compatible reference apart from P→L. The fresh check finds alternate bases in both orientations with median quality 37. Rare QC-passing population alleles occur in both queried cohorts, with no South Asian alternate observed in those records and a largest recorded group frequency below 0.003%. Segmented GERMQ remains 93 in the higher-prior diagnostic, and the no-segmentation diagnostic passes. Best predicted affinity is 243.74 nM. Its low weighted score does not create a new eligibility failure, but the limited RNA evidence remains explicit.

## AG18 · F5

Sequence: `ALGQMPIFPDLSHTTLSLDFSQTNL`. Maximum KD19: 0.363.

ND in both searched normal catalogs. Two complete-window fragment names agree with the nominated direction. The selected sequence has an additional S→F difference whose targeted check found 0 alternate / 1,300 reference DNA fragments and 5 / 1,441 RNA fragments; its cause is unresolved. This is retained under the stated policy, which adds no separate substitution gate.

## AG19 · MUC5B

Sequence: `LTTTATTPAATISTVTPSSALGTTH`. Maximum KD19: 0.500.

ND in both searched normal catalogs. Two complete-window fragment names agree with the nominated direction. The sequence is supported in the annotated upstream frame. This event has not received an additional per-substitution audit or an alignment-pattern independence requirement.

## AG20 · MUC6

Sequence: `PSSFSMAKTSTSLPYHTSSTHHPEV`. Maximum KD19: -0.426.

ND in both searched normal catalogs. Three complete-window fragment names support the sequence; two agree with the nominated direction and one opposes it. The selected window also contains RNA differences beyond the splice join; its T→M-associated base has 8 alternate / 5,583 reference DNA fragments and 9 / 269 RNA fragments. These are observations, not a validated somatic classification.

## AG21 · MUC3A

Sequence: `TSQPTTITPSSVGFSGSLPMMTDLT`. Maximum KD19: 0.532.

ND in both searched normal catalogs. Three complete-window fragment names agree with the nominated direction. The sequence is supported in the annotated upstream frame. Repetitive context, individual substitutions and original-molecule independence remain separate review questions.

## AG22 · HLA-B

Sequence: `ITQRKWEAAREAEQRSLPGGHVRGV`. Maximum KD19: -0.879.

ND in both searched normal catalogs. Two complete-window fragment names agree with the nominated direction. The sequence reconstructs in the verified annotated upstream frame. The motif points opposite the nominated gene; BAM orientation, rather than motif-derived XS, determines the limited strand decision. Normal allele-specific or unrepresented transcripts are not ruled out.

## AG23 · PRAME

Sequence: `SYEDIHGTLHLERLAYLHARLRELL`. Maximum KD19: 0.047.

The selected PRAME reference region has positive supplied expression. Every selected amino-acid codon has at least 2 reference-matching RNA fragment names under the recorded base-quality filter; 2 fragment names cover the best predicted peptide YLHARLREL. No single merged fragment pair spans the entire 25-aa sequence. A small minority of observations differs at one other codon, so this is not a fully reconstructed patient protein haplotype.

The native 85-fragment figure is transcript-level, not support for this particular peptide. PRAME is testis-enriched in the current Human Protein Atlas, which also reports RNA detection in some normal tissues. I found no PRAME IHC or tumour-cell peptide measurement in the reviewed candidate evidence. This is a reasonable CTA research hypothesis, with regional/protein confirmation still needed. [Human Protein Atlas PRAME](https://www.proteinatlas.org/ENSG00000185686-PRAME/tissue).

## AG24 · SPATA31A6

Sequence: `SLSASQPPEPSLPLEHPSPEPPALF`. Maximum KD19: -0.637.

Positive gene expression and 98 transcript-level fragments do not support this particular reference window. In the recorded regional check, 24 of its 25 codons have zero reference-matching RNA fragments; the last codon has one. Neither the complete region nor the strongest predicted peptide ASQPPEPSL is supported by a qualifying fragment.

That strongest peptide also occurs in SPATA31A1, A3, A5 and A7 reference proteins. The existing CTA policy allows matches within the curated CTA family, so this is not a pipeline violation; it is an important limit on gene-specific attribution. The family has segmental duplications. The source normal-tissue field is missing, while current HPA labels SPATA31A6 testis-enriched. I would deprioritize this window until its local sequence and expression attribution are established. [HPA SPATA31A6](https://www.proteinatlas.org/ENSG00000185775-SPATA31A6/tissue), [SPATA31 duplication study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5338094/).

## AG25 · LRRC36

Sequence: `GSYAGKIHSIGDAFRNFKNLRSLDL`. Maximum KD19: -0.226.

The native 39-fragment support is transcript-level. The selected region has just one reference-matching fragment at each of 24 codons and none at the remaining codon under the new quality filter. One fragment matches SYAGKIHSI; no fragment pair establishes the entire 25-aa window. This is substantially thinner evidence than the gene-level number suggests.

The saved normal-tissue list says testis, but current HPA calls RNA testis-enriched and detected in some tissues, rather than testis-exclusive. Its normal protein datasets also warrant inspection before making a restriction claim. I would retain this as provisional, with local sequence and normal-tissue review ahead of claims about tumour specificity. [HPA LRRC36](https://www.proteinatlas.org/ENSG00000159708-LRRC36/tissue).

## AG26 · YBX2

Sequence: `PNRRKSRRFIPRPPSVAPPPMVAEI`. Maximum KD19: -0.089.

One fragment pair matches the full selected 25-aa region; two match VAPPPMVAEI. All codons have at least two reference-matching fragment names. This is direct regional RNA evidence, although it remains a small sample and does not establish protein presentation.

The saved mTEC mean/median are approximately 0.219/0.208 TPM, higher than the PRAME values but not an immunogenicity assay. Current HPA describes testis-enriched RNA with detection in many tissues. I would keep YBX2 as a reference-antigen hypothesis and inspect normal cell types rather than rely on the source testis-only label. [HPA YBX2](https://www.proteinatlas.org/ENSG00000006047-YBX2/tissue).

## AG27 · TEKT5

Sequence: `FTIDDTLQTLKLRLRETQDTLQLLV`. Maximum KD19: -0.089.

Two fragment pairs match the full selected 25-aa sequence and the best predicted peptide TIDDTLQTL. Every codon has at least two matching fragment names. This provides direct local RNA support despite a relatively small supplied gene value.

The source normal-tissue name Actin is not an organ and must not be interpreted as a normal-tissue safety annotation. Current HPA classifies TEKT5 RNA as testis-enriched and detected in a single tissue category; normal protein annotation is pending. Its recorded mTEC mean/median are approximately 0.212/0.055 TPM. Keep the original field for provenance, use a corrected sourced interpretation, and retain protein/normal-tissue validation as unresolved. [HPA TEKT5](https://www.proteinatlas.org/ENSG00000153060-TEKT5/tissue).

## AG28 · CDH1--NRG1

Sequence: `EPDTFMEQKITYIYIHHWDKPSCKM`. Maximum KD19: -0.558.

The selected region belongs to the predicted frameshift context at CDH1 chr16:68815759 → NRG1 chr8:32727949. Exact-breakpoint FusionInspector evidence reports 62 junction-plus-spanning support units, versus 473 unfused left-partner units. Ranking uses 62/(62+473) times supplied CDH1 expression. This context is distinct from the explicitly included in-frame fusion and is not evidence of translation of the functional driver protein. Detailed saved artifact fields remain uncomputed.

## AG29 · CDH1::NRG1

Sequence: `STYTALIIATDNATSTSTTGTSHLV`. Maximum KD19: 0.184.

The source call reports 200 junction reads and one spanning fragment at the exact in-frame breakpoint pair. A 38-peptide / six-allele screen excludes eight reference-matching peptides. The strongest remaining prediction is 508.95 nM, included as a documented exception to the general 500 nM gate. No expression-weighted score or support fraction is borrowed from the other fusion.

The [hydropathy audit](selection/hydrophobicity-screen.tsv) retains all 223 nominated regions, including excluded PTCH1. Its DNA/RNA observations remain in the shared evidence; exclusion is a construct-compatibility decision.
