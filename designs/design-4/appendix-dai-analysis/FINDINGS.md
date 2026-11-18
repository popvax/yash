# Findings

All binding values below are model predictions. **DAI = WT nM / mutant nM on the same HLA allele.** The tables describe predicted binding differences; they do not measure T-cell discrimination or clinical outcomes. “WT” means the supplied single-residue reversion of the somatic change; other residues in the source context are retained. It need not equal an unmodified canonical-reference peptide.

## 1. Ten regions have a higher-DAI peptide–HLA pair

Design 4 chooses the mutation-region lead from eligible candidates by predicted Class-I affinity, without a DAI ranking. Inclusion also uses DNA, RNA and reference-sequence evidence; the fusion is an explicit inclusion exception. Its recorded best target affinity is **508.95 nM**. See [curation code](../pipeline/curate.py) and [selected regions](../selection/selected.tsv).

Among the 14 missense regions with applicable WT-reversion comparisons, 10 have a pair with higher DAI and predicted mutant affinity ≤500 nM:

| Region | Existing lead | Lead DAI | Maximum eligible DAI | Peptide at that maximum | Difference |
|---|---|---:|---:|---|---|
| AQP3 | SRYGEMLHI | 3.76 | 133.67 | RQKELVSRY | Peptide |
| H2BC12 | AVSKGTKAV | 0.44 | 21.45 | ELAKHAVSK | Peptide |
| SF1 | IYNSEGKQL | 1.76 | 12.07 | KQLNTREFRTR | Peptide |
| TJP2 | EAVWVPEGK | 1.58 | 7.60 | WVPEGKMEGM | Peptide |
| PPIC | ATDGHDRSL | 1.37 | 5.35 | Same peptide, C*07:02 | Allele |
| PLEC | ERAEQESAWQL | 0.59 | 4.33 | ERAEQESAW | Peptide |
| RXRA | SVEPKTETY | 1.62 | 3.98 | LSVEPKTETY | Peptide |
| KLF5 | WRFARSNEL | 1.40 | 3.37 | Same peptide, B*08:01 | Allele |
| VPS11 | FTAGSVTL | 0.44 | 2.39 | FTAGSVTLNK | Peptide |
| ARHGAP40 | LHQGWPPKL | 0.65 | 1.76 | HQGWPPKLPK | Peptide |

Five existing lead pairs have DAI ≤1: VPS11, H2BC12, PLEC, ARHGAP40 and FES. These are comparisons between predictions, not evidence that a lead is ineffective or harmful. They can supply an additional variable when planning research comparisons; no superior experimental or clinical performance has been established.

Sources: [all DAI rows](data/encoded-peptide-dai.tsv.gz), [region summary, including both HLA assignments](data/region-dai-summary.tsv), and the [original WT-reversion contrasts](../immunogenicity/class1-counterfactual-contrasts.tsv.gz).

## 2. Changing a lead label does not change the encoded sequence

The encoded windows are 25 aa except KCTD19 at 23 aa and FES at 15 aa. Across the assessed missense comparisons, **18 peptide–HLA pairs covering 17 distinct peptides** have DAI ≤1 and predicted mutant affinity ≤500 nM. KCTD19 and HSPG2, the two regions with the highest-DAI leads, each also contain a peptide with predicted binding decreased by the mutation.

All these sequences remain in their windows when a lead label changes. This establishes sequence membership, not processing or cell-surface presentation. Reducing their encoded content would require a sequence or selection change, with consequences outside this appendix's analysis. No such change is proposed.

Sources: [full residual list](data/residual-low-dai-peptides.tsv) and [selected sequences](../selection/selected.tsv).

## 3. Two C-terminal substitutions have large predicted ratios

| Variant | Mutant peptide | Mutant nM | WT peptide | WT nM | DAI |
|---|---|---:|---|---:|---:|
| KCTD19 p.G920R | RSVSYSILR | 36.9 | RSVSYSILG | 11,691.6 | 316.5 |
| H2BC12 p.E114K | ELAKHAVSK | 139.8 | ELAKHAVSE | 2,998.1 | 21.5 |

Both comparisons use HLA-A*33:03. The changed residue lies at the peptide C terminus. The ratios do not establish actual HLA loading or T-cell recognition, and weak predicted WT binding does not prove that the WT peptide is unpresented or that tolerance is absent.

KCTD19 has **0 coding-strand mutant fragments and 4 coding-strand reference fragments** in the supplied historical tumour-RNA audit, plus **2 opposite-strand mutant fragments**. That shallow coverage does not prove absence of expression. The predicted binding ratio does not establish coding-RNA expression in this specimen or current tumour lesions.

Sources: [WT-reversion contrasts](../immunogenicity/class1-counterfactual-contrasts.tsv.gz), [strand-audit counts](../exports/rna-strand-audit.tsv), and [current interpretation](../review/results-summary.json.gz).

## 4. The original TIER1 flags remain valid under their stated rule

In the original homolog table, `query_best_nm` is the query's minimum across alleles, while `shared_allele_homolog_nm` refers to the recorded shared allele. Those two columns need not describe the same HLA restriction. This appendix compares both peptides on the same recorded allele.

| Model-derived comparison | Pairs |
|---|---:|
| Homolog has lower or equal MHCflurry nM | 16 |
| Homolog has higher MHCflurry nM | 7 |
| Both MHCflurry affinities exceed 500 nM | 8 |
| Of those 8, both MixMHCpred ranks are ≤2% | 8 |
| Of those 8, both MHCflurry percentiles are ≤2% | 1 |
| Both MixMHCpred values are available | 21 |
| Of those 21, relative ordering reverses between models | 2 |

The eight over-500-nM pairs overlap the 16/7 split: six fall in the first group and two in the second. They are not a third disjoint category.

The original shared-HLA predicate requires **each peptide on the same allele** to satisfy at least one branch:

> MHCflurry nM ≤500 **OR** MHCflurry percentile ≤2 **OR** MixMHCpred rank ≤2.

All eight over-500-nM pairs satisfy the MixMHCpred branch. Applying the absolute-nM threshold alone therefore does not invalidate a flag. PRIME is not part of this predicate.

MHCflurry affinity and its calibrated percentile are two representations of one model output. Their different cutoffs can classify the same pair differently. This is not a contradiction between models or evidence that an allele binds weakly in general.

The two genuine cross-model ordering reversals are:

| Pair / allele | MHCflurry query / homolog nM | MixMHCpred query / homolog rank |
|---|---:|---:|
| RFARSNELTR / RFARSDELTR, A*33:03 | 243.91 / 234.88 | 8.284041 / 10.628100 |
| RILEAELSV / RILEAELAV, C*01:02 | 394.41 / 511.48 | 1.105579 / 0.948897 |

This comparison establishes neither occupancy nor cross-reactivity, harm or safety. TIER1 remains a computational prioritisation category.

Sources: [paired values](data/tier1-direction.tsv), [summary](data/tier1-direction-summary.json), [original homolog rows](../safety/offtarget-homologs.tsv.gz), [shared-HLA rule](../pipeline/safety_audit.py), and [MixMHCpred values](../immunogenicity/mixmhcpred3.tsv.gz). The optional verifier reruns MHCflurry inference; MixMHCpred values were inherited, not rerun.
