# Candidate assessment

The main protein encodes 15 mutation windows and one in-frame fusion window. Eight candidates remain in reserve. All are research candidates; presentation, specificity and benefit have not been established.

## Encoded regions

| Gene / change | AA | Coding / opposite mutant fragments | Predicted lead / HLA / nM | Limitation |
| --- | --- | --- | --- | --- |
| KLF5 p.D418N | 25 | 30 / 0 | WRFARSNEL / HLA-C*07:02 / 37.04 | Normal KLF-family homologs warrant specificity testing; cross-recognition is unmeasured. |
| SF1 p.R218Q | 25 | 19 / 2 | IYNSEGKQL / HLA-C*07:02 / 41.47 | The best recorded Class-II core is unchanged; mutation-specific CD4 recognition is unproven. |
| AQP3 p.C11Y | 25 | 2 / 0 | SRYGEMLHI / HLA-C*07:02 / 29.42 | Two coding-strand mutant fragments provide limited RNA sampling. |
| RXRA p.A241S | 25 | 22 / 1 | SVEPKTETY / HLA-B*15:08 / 58.08 | Wild-type and RXR-family homolog comparisons remain necessary; clinical isoform labels differ. |
| HSPG2 p.R3952L | 25 | 5 / 0 | LTNTHHELL / HLA-C*01:02 / 81.11 | DNA/RNA support does not establish malignant-cell presentation. |
| TGFBR2 p.N515fs | 25 | 85 / 0 | ESMKDKRVER / HLA-A*33:03 / 38.03 | Frameshift numbering differs by isoform; clonality is unestablished. |
| PPIC p.P187S | 25 | 9 / 0 | ATDGHDRSL / HLA-C*01:02 / 36.51 | Class-II support is weak and the best recorded core is unchanged. |
| KCTD19 p.G920R | 23 | 0 / 2 | RSVSYSILR / HLA-A*33:03 / 36.94 | Somatic DNA-supported; clonality and coding-RNA expression remain unestablished. The separate protein alternative excludes it. |
| VPS11 p.D151A | 25 | 7 / 0 | FTAGSVTL / HLA-C*01:02 / 115.60 | Predicted presentation and mutation-specific T-cell recognition remain untested. |
| ARHGAP40 p.R438W | 25 | 9 / 0 | LHQGWPPKL / HLA-C*07:02 / 138.22 | Model affinity does not establish presentation or recognition. |
| PLEC p.R2171W | 25 | 7 / 1 | ERAEQESAWQL / HLA-C*07:02 / 413.26 | Mutant-versus-wild-type recognition remains untested. |
| SYT12 p.P391L | 25 | 1 / 1 | NVGHVIIGL / HLA-C*01:02 / 243.74 | One coding-strand and one opposite-strand mutant fragment; RNA support is shallow. |
| H2BC12 p.E114K | 25 | 10 / 1 | AVSKGTKAV / HLA-C*01:02 / 136.46 | Supplied gene-level expression value is missing; histone homologs warrant specificity testing. |
| TJP2 p.S886P | 25 | 29 / 0 | EAVWVPEGK / HLA-A*33:03 / 212.51 | Readthrough/isoform annotation remains ambiguous; the observed peptide window is the stable identity. |
| FES p.R804W | 15 | 1 / 1 | AYEPGQWPSF / HLA-C*07:02 / 94.54 | Only the observed 15-aa window is encoded; one coding-strand mutant fragment. Regional spliced reads are not mutant support. |
| CDH1::NRG1 | 25 | n/a (fusion) | ATDNATSTST / HLA-A*01:01 / 508.95 | Approximately 509 nM Class-I prediction is an explicit inclusion exception. The source records 200 junction reads and one spanning fragment in distinct, non-UMI units; peptide presentation is untested. |

Affinities and HLA assignments in this table come from `immunogenicity/antigen-model-summary.tsv`. They describe the recorded lead in the engineered context, not measured binding. Per-region target and context scores are retained in [the TSV](antigen-review.tsv) and [prediction files](../immunogenicity/). Fragment counts are read-name-based, not UMI molecules. Fusion support is reported separately: 200 recorded junction reads and one spanning fragment in distinct, non-UMI units; these are not 201 independently counted molecules.

The 12 mutation windows with at least two coding-strand mutant fragments sum to 234; FES and SYT12 bring the coding total to 236. The nine opposite-strand fragments include KCTD19's two. These sums do not measure tumour coverage or clonality. The KCTD19-excluded [481-aa protein](../strand-strict-15/construct/protein/cassette.fasta) and its [1,586-nt known RNA segment](../strand-strict-15/construct/mrna/mrna-segments.fasta) remain a computational alternative.

## Reserves

| Gene | Source | AA | Reason for reserve |
| --- | --- | --- | --- |
| SERPINA1 | MUTATION | 25 | S22Y lies in the annotated signal peptide; Class-I nomination exceeds 500 nM. Signal-derived presentation remains unproven. |
| POLE | SPLICE | 25 | RNA-only splice hypothesis with a noncanonical CT-GC motif; translation and tumour specificity remain unconfirmed. |
| PRAME | CTA/SELF | 25 | Unmutated self/CTA candidate; no full-window merged fragment in the recorded audit. Tumour specificity and normal-tissue safety are unestablished. |
| NUTM1 | CTA/SELF | 25 | Unmutated self/CTA candidate; full-window RNA support does not establish malignant-cell origin. |
| MAGEA6 | CTA/SELF | 25 | Unmutated self/CTA candidate; homologous mapping and T-cell specificity require validation. |
| TSPYL6 | CTA/SELF | 25 | Unmutated self/CTA candidate; normal-tissue presentation and useful T-cell recognition are unestablished. |
| MAGEA1 | CTA/SELF | 25 | Unmutated self/CTA candidate; normal-tissue presentation and useful T-cell recognition are unestablished. |
| CDH1--NRG1 | FUSION | 25 | Alternative out-of-frame fusion hypothesis, not a second independently established driver. |

Reserves are not encoded. Their earlier expression and model scores do not establish malignant-cell origin or normal-organ safety. See [the comparison](comparison.md), [methods](../methods.md), [safety notes](../SAFETY.md) and [disclaimer](../DISCLAIMER.md).
