# Matched-normal reassessment

This analysis uses the new white-blood-cell DNA with the original tumour DNA and RNA. The recorded initial design remains unchanged for comparison.

The normal and tumour fingerprints matched. Within the explicitly requested normal coding-plus-flank intervals, 73,917,616 of 76,021,670 bases (97.2%) had native genotype quality at least 20 and depth at least 20. This is not a whole-genome callability percentage.

## The 17 previously selected mutations

| Gene | Protein change | Normal REF / ALT fragments | Reassessment | Selected now |
| --- | --- | --- | --- | --- |
| KLF5 | D418N | 50 / 0 | SUPPORTED_PAIRED_SOMATIC | True |
| SF1 | R218Q | 57 / 0 | SUPPORTED_PAIRED_SOMATIC | True |
| KMT2D | G3189R | 35 / 26 | PRESENT_IN_NORMAL | False |
| RXRA | A241S | 42 / 0 | SUPPORTED_PAIRED_SOMATIC | True |
| HSPG2 | R3952L | 46 / 0 | SUPPORTED_PAIRED_SOMATIC | True |
| PARP3 | K37Q | 29 / 23 | PRESENT_IN_NORMAL | False |
| CAST | I309V | 15 / 22 | PRESENT_IN_NORMAL | False |
| CCNI | P285S | 23 / 17 | PRESENT_IN_NORMAL | False |
| PPIC | P187S | 47 / 0 | SUPPORTED_PAIRED_SOMATIC | True |
| PPP2R1A | K280N | 22 / 26 | PRESENT_IN_NORMAL | False |
| KCTD19 | G920R | 49 / 0 | SUPPORTED_PAIRED_SOMATIC | True |
| VPS11 | D151A | 32 / 0 | SUPPORTED_PAIRED_SOMATIC | True |
| ARHGAP40 | R438W | 49 / 0 | SUPPORTED_PAIRED_SOMATIC | True |
| NRXN2 | D203H | 33 / 23 | PRESENT_IN_NORMAL | False |
| PLEC | R2171W | 26 / 0 | SUPPORTED_PAIRED_SOMATIC | True |
| DOCK8 | G1251R | 18 / 17 | PRESENT_IN_NORMAL | False |
| SYT12 | P391L | 56 / 0 | SUPPORTED_PAIRED_SOMATIC | True |

Missing normal coverage is unresolved. A missing exact allele in independent discovery is not labelled a native filtering failure. Blood-cell presence is not automatically labelled inherited.

## Independent discovery and RNA reconstruction

New paired Mutect2 discovery covered the planned primary nuclear chromosome intervals without using Caris calls or the old 17 mutations as a nomination list. Native sample-wide filtering used the matched contamination/segmentation estimates, read-orientation evidence and a 5% model-FDR target; this is not an empirically measured false-discovery rate.

Annotation retained 99 primary DNA candidates. Isovar reconstructed 17 usable leading protein contexts; 12 mutation regions passed the explicit primary RNA/expression/class-I ranking gates before construct hydropathy selection.

## Selected design and checks

The new selection contains 24 regions: 12 MUTATION, 5 SPLICE, 5 CTA/SELF, 2 FUSION.
All passing primary mutation regions were considered; the mutation branch has no fixed top-N quota. The recorded RNA-source discovery and normal-tissue splice screen were reused with their provenance. Altered-target scoring also excludes possible matches in the new normal coding/HLA search models. CTA regions remain intentional unmutated self targets.

Every internal and terminal linker/order option was freshly predicted with the recorded class-I affinity model. The joint optimizer, independent full-chain junction scan, separate fresh selected-junction predictions and hydropathy checks passed. The fixed protein was re-encoded with the recorded PopRNA settings; translation and GenBank annotations were verified. PolyA and the downstream extension remain unspecified, so the exported nucleotide segments are partial.

## Practical limits

Supported paired somatic status is evidence from these samples and models, not a guarantee of healthy-tissue absence. Normal coding models include possible local unphased combinations and explicitly retained unresolved cases; no search establishes a complete personal proteome. Normal DNA does not measure normal RNA splicing, protein translation or HLA presentation. Affinity predictions and synonymous encoding heuristics do not establish immunogenicity, efficacy, topology, expression or manufacturing suitability.

The normal sequencing instrument is unconfirmed. BaseRecalibrator uses its native --default-platform UNKNOWN setting, preserving the original alignment header. The original Caris reference with its custom masks is unavailable; the normal was aligned to the same compatible GRCh38 analysis reference used for this reassessment.

See summary.json for native receipts and callability scope, current-17-final-audit.json for full focal evidence, and the integrate comparison/selection tables for every added, removed and retained region. No website or public repository was changed by this analysis.
