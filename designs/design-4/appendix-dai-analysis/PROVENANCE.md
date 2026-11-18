# Source provenance

The source snapshot is PopVax's [`yash` repository at commit `7ece8f7`](https://github.com/popvax/yash/tree/7ece8f7a77c8db935422f6502d2c4db2352d7bf0). The project attributes case `TN25-112103` to Yash. Repository attribution and molecular concordance do not independently authenticate specimen custody.

## Patient evidence and derived analysis

The five files in [the input-digest receipt](data/design4-input-digests.json) are **derived Design 4 research outputs**, not raw patient measurements or byte-identical files published in that upstream snapshot. They contain selection decisions, sequence contexts and model predictions. Binding predictions are not findings printed in the clinical PDFs.

The separate source audit matched all 20 upstream foundation exports recorded in Design 4's [source verification](../review/source-verification.json.gz) to the pinned Git objects. All four original [clinical PDFs](https://github.com/popvax/yash/tree/7ece8f7a77c8db935422f6502d2c4db2352d7bf0/data/reports) matched those Git objects and the SHA-256 values in the pinned [file catalog](https://github.com/popvax/yash/blob/7ece8f7a77c8db935422f6502d2c4db2352d7bf0/data/files.tsv). The downloaded vendor VCF, paired variant VCF, normal coding gVCF and complete RNA BAM also matched their catalog checksums.

All 15 selected somatic genomic alleles matched the supplied paired VCF and upstream evidence exports. The fusion is supported separately by the supplied fusion evidence and clinical report. These comparisons do not repeat alignment or variant calling from the full DNA datasets.

Reference-context qualifications matter:

* **ARHGAP40:** the inherited `p.R438W` label uses the older `ENST00000373345.8` transcript. The same local WT window occurs at R491 in `ENST00000373345.9`. The genomic allele is `chr20:38643812 C>T`; the local mutant/WT peptide comparison is unchanged by that renumbering.
* **PPIC:** the single-residue somatic reversion preserves a nearby N190S background change. The matched-normal gVCF calls `chr5:123023945 T>C` homozygous alternate. Thus this counterfactual need not equal the canonical-reference sequence.
* **RXRA:** the local WT window occurs at A214 in the shorter `ENST00000672570` protein and A241 in `ENST00000481739.2`. The audited baseline used the incorrect `NM_002957.6:c.640G>T` citation; the correction is `c.721G>T` (GCC to TCC, A241S). This reporting-label correction changes no genomic allele or peptide in this appendix. The clinical report itself prints A214S; it does not provide a nucleotide-level plasma sequence for independent allele authentication.

The companion baseline correction also separates Guardant's DNMT3A C351Afs*56 label from the supplied deletion's versioned-reference consequence, `NM_022552.5:c.1048_1051del`, p.Phe350AlafsTer56. DNMT3A is outside this appendix's selected regions and remains excluded from targeting.

## Clinical annotations

The `guardant360_covered` column is a curated historical annotation: RXRA, TGFBR2 and CDH1::NRG1 were detected in the 14 October 2025 sample. Its sources are the [Guardant fusion results, physical page 3](../../../data/reports/guardant360-2025-10-21.pdf#page=3), [RXRA/TGFBR2 results, physical page 4](../../../data/reports/guardant360-2025-10-21.pdf#page=4), and [panel list, physical page 9](../../../data/reports/guardant360-2025-10-21.pdf#page=9). It is not calculated from the five prediction inputs and says nothing about present-day detection. The HLA predictions and their clinical-typing limitation appear on [Caris physical page 3](../../../data/reports/caris-tissue-2025-02-22.pdf#page=3).

False values in this annotation do not establish that an alteration is absent; most selected mutation genes are outside that report's panel. Protein-label correspondence across isoforms is distinct from direct comparison of raw plasma sequencing reads, which were not supplied.

A fresh bounded comparison of the supplied RNA BAM and normal coding gVCF reproduced 93 exact threshold-genotype matches among 110 selected markers, ten uncertain RNA calls and seven normal-heterozygote/RNA-homozygote calls. All eight normal homozygous-alternate markers matched in RNA; there were no opposite-homozygote conflicts. This supports concordance, not a calibrated identity probability or exhaustive contamination testing. See Design 4's [provenance limits](../PROVENANCE.md).

The full tumour and normal DNA BAM payloads were not independently rehashed in this review. The upstream DNA fingerprint is a recorded source result, not a newly repeated fingerprint test. None of these provenance checks establishes tumour presentation, safe administration or clinical benefit.
