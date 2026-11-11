# Tumour-only germline-model sensitivity

These diagnostics assess the 18 upstream mutation candidates, before the separate hydropathy screen. They do not
replace the segmented selection policy or establish inherited/somatic status.

[Mutation measurements](mutation-review.tsv) retain exact DNA/RNA support,
population counts and frequencies, the prior actually supplied to GATK and the
regional allele balance. [Sensitivity results](mutation-sensitivity.tsv) separate
native FILTER from the additional GERMQ ≥20 research gate.

Both native GATK 4.6.2.0 FilterMutectCalls runs used all 879,912 original raw
records, original contamination/orientation resources and the same 5% model-FDR
strategy. They did not recall variants or nominate calls from the Caris VCF.

1. Retain segmentation and raise eligible loci's population priors, where higher,
   to the largest recorded QC-PASS global/group frequency. This changes 35 priors
   among 105 matched loci. No patient ancestry is inferred from location.
2. Retain the original global priors and omit optional tumour segmentation.
   Native GATK then uses a 0.5 heterozygote allele-balance expectation. This is
   a diagnostic assumption, not a decision to replace the selected model.

All 18 pre-screen calls retain native PASS with higher priors, but CAST falls
below the additional GERMQ gate. Without segmentation, KMT2D, CCNI and PTCH1
receive native germline filters; HSPG2, PARP3 and CAST also fail the extra GERMQ
gate. KLF5, SF1, RXRA, PPIC, PPP2R1A, KCTD19, VPS11, ARHGAP40, PLEC, DOCK8 and SYT12 pass the research gates in both diagnostics. NRXN2 is native PASS without segmentation but has GERMQ 19, below the additional gate.

The selected global population priors agree with the recorded population resource.
Rare inherited alternatives nevertheless remain possible. Regional allele-balance
estimation can make the germline likelihood sensitive to the chosen segmentation,
especially at very high depth. Neither diagnostic is ground truth for a locus.

Exact commands, edited priors, per-call results and native filtering statistics
are retained alongside this README. The [normal confirmation list](normal-confirmation-targets.tsv)
provides exact GRCh38 alleles for future matched-normal checks. The current
[selection policy](../policy.json) continues to use segmented filtering.

The native per-call JSON retains its original run-time `selected` flags; those flags are not admission decisions for the current design. Current tables join all 18 by stable variant ID.

The construct-compatibility policy excludes PTCH1 from assembly. Its somatic-model
measurements remain here; 17 mutation regions enter the initial design.
