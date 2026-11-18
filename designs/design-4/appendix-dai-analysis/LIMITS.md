# Interpretation and verification limits

DAI compares two predicted nM values for a specified peptide pair and HLA allele. **DAI =1 is the mathematical boundary between predicted increases and decreases in binding.** Neither that boundary nor any other DAI threshold is established here as an antigen-selection rule or clinical decision boundary. The DAI ≤1 and mutant ≤500 nM filters describe table subsets. Both filters use unrounded source-derived values; displayed DAI values are rounded to four decimal places.

## Biological scope

* DAI does not measure T-cell recognition, tolerance, processing, presentation, safety or benefit. A higher ratio does not establish a better target; a lower ratio does not establish a harmful or ineffective one.
* These are **Class-I, mutation-spanning missense comparisons**. The 3,024 rows do not cover all peptides encoded by the construct. Unchanged flanks have identical WT counterparts and a ratio of 1 under the same model and allele.
* TGFBR2 frameshift and CDH1::NRG1 fusion regions have no applicable single-amino-acid WT-reversion comparison in this analysis. Their DAI fields are blank and `dai_defined=False`. This does not mean every peptide in those windows is absent from normal proteins.
* Class-II predictions are not recomputed by this appendix. It makes no new Class-II comparison or processing claim; the original Design 4 analyses remain outside this appendix's verification scope.
* The six Class-I alleles are A*01:01, A*33:03, B*08:01, B*15:08, C*01:02 and C*07:02. The [Caris report, physical page 3](../../../data/reports/caris-tissue-2025-02-22.pdf#page=3), calls these predicted results and states that tumour sequencing is not a substitute for peripheral-blood histocompatibility testing. Agreement with research predictions does not replace that confirmation.
* The KCTD19 RNA counts describe the analysed historical specimen, not every tumour cell or present-day lesion. Their source and limitations are given in [FINDINGS.md](FINDINGS.md).

Assessing recognition requires experimental evidence under stated conditions and detection limits. A peptide-response assay could address recognition under its assay conditions; it would not establish natural tumour presentation or clinical safety. A suitably validated strand-specific RNA assay could address KCTD19 coding-transcript detection in a tested sample. No particular assay outcome is predicted here.

## Computational scope

Ordinary verification hashes these five Design 4 inputs:

1. `selection/selected.tsv`
2. `immunogenicity/class1-counterfactual-contrasts.tsv`
3. `immunogenicity/mixmhcpred3.tsv`
4. `safety/offtarget-homologs.tsv.gz`
5. `construct/all-junction-affinity.tsv.gz`

It compares derived tables against those inputs, checks finite numeric values, required rows and keys, and derives summary counts from numeric cells.

The 46 paired MHCflurry affinity cells are corroborated on exact `(peptide, allele)` keys: 26 use the junction cache and 20 use `shared_allele_homolog_nm` from the homolog table. For the frozen submitted table, the largest absolute difference is **0.00548733 nM**, within the declared 0.05 nM tolerance. The verifier reruns this lookup and checks the saved [corroboration receipt](data/cell-corroboration.json). This is agreement between computational outputs, not independent biological validation. Stored-source corroboration of nM cells does not by itself verify percentile cells.

`--recompute` additionally loads the explicitly selected affinity-model directory, checks its inference-file digests against the recorded assets, and repeats affinity and percentile inference. `--models-dir` allows a local path while retaining the digest requirement. Training archives are excluded from the asset hash, as recorded. MixMHCpred is not rerun; PRIME is not consulted for the shared-HLA comparison.

A passing verifier establishes only the checks it performs. It does not authenticate other Design 4 files, prove the original predictor's biological accuracy, verify every historical review narrative, or demonstrate that a manufactured product is safe or effective.

Produced by an AI agent under non-expert human guidance. This research analysis is not medical advice or a recommendation to manufacture or administer a product. It establishes no safety, immunogenicity or clinical benefit.
