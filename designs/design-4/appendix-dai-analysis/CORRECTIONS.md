# Correction record

The current submission is an arithmetic research appendix. An earlier proposal to replace Design 4 with a 14-region construct under a mandatory DAI ≥3 rule was withdrawn. It is not part of the current proposal.

## Withdrawn redesign

The withdrawn package had inconsistent DNA and RNA stop-codon specifications, omitted the fusion from lead-specific homolog/control analyses, and failed to inherit its intended MixMHCpred/PRIME values. Its passing checks did not detect those defects. The artifact inconsistency meant the package was not a verified product specification; it did not demonstrate biological nonfunction.

Its manifest also had two distinct write-order defects: a self-entry describing the previous manifest, and a software receipt written after hashing. The latter becomes stale only when the rewritten receipt's contents change; an identical-byte rerun retains the correct digest.

## Interpretation corrections

* DAI describes a ratio of predictions, not TCR discrimination, immune tolerance or response probability. The proposed DAI ≥3 selection rule and resulting redesign were withdrawn.
* Unchanged flanking peptides have identical WT counterparts. The DAI table covers mutation-spanning single-amino-acid reversion comparisons, not all encoded peptides. KCTD19 and FES windows are 23 and 15 aa respectively.
* The original TIER1 rule includes absolute MHCflurry affinity, MHCflurry percentile and MixMHCpred rank as alternative branches. Applying only the absolute-nM branch does not invalidate a flag. Earlier “unoccupied,” “adverse” and “favourable” interpretations were unsupported and are withdrawn.
* Affinity and its calibrated percentile are not independent models. The current findings distinguish their different cutoffs from the two relative-order reversals between MHCflurry and MixMHCpred.

## Verification corrections

Earlier checks incompletely validated displayed values, non-finite numbers, duplicate or missing rows, undefined comparisons, source changes and derived receipts. Model pinning initially hashed a different directory from the one implicitly loaded. Current verification derives table expectations from the five named inputs, reruns stored-cell corroboration, and explicitly loads the affinity directory whose asset digests it checks. See [LIMITS.md](LIMITS.md) for the exact scope.

## Historical source facts

* Reference checks corrected the baseline's RXRA `NM_002957.6:c.640G>T` citation to `c.721G>T` and separated Guardant's DNMT3A C351Afs*56 label from the versioned-reference p.Phe350AlafsTer56 consequence. The companion annotation correction preserves genomic alleles, peptides and all numerical predictions. See [source provenance](PROVENANCE.md).
* Guardant reports **31.9% tumour fraction for 14 October 2025**, not present-day disease burden. It reports CDH1–NRG1 fusion values of **0.4% and 9.3%** for 29 March and 14 October 2025. Source: [Guardant report, physical page 3](../../../data/reports/guardant360-2025-10-21.pdf#page=3).
* The repository's project timeline records afatinib starting in November 2025 and rash limiting its use, citing the patient's account on the project website. It does not establish continuing treatment or suppression of an ERBB4 clone. Source: [case timeline](../../../data/clinical-timeline.md).
* The prior “10 of 12 off-panel mutation genes” count concerned the **withdrawn 14-region redesign**. It is not a count for the unchanged 16-region Design 4. Panel membership does not establish mutation presence or absence. Source: [Guardant panel, physical page 9](../../../data/reports/guardant360-2025-10-21.pdf#page=9).
* H2BC12's omission from the vendor expression table does not establish absent expression or explain assay chemistry. Other histone genes occur in that table. Source: [vendor expression export](../exports/expression.tsv.gz).

These corrections preserve the distinction between repository evidence, computational predictions and unmeasured biology. No clinical improvement, safety or cure has been demonstrated by the appendix.

Reference reporting was refreshed to retain complete homolog-gene lists and distinguish selected-target exclusion from junction reference-match review. The inherited gene annotations and corresponding digests were updated. All numerical prediction values and DAI findings are unchanged; this metadata refresh does not claim new predictor inference.

The source homolog table now distinguishes unavailable biopsy-expression values from recorded zero and retains explicit gene coverage. Its pinned digest was refreshed after checking that all peptide, HLA, model-score and priority fields were unchanged. This appendix's numerical tables and findings are unchanged.
