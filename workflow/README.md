# Upstream execution

Start with [usage](../docs/usage.md) for local reproduction and independent designs.
The lead matched-normal controllers are in [Design #2 pipeline](../designs/design-2/pipeline/README.md).
They implement paired discovery, normal/self review and fresh downstream inference.

The source-specific stages below record Design #1 discovery and shared RNA-source
evidence reused by the lead. The CLI `upstream dna/rna` commands address these
earlier tumour-only adapters, not the paired branch.

| Stage | Code | Output |
|---|---|---|
| Independent DNA discovery and segmented filtering | stages/01_mutation_calling/discovery/ | Native calls, population evidence and DNA eligibility |
| RNA-supported mutation reconstruction and scoring | stages/01_mutation_calling/discovery/review/ | Isovar contexts and Vaxrank affinity-ranked regions |
| RNA-only LENS discovery and review | stages/02_lens/ | CTA/self, fusion and pre-HLA splice evidence |
| CTA/self and fusion region scoring | stages/04_region_integration/ | Reference-filtered windows and affinity contributions |
| Actual patient-RNA splice reconstruction | stages/04_region_integration/patient_splice/ | Observed haplotypes, upstream frames, QC, affinity and gene ranking |

[Stage requirements](../docs/upstream.md) describe inputs and runtime contracts.
Controllers print plans unless --execute is given. Some scientific workers run
at import; inspect their source without bulk-importing them.

The complete initial design is under designs/initial/. Its local reproduction
uses the bundled measurements and predictions, without cloud jobs or downloads.
