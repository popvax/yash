# Reference matches and screening scope

These notes record sequence and model findings. They do not demonstrate peptide production, cell-surface presentation, T-cell cross-recognition, harm or safety. No sequence or selection change is proposed.

## Junction reference matches

* **15-region alternative: `STVMAPNLF`.** This 9-mer occurs at precursor coordinates `[25,34)` (zero-based), across the annotated SEC/ARHGAP40 boundary. It matches ARHGAP40 and ARHGAP28 reference translations. Its recorded MHCflurry prediction is 753.18 nM on HLA-B*15:08, percentile 1.065%. The row is retained in the [alternative homolog table](../strand-strict-15/safety/offtarget-homologs.tsv.gz). The normal ARHGAP40 flank contains the same sequence upstream of the mutation. Processing across the signal-peptide boundary remains unmeasured; the match proves neither presentation nor its absence.
* **15-region alternative: `STVMAPNLFL`.** The overlapping 10-mer at `[25,35)` also matches ARHGAP40 and is retained as `EXACT_REFERENCE_MATCH_JUNCTION_REVIEW` in the alternative homolog table. Its shared-binding-allele field is blank: it does not meet the defined shared-HLA predicate. This does not establish absence of presentation or cross-recognition.
* **Main construct: `GAGGGRGD`.** This 8-mer occurs at `[75,83)` across the KLF5/linker/SYT12 join and matches a SEPTIN12 reference translation. Its lowest recorded MHCflurry affinity is 22,921.53 nM; all six affinity percentiles are at least 43.69%. It does not enter the model-admitted homolog query set. This differs from the alternative's model-active query, while weak prediction still does not prove non-presentation. Sources: [junction scores](../construct/all-junction-affinity.tsv.gz), [protein components](../construct/protein/components.tsv), and the [v46](../references/reference-receipt.json)/[v50](../references/v50-audit/reference-receipt.json) GENCODE translation receipts.

The fatal exact-reference exclusion is limited to nominated `TARGET` queries. Junction exact matches remain recorded for review rather than automatically causing that selection exclusion. The categories `EXACT_REFERENCE_MATCH_SELECTED_TARGET_EXCLUSION` and `EXACT_REFERENCE_MATCH_JUNCTION_REVIEW` identify those different computational scopes, not measured clinical risk. Zero nominated-target matches replay the prior exclusion rule; they are not independent proof of safety.

## Eleven UTR-frame peptide matches

The 71-nt 3′ UTR is inherited from the supplied HBA1-labelled backbone. Its hypothetical alternative frames contain eleven sense 8–11-mers with stored MHCflurry predictions ≤500 nM, including `GLPPAPPPL` at 14.95 nM. The candidates are recorded in both constructs' cryptic-ORF tables ([main](../safety/cryptic-orf-peptides.tsv), [alternative](../strand-strict-15/safety/cryptic-orf-peptides.tsv)).

All eleven match alternative annotated HBA1 translations in GENCODE v50; none matches a v46 translation. For example `GLPPAPPPL` matches HBA1-216 (`ENSP00000767388.1`, `ENST00001097582.1`), while `LPLGPPPSP` and `LGPPPSPSSP` match HBA1-219 (`ENSP00000767391.1`). These are alternative-reference sequence matches, not an observation that the canonical haemoglobin chain contains these peptides or that the alternative proteins are expressed in relevant normal tissue.

The stop-delimited enumeration does not require demonstrated translation initiation. All eleven entries lack an AUG/Methionine in their enumerated run. Non-AUG initiation or frameshifting is not thereby excluded, but neither is demonstrated. The three primary-frame stops do not specify every alternative frame, and readthrough alone does not explain entry into a different frame. No experiment here establishes production or presentation of these UTR peptides. Sources: [backbone](../construct/mrna/backbone.json), [enumeration code](../pipeline/safety_audit.py) and the recorded reference translation receipts.

## Declared classification rules

Broad homolog-query admission includes PRIME percentile ≤2%. Downstream homolog classification uses a stricter PRIME ≤0.5% branch alongside the other stated model criteria. Consequently, admission alone need not yield a downstream priority label. Neither cutoff is a measured response probability.

Shared-HLA restriction is evaluated separately: **each peptide on the same allele** must satisfy MHCflurry nM ≤500 **or** MHCflurry percentile ≤2 **or** MixMHCpred rank ≤2. PRIME does not enter that predicate. The complete policy is recorded in [the main summary](../safety/safety-summary.json) and [alternative summary](../strand-strict-15/safety/safety-summary.json).

`homolog_genes` retains all matching genes, including the query's own parent where applicable. `offtarget_gene_count` excludes the own-parent/counterpart set. The specificity plate's `human_genes` lists off-target genes when present, otherwise the own-WT counterpart genes. These fields need not have equal counts and do not measure tissue expression or immune recognition.

## Product-context limits

Internal restriction motifs matter for a cloning or linearisation workflow that uses the corresponding enzyme. An isolated `AATAAA` does not establish premature polyadenylation of a T7 in-vitro transcript; nuclear expression is a different, context-dependent question. No synonymous substitutions or cloning strategy are specified by these observations.

The main protein contains six cysteines and one N-X-S/T sequon (`NAT`, residue 366 one-based). Their presence does not establish which cysteines are free or paired, aggregation, glycosylation occupancy or epitope masking. Cap chemistry, nucleotide modifications, PolyA, delivery and manufactured-product behaviour remain unspecified or untested as described in [SAFETY.md](../SAFETY.md).
