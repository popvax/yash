# Design 4: a computational vaccine research candidate

Design 4 is an alternative construct based on the matched-normal patient data used in [Design 2](../design-2/README.md) and PopVax's lead [Design 3](../design-3/README.md). It retains DNA-supported mutation candidates and an in-frame CDH1::NRG1 fusion window, checks strand-specific RNA evidence, and optimizes antigen order and linkers using recorded prediction scores. Tumour presentation, safety and clinical benefit have not been demonstrated for this construct.

Read the [AI disclosure and contributor disclaimer](DISCLAIMER.md).

The main protein contains 16 regions: 15 mutation windows and one fusion window. It is 507 amino acids long (`1,521 nt` sense coding; `1,664 nt` known RNA segment before unspecified PolyA), including the routing domains and linkers. Twelve mutation windows have at least two coding-strand mutant RNA fragments; FES and SYT12 have one each. KCTD19 has somatic DNA support (`AF 0.318` in liver biopsy) but no coding-strand mutant RNA fragments in the audited sample (`0` coding / `2` opposite-strand); its tumour expression and clonality remain unestablished. A separate 15-region alternative (`strand-strict-15/`) removes KCTD19, reoptimizes the remaining order, linkers and known-segment PopRNA/GenBank nucleotide encoding (`481 aa`; `1,443 nt` sense coding; `1,586 nt` known RNA segment before unspecified PolyA), and includes its own fresh junction/presentation recheck and homolog/cryptic-ORF/control-panel screens (`strand-strict-15/safety/`). Neither option is established as clinically safer or more effective.

Relative to Design 2 (`24` regions, `724 aa`) and Design 3 (`30` regions, `885 aa`), the main construct (`16` regions, `507 aa`) retains the 12 paired-DNA mutation windows from Design 2 plus `H2BC12` (also included in Design 3), `TJP2` and the observed 15-aa `FES` window, alongside the in-frame `CDH1::NRG1` fusion window. It omits the five reference-derived Guardant360 VUS hypotheses added in Design 3 (`AR P555S`, `DHX36 R160K`, `ERBB4 D113V`, `KMT2D R2789Q`, `RAF1 G535S`), the five RNA-only splice windows, the five intentional-self/CTA windows and the second `CDH1::NRG1` fusion hypothesis. Eight additional candidates remain in reserve. These choices trade target breadth against evidential uncertainty; they do not establish a better treatment.

## Reading guide

**How this contribution fits together:** Design 4 proposes an alternative antigen selection and construct using the supplied patient evidence. The accompanying [DAI appendix](appendix-dai-analysis/README.md) examines predicted mutant-versus-wild-type binding within that design. It proposes no selection or sequence change and provides no new evidence of clinical safety or benefit. PopVax designates Design 3 as its lead candidate (`designs/design-3/`), succeeding the earlier matched-normal Design 2 (`designs/design-2/`). Design 4 is an independent computational alternative; its relative clinical safety and effectiveness are unestablished.

The main Design 4 construct contains 16 regions. The separate 15-region alternative excludes KCTD19; it is part of Design 4, rather than a proposal from the DAI appendix.

Read in this order:

1. **Overview:** the opening paragraphs above explain the proposal and its two constructs.
2. **[Comparison with Design 3 and Design 2](review/comparison.md):** what changes and the limits of those comparisons. Exclusion means a candidate was left out under the stated selection rules; it does not establish harm or lack of benefit.
3. **[Candidate evidence](review/antigen-review.md) and [safety questions](SAFETY.md):** the support and unresolved concerns for the selected regions.
4. **[Optional DAI appendix](appendix-dai-analysis/README.md):** an additional analysis of existing predictions. A higher DAI does not establish a better vaccine target.
5. **[Verification and replay](#verify-and-replay):** what the computer checks reproduce and what they do not establish.

## Recorded results

| Junction prediction | Main construct (`507 aa` / `1,521 nt` CDS) | Design 4 without KCTD19 (`481 aa` / `1,443 nt` CDS) |
|---|---:|---:|
| Class-I peptide–HLA occurrences | 4,380 | 4,104 |
| MHCflurry affinity percentile ≤0.5% | 0 | 0 |
| MixMHCpred rank ≤0.5% | 0 | 0 |
| Conditional Class-II synthetic-core rank ≤0.5% | 0 | 0 |
| PRIME rank ≤0.5% | 1 | 1 |
| MHCflurry affinity percentile ≤2% | 20 | 21 |
| MHCflurry predicted affinity ≤500 nM | 16 | 14 |
| Conditional Class-II synthetic-core rank ≤2% | 12 | 11 |

These are prediction counts under defined thresholds, not measured epitopes or immune responses. Target queries are retained after exact-reference exclusion against GENCODE v46 and v50 provider translations (recorded with EBI download URLs, SHA-256 digests and exact-match outputs in [`references/`](references/)) and the bundled constructed patient-normal coding/HLA reference (`references/patient-normal-and-hla-contexts.fasta.gz`). Reproducing zero matches checks that selection condition; it is not an independent biological safety endpoint. The encoded normal flanks and routing sequences are not asserted to be absent from normal proteins. The proposed specificity control panels (`197` Class-I and `39` Class-II peptides in `safety/` for the 16-region construct; `184` Class-I and `37` Class-II peptides in `strand-strict-15/safety/` for the 15-region alternative) cover specified selection rules; they do not exhaust possible T-cell cross-recognition.

See the [candidate assessment](review/antigen-review.md), [Design 3 and Design 2 comparison](review/comparison.md), [open evidence questions](review/missing-data-priorities.tsv), [methods](methods.md) and [unresolved safety questions](SAFETY.md).

## Verify and replay

Use Python 3.11 in a virtual environment. From this directory:

```sh
python -m pip install -r requirements.txt
python verify_results.py --rescan-reference
python -m unittest discover -s tests -p 'test_*.py'
python run.py --output /absolute/path/to/a/new/replay-directory
```

This revision passed all unit tests, the three reference rescans and full replay of both the 16-region main construct and the 15-region KCTD19-excluded alternative in Python 3.11 (recorded receipt environment Python 3.11.18 in review/verification.json; verified compatible on Python 3.11.2). All ten optimization objectives reproduced for both constructs. Protein and nucleotide sequences and the checked annotation files matched byte for byte; diagnostic floats stayed within the `1e-12` tolerance. Replay uses stored predictions and does not repeat model inference, raw-read processing or experiments. Optional live AlphaGenome 2 (`alphagenome==2.0.0`) re-execution is documented in `pipeline/alphagenome2_audit.py --execute` (requiring `ALPHAGENOME_API_KEY` in the environment; no credentials are stored), and optional raw-BAM/MHCflurry locus recomputation for the 13 Guardant360 variants is provided in `tools/audit_ctdna_loci.py --execute`. The [verification receipt](review/verification.json) separates historical upstream execution notes from the current executed verification logs (`review/logs/`); later edits require their own verification.

The tumour specimen dates to 15 January 2025 (liver metastasis Block A2, ~60–70% cellularity). Longitudinal liquid biopsy on 14 October 2025 reported 31.9% epigenomic ctDNA tumour fraction alongside clinical progression in lung and peritoneum while liver and pancreatic lesions remained controlled on zenocutuzumab; 13 of the 15 Design 4 somatic mutation genes are absent from that targeted cfDNA panel (`review/clinical-reconciliation.json`), so persistence in unbiopsied lung/peritoneal lesions remains unmeasured (`review/missing-data-priorities.tsv`). The project reports matched-normal blood collection on 30 October 2026. Identity checks support donor concordance but do not authenticate named-patient chain of custody; see [data provenance](PROVENANCE.md). Cap chemistry, nucleotide modifications, PolyA, delivery, dose and manufacturing remain unspecified; the nucleotide files describe known segments of an incomplete product.

Source permissions and attribution are described in the [data and attribution notice](review/data-and-license-notice.md).
