# Evidence review

- [Candidate assessment](antigen-review.md) and [TSV](antigen-review.tsv): the 16 encoded regions and eight reserves, including RNA limitations and model predictions.
- [Design 3 and Design 2 comparison](comparison.md) and [TSV](comparison.tsv): recorded computational metrics with their denominators.
- `results-summary.json.gz`: structured selection, protein, RNA, safety and AlphaGenome 2 results for both the 16-region main construct and the 15-region KCTD19-excluded alternative.
- `clinical-reconciliation.json`: structured reconciliation separating the four repository clinical PDFs from project narrative/sample metadata, longitudinal cfDNA context, pharmacogenomics and explicit limits on clinical interpretation.
- [Open evidence questions](missing-data-priorities.tsv): six bounded molecular/clinical evidence gaps (`MD-01`–`MD-06`) and what measurements can or cannot resolve, without clinical or assay prescriptions.
- [Verification receipt](verification.json) and `logs/`: current cross-platform verification commands, exit codes and logs separated from historical upstream execution notes.
- `source-verification.json.gz`: per-file comparison of bundled artifacts against upstream and earlier local research. Each record's `hash_target` states what its recorded SHA-256 covers — `file_bytes` for the stored bytes of `bundle_file`, `decompressed_bytes` for the gzipped payload actually present — because both conventions appear in this receipt; read that field rather than assuming one. [`verify_results.py`](../verify_results.py) re-checks every record under its declared convention.
- [Identity concordance](identity-concordance.json): a bounded RNA/normal comparison and the recorded upstream DNA fingerprint; neither authenticates named-patient custody.
- [Literature](literature.tsv) and `current-trial-registry.json`: references and the recorded registry snapshot. These do not establish benefit for this construct or current eligibility.
- [Data and attribution notice](data-and-license-notice.md): source identities and third-party terms.

Offline verification re-derives the prediction-dependent columns of `construct/protein/edge-costs.tsv.gz` for every edge in the table, not only the edges on the chosen path, from the recorded per-peptide prediction tables, the selected regions and the linker set. Boundary processing scores are re-derived across all 2,992 edges offline by summing the 128 recorded per-query MHCflurry 2.3.0 processing predictions in [`construct/protein/boundary-processing-scores.tsv.gz`](../construct/protein/boundary-processing-scores.tsv.gz) (814 non-zero edges and 2,178 zero-by-construction edges; see [`construct/protein/boundary-processing-receipt.json`](../construct/protein/boundary-processing-receipt.json)), and `verify_results.py` reports those exact counts under `boundary_processing_score_derivation`.

The sequences and predictions remain research candidates. Read [the safety notes](../SAFETY.md) and [contributor disclaimer](../DISCLAIMER.md) with these results.
