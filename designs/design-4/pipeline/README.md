# Pipeline guide

Start with the [design overview](../README.md) and [methods](../methods.md). Run the top-level `verify_results.py` and `run.py` entry points from the design directory after installing its pinned requirements. They use local files and stored predictions.

| Module | Purpose |
|---|---|
| `curate.py` | Re-derive the candidate and core/reserve selections from frozen evidence |
| `strand_audit.py` | Inspect library orientation and mutation-supporting RNA fragments when an alignment is supplied |
| `reference_screen.py` | Search the recorded peptide set in the supplied coding references |
| `hydrophobicity.py` | Apply the 19-aa construct heuristic and report routing overlaps |
| `protein.py`, `consensus.py` | Score cached edge options and solve the fixed ten-objective ordering/linker problem |
| `build_strand_strict_contingency.py` | Assemble and check the KCTD19-excluded 481-aa protein and its known-segment PopRNA/GenBank nucleotide encoding |
| `alphagenome2_audit.py` | Verify recorded AlphaGenome 2 / Atlas tables and retained per-track `imputation_mode` rows offline, or unconditionally re-query live with `--execute` (`ALPHAGENOME_API_KEY`) |
| `immunogenicity.py`, `predict_mix_alleles.py` | Prepare or evaluate external predictor requests |
| `safety_audit.py` | Generate or verify homolog, control-panel and RNA cryptic-ORF/MFE screens for both the 16-region main construct and 15-region alternative |
| `research_qualifications.py` | Attach bounded interpretation to cached screening metadata without running new predictions |
| `prepare_mrna.py`, `recheck_construct.py`, `benchmark_mrna.py` | Prepare RNA inputs, run fresh MHCflurry junction/presentation rechecks (`--target {main,strand-strict-15,both}`) and run optional local-tool comparisons |
| `build_review.py` | Generate the candidate review, comparison, clinical reconciliation and open evidence questions from the supplied artifacts |
| `record_verification.py` | Execute cross-platform verification, unit tests, reference rescans and full dual-construct replay into a temporary directory and write `review/verification.json` + `review/logs/` |
| `package_review.py` | Check for credential-shaped content and symlinks, refresh `manifest.json` and build clean local ZIP archives |
| `replay_helpers.py`, `dataio.py` | Handle deterministic file formats, canonical paths and bounded floating-point comparisons |

Fresh predictor runs require compatible tools, weights and their applicable licences. Raw-read auditing requires the relevant BAM and index. These operations are distinct from offline verification and replay. Successful execution does not establish presentation, T-cell specificity or clinical performance.
