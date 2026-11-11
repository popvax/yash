# Shared data and evidence

This directory contains evidence that different vaccine designs can reuse.
A design's choices and cassette live under designs/, alongside its recipe.

## Distribution and acknowledgements

Large sequencing files and case-analysis outputs are stored in a private R2
bucket. Public read-only downloads need no credentials; check the
[current status and endpoint](catalog/access.json).
Third-party software, trained weights and general reference libraries are
[obtained directly upstream](../docs/external-resources.md). Full RNA-tool source
modules are not included in the checkout or snapshot; only minimal patches and
hashes are retained. Prediction-result caches are not model weights.

[Required credits](../ACKNOWLEDGEMENTS.md) and
[retained-source notices](../THIRD_PARTY_NOTICES.md) identify the applicable sources.

## Included with the repository

The bundled mutation-evidence/checkpoint tables record the earlier **Design #1**
tumour-only workflow. The lead **Design #2** uses paired tumour–normal evidence in
[its exports](../designs/design-2/exports/) and the
[matched-normal data catalog](matched-normal/README.md). Shared RNA-source
measurements retain their original provenance.

[Clinical reports](reports/) contains the Stanford liver pathology/IHC
and STAMP genetic panel, Caris tumour profiling and Guardant360 liquid biopsy
PDFs. They are included when you clone the repository.

[data/checkpoints](checkpoints/) contains the compact evidence needed for local
selection and construction. Its [manifest](checkpoints/manifest.json) records
file sizes, SHA-256 hashes and the derivation of extracted inputs.

| Evidence | Used for |
|---|---|
| mutation-regions.tsv | Mutation-region sequences, support, annotations and scores; common pool rows are derived directly from this table |
| lens-scored-windows.tsv.gz | LENS CTA/self and fusion windows and scores |
| target-occurrences.tsv.gz | Independent recount of window contributions |
| splice-evidence/ | Actual patient RNA haplotypes, annotated upstream frames, native review, fresh affinities, reference checks and rankings |
| construct-predictions.npz | Freshly computed internal-junction affinity predictions for the selected regions |
| splice-evidence/normal-screen/ | Coordinate-aware normal measurements, ND-only admission, strand checks and rankings |
| fusion-evidence/ | In-frame driver-junction sequence, binding predictions and explicit admission policy |

Tables are tab-separated text; .gz files are compressed text. The NPZ is the
large numerical prediction cache. Its model and environment records are in
[checkpoint provenance](checkpoints/provenance/). Recorded inputs for the
independent mutation stages are in [mutation-evidence](mutation-evidence/README.md).

[Patient-RNA splice evidence](splice-evidence/README.md) includes all normal-filtered
input events, outcomes, full-window sequence/path witnesses, stop-free upstream
coding anchors, native QC, exact-self matches and fresh six-allele predictions.
Local replay independently recounts all recorded scores and checks the separate
BAM-verification receipt. It needs no raw-read access or model download.

The [mutation evidence](mutation-evidence/README.md) includes segmented DNA decisions,
RNA reconstruction, all-window scores, exact-reference checks and native HLA-I
predictions. The Design #1 selection replay independently rebuilds its 18-row
mutation checkpoint before construct-compatibility filtering. The default replay
targets the lead Design #2 and uses its matched-normal evidence and
[integrated candidate pool](../designs/design-2/exports/integrated-candidate-pool.tsv).

## Matched-normal data for the lead design

[Matched-normal DNA and Design #2 outputs](matched-normal/README.md) provides the
final normal BAM/index, coding-plus-flank GVCF, paired tumour–normal calls,
QC and supporting analysis files. The original normal FASTQs are excluded.
[Design #2](../designs/design-2/README.md) is the lead and uses matched-normal DNA;
[Design #1](../designs/initial/README.md) was generated without matched-normal DNA.

## Stored externally

The [file catalog](catalog/README.md) lists the externally stored patient/results data.
The [retained manifest](catalog/public-object-manifest.json) records locations and checksums;
[access.json](catalog/access.json) configures the read-only downloader.

```sh
yash-vaccine data list
yash-vaccine data restore --category raw --destination ../yash-data
```

The first command lists sequencing inputs. The second prints a download plan;
add --execute to download. Neither command starts a compute job.

Downloads need no credentials. Authorized maintainers can also use private
S3 access. See the
[download instructions](../tools/r2/README.md). Offline replay needs no external data access.

Upstream outputs become input evidence for downstream design. Selecting a region
or choosing a linker is recorded in a named design, rather than changing this
shared evidence. New upstream analyses need their own provenance and checksums.
