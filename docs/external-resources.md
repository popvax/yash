# Obtain external tools and references from their providers

## Inspect or replay our designs

No third-party source archive, trained weights or general reference library is
needed for the recorded offline replay. Clone this repository, install the
pinned dependencies, and use `yash-vaccine verify` or `yash-vaccine reproduce`.
Included NPZ/TSV caches are **our recorded prediction results, not trained model
weights**. Scientific measurements and final design outputs remain unchanged.

The [retained file catalog](../data/catalog/public-object-manifest.json) identifies
patient-data and case-analysis files in private R2 storage; see the
[download guide](../tools/r2/README.md) for current public/private availability. Third-party
software, weights and reference libraries come from their providers.
Original normal FASTQs remain excluded. Clinical reports remain in data/reports/.
Case-specific GATK orientation-model outputs are retained as patient-derived QC
evidence; they are not third-party pretrained models. Duplicate project scripts
are obtained from GitHub rather than mirrored in the data endpoint.

## New inference or raw-read discovery

Obtain the tools and resources directly from their providers under their terms.
Use the exact versions and input identities in the design methods/receipts;
another file called "hg38" or another model release is not an equivalent input.
[External file identities](../data/catalog/external-resources.json) record the
historical hashes, intended local paths and provider pointers. Only entries
marked `needed_for_recorded_upstream_stages` are planning prerequisites; unused
historical resources are not prerequisites.

| Needed for | Resource and provider | Version/scope |
|---|---|---|
| Fresh HLA-I prediction | [MHCflurry](https://github.com/openvax/mhcflurry), [native model downloads](https://openvax.github.io/mhcflurry/model_downloads.html) | Code 2.3.14; recorded affinity models 2.3.0. Processing models and training datasets are not required for our affinity-only objective. |
| DNA alignment/calling | [GATK resources](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle) and the recorded [reference preparation](../designs/design-2/pipeline/resources.py) | Assembly/reference bases, dbSNP, Mills/1000G and fingerprint resources must match the recorded workflow. Derived indices can be rebuilt; the original vendor masks remain unverified. |
| Population priors/filtering | [gnomAD downloads](https://gnomad.broadinstitute.org/downloads) | v4.1 genomes/exomes as specified in the methods; primary frequencies, not separately licensed third-party annotations. |
| Splice frames/annotations | [GENCODE 37](https://www.gencodegenes.org/human/release_37.html) and [Ensembl release 103](https://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/pep/) | Recorded GTF and protein translations. Do not silently substitute a current annotation. |
| Normal-junction filtering | [SNAF reference installation](https://github.com/frankligy/SNAF/blob/main/docs/installation.rst) | Historical GTEx/TCGA-derived processed controls; these do not reflect current/accurate GTEx Portal data. Match the recorded source hashes. |
| Fusion rediscovery | [CTAT genome libraries](https://github.com/TrinityCTAT/ctat-genome-lib-builder/wiki), [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion) | GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play. |
| CTA/normal-tissue prioritization | [LENS reference provider](https://uselens.io/) | Recorded curated CTA list dated 06APR2026 and mTEC/normal annotation dated 21APR2026. |
| HLA reference checks | [IPD-IMGT/HLA](https://www.ebi.ac.uk/ipd/imgt/hla/licence/), [recorded allele release](https://github.com/ANHIG/IMGTHLA/tree/d972ae287948dc2e5468a40ff3b3d27227541189) | Original allele references and terms; our patient read evidence is separate. |
| Upstream software | [OpenVax](https://github.com/openvax), [GATK](https://github.com/broadinstitute/gatk), [LENS](https://uselens.io/), [MSNAF](https://github.com/spvensko/msnaf) | Install from their providers. Current optional dependency pins are in pyproject.toml; execution receipts identify actual historical versions. |

For MHCflurry, use its own download manager rather than a project mirror:

```sh
python -m pip install -e '.[inference]'
mhcflurry downloads fetch models_class1_pan --release 2.3.0
mhcflurry downloads info models_class1_pan --release 2.3.0
```

Check model hashes against the recorded model-run evidence before claiming a
historical rerun. A different compatible model is a new analysis, not a replay.
Obtain these resources under their providers’ terms.

Unused training datasets, processing models and viral/HERV resources are not
prerequisites for the recorded design.

## Reconstruct the three patched RNA modules

Only minimal zero-context patches are retained in
[patches/](../workflow/stages/02_lens/patches/README.md). Full LENSTools and MSNAF
modules are absent from the checkout and source snapshot. Provider image names,
original/target hashes and extraction paths are recorded in
[upstream-sources.json](../workflow/stages/02_lens/upstream-sources.json).

After obtaining the provider images under their terms, an Apptainer/Singularity
runtime can extract the originals locally. Run from the repository root:

```sh
mkdir -p work/upstream-originals
apptainer pull work/lenstools.sif docker://spvensko/lenstools:1.9.2
apptainer pull work/msnaf.sif docker://spvensko/msnaf:1.2.1
apptainer exec --cleanenv work/lenstools.sif cat /opt/lenstools/src/metadata.py > work/upstream-originals/lenstools-metadata.py
apptainer exec --cleanenv work/msnaf.sif cat /usr/local/lib/python3.10/site-packages/msnaf/core.py > work/upstream-originals/msnaf-core.py
python workflow/stages/02_lens/prepare_upstream.py --sources-dir work/upstream-originals --output work/upstream-patched --execute
```

Preparation refuses a changed base, changed patch, existing output directory or
incorrect target hash. It does not run inference or publish the reconstructed
files. Keep those files under ignored work/, not in Git. Set
`YASH_UPSTREAM_SOURCES_DIR` to another prepared directory if desired.

See [retained-patch notices](../THIRD_PARTY_NOTICES.md) before using the provider source.
The saved-volume LENS runner extracts originals from its provider images and
applies the same verified patches before starting new work. Its dependencies and
saved-image/volume assumptions remain documented in [upstream](upstream.md).

## Verification scope

The bundled source snapshot contains project source and documentation, with
three upstream modules replaced by their original expected hashes and provider
instructions. Verification checks every bundled byte and scientific artifact;
it explicitly reports when those three upstream source files have not been
obtained and byte-checked locally. It never treats a missing source file as
verified code. The default offline design replay does not execute those modules.
