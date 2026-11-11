# Raw-read processing and execution tools

The lead design is complete. Start with [its guide](../README.md) for verification,
offline replay and current results. These tools are for new raw-read processing,
new model inference or maintaining the recorded Modal execution.

They require an authenticated Modal account, the recorded volume, staged inputs
and compatible images. Image IDs and dated run roots identify the original
PopVax execution; they are not a verified fresh-account deployment recipe.
Patient data can be restored using [the data guide](../../../data/README.md);
current public/private availability is recorded in data/catalog/access.json.
Original normal FASTQs are excluded; alignment from those originals requires
separate access. Offline lead replay needs no reads, Modal or model download.

Controllers print a plan unless `--execute` is supplied. Workers are not importable
libraries: some execute at import. Use the controllers, not bulk imports.

## Controllers

Run with a Python environment containing the authenticated Modal SDK. The
scientific packages and reference paths are in the recorded remote images and
volume. Each command defaults to a plan unless `--execute` is supplied.

```sh
python designs/design-2/pipeline/launch_alignment.py --execute
python designs/design-2/pipeline/run_bqsr.py --execute --wait-alignment
python designs/design-2/pipeline/continue_after_normal.py --execute
python designs/design-2/pipeline/run_normal_coding.py --execute --wait-identity
python designs/design-2/pipeline/run_downstream.py --execute
python designs/design-2/pipeline/run_construct.py --execute
```

The controllers after alignment can run concurrently on the laptop: they wait for
verified prerequisites and stage receipts. BQSR trains on primary-nuclear
20-Mb intervals and gathers native reports. Its application uses disjoint whole-
contig groups and a separate no-coordinate tail, followed by native
GatherBamFiles. A native comparison on 16,111 actual normal records produced
identical recalibrated records to whole-file application, including placed and
unplaced unmapped reads, duplicates and supplementary records. Final gathering
must preserve all per-contig mapped/unmapped and no-coordinate record counts.
The initial whole-file application was deliberately superseded after its
observed rate implied about 4.7 hours; its source and stop record are retained.
Training and the gathered recalibration table were not repeated.

`run_bqsr.py` defaults to this checked parallel application. To resume after
gathering without repeating training:

```sh
python designs/design-2/pipeline/run_bqsr_scatter.py --execute
```

The controller reattaches matching successful/live workers. A stopped failed
worker requires explicit inspection; an observation timeout is not a reason to
restart it. The optional `run_bqsr.py --application whole` preserves the original
unscattered route for a new run; do not run both writers simultaneously.

An optional early identity check can run on the indexed, duplicate-marked BAM
while calibration trains:

```sh
python designs/design-2/pipeline/run_stage.py identity_precal --execute
```

Its result is labelled pre-BQSR and does not classify mutations. The required
identity/contamination stage repeats fingerprinting on the calibrated BAM.

The Broad fingerprint map has a 3,366-contig dictionary starting with chr1;
both case BAMs have the same 152-contig dictionary starting with chrM.
`fingerprint_panel.py` verifies the two BAM dictionaries agree, checks each
marker coordinate and reference allele, and writes a separate map with their
shared dictionary. Any incompatible marker excludes its entire haplotype block.
Marker alleles, positions, frequencies and anchor relationships are preserved;
the original map and BAM headers stay intact. Native dictionary validation
remains enabled. The coordinate audit records every retained/excluded block.

After BQSR, `continue_after_normal.py` performs native fingerprinting and
contamination estimation, launches the explicit 17-target audit, and independently
runs full paired discovery and native filtering. The audit does not nominate
variants for that discovery. Every newly emitted ALT is reviewed against the
cached full-chromosome gnomAD genome/exome references before native filtering;
sample fields always name the tumour and normal explicitly.

`run_normal_coding.py` generates normal coding-region GVCFs for later normal
callability and patient-specific self-sequence review. Its scope is primary
nuclear CDS plus 100-bp flanks. It is not a complete normal transcriptome or an
inherited-disease interpretation. Conflicting scatter measurements are held for
review, not selected by their most favorable score.

`run_downstream.py` waits for actual successful native normal-genotype and paired
filtering jobs, then runs the normal-sequence library, annotation, fresh Isovar
reconstruction and native Vaxrank/MHCflurry affinity ranking. A source file being
present is not a result; each stage requires a successful remote exit and hashed
completion receipt. It then rechecks the recorded RNA-source pool and combines
it with the new mutation regions. This step does not relabel reused LENS
discovery or splice affinities as a fresh run. Fusion affinities are re-inferred
and checked against the original contributions; unchanged splice and CTA
measurements keep their recorded provenance.

`run_construct.py` waits for that completed selection. It submits fresh A10
affinity inference for all internal and SEC/MITD junction/linker options, runs
the recorded joint path/linker optimizer and hydropathy gate, scans the full
protein independently, and submits a separate fresh prediction of the chosen
junctions. It then re-encodes the fixed protein with the recorded PopRNA settings
and runs the maintained sequence/GenBank verifier. PolyA and any downstream
extension stay unspecified. The selection, final junction verification and RNA receipts have all completed and passed their recorded hard gates.


## Tool map

- `launch_alignment.py`, `run_bqsr.py`, `run_bqsr_scatter.py`: normal alignment/calibration.
- `continue_after_normal.py`, `run_paired.py`, `run_normal_coding.py`: identity, paired discovery and normal genotyping.
- `run_downstream.py`, `run_construct.py`: RNA/self review, selection, new inference and construct checks.
- `run_stage.py`: an explicit single-stage controller.
- `collect_results.py`: collect completed results, preserving canonical local paths.
- `review_selected.py`: regenerate the human-readable review from current evidence.
- `status.py`: inspect locally recorded Modal handles. A fresh clone has no live-job state; use `yash-vaccine verify` for the completed design.
- `inputs/`: recorded reference dictionaries, HLA models and nominated-target audit inputs.

Successful receipts remain under [execution/](../execution/README.md). Historical
source bytes are preserved in [the source snapshot](../../../provenance/snapshots/README.md).
