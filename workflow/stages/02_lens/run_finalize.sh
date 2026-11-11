#!/usr/bin/env bash
set -euo pipefail
ROOT=${LENS_RUN_ROOT:-/data/results/TN25-112103/lens-expression-aware-20261030}
WORK=/work/lens-expression-aware
export APPTAINER_BINDPATH=/data,/work
singularity exec --cleanenv "$WORK/containers/lenstools.sif" python3 /opt/lenstools/lenstools.py filter-high-confidence \
 -i "$ROOT/prediction/lens-self-filtered.tsv" -o "$ROOT/prediction/native-source-filter.tsv"
singularity exec --cleanenv "$WORK/containers/lenstools.sif" python3 "$WORK/finalize_candidates.py" --root "$ROOT"
# Preserve the native, stricter source tier separately: it requires DNA linkage
# for splice events, which this RNA-only run cannot establish.
singularity exec --cleanenv "$WORK/containers/lenstools.sif" python3 /opt/lenstools/lenstools.py \
 filter-high-confidence -i "$ROOT/final/pass-candidates.tsv" -o "$ROOT/final/native-high-confidence.tsv"
