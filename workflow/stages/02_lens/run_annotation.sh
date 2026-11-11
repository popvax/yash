#!/usr/bin/env bash
set -euo pipefail
ROOT=${LENS_RUN_ROOT:-/data/results/TN25-112103/lens-expression-aware-20261030}
WORK=/work/lens-expression-aware
export APPTAINER_BINDPATH=/data,/work
python3 "$WORK/format_predictions.py" --root "$ROOT" --max-affinity 1000 --affinity-file "$WORK/affinity-all-unique.tsv"
lens() { singularity exec --cleanenv --bind "$WORK/metadata.annotation-header-fix.py:/opt/lenstools/src/metadata.py:ro" "$WORK/containers/lenstools.sif" python3 /opt/lenstools/lenstools.py "$@"; }
lens aggregate-pmhc-outputs --mhcflurry "$ROOT/prediction/mhcflurry.working.tsv" -o "$ROOT/prediction/aggregated.tsv"
lens annotate-pmhcs -s "$ROOT/prediction/aggregated.tsv" -f "$ROOT/prediction/all-sources.pep.fa" -o "$ROOT/prediction/lens-annotated.tsv"
