#!/usr/bin/env bash
set -euo pipefail
ROOT=${LENS_RUN_ROOT:-/data/results/TN25-112103/lens-expression-aware-20261030}
WORK=/work/lens-expression-aware
export APPTAINER_BINDPATH=/data,/work
python3 "$WORK/prepare_fusion_review.py" --root "$ROOT"
bash "$WORK/run_fusion_review_exact.sh"
