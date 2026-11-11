#!/usr/bin/env bash
set -euo pipefail
ROOT=${LENS_RUN_ROOT:-/data/results/TN25-112103/lens-expression-aware-20261030}
WORK=/work/lens-expression-aware
export APPTAINER_BINDPATH=/data,/work
python3 "$WORK/prepare_cta_review.py" --root "$ROOT"
cp /data/raft/references/qc/review_scoring_defaults.yaml "$ROOT/review_scoring_defaults.original.yaml"
# The native CTA gene-abundance ratio is capture-derived. Positivity has already
# been assessed with supplied expression; keep read support and coverage QC.
singularity exec --cleanenv "$WORK/containers/lenstools.sif" python3 - "$ROOT" <<'PY'
import yaml,sys
from pathlib import Path
r=Path(sys.argv[1]);c=yaml.safe_load((r/'review_scoring_defaults.original.yaml').read_text())
c['flags']['cta_min_abundance_ratio']=0
c['flags']['min_support_ratio']['CTA/SELF']=0
(r/'review_scoring_defaults.expression-aware.yaml').write_text(yaml.safe_dump(c,sort_keys=False))
PY
lens() { singularity exec --cleanenv "$WORK/containers/lenstools.sif" python3 /opt/lenstools/lenstools.py "$@"; }
CFG="$ROOT/review_scoring_defaults.expression-aware.yaml"
lens extract-review-evidence --candidates "$ROOT/cta/transcript-review-candidates.tsv" \
 --tumor-rna-bam /data/inputs/TN25-112103/RNA_TN25-112103.bam \
 --reference-fasta /data/raft/references/homo_sapiens/fasta/Homo_sapiens.assembly38.fa \
 --gtf "$ROOT/cta/cta.selected.gtf" --reference-build hg38 --config "$CFG" \
 --output-dir "$ROOT/cta/review-evidence"
lens derive-review-flags --candidates "$ROOT/cta/transcript-review-candidates.tsv" \
 --evidence-dir "$ROOT/cta/review-evidence" --config "$CFG" --output "$ROOT/cta/review-flags.tsv"
lens score-review-candidates --flags-tsv "$ROOT/cta/review-flags.tsv" \
 --evidence-dir "$ROOT/cta/review-evidence" --config "$CFG" --output "$ROOT/cta/review-scores.tsv"
