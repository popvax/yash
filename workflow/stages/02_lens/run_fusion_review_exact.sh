#!/usr/bin/env bash
set -euo pipefail
ROOT=${LENS_RUN_ROOT:-/data/results/TN25-112103/lens-expression-aware-20261030}
WORK=/work/lens-expression-aware
export APPTAINER_BINDPATH=/data,/work
python3 "$WORK/correct_fusion_review.py" "$ROOT"
while read -r cid; do
 singularity exec --cleanenv "$WORK/containers/lenstools.sif" python3 /opt/lenstools/lenstools.py \
 extract-review-evidence --candidates "$ROOT/fusion/per-breakpoint-review/$cid.candidate.tsv" \
 --tumor-rna-bam /data/inputs/TN25-112103/RNA_TN25-112103.bam \
 --reference-fasta /data/raft/references/homo_sapiens/fasta/Homo_sapiens.assembly38.fa \
 --starfusion-predictions "$ROOT/fusion/per-breakpoint-review/$cid.fusioninspector.tsv" \
 --reference-build hg38 --config "$ROOT/review_scoring_defaults.expression-aware.yaml" \
 --output-dir "$ROOT/fusion/event-review-evidence"
done < "$ROOT/fusion/per-breakpoint-review/candidate-ids.txt"
