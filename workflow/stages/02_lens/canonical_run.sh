#!/usr/bin/env bash
# Clean rerun phases on the saved references/inputs. Use a NEW output directory.
set -euo pipefail
: "${LENS_RUN_ROOT:?Set a new output directory on /data}"
ROOT="$LENS_RUN_ROOT"
WORK=/work/lens-expression-aware
SOURCE=/data/results/TN25-112103/lens-expression-aware-20261030
export APPTAINER_BINDPATH=/data,/work
mkdir -p "$WORK/containers" "$ROOT/logs" "$ROOT/source" "$ROOT/annotations"
for image in lenstools starfusion msnaf1.2.1; do
 cp "$SOURCE/containers/$image.sif" "$WORK/containers/$image.sif"
 expected=$(cut -d ' ' -f1 "$SOURCE/containers/$image.sha256")
 actual=$(sha256sum "$WORK/containers/$image.sif" | cut -d ' ' -f1)
 test "$actual" = "$expected"
done
case "${1:?pre or post}" in
pre)
 test "$ROOT" != "$SOURCE"
 test ! -e "$ROOT/cta/cta.pep.fa"
 mkdir -p "$ROOT/upstream-originals"
 singularity exec --cleanenv "$WORK/containers/lenstools.sif" cat /opt/lenstools/src/metadata.py > "$ROOT/upstream-originals/lenstools-metadata.py"
 singularity exec --cleanenv "$WORK/containers/msnaf1.2.1.sif" cat /usr/local/lib/python3.10/site-packages/msnaf/core.py > "$ROOT/upstream-originals/msnaf-core.py"
 python3 "$WORK/prepare_upstream.py" --sources-dir "$ROOT/upstream-originals" --output "$ROOT/upstream-patched" --execute
 for name in metadata.annotation-header-fix.py metadata.splice-review.py msnaf-core.8-11.py; do
  ln -sf "$ROOT/upstream-patched/$name" "$WORK/$name"
 done
 cp "$WORK"/*.py "$WORK"/*.sh "$ROOT/source/"
 cp "$WORK/msnaf-core.8-11.py" "$ROOT/source/"
 (cd "$ROOT/annotations"; python3 "$WORK/prepare_cta_annotation_native.py" /data/raft/references/homo_sapiens/cta_self/canonical.annotated.mtec.norm_tiss.subcell.21APR2026.tsv)
 bash "$WORK/run_cta.sh" > "$ROOT/logs/discovery-cta.log" 2>&1 & p1=$!
 bash "$WORK/run_splice_sharded.sh" > "$ROOT/logs/discovery-splice.log" 2>&1 & p2=$!
 bash "$WORK/run_discovery.sh" fusion > "$ROOT/logs/discovery-fusion.log" 2>&1 & p3=$!
 wait "$p1";echo 0 > "$ROOT/logs/cta.exit"
 wait "$p2";echo 0 > "$ROOT/logs/splice.exit"
 wait "$p3";echo 0 > "$ROOT/logs/fusion.exit"
 bash "$WORK/run_fusion_sequences.sh" > "$ROOT/logs/fusion-sequences.log" 2>&1
 bash "$WORK/run_cta_review.sh" > "$ROOT/logs/cta-review.log" 2>&1
 bash "$WORK/run_fusion_review.sh" > "$ROOT/logs/fusion-review.log" 2>&1
 python3 "$WORK/prepare_splice_review.py" --root "$ROOT"
 python3 "$WORK/prepare_pmhc_inputs.py" --root "$ROOT"
 ;;
post)
 # GPU output has been explicitly materialized in WORK by the controller.
 bash "$WORK/run_annotation.sh" > "$ROOT/logs/annotation.log" 2>&1
 bash "$WORK/run_self_filter.sh" > "$ROOT/logs/self-filter.log" 2>&1
 python3 "$WORK/prepare_peptide_review.py" --root "$ROOT"
 bash "$WORK/rerun_splice_contract.sh" > "$ROOT/logs/peptide-review.log" 2>&1
 bash "$WORK/run_finalize.sh" > "$ROOT/logs/finalize.log" 2>&1
 python3 -m pip install --target "$WORK/verification-deps" pyahocorasick==2.2.0 > "$ROOT/logs/verification-dependencies.log" 2>&1
 python3 "$WORK/audit_reference_filter.py" "$ROOT" > "$ROOT/logs/reference-filter-audit.log" 2>&1
 python3 "$WORK/hash_inputs.py" --root "$ROOT" > "$ROOT/logs/input-hashes.log" 2>&1
 singularity exec --cleanenv "$WORK/containers/lenstools.sif" python3 "$WORK/audit_all_splice_motifs.py" "$ROOT" > "$ROOT/logs/all-splice-motifs-audit.log" 2>&1
 python3 "$WORK/verify_final.py" --root "$ROOT" > "$ROOT/logs/completion-audit.log" 2>&1
 ;;
*) exit 2;;
esac
