#!/usr/bin/env bash
# Equivalent to native lenstools_make_fusion_peptides_context (reference CDS path).
set -euo pipefail
ROOT=${LENS_RUN_ROOT:-/data/results/TN25-112103/lens-expression-aware-20261030}
WORK=/work/lens-expression-aware
REF=/data/raft/references/homo_sapiens
OLD=/data/results/TN25-112103/original-fastqs/outputs/samples/TN25_research/TN25-112103/ar-TN25-112103
FUSIONS="$OLD/starfusion/TN25_research-TN25-112103-ar-TN25-112103.star-fusion.fusion_predictions.abridged.coding_effect.tsv"
export APPTAINER_BINDPATH=/data,/work
mkdir -p "$ROOT/fusion"
# Select by header, avoiding dependence on the physical TSV column order.
python3 - "$FUSIONS" "$REF/annot/gencode.v37.annotation.with.hervs.gtf" "$ROOT/fusion" <<'PY'
import csv,re,sys
from pathlib import Path
calls,gtf,out=sys.argv[1:];out=Path(out)
with open(calls) as f:
 rows=list(csv.DictReader(f,delimiter='\t'))
txs={r[k] for r in rows for k in ('CDS_LEFT_ID','CDS_RIGHT_ID') if r[k].startswith('ENS')}
(out/'fusion-transcripts.txt').write_text('\n'.join(sorted(txs))+'\n')
pattern=re.compile(r'transcript_id "([^"]+)"')
with open(gtf) as fi,(out/'fusion.selected.gtf').open('w') as fo:
 for line in fi:
  m=pattern.search(line)
  if m and m[1] in txs:fo.write(line)
PY
singularity exec --cleanenv "$WORK/containers/lenstools.sif" \
 python3 /opt/lenstools/lenstools.py make-fusion-peptides-context \
 -f "$FUSIONS" -g "$ROOT/fusion/fusion.selected.gtf" \
 -t "$ROOT/fusion/fusion-transcripts.txt" \
 -o "$ROOT/fusion/fusion.pep.fa" -n "$ROOT/fusion/fusion.nt.fa"
