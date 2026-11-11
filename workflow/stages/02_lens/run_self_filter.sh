#!/usr/bin/env bash
# Native LENS self-peptidome filter, parallelized by peptide to retain cache reuse.
set -euo pipefail
ROOT=${LENS_RUN_ROOT:-/data/results/TN25-112103/lens-expression-aware-20261030}
WORK=/work/lens-expression-aware
export APPTAINER_BINDPATH=/data,/work
mkdir -p "$ROOT/prediction/self-filter-shards"
python3 - "$ROOT" <<'PY'
import csv,hashlib,sys
from pathlib import Path
r=Path(sys.argv[1]);d=r/'prediction/self-filter-shards'
handles=[(d/f'input.{i}.tsv').open('w') for i in range(4)];writers=[]
with (r/'prediction/lens-annotated.tsv').open() as f:
 reader=csv.DictReader(f,delimiter='\t')
 for h in handles:
  w=csv.DictWriter(h,fieldnames=reader.fieldnames,delimiter='\t');w.writeheader();writers.append(w)
 for row in reader:
  i=int(hashlib.md5(row['peptide'].encode()).hexdigest(),16)%4;writers[i].writerow(row)
for h in handles:h.close()
PY
pids=()
for i in 0 1 2 3; do
 singularity exec --cleanenv "$WORK/containers/lenstools.sif" \
  python3 /opt/lenstools/lenstools.py filter-peptides \
  -i "$ROOT/prediction/self-filter-shards/input.$i.tsv" \
  -p /data/raft/references/homo_sapiens/protein/gencode.v37.pc_translations.fa \
  -c /data/raft/references/homo_sapiens/cta_self/cta.pirlygenes.06APR2026.gene_list \
  -o "$ROOT/prediction/self-filter-shards/output.$i.tsv" > "$ROOT/logs/self-filter-$i.log" 2>&1 &
 pids+=("$!")
done
failed=0;for p in "${pids[@]}"; do wait "$p" || failed=1; done
test "$failed" = 0
python3 - "$ROOT" <<'PY'
import sys,csv,json,collections
from pathlib import Path
r=Path(sys.argv[1]);counts=collections.Counter()
with (r/'prediction/lens-self-filtered.tsv').open('w') as out:
 writer=None
 for i in range(4):
  with (r/f'prediction/self-filter-shards/output.{i}.tsv').open() as f:
   reader=csv.DictReader(f,delimiter='\t')
   if writer is None:writer=csv.DictWriter(out,fieldnames=reader.fieldnames,delimiter='\t');writer.writeheader()
   assert reader.fieldnames==writer.fieldnames
   for row in reader:writer.writerow(row);counts[row['antigen_source']]+=1
(r/'prediction/self-filter-summary.json').write_text(json.dumps(dict(counts),indent=2)+'\n');print(json.dumps(dict(counts)))
PY
