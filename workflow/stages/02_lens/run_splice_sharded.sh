#!/usr/bin/env bash
# Independent native maxmin junction filtering/translation on disjoint UID shards.
# No cohort percentile or cross-junction abundance filter is applied.
set -euo pipefail
ROOT=${LENS_RUN_ROOT:-/data/results/TN25-112103/lens-expression-aware-20261030}
WORK=/work/lens-expression-aware
REF=/data/raft/references/homo_sapiens
export APPTAINER_BINDPATH=/data,/work
mkdir -p "$ROOT/splice/shards"
python3 - "$ROOT" <<'PY'
from pathlib import Path
import sys,json,hashlib
r=Path(sys.argv[1]);source=Path('/data/results/TN25-112103/original-fastqs/outputs/samples/multisample/msnaf/counts.msnaf.input.txt')
handles=[(r/f'splice/shards/counts.{i}.txt').open('w') for i in range(4)]
seen=set();counts=[0]*4
with source.open() as f:
 header=next(f)
 for h in handles:h.write(header)
 for i,line in enumerate(f):
  uid=line.split('\t')[0];assert uid not in seen;seen.add(uid)
  handles[i%4].write(line);counts[i%4]+=1
for h in handles:h.close()
(r/'splice/shards/partition.json').write_text(json.dumps({'rows':counts,'input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'method':'round-robin disjoint UID rows; native per-UID normal filtering and translation'},indent=2)+'\n')
PY
pids=()
for i in 0 1 2 3; do
 mkdir -p "$ROOT/splice/shards/$i"
 singularity exec --cleanenv --pwd "$WORK" \
  --bind "$WORK/msnaf-core.8-11.py:/usr/local/lib/python3.10/site-packages/msnaf/core.py:ro" \
  "$WORK/containers/msnaf1.2.1.sif" msnaf \
  --counts "$ROOT/splice/shards/counts.$i.txt" --refs "$REF/splice_variants/snaf/snaf-data" \
  --genome-fasta "$REF/fasta/Homo_sapiens.assembly38.fa" \
  --output "$ROOT/splice/shards/$i/msnaf_peptides.csv" --species human \
  --filter-mode maxmin --t-min 20 --n-max 3 > "$ROOT/logs/splice-shard-$i.log" 2>&1 &
 pids+=("$!")
done
failed=0
for p in "${pids[@]}"; do wait "$p" || failed=1; done
test "$failed" = 0
python3 - "$ROOT" <<'PY'
from pathlib import Path
import sys,csv,json
r=Path(sys.argv[1]);seen=set();n=0
with (r/'splice/msnaf_peptides.csv').open('w') as out:
 writer=None
 for i in range(4):
  with (r/f'splice/shards/{i}/msnaf_peptides.csv').open() as f:
   reader=csv.DictReader(f)
   if writer is None:writer=csv.DictWriter(out,fieldnames=reader.fieldnames);writer.writeheader()
   assert writer.fieldnames==reader.fieldnames
   shard_uids=set()
   for row in reader:
    assert len(row['peptide']) in (8,9,10,11)
    assert row['uid'] not in seen
    writer.writerow(row);shard_uids.add(row['uid']);n+=1
   seen.update(shard_uids)
with (r/'splice/snaf_intermediates.tsv').open('w') as out:
 for i in range(4):out.write((r/f'splice/shards/{i}/snaf_intermediates.tsv').read_text())
with (r/'splice/NeoJunction_statistics_maxmin.txt').open('w') as out:
 header=None
 for i in range(4):
  with (r/f'splice/shards/{i}/NeoJunction_statistics_maxmin.txt').open() as f:
   h=next(f)
   if header is None:header=h;out.write(h)
   assert h==header
   for line in f:out.write(line)
summary={'peptide_rows':n,'translated_junctions':len(seen),'msnaf':'1.2.1','lengths':[8,9,10,11],'shards':4}
(r/'splice/translation-summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))
PY
