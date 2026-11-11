#!/usr/bin/env bash
# Recompute every primary splice peptide with the verified MSNAF contract.
set -euo pipefail
ROOT=${LENS_RUN_ROOT:-/data/results/TN25-112103/lens-expression-aware-20261030}
WORK=/work/lens-expression-aware
export APPTAINER_BINDPATH=/data,/work
python3 - "$ROOT" <<'PY'
import csv,json,sys,hashlib
from pathlib import Path
r=Path(sys.argv[1]);p=r/'peptide-review';bak=p/'evidence.automatic-coordinate-withdrawn';bak.mkdir(exist_ok=True)
with (p/'candidates.tsv').open() as f:
 reader=csv.DictReader(f,delimiter='\t');columns=reader.fieldnames;rows=[x for x in reader if x['source_type']=='SPLICE']
shards=[[] for _ in range(8)]
for x in rows:
 source=p/'evidence'/(x['candidate_id']+'.json');target=bak/source.name
 if source.exists() and not target.exists():source.rename(target)
 i=int(hashlib.sha256(x['splice_description'].encode()).hexdigest(),16)%8;shards[i].append(x)
for i,xs in enumerate(shards):
 xs.sort(key=lambda x:(x['splice_description'],x['identity']))
 with (p/f'splice-contract.{i}.tsv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=columns,delimiter='\t');w.writeheader();w.writerows(xs)
(p/'coordinate-rerun-summary.json').write_text(json.dumps({'splice_peptide_reviews':len(rows),'rows_per_shard':[len(x) for x in shards],'partition':'entire junction UID assigned to one worker; native bounded cache reuses its read metrics','coordinate_contract':'MSNAF exon-boundary coordinates -> intron start+1,end-1; exact CIGAR intron match'},indent=2)+'\n');print('Contract-correct splice reviews:',len(rows),flush=True)
PY
pids=()
for i in 0 1 2 3 4 5 6 7; do
 singularity exec --cleanenv --env LENS_SPLICE_EVIDENCE_CACHE_SIZE=16 --env LENS_MSNAF_EXON_BOUNDS=1 \
 --bind "$WORK/metadata.splice-review.py:/opt/lenstools/src/metadata.py:ro" \
 "$WORK/containers/lenstools.sif" python3 /opt/lenstools/lenstools.py \
 extract-review-evidence --candidates "$ROOT/peptide-review/splice-contract.$i.tsv" \
 --tumor-rna-bam /data/inputs/TN25-112103/RNA_TN25-112103.bam \
 --reference-fasta /data/raft/references/homo_sapiens/fasta/Homo_sapiens.assembly38.fa \
 --gtf "$ROOT/splice/review.exons.complete.gtf" --reference-build hg38 \
 --config "$ROOT/review_scoring_defaults.expression-aware.yaml" \
 --output-dir "$ROOT/peptide-review/evidence" > "$ROOT/logs/peptide-review-contract-$i.log" 2>&1 &
 pids+=("$!")
done
failed=0;for p in "${pids[@]}"; do wait "$p" || failed=1; done
test "$failed" = 0
lens() { singularity exec --cleanenv "$WORK/containers/lenstools.sif" python3 /opt/lenstools/lenstools.py "$@"; }
lens derive-review-flags --candidates "$ROOT/peptide-review/candidates.tsv" \
 --evidence-dir "$ROOT/peptide-review/evidence" --config "$ROOT/review_scoring_defaults.expression-aware.yaml" \
 --output "$ROOT/peptide-review/flags.tsv"
lens score-review-candidates --flags-tsv "$ROOT/peptide-review/flags.tsv" \
 --evidence-dir "$ROOT/peptide-review/evidence" --config "$ROOT/review_scoring_defaults.expression-aware.yaml" \
 --output "$ROOT/peptide-review/scores.tsv"
