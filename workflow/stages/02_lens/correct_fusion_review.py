"""Use one exact breakpoint-matched FI row per native fusion-review invocation.
Native LENSTools indexes its external table by FusionName, losing multiple
breakpoints for the same pair. Isolating one row avoids changing native biology.
"""
import csv,json,sys,shutil
from pathlib import Path
r=Path(sys.argv[1]);out=r/'fusion/per-breakpoint-review';out.mkdir(exist_ok=True)
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows,columns):
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=columns,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
fi=read(r/'fusion/fusioninspector/finspector.FusionInspector.fusions.abridged.tsv');lookup={(x['#FusionName'],x['LeftBreakpoint'],x['RightBreakpoint']):x for x in fi};assert len(lookup)==len(fi)
rows=read(r/'fusion/event-review-candidates.tsv')
for x in rows:
 key=x['fusion_id'],x['fusion_left_breakpoint'],x['fusion_right_breakpoint'];assert key in lookup
 cid=x['candidate_id'];write(out/(cid+'.candidate.tsv'),[x],list(rows[0]));write(out/(cid+'.fusioninspector.tsv'),[lookup[key]],list(fi[0]))
(out/'candidate-ids.txt').write_text('\n'.join(x['candidate_id'] for x in rows)+'\n')
bak=r/'fusion/event-review-evidence.name-keyed-withdrawn'
if (r/'fusion/event-review-evidence').exists() and not bak.exists():shutil.copytree(r/'fusion/event-review-evidence',bak)
print('Exact-breakpoint native review inputs:',len(rows))
