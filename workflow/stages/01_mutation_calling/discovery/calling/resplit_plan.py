"""Subdivide a slow interval without changing scope, nomination or caller settings."""
from pathlib import Path
import argparse,hashlib,json

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--shard',required=True);p.add_argument('--parts',type=int,required=True);a=p.parse_args()
 assert a.input.resolve()!=a.output.resolve(),'Preserve the original plan separately'
 assert 2<=a.parts<=32
 data=a.input.read_bytes();plan=json.loads(data);parent=next(s for s in plan['shards'] if s['id']==a.shard)
 length=parent['end']-parent['start']+1;assert length>=a.parts
 children=[]
 for i in range(a.parts):
  lo=parent['start']+length*i//a.parts;hi=parent['start']+length*(i+1)//a.parts-1
  children.append(dict(id=parent['id']+f'_part{i:02d}',chrom=parent['chrom'],start=lo,end=hi,interval=f"{parent['chrom']}:{lo}-{hi}",estimated_read_weight=parent['estimated_read_weight']/a.parts))
 assert children[0]['start']==parent['start'] and children[-1]['end']==parent['end']
 assert all(a['end']+1==b['start'] for a,b in zip(children,children[1:]))
 before=list(plan['shards']);plan['shards']=[child for s in before for child in (children if s['id']==parent['id'] else [s])]
 plan.setdefault('scatter_revisions',[]).append(dict(parent_plan_sha256=hashlib.sha256(data).hexdigest(),superseded_shard=parent,children=[c['id'] for c in children],reason='Runtime imbalance; unchanged full genomic scope and caller/filtering settings; completed other intervals reused'))
 assert len(plan['shards'])==len(before)-1+a.parts
 a.output.write_text(json.dumps(plan,indent=2)+'\n');print(json.dumps(dict(original_intervals=len(before),new_intervals=len(plan['shards']),subdivided=parent['interval'],parts=a.parts,coverage_unchanged=True),indent=2))
if __name__=='__main__':main()
