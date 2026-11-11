"""Attach exact-breakpoint FusionInspector counts and two expression proxies."""
import csv,json,argparse
from pathlib import Path
from prepare_pmhc_inputs import fasta
from expression_weights import fusion_proxies
p=argparse.ArgumentParser();p.add_argument('--root',required=True);a=p.parse_args();r=Path(a.root)
expr={}
for line in Path('/data/inputs/TN25-112103/vendor-transformed-expression.tsv').read_text().splitlines():
 if line and not line.startswith('#'):
  gene,v=line.split('\t')[:2];expr[gene]=float(v)
with (r/'fusion/fusioninspector/finspector.FusionInspector.fusions.abridged.tsv').open() as f:fi=list(csv.DictReader(f,delimiter='\t'))
lookup={}
for x in fi:
 key=x['#FusionName'],x['LeftBreakpoint'],x['RightBreakpoint'];assert key not in lookup;lookup[key]=x
rows=[]
for identity,seq,m in fasta(r/'fusion/fusion.pep.fa'):
 left=m['fusion_left_gene'].split('^')[0];right=m['fusion_right_gene'].split('^')[0]
 x=lookup.get((m['fusion_id'],m['fusion_left_breakpoint'],m['fusion_right_breakpoint']))
 row=dict(candidate_id=identity,source_type='FUSION',peptide='',hla_allele='',peptide_length=0,gene=m['fusion_id'],transcript_id='',lens_row_id=identity,**m)
 row.update(left_supplied_expression=expr.get(left,''),right_supplied_expression=expr.get(right,''),fusion_inspector_match='exact_breakpoints' if x else 'not_exactly_reconfirmed')
 if x:
  F=float(x['JunctionReadCount'])+float(x['SpanningFragCount']);L=float(x['NumCounterFusionLeft']);R=float(x['NumCounterFusionRight'])
  row.update(fusion_fragments=F,left_unfused_fragments=L,right_unfused_fragments=R,FAR_left_reported=x['FAR_left'],FAR_right_reported=x['FAR_right'],**fusion_proxies(expr.get(left),expr.get(right),F,L,R))
 rows.append(row)
cols=list(dict.fromkeys(k for row in rows for k in row))
with (r/'fusion/event-review-candidates.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows(rows)
summary={'native_sequence_contexts':len(rows),'fusion_inspector_calls':len(fi),'exact_reconfirmed_sequence_contexts':sum(x['fusion_inspector_match']=='exact_breakpoints' for x in rows),'count_definition':'JunctionReadCount + SpanningFragCount; bounded fractions use raw NumCounterFusionLeft/Right without FAR pseudocounts','reference_discovery_pool':'original full STAR-Fusion coding-effect calls, before HLA prediction; FusionInspector provides confirmation and counter-fusion evidence'}
(r/'fusion/expression-proxy-summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))
