"""Include weaker event-overlapping peptides without admitting weak-only windows."""
import bisect
import collections
import json
from pathlib import Path
import sys

R=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
original=R/'seed-sequence-contexts.jsonl'
if not original.exists():original.write_bytes((R/'sequence-contexts.jsonl').read_bytes())
contexts=[json.loads(l) for l in original.read_text().splitlines()]
queries={t['peptide'] for c in contexts for t in c['targets']};counts=collections.Counter()
for c in contexts:
    seq=c['sequence'];c['seed_targets']=c['targets']
    starts=sorted({s for t in c['seed_targets'] for s in range(max(0,t['end']-25),min(t['start'],len(seq)-25)+1)})
    c['window_starts']=starts
    targets=[];seeds={(t['peptide'],t['start']) for t in c['seed_targets']}
    for k in (8,9,10,11):
        for i in range(len(seq)-k+1):
            if not (i<c['target_end'] and i+k>c['target_start']):continue
            j=bisect.bisect_left(starts,max(0,i+k-25))
            if j==len(starts) or starts[j]>min(i,len(seq)-25):continue
            pep=seq[i:i+k]
            if not set(pep)<=set('ACDEFGHIKLMNPQRSTVWY'):continue
            targets.append(dict(peptide=pep,start=i,end=i+k,individually_admitted_seed=(pep,i) in seeds))
            queries.add(pep)
    c['targets']=targets;counts[c['source']]+=len(targets)
with (R/'sequence-contexts.jsonl').open('w') as f:
    for c in contexts:f.write(json.dumps(c)+'\n')
(R/'expanded-target-peptides.txt').write_text('\n'.join(sorted(queries))+'\n')
summary=dict(contexts=len(contexts),scoring_occurrences_by_source=dict(counts),unique_scoring_peptides=len(queries),
             seed_gate='Every candidate window contains a previously admitted LENS target; additional weaker peptides contribute only within these windows',
             RNA_scope='Additional peptides inherit source/frame context support; individual LENS peptide review is not asserted')
(R/'expanded-scoring-summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
