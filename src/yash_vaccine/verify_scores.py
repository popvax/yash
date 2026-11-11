"""Independently recount every saved LENS window from native peptide scores."""
import collections,math
from .io import rows,write_json

def verify(checkpoints,out):
    per_context=collections.defaultdict(list)
    occurrences=0
    for row in rows(checkpoints/'target-occurrences.tsv.gz'):
        occurrences+=1
        if row['excluded_exact_self']=='True':continue
        per_context[row['context_id']].append((int(row['start']),int(row['end']),float(row['affinity_contribution'])))
    score_at=collections.defaultdict(lambda:collections.defaultdict(float))
    # Window length is fixed: a target at [a,b) fits starts [b-25,a].
    for context,targets in per_context.items():
        for begin,end,score in targets:
            for start in range(max(0,end-25),begin+1):score_at[context][start]+=score
    checked=0;maxdiff=0.
    for row in rows(checkpoints/'lens-scored-windows.tsv.gz'):
        value=score_at[row['context_id']][int(row['window_start'])]
        original=float(row['native_affinity_sum']);maxdiff=max(maxdiff,abs(value-original))
        assert math.isclose(value,original,rel_tol=1e-9,abs_tol=1e-9),row['window_id']
        checked+=1
    assert checked==139677
    write_json(out/'window-score-verification.json',dict(status='VERIFIED',windows=checked,target_occurrences=occurrences,
                max_absolute_score_difference=maxdiff,scope='Independent window sum from frozen native per-peptide contributions; model inference not rerun'))
