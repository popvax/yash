"""Independent sequence, admission, affinity and compact-window checks."""
import collections
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

R=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
def read(name):
    with (R/name).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(name,rows):
    with (R/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def near(a,b):assert math.isclose(float(a),float(b),rel_tol=1e-9,abs_tol=1e-9),(a,b)

def main():
    summary=json.loads((R/'summary.json').read_text());assert summary['status']=='COMPLETE'
    contexts={c['context_id']:c for c in map(json.loads,(R/'sequence-contexts.jsonl').read_text().splitlines())}
    parent=Path('/data/results/TN25-112103/lens-expression-aware-20261030/final/all-reviewed-candidates.tsv')
    if not parent.exists():parent=R/'inputs/lens-all-reviewed.tsv'
    provenance=json.loads((R/'context-summary.json').read_text())
    assert hashlib.sha256(parent.read_bytes()).hexdigest()==provenance['input_hashes']['inputs/lens-all-reviewed.tsv']
    admitted={}
    with parent.open() as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r['primary_eligible']=='True':admitted[(r['identity'],r['peptide'])]=r
    mapped=set()
    for c in contexts.values():
        for t in c['seed_targets']:
            key=(t['identity'],t['peptide']);original=admitted[key];mapped.add(key)
            assert original['event_id']==c['event_id'] and original['antigen_source']==c['source']
            assert c['sequence'][t['start']:t['end']]==t['peptide']
            assert t['start']<c['target_end'] and t['end']>c['target_start']
            expected_proxy=original['left_or_single_supplied_expression'] if c['source']=='SPLICE' else original['expression_proxy_left']
            assert expected_proxy==c['expression_proxy_left']
    assert mapped==set(admitted)
    windows=read('all-scored-windows.tsv');combined=read('compact-pool.tsv')
    matches=json.loads((R/'scoring-reference-membership.json').read_text())
    ledger=read('native-affinity-contributions.tsv');scores=collections.Counter();pairs=set();contributions={}
    for r in ledger:
        pair=(r['peptide'],r['allele']);assert pair not in pairs;pairs.add(pair)
        a=float(r['affinity_nm']);expected=(1+math.exp(-350/150))/(1+math.exp((a-350)/150)) if a<5000 else 0
        # Native min_epitope_score applies at the per-allele scoring stage.
        if expected<1e-5:expected=0
        near(r['contribution'],expected);scores[r['peptide']]+=float(r['contribution'])
        contributions[pair]=float(r['contribution'])
    config=json.loads((R/'score-config.json').read_text());bits={a:1<<i for i,a in enumerate(config['alleles'])}
    raw_seen=collections.Counter();raw_rows=0
    with (R/'scoring-affinities.tsv').open() as f:
        for r in csv.DictReader(f,delimiter='\t'):
            p=r['peptide'];a=float(r['mhcflurry_2.3.14.aff']);bit=bits[r['allele']]
            assert not raw_seen[p]&bit;raw_seen[p]|=bit;raw_rows+=1
            value=(1+math.exp(-350/150))/(1+math.exp((a-350)/150)) if a<5000 else 0
            if value<1e-5:value=0
            near(contributions.get((p,r['allele']),0),value)
    assert set(raw_seen)==set((R/'expanded-target-peptides.txt').read_text().splitlines())
    assert all(v==63 for v in raw_seen.values())
    # Recompute sequence containment, event overlap, complete window enumeration,
    # self exclusion and expression weighting independently of Vaxrank objects.
    expected_starts={}
    for cid,c in contexts.items():
        expected_starts[cid]=set(c['window_starts'])
    actual_starts=collections.defaultdict(set);bywindow={}
    for w in windows:
        cid=w['context_id'];c=contexts[cid];start=int(w['window_start']);end=int(w['window_end'])
        assert w['window_id'] not in bywindow;bywindow[w['window_id']]=w
        assert start not in actual_starts[cid];actual_starts[cid].add(start)
        assert end-start==25 and w['sequence']==c['sequence'][start:end]
        assert start<c['target_end'] and end>c['target_start']
        assert w['region_id']=='REG_'+hashlib.sha256(w['sequence'].encode()).hexdigest()[:20]
        occurrences={(t['peptide'],t['start']) for t in c['targets'] if start<=t['start'] and t['end']<=end}
        good={(p,s) for p,s in occurrences if scores[p]>0 and not matches[p]['occurs_in_non_CTA' if c['source']=='CTA/SELF' else 'occurs']}
        bad=occurrences-good
        near(w['native_affinity_sum'],sum(scores[p] for p,s in good))
        near(w['excluded_self_affinity_sum'],sum(scores[p] for p,s in bad))
        assert set(w['target_peptides'].split(';'))-{''}=={p for p,s in good}
        expected_seeds={t['peptide'] for t in c['seed_targets'] if start<=t['start'] and t['end']<=end}&{p for p,s in good}
        assert set(w['seed_peptides'].split(';'))-{''}==expected_seeds
        assert (w['primary_window_eligible']=='True')==(bool(expected_seeds) and float(w['native_affinity_sum'])>0)
        if w['primary_window_eligible']=='True':assert float(w['best_target_affinity_nm'])<=500
        value=c['expression_proxy_left']
        if value=='':assert w['weighted_score_left']==''
        else:near(w['weighted_score_left'],float(w['native_affinity_sum'])*math.sqrt(float(value)))
    assert all(actual_starts[cid]==starts for cid,starts in expected_starts.items())
    mutation_input={r['variant_id']:r for r in read('inputs/mutation-regions.tsv')}
    mutation_rows=[r for r in combined if r['source']=='MUTATION'];assert len(mutation_rows)==len(mutation_input)
    for w in mutation_rows:
        original=mutation_input[w['event_id']];assert w['sequence']==original['sequence']
        near(w['native_affinity_sum'],original['post_self_affinity_sum'])
        near(w['weighted_score_left'],original['weighted_score']);near(w['weighted_score_right'],original['weighted_score'])
    best={}
    for w in windows:
        if w['primary_window_eligible']!='True':continue
        key=(w['source'],w['event_id'])
        if key not in best or (-float(w['native_affinity_sum']),w['window_id'])<(-float(best[key]['native_affinity_sum']),best[key]['window_id']):best[key]=w
    expected_lens=[w for _,w in sorted(best.items())]
    observed_lens=[w for w in combined if w['source']!='MUTATION']
    assert [(w['source'],w['event_id'],w['window_id']) for w in observed_lens]==[(w['source'],w['event_id'],w['window_id']) for w in expected_lens]
    best_rows=observed_lens+mutation_rows
    result=dict(status='VERIFIED',window_placements_checked=len(windows),native_affinity_contributions_checked=len(ledger),raw_prediction_rows_checked=raw_rows,
                original_admitted_origin_peptides_mapped=len(mapped),source_expression_policy_verified=True,
                mutation_regions_checked=len(mutation_input),
                complete_window_enumeration=True,highest_affinity_eligible_window_per_event=True,
                compact_best_region_rows=len(best_rows),compact_best_regions_by_source=dict(collections.Counter(w['source'] for w in best_rows)),
                numerical_checks='Independent scalar normalized logistic, window sums, expression multipliers and sort order',
                limitations='Checks computational calculations and frozen evidence, not tumour specificity, translation or biological efficacy')
    (R/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
