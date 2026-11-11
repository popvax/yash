"""Create a common event/context/target registry from frozen source results."""
import collections
import csv
import hashlib
import json
from pathlib import Path
from Bio.Seq import Seq

R = Path(__file__).resolve().parent
def load_tsv(p):
    with p.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))
def offsets(sequence, peptide):
    start = sequence.find(peptide)
    while start >= 0:
        yield start
        start = sequence.find(peptide, start+1)

def main():
    rows = [r for r in load_tsv(R/'inputs/lens-all-reviewed.tsv') if r['primary_eligible']=='True' and r['antigen_source'] in ['CTA/SELF','FUSION']]
    byid = collections.defaultdict(list)
    for r in rows:
        byid[r['identity']].append(r)
    metadata = {}
    with (R/'inputs/prediction/source-metadata.jsonl').open() as f:
        for line in f:
            m=json.loads(line)
            if m['identity'] in byid and m['source']!='SPLICE':
                metadata[m['identity']]=m
    fusion_calls = load_tsv(R/'inputs/inputs-checkpoints/original-starfusion-coding-effect.tsv')
    contexts=[]; mappings=[]; failures=[]
    for oid,m in metadata.items():
        source=m['source'];rs=byid[oid];first=rs[0]
        if source=='CTA/SELF':
            seq=m['sequence'];lo=0;hi=len(seq);mode='intentional_reference_CTA_region'
        else:
            calls=[f for f in fusion_calls if f['#FusionName']==m['fusion_id'] and
                   f['LeftBreakpoint']==m['fusion_left_breakpoint'] and f['RightBreakpoint']==m['fusion_right_breakpoint'] and
                   f['CDS_LEFT_ID']==m['fusion_left_transcript'] and f['CDS_RIGHT_ID']==m['fusion_right_transcript']]
            assert len(calls)==1
            call=calls[0];seq=str(Seq(call['FUSION_CDS']).translate(to_stop=True))
            assert seq==call['FUSION_TRANSL'].split('*')[0]
            boundary_nt=int(call['CDS_LEFT_RANGE'].split('-')[1]);lo=boundary_nt//3
            hi=len(seq) if call['PROT_FUSION_TYPE']=='FRAMESHIFT' else lo+(boundary_nt%3 != 0)
            mode='frameshift_tail' if call['PROT_FUSION_TYPE']=='FRAMESHIFT' else 'junction_codon' if hi>lo else 'junction_boundary'
            expected_context_start=max(0,lo-16)
            assert seq[expected_context_start:expected_context_start+len(m['pep_context'])]==m['pep_context']
            assert m['sequence'] in seq
        targets=[]
        for pep in sorted({r['peptide'] for r in rs}):
            found=False
            for start in offsets(seq,pep):
                end=start+len(pep)
                if not (start<hi and end>lo):
                    continue
                targets.append(dict(identity=oid,peptide=pep,start=start,end=end));found=True
            if not found:failures.append(dict(identity=oid,peptide=pep,reason='NO_EVENT_OVERLAPPING_OCCURRENCE',source=source))
        contexts.append(dict(context_id=source.replace('/','_')+'_'+oid,source=source,event_id=first['event_id'],gene=first['gene'],
                             sequence=seq,target_start=lo,target_end=hi,target_mode=mode,phase=None,
                             transcript_id=m.get('transcript_id'),targets=targets))
    admitted={(r['identity'],r['peptide']) for r in rows}
    covered=set()
    for c in contexts:
        for t in c['targets']:
            assert c['sequence'][t['start']:t['end']]==t['peptide']
            assert t['start']<c['target_end'] and t['end']>c['target_start'],(c,t)
            assert (t['identity'],t['peptide']) in admitted
            covered.add((t['identity'],t['peptide']))
        by_origin=collections.defaultdict(list)
        for t in c['targets']:by_origin[t['identity']].append(t)
        # Canonical expression metadata is event-level; all origins must agree.
        evidence_rows=[byid[oid][0] for oid in by_origin]
        for field in ['expression_proxy_left','expression_proxy_right','left_or_single_event_fraction','right_event_fraction','expression_proxy_notes']:
            vals={r[field] for r in evidence_rows};assert len(vals)<=1,(c['context_id'],field,vals)
            c[field]=next(iter(vals),'')
        c['source_identity_count']=len(by_origin)
        if evidence_rows:
            c['normal_reference_annotations']={k:v for k,v in evidence_rows[0].items()
                if k.startswith('cta_') or k in ['fraction_baseline_quality','native_splice_psi_like','rna_support_count','rna_baseline_count']}
    missing=admitted-covered
    assert missing=={(x['identity'],x['peptide']) for x in failures},(len(missing),len(failures))
    with (R/'sequence-contexts.jsonl').open('w') as f:
        for c in contexts:f.write(json.dumps(c)+'\n')
    peptides=sorted({t['peptide'] for c in contexts for t in c['targets']})
    (R/'target-peptides.txt').write_text('\n'.join(peptides)+'\n')
    (R/'unmapped-admitted-targets.json').write_text(json.dumps(failures,indent=2)+'\n')
    summary=dict(contexts_by_source=dict(collections.Counter(c['source'] for c in contexts)),
                 input_primary_rows=len(rows),input_origin_peptides=len(admitted),mapped_origin_peptides=len(covered),
                 unmapped_origin_peptides=len(missing),unique_target_peptides=len(peptides),
                 short_contexts=sum(len(c['sequence'])<25 for c in contexts),
                 input_hashes={str(p.relative_to(R)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                               [R/'inputs/lens-all-reviewed.tsv',R/'inputs/mutation-regions.tsv']})
    (R/'context-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
