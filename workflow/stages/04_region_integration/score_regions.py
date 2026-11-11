"""Native Vaxrank affinity scoring of all eligible fixed-length windows."""
import collections
import csv
import dataclasses
import hashlib
import importlib.metadata
import inspect
import json
import math
from pathlib import Path
import shutil
import sys
import msgspec

from mhctools.pred import Prediction, Kind
from vaxrank.candidate_epitope import CandidateEpitope
from vaxrank.epitope_config import EpitopeConfig
from vaxrank.epitope_dsl import attach_per_allele_scores, epitopes_for_ranking, default_score_expr
from vaxrank.epitope_logic import slice_epitopes
from vaxrank.vaccine_antigen import VaccineAntigen, TargetableMask, AminoAcidInterval, TumorSpecificityAttestation, SelfReferenceMatch
from vaxrank.vaccine_peptide import VaccinePeptide

R=Path(sys.argv[1]);S=Path('/data/results/TN25-112103');L=S/'lens-expression-aware-20261030'
def save(name,data):
    (R/name).write_text(json.dumps(data,indent=2,default=str)+'\n')
def stage(name,**kw):
    save('status.json',dict(stage=name,**kw));print(name,kw,flush=True)
def read(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(name,rows):
    if not rows:return
    columns=list(dict.fromkeys(k for r in rows for k in r))
    with (R/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=columns,delimiter='\t');w.writeheader();w.writerows(rows)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def close(a,b):return math.isclose(a,b,rel_tol=1e-9,abs_tol=1e-9)

def main():
    contexts=[json.loads(l) for l in (R/'sequence-contexts.jsonl').read_text().splitlines()]
    wanted=set((R/'expanded-target-peptides.txt').read_text().splitlines())
    raw=collections.defaultdict(list)
    stage('loading_cached_predictions',target_peptides=len(wanted))
    with (R/'scoring-affinities.tsv').open() as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r['peptide'] in wanted:raw[r['peptide']].append(r)
    assert set(raw)==wanted and all(len(rs)==6 and len({x['allele'] for x in rs})==6 for rs in raw.values())
    alleles=tuple(sorted(r['allele'] for r in next(iter(raw.values()))))
    ec=EpitopeConfig();bases={};ledger=[];unscored=[]
    for pep,rs in sorted(raw.items()):
        preds=tuple(Prediction(kind=Kind.pMHC_affinity,peptide=pep,allele=r['allele'],
                               value=float(r['mhcflurry_2.3.14.aff']),score=float(r['mhcflurry_2.3.14.aff']),
                               percentile_rank=float(r['mhcflurry_2.3.14.aff_perc']),predictor_name='mhcflurry',
                               predictor_version='2.3.14+models2.3.0') for r in rs)
        e=CandidateEpitope(sequence=pep,source_sequence=pep,offset=0,predictions=preds,
                           overlaps_targetable=True,patient_alleles=alleles)
        unscored.append(e)
    scored=[]
    for i in range(0,len(unscored),5000):
        scored.extend(epitopes_for_ranking(attach_per_allele_scores(unscored[i:i+5000],ec),ec))
        print("SCORED_PEPTIDES",min(i+5000,len(unscored)),"/",len(unscored),flush=True)
    for e in scored:
        pep=e.sequence;bases[pep]=e
        for a,v in e.per_allele_scores.items():
            ledger.append(dict(peptide=pep,allele=a,affinity_nm=next(float(r['mhcflurry_2.3.14.aff']) for r in raw[pep] if r['allele']==a),contribution=v))
    write('native-affinity-contributions.tsv',ledger)
    save('score-config.json',dict(version=importlib.metadata.version('vaxrank'),
         epitope_config=msgspec.to_builtins(ec),resolved_expression=default_score_expr(ec),
         weighted_expression='target_epitope_score * sqrt(expression_proxy)',alleles=alleles,
         admission='Every region must contain an admitted LENS RNA-PASS/source-filtered seed (<=500 nM); score all overlapping target peptides, including weaker binders, with native defaults',
         CTA_self_policy='Preserve frozen LENS curated-CTA source exception; exact non-CTA-reference exclusion already verified',
         region_length=25,window_selection='Highest native affinity sum per event/gene; window ID resolves ties'))

    membership=json.loads((R/'scoring-reference-membership.json').read_text())
    evidence_summary=json.loads((R/'scoring-evidence-summary.json').read_text())
    assert evidence_summary['native_membership_parity'] and set(membership)==wanted
    from vaxrank.vaccine_antigen import SelfReferenceSource
    att=TumorSpecificityAttestation(status='overridden',evidence_kind='RNA_only_research_region_catalogue',
        evidence_source='Frozen LENS RNA review',requires_review=True,
        override_reason='Computational region comparison only; not validated tumour specificity or construct approval')
    stage('scoring_25aa_windows',contexts=len(contexts),scoring_peptides=len(wanted),positive_scoring_peptides=len(bases))
    allwindows=[];occurrences=[];context_outcomes=[]
    for ci,c in enumerate(contexts):
        kind={'CTA/SELF':'CTA','SPLICE':'splice','FUSION':'fusion'}[c['source']]
        seq=c['sequence'];lo=c['target_start'];hi=c['target_end']
        antigen=VaccineAntigen(kind=kind,amino_acids=seq,targetable_mask=TargetableMask((AminoAcidInterval(lo,hi),)),
                               tumor_specificity=att,gene_name=c['gene'],source_identifier=c['context_id'])
        seen=set();eps=[];ep_targets={}
        for t in c['targets']:
            key=(t['peptide'],t['start'])
            if key in seen:continue
            seen.add(key)
            if t['peptide'] not in bases:continue
            member=membership[t['peptide']]
            is_cta=c['source']=='CTA/SELF'
            excluded=member['occurs_in_non_CTA'] if is_cta else member['occurs']
            witness=member['non_CTA_witness'] if is_cta else member['witness']
            match=SelfReferenceMatch(peptide=t['peptide'],occurs=excluded,antigen_kind=kind,
                sources=() if not excluded else (SelfReferenceSource(**witness),),
                excluded_gene_ids=tuple(evidence_summary['cta_excluded_gene_ids']) if is_cta else (),
                source_provenance_complete=False,genome_release='103',reference_sha256=evidence_summary['reference_sha256'])
            e=dataclasses.replace(bases[t['peptide']],source_sequence=seq,source_name=c['context_id'],
                                  offset=t['start'],self_reference_match=match,occurs_in_reference=member['occurs'],
                                  overlaps_targetable=True,overlaps_mutation=False)
            eps.append(e)
            occurrences.append(dict(context_id=c['context_id'],source=c['source'],event_id=c['event_id'],
                                    peptide=t['peptide'],start=t['start'],end=t['end'],
                                    individually_admitted_LENS_seed=t.get('individually_admitted_seed',False),
                                    affinity_contribution=e.epitope_score,excluded_exact_self=match.occurs))
        starts=c['window_starts']
        positive=0
        for start in starts:
            end=start+25;sliced=slice_epitopes(eps,start,end)
            if not sliced:continue
            vp=VaccinePeptide(antigen=antigen.sliced(start,end),epitopes=sliced,num_target_epitopes_to_keep=0,
                              combined_score_expr='target_epitope_score',ranking_rules=('target_epitope_score',))
            expected=sum(e.epitope_score for e in eps if not e.occurs_in_self_reference and e.offset>=start and e.offset+len(e.sequence)<=end)
            assert close(expected,vp.target_epitope_score)
            rid='REG_'+hashlib.sha256(vp.amino_acids.encode()).hexdigest()[:20]
            wid=c['context_id']+':'+str(start)
            target_peptides=sorted({e.sequence for e in vp.target_epitopes})
            seeds={t['peptide'] for t in c['seed_targets'] if start<=t['start'] and t['end']<=end}
            retained_seeds=sorted(seeds & set(target_peptides))
            rec=dict(window_id=wid,region_id=rid,source=c['source'],event_id=c['event_id'],gene=c['gene'],context_id=c['context_id'],
                     context_phase=c.get('phase'),window_start=start,window_end=end,sequence=vp.amino_acids,
                     target_start_in_window=max(0,lo-start),target_end_in_window=min(25,hi-start),target_mode=c['target_mode'],
                     native_affinity_sum=vp.target_epitope_score,excluded_self_affinity_sum=vp.self_epitope_score,
                     target_peptides=';'.join(target_peptides),distinct_target_peptides=len(target_peptides),
                     seed_peptides=';'.join(retained_seeds),distinct_seed_peptides=len(retained_seeds),
                     primary_window_eligible=vp.target_epitope_score>0 and bool(retained_seeds),
                     best_target_affinity_nm=min((p.value for e in vp.target_epitopes for p in e.predictions_flat() if p.kind=='pMHC_affinity'),default=None),
                     expression_proxy_left=c['expression_proxy_left'],expression_proxy_right=c['expression_proxy_right'],
                     expression_notes=c['expression_proxy_notes'],research_only=True)
            rec.update(c.get('normal_reference_annotations',{}))
            value=c['expression_proxy_left']
            rec['weighted_score_left']='' if value=='' else vp.target_epitope_score*math.sqrt(float(value))
            # Compatibility annotation only; no alternative partner ranking.
            rec['weighted_score_right']=rec['weighted_score_left']
            assert len(rec['sequence'])==25 and start<hi and end>lo
            allwindows.append(rec)
            if rec['primary_window_eligible']:
                positive+=1
        context_outcomes.append(dict(context_id=c['context_id'],source=c['source'],event_id=c['event_id'],gene=c['gene'],
                                     context_length=len(seq),target_occurrences=len(eps),candidate_windows=len(starts),positive_windows=positive,
                                     status='HAS_POSITIVE_25AA_WINDOWS' if positive else 'CONTEXT_SHORTER_THAN_25AA' if len(seq)<25 else 'NO_POSITIVE_TARGET_WINDOWS'))
        if ci%100==0:print('WINDOW_CONTEXT',ci,'/',len(contexts),'windows',len(allwindows),flush=True)
    write('all-scored-windows.tsv',allwindows);write('target-occurrences.tsv',occurrences);write('context-outcomes.tsv',context_outcomes)
    # Mutation windows are independently verified by the mutation stage.
    mutations=read(R/'inputs/mutation-regions.tsv')
    combined=[]
    for m in mutations:
        score=float(m['post_self_affinity_sum'])
        proxy=float(m['supplied_transformed_tpm'])*float(m['rna_alternate_fraction'])
        assert close(score*math.sqrt(proxy),float(m['weighted_score'])),m['variant_id']
        combined.append(dict(source='MUTATION',event_id=m['variant_id'],gene=m['gene'],
            region_id='REG_'+hashlib.sha256(m['sequence'].encode()).hexdigest()[:20],
            sequence=m['sequence'],native_affinity_sum=score,weighted_score_left=m['weighted_score'],
            weighted_score_right=m['weighted_score'],current_construct_region=True,
            prior_rank=m['rank'],research_only=True))
    write('mutation-regions.tsv',combined)
    from consolidate_regions import main as consolidate
    summary=consolidate(R)
    (R/'environment.freeze.txt').write_text(__import__('subprocess').check_output(['python3','-m','pip','freeze'],text=True))
    source_hashes={}
    for modname in ['vaxrank.epitope_dsl','vaxrank.candidate_epitope','vaxrank.vaccine_peptide','vaxrank.reference_proteome']:
        mod=__import__(modname,fromlist=['']);source_hashes[modname]=digest(Path(mod.__file__))
    save('native-source-hashes.json',source_hashes)
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':
    try:main()
    except Exception:
        import traceback
        save('status.json',dict(status='FAILED',error=traceback.format_exc()));raise
