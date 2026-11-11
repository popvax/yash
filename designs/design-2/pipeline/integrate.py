"""Combine fresh paired mutation regions with rechecked recorded RNA sources."""
from pathlib import Path
import collections
import csv
import gzip
import hashlib
import json
import math
import sys
import tarfile

import numpy as np
from common import save,receipt,sha
from self_search import patient_self_matches


def rows(path):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as stream:return list(csv.DictReader(stream,delimiter='\t'))


def table(path,values):
    fields=list(dict.fromkeys(key for row in values for key in row)) or ['source','status']
    with path.open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,delimiter='\t');writer.writeheader()
        for row in values:writer.writerow({k:json.dumps(v) if isinstance(v,(list,dict,tuple)) else v for k,v in row.items()})


def contribution(nm):
    value=0 if nm>=5000 else (1+math.exp(-350/150))/(1+math.exp((nm-350)/150))
    return value if value>=1e-5 else 0.


def restore(root):
    manifest=json.loads((root/'inputs/frozen-inputs.json').read_text());out=root/'frozen'
    with tarfile.open(root/'inputs/frozen-inputs.tar.gz') as archive:
        assert {m.name for m in archive.getmembers()}==set(manifest['files'])
        for member in archive.getmembers():
            assert member.isfile() and not Path(member.name).is_absolute() and '..' not in Path(member.name).parts
            path=out/member.name;path.parent.mkdir(parents=True,exist_ok=True)
            with archive.extractfile(member) as source:path.write_bytes(source.read())
            assert path.stat().st_size==manifest['files'][member.name]['bytes']
            assert sha(path)==manifest['files'][member.name]['sha256']
    return out,manifest


def main():
    root=Path(sys.argv[1]);out=root/'integrate';out.mkdir(parents=True,exist_ok=True)
    assert json.loads((root/'rank/complete.json').read_text())['status']=='COMPLETED'
    frozen,manifest=restore(root);sys.path.insert(0,str(frozen/'src'))
    from yash_vaccine.hydrophobicity import summarize
    from yash_vaccine.io import ranked as original_ranked
    alleles=json.loads((root/'case.json').read_text())['hla_class_I']
    native=rows(frozen/'data/checkpoints/lens-scored-windows.tsv.gz')
    lens=[row for row in native if row['source']!='SPLICE' and row['primary_window_eligible']=='True']
    splice=[row for row in rows(frozen/'data/splice-evidence/normal-screen/all-windows.tsv.gz') if row['primary_window_eligible']=='True']
    assert all(r['normal_detection']=='NOT_DETECTED_IN_SEARCHED_REFERENCE' for r in splice)
    occurrences=collections.defaultdict(list)
    for row in rows(frozen/'data/checkpoints/target-occurrences.tsv.gz'):
        if row['source']=='FUSION' and row['excluded_exact_self']!='True':occurrences[row['context_id']].append(row)
    fusion_peptides=sorted({r['peptide'] for group in occurrences.values() for r in group})
    from mhctools import MHCflurry
    models='/data/results/TN25-112103/openvax-mutations-20261018/references/mhcflurry/2.3.0/models_class1_presentation/models'
    fresh=MHCflurry(alleles=alleles,default_peptide_lengths=[8,9,10,11],models_path=models).predict_peptides(fusion_peptides).to_dataframe()
    assert len(fresh)==len(fusion_peptides)*6 and not fresh.duplicated(['peptide','allele']).any()
    assert np.isfinite(fresh.affinity).all() and (fresh.affinity>0).all()
    fusion_aff={(r['peptide'],r['allele']):float(r['affinity']) for r in fresh.to_dict('records')}
    fresh.to_csv(out/'fresh-fusion-affinities.tsv',sep='\t',index=False)
    splice_aff={(r['peptide'],r['allele']):float(r['affinity_nm']) for r in rows(frozen/'data/splice-evidence/scoring/affinities.tsv.gz')}
    request=json.loads((frozen/'data/fusion-evidence/inframe-request.json').read_text())
    admission=json.loads((frozen/'data/fusion-evidence/inframe-selection.json').read_text())
    review=json.loads((frozen/'data/fusion-evidence/inframe-binding-review.json').read_text())
    inframe_nonself={r['peptide'] for r in review['predictions'] if not r['exact_reference_match']}
    wanted=set(fusion_peptides)|inframe_nonself|{p for r in splice for p in r['target_peptides'].split(';') if p}
    self_matches=patient_self_matches(root,wanted);save(out/'RNA-source-patient-self-matches.json',self_matches)
    def retained(peptide):return not self_matches[peptide]['occurs_in_searched_models']
    rescored=[]
    for row in lens:
        value=dict(row);source=row['source']
        if source=='CTA/SELF':
            value.update(patient_self_exclusion_applied=False,normal_sequence_interpretation='Intentional unmutated CTA target; normal DNA is not expression evidence')
        elif source=='FUSION':
            a,b=int(row['window_start']),int(row['window_end'])
            before=[r for r in occurrences[row['context_id']] if int(r['start'])>=a and int(r['end'])<=b]
            assert math.isclose(sum(float(r['affinity_contribution']) for r in before),float(row['native_affinity_sum']),rel_tol=1e-8,abs_tol=1e-8)
            for r in before:
                expected=float(r['affinity_contribution']);actual=sum(contribution(fusion_aff[r['peptide'],a]) for a in alleles)
                assert math.isclose(actual,expected,rel_tol=1e-4,abs_tol=1e-6),(r['peptide'],actual,expected)
            after=[r for r in before if retained(r['peptide'])]
            score=sum(sum(contribution(fusion_aff[r['peptide'],a]) for a in alleles) for r in after)
            seeds=set(row['seed_peptides'].split(';'))&{r['peptide'] for r in after}
            best=min((fusion_aff[p,a] for p in seeds for a in alleles),default=None)
            value.update(native_affinity_sum=score,best_target_affinity_nm=best,
                primary_window_eligible=best is not None and best<=500,
                weighted_score_left=score*math.sqrt(float(row['expression_proxy_left'])) if row['expression_proxy_left'] else '',
                patient_self_exclusion_applied=True,patient_self_excluded_occurrences=len(before)-len(after),
                original_affinity_sum=row['native_affinity_sum'],fusion_affinity_source='Fresh matched-normal branch reinference with original-model parity check')
        else:raise ValueError(source)
        rescored.append(value)
    for row in splice:
        sequence=row['sequence'];lo,hi=int(row['target_start_in_window']),int(row['target_end_in_window'])
        targets=set(row['target_peptides'].split(';'))-{''}
        before=[sequence[i:i+k] for k in [8,9,10,11] for i in range(len(sequence)-k+1)
                if (i<hi and i+k>lo if hi>lo else i<lo<i+k) and sequence[i:i+k] in targets]
        original=sum(contribution(splice_aff[p,a]) for p in before for a in alleles)
        assert math.isclose(original,float(row['native_affinity_sum']),rel_tol=1e-8,abs_tol=1e-8)
        after=[p for p in before if retained(p)];score=sum(contribution(splice_aff[p,a]) for p in after for a in alleles)
        best=min((splice_aff[p,a] for p in after for a in alleles),default=None);expression=float(row['supplied_transformed_gene_tpm'])
        rescored.append(dict(row,context_id=row['window_id'],window_start=0,window_end=25,
            native_affinity_sum=score,weighted_score_left=score*math.sqrt(expression),expression_proxy_left=expression,
            best_target_affinity_nm=best,primary_window_eligible=best is not None and best<=500,
            original_affinity_sum=original,patient_self_exclusion_applied=True,
            patient_self_excluded_occurrences=len(before)-len(after),affinity_source='Recorded fresh patient-RNA splice inference; unchanged model and HLA alleles'))
        rescored[-1].update(patient_RNA_sequence_confirmed=True,upstream_observed_stop_codons_excluded=True,
                            annotated_upstream_frame=True,splice_fraction_used_for_weighting=False)
    a,b=admission['source_protein_window_0based'];sequence=request['source_sequence'][a:b]
    possibilities=[r for r in review['predictions'] if r['peptide'] in sequence and not r['exact_reference_match'] and retained(r['peptide'])]
    best=min(possibilities,key=lambda r:r['affinity_nm']) if possibilities else None
    event='CDH1::NRG1:chr16:68819425:+--chr8:32727949:+:INFRAME'
    rescored.append(dict(source='FUSION',event_id=event,gene='CDH1::NRG1',sequence=sequence,
        region_id='REG_'+hashlib.sha256((event+':'+str(a)+':'+sequence).encode()).hexdigest()[:20],
        window_id='CDH1_NRG1_INFRAME',context_id='CDH1_NRG1_INFRAME',window_start=a,window_end=b,
        manual_inclusion=True,admission_exception='DRIVER_JUNCTION_ABOVE_500NM',native_affinity_sum='',weighted_score_left='',
        expression_proxy_left='',best_target_affinity_nm=best['affinity_nm'] if best else None,primary_window_eligible=best is not None,
        patient_self_exclusion_applied=True,original_junction_read_count=admission['original_junction_read_count']))
    # Select best window per RNA event before the existing 19-aa hydropathy gate.
    event_bests={}
    for row in rescored:
        if str(row['primary_window_eligible'])!='True':continue
        key=row['source'],row['event_id'];current=event_bests.get(key)
        order=lambda r:(-float(r['native_affinity_sum'] or 0),r.get('window_id',''))
        if current is None or order(row)<order(current):event_bests[key]=row
    mutation_evidence={r['variant_id']:r for r in json.loads((root/'annotate/variant-evidence.json').read_text())}
    pool=list(event_bests.values())
    for row in json.loads((root/'rank/primary-ranked-regions.json').read_text()):
        evidence=mutation_evidence[row['variant_id']]
        pool.append(dict(row,source='MUTATION',event_id=row['variant_id'],
            region_id='REG_'+hashlib.sha256(row['sequence'].encode()).hexdigest()[:20],
            protein_change=row['DNA_annotation_protein_change'],normalized_allele=evidence['normalized_allele'],
            native_affinity_sum=row['affinity_contribution_sum'],weighted_score_left=row['weighted_score'],
            best_target_affinity_nm=row['best_target_affinity_nm'],primary_window_eligible=True))
    selected=[];audit=[];used=set()
    for source,quota in [('MUTATION',None),('SPLICE',5),('CTA/SELF',5),('FUSION',2)]:
        # Preserve the recorded recipe's event-id tie break. A different tie
        # order can swap equally scoring RNA events without new evidence.
        candidates=original_ranked([r for r in pool if r['source']==source])
        genes=set();count=0
        for row in candidates:
            hydro=summarize(row['sequence']);reason='SELECTED'
            if hydro['status']!='PASS':reason='EXCLUDE_HYDROPHOBIC_REGION'
            elif source!='FUSION' and row['weighted_score_left']=='':reason='MISSING_EXPRESSION'
            elif source=='SPLICE' and row['gene'] in genes:reason='GENE_ALREADY_SELECTED'
            elif row['sequence'] in used:reason='DUPLICATE_EXACT_SEQUENCE'
            elif quota is not None and count>=quota:reason='QUOTA_FILLED'
            audit.append(dict(source=source,event_id=row['event_id'],region_id=row['region_id'],gene=row['gene'],action=reason,hydrophobicity=hydro))
            if reason!='SELECTED':continue
            count+=1;genes.add(row['gene']);used.add(row['sequence'])
            selected.append(dict(row,antigen_id=f'AG{len(selected)+1:02d}',hydrophobicity=hydro,
                                 protein_change=row.get('protein_change','')))
    assert selected and all(len(r['sequence'])==25 for r in selected)
    original=rows(frozen/'designs/initial/selection/selected.tsv')
    before={(r['source'],r['gene'],r['sequence']) for r in original};after={(r['source'],r['gene'],r['sequence']) for r in selected}
    comparison=[dict(source=s,gene=g,sequence=q,change='RETAINED' if (s,g,q) in before&after else 'REMOVED' if (s,g,q) in before else 'ADDED') for s,g,q in sorted(before|after)]
    for name,values in [('RNA-all-rescored-windows.tsv',rescored),('candidate-pool.tsv',pool),('selected.tsv',selected),('selection-audit.tsv',audit),('comparison.tsv',comparison)]:table(out/name,values)
    save(out/'selected.json',selected);save(out/'reused-evidence-manifest.json',manifest)
    (out/'selected.fasta').write_text(''.join('>'+r['antigen_id']+'|'+r['source']+'|'+r['gene']+'\n'+r['sequence']+'\n' for r in selected))
    receipt(out,[out/name for name in ['RNA-all-rescored-windows.tsv','candidate-pool.tsv','selected.tsv','selected.json','selected.fasta','selection-audit.tsv','comparison.tsv','reused-evidence-manifest.json','RNA-source-patient-self-matches.json','fresh-fusion-affinities.tsv']],
        regions=len(selected),source_counts=dict(collections.Counter(r['source'] for r in selected)),
        mutation_quota=None,unchanged_sequences=len(before&after),removed=len(before-after),added=len(after-before),
        reused_RNA_source_discovery=True,new_matched_normal_self_search=True,CTA_self_targets_intentionally_retained=True,
        selection_order_source='Frozen initial src/yash_vaccine/io.py ranked',
        selection_order_source_sha256=sha(frozen/'src/yash_vaccine/io.py'),
        normal_DNA_does_not_measure_normal_splicing=True,construct_rebuilt=False,worker_sha256=sha(__file__))


if __name__=='__main__':main()
