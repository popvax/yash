"""Offline verification/replay of the recorded matched-normal lead design.

Replays selection from the recorded candidate pool, solves the recorded edge
costs again, checks chosen junctions against included independent predictions,
and reruns synonymous encoding. No variant calling, new inference or cloud jobs.
"""
from yash_vaccine.paths import resolve_input
from pathlib import Path
import csv
import importlib.util
import json
import math
import runpy
import shutil
import subprocess
import sys
import tempfile
from .io import rows, sha256, write_json
from .provenance import RecordedSource
from .hydrophobicity import summarize


def verify(repo):
    result = subprocess.run([sys.executable, str(Path(repo)/'designs/design-2/verify_results.py')],
                            cwd=repo, capture_output=True, text=True)
    if result.returncode:
        raise ValueError('Lead-design verification failed:\n'+result.stdout+result.stderr)
    return json.loads(result.stdout)


def _helpers(path):
    before = list(sys.path)
    prior = sys.modules.pop('common', None)
    try:
        sys.path.insert(0, str(path.parent))
        return runpy.run_path(str(path), run_name='recorded_replay_helpers')
    finally:
        sys.path[:] = before
        sys.modules.pop('common', None)
        if prior is not None:sys.modules['common'] = prior


def _table(path, values, reference):
    fields = next(csv.reader(resolve_input(reference).open(), delimiter='\t'))
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows({k:json.dumps(v) if isinstance(v,(list,dict)) else v for k,v in row.items()} for row in values)


def _selection(design, destination):
    pool = list(rows(design/'exports/integrated-candidate-pool.tsv'))
    selected=[];used=set()
    for source,quota in [('MUTATION',None),('SPLICE',5),('CTA/SELF',5),('FUSION',2)]:
        candidates=sorted([r for r in pool if r['source']==source],key=lambda r:(r['weighted_score_left']=='',-float(r['weighted_score_left'] or 0),r['event_id'],r['region_id']))
        genes=set();count=0
        for row in candidates:
            hydro=summarize(row['sequence'])
            if hydro['status']!='PASS':continue
            if source!='FUSION' and row['weighted_score_left']=='':continue
            if source=='SPLICE' and row['gene'] in genes:continue
            if row['sequence'] in used:continue
            if quota is not None and count>=quota:continue
            count+=1;genes.add(row['gene']);used.add(row['sequence'])
            selected.append(dict(row,antigen_id=f'AG{len(selected)+1:02d}',hydrophobicity=hydro))
    expected=list(rows(design/'selection/selected.tsv'))
    if [(r['region_id'],r['antigen_id']) for r in selected]!=[(r['region_id'],r['antigen_id']) for r in expected]:
        raise ValueError('Replayed lead selection differs')
    destination.mkdir()
    _table(destination/'selected.tsv',selected,design/'selection/selected.tsv')
    (destination/'selected.fasta').write_text(''.join('>'+r['antigen_id']+'|'+r['source']+'|'+r['gene']+'\n'+r['sequence']+'\n' for r in selected))
    return selected


def reproduce(repo, out=None):
    repo=Path(repo).resolve();design=repo/'designs/design-2'
    recorded_check=verify(repo)
    if out is None:
        with tempfile.TemporaryDirectory(prefix='yash-lead-replay-') as temporary:
            result=reproduce(repo,Path(temporary)/'design-2');result['output']=None;return result
    out=Path(out).resolve()
    if resolve_input(out).exists():raise ValueError('Choose a new reproduction output directory: '+str(out))
    out.mkdir(parents=True)
    selected=_selection(design,out/'selection')
    frozen=json.loads((resolve_input(design/'execution/frozen-inputs.json')).read_text())
    with tempfile.TemporaryDirectory(prefix='yash-recorded-source-') as temporary:
        source=Path(temporary)
        recorded_source=RecordedSource(repo)
        recorded_source.materialize(source,frozen['files'],allow_external_references=True)
        solve=_helpers(source/'designs/initial/construct/optimize.py')['solve']
        helpers=_helpers(design/'pipeline/cassette.py')
        edges=list(rows(design/'construct/protein/edge-costs.tsv'))
        for e in edges:
            for k in ['left_index','right_index','strong_rank','mild_rank','affinity_50','affinity_500','added_aa','peptide_occurrences']:e[k]=int(e[k])
            for k in ['min_rank','min_affinity_nm','max_pair_mean_kd']:e[k]=float(e[k])
            e['hydrophobicity_pass']=e['hydrophobicity_pass']=='True'
            if e['category']=='internal' and e['hydrophobicity_pass']!=(summarize(e['left']+e['linker_sequence']+e['right'])['status']=='PASS'):
                raise ValueError('Recorded internal-edge hydropathy mismatch')
        path,solver=solve(len(selected),[e for e in edges if e['hydrophobicity_pass']],seconds=600)
        prior=json.loads((resolve_input(design/'construct/protein/construct.json')).read_text())
        actual_objective=solver['objective_tuple'];recorded_objective=prior['solver']['objective_tuple']
        if (not solver['all_phases_proven_optimal']
                or list(actual_objective[:2])!=recorded_objective[:2]
                or actual_objective[3]!=recorded_objective[3]
                or not math.isclose(actual_objective[2],recorded_objective[2],rel_tol=1e-12,abs_tol=1e-12)):
            raise ValueError('Replayed lead solver objective differs')
        protein,components,joins=helpers['assemble'](path,selected,prior['SEC'],prior['MITD'])
    output=out/'construct/protein';output.mkdir(parents=True)
    (output/'cassette.fasta').write_text('>Yash_matched_normal_SEC_antigens_MITD\n'+protein+'\n')
    _table(output/'components.tsv',components,design/'construct/protein/components.tsv')
    _table(output/'junctions.tsv',joins,design/'construct/protein/junctions.tsv')
    by_id={r['region_id']:r for r in selected}
    ordered=[by_id[c['region_id']] for c in components if c['kind']=='antigen']
    _table(output/'antigen-order.tsv',[dict(position=i+1,antigen_id=r['antigen_id'],region_id=r['region_id'],gene=r['gene'],sequence=r['sequence']) for i,r in enumerate(ordered)],design/'construct/protein/antigen-order.tsv')
    # Recount full-chain queries against independently recorded chosen-junction
    # predictions. Values remain recorded model measurements, not new inference.
    fresh={(r['peptide'],r['allele']):(float(r['affinity_nm']),float(r['affinity_percentile'])) for r in rows(design/'construct/selected-inference/predictions.tsv')}
    alleles=json.loads((resolve_input(repo/'config/case.json')).read_text())['hla_class_I']
    scan=helpers['full_scan'](protein,joins)
    measured=[dict(**o,allele=a,affinity_nm=fresh[o['peptide'],a][0],affinity_percentile=fresh[o['peptide'],a][1]) for o in scan for a in alleles]
    totals=helpers['summarize_predictions'](measured)
    if totals!=json.loads((resolve_input(design/'construct/final-check.json')).read_text())['totals']:
        raise ValueError('Independent recorded junction recount differs')
    for name in ['construct/protein/construct.json','construct/architecture.json','construct/final-check.json']:
        target=out/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(resolve_input(design/name),target)
    mrna=out/'construct/mrna';mrna.mkdir()
    for name in ['run.py','optimizer.py','genbank.py','backbone.json','annotation-sources.json','human-codon-usage.json','requirements.txt']:
        shutil.copy2(resolve_input(design/'construct/mrna'/name),mrna/name)
    encoded=subprocess.run([sys.executable,str(mrna/'run.py'),'--output',str(mrna)],cwd=repo,capture_output=True,text=True)
    if encoded.returncode:raise ValueError('Lead RNA replay failed:\n'+encoded.stdout+encoded.stderr)
    checked=['selection/selected.tsv','selection/selected.fasta','construct/protein/cassette.fasta','construct/protein/components.tsv','construct/protein/junctions.tsv','construct/protein/antigen-order.tsv']
    checked += ['construct/mrna/'+name for name in json.loads((resolve_input(design/'construct/mrna/verification.json')).read_text())['output_sha256']]
    for name in checked:
        if sha256(out/name)!=sha256(design/name):raise ValueError('Replayed lead artifact differs: '+name)
    result=dict(status='VERIFIED',name='design-2',selected_regions=24,protein_length_aa=724,
                all_recorded_artifacts_identical=True,checked_artifacts=checked,
                all_solver_phases_optimal=True,polyA_status='unspecified',rna_length_nt=None,
                model_inference_rerun=False,cloud_jobs_started=False,
                source_snapshot_verification=recorded_check['source_snapshot_verification'],
                stages=['selection_from_recorded_pool','joint_order_linkers_from_recorded_costs','independent_recorded_junction_recount','RNA_encoding'],
                input_sha256={name:sha256(design/name) for name in ['exports/integrated-candidate-pool.tsv','construct/protein/edge-costs.tsv','construct/selected-inference/predictions.tsv']},output=str(out))
    write_json(out/'verification.json',result)
    return result
