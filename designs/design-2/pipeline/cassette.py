"""Joint routed-cassette rebuild using the recorded recipe and new selected regions."""
from pathlib import Path
from collections import Counter
import csv
import importlib.util
import hashlib
import json
import math
import shutil
import subprocess
import sys
from common import save,receipt,sha,run
LENGTHS=[8,9,10,11]

def junction_peptides(left, linker, right):
    sequence = left[-10:] + linker + right[:10]
    begin, end = len(left[-10:]), len(left[-10:]) + len(linker)
    return [sequence[start:start+k] for k in LENGTHS for start in range(len(sequence)-k+1)
            if (start < end and start+k > begin if linker else start < begin < start+k)]

def full_scan(sequence, joins):
    """Independently scan full protein against placed linker intervals/direct joins."""
    found = []
    for length in LENGTHS:
        for start in range(len(sequence)-length+1):
            end = start+length
            touched = [j for j in joins if
                (start < j["linker_end_aa"] and end > j["linker_start_aa"]
                 if j["linker_end_aa"] > j["linker_start_aa"]
                 else start < j["linker_start_aa"] < end)]
            if len(touched) > 1:
                raise ValueError("One peptide spans multiple joins; pairwise objective is invalid")
            if touched:
                j = touched[0]
                found.append(dict(junction_index=j["junction_index"], category=j["category"],
                    start_aa=start, end_aa=end, peptide=sequence[start:end]))
    return found

def assemble(path, regions, sec, mitd):
    sequence, components, joins = "", [], []
    def add(kind, name, aa, source=None):
        nonlocal sequence
        source = source or {}
        components.append(dict(kind=kind, name=name, aa=aa, start_aa=len(sequence),
            end_aa=len(sequence)+len(aa), region_id=source.get("region_id", ""),
            event_id=source.get("event_id", ""), gene=source.get("gene", "")))
        sequence += aa
    add("signal_peptide", "HLA_C_01_02_SEC", sec)
    for i, edge in enumerate(path):
        begin = len(sequence)
        add("linker", edge["linker"], edge["linker_sequence"])
        end = len(sequence)
        joins.append(dict(junction_index=i, category=edge["category"],
            left_id=edge["left_id"], right_id=edge["right_id"], linker=edge["linker"],
            linker_sequence=edge["linker_sequence"], linker_start_aa=begin, linker_end_aa=end,
            left=edge["left"], right=edge["right"]))
        if edge["right_index"] == len(regions):
            add("mitd", "HLA_B_15_08_MITD", mitd)
        else:
            r = regions[edge["right_index"]]
            add("antigen", r["antigen_id"], r["sequence"], r)
    assert sequence.startswith(sec) and sequence.endswith(mitd)
    assert len(sequence) == sum(len(r["sequence"]) for r in regions)+len(sec)+len(mitd)+sum(e["added_aa"] for e in path)
    assert sorted((p["region_id"],p["aa"]) for p in components if p["kind"]=="antigen") == \
           sorted((r["region_id"],r["sequence"]) for r in regions)
    assert all(sequence[p["start_aa"]:p["end_aa"]] == p["aa"] for p in components)
    return sequence, components, joins

def summarize_predictions(predictions):
    return dict(peptide_occurrences=len(predictions)//6, peptide_hla_rows=len(predictions),
        strong_rank=sum(r["affinity_percentile"]<=.5 for r in predictions),
        mild_rank=sum(r["affinity_percentile"]<=2 for r in predictions),
        affinity_50=sum(r["affinity_nm"]<=50 for r in predictions),
        affinity_500=sum(r["affinity_nm"]<=500 for r in predictions),
        min_rank=min(r["affinity_percentile"] for r in predictions),
        min_affinity_nm=min(r["affinity_nm"] for r in predictions))

def assess(sequence, joins, path, alleles, lookup):
    occurrences = full_scan(sequence, joins)
    predicted = [dict(**row, allele=allele, affinity_nm=lookup(row["peptide"],h)[0],
                      affinity_percentile=lookup(row["peptide"],h)[1])
                 for row in occurrences for h,allele in enumerate(alleles)]
    for i, edge in enumerate(path):
        actual = Counter(r["peptide"] for r in occurrences if r["junction_index"] == i)
        expected = Counter(junction_peptides(edge["left"],edge["linker_sequence"],edge["right"]))
        if actual != expected:
            raise ValueError("Full-protein scan disagrees with pairwise junction enumeration")
    totals = summarize_predictions(predicted)
    for key in ["strong_rank", "mild_rank", "affinity_50", "affinity_500", "peptide_occurrences"]:
        if totals[key] != sum(e[key] for e in path):
            raise ValueError("Full-chain cost mismatch: " + key)
    by_category = [dict(category=category, **summarize_predictions([r for r in predicted if r["category"]==category]))
                   for category in ["SEC", "internal", "MITD"]]
    return predicted, totals, by_category


def table(path,values):
    with path.open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(values[0]),delimiter='\t');writer.writeheader();writer.writerows(values)


def data(root):
    selected=json.loads((root/'integrate/selected.json').read_text())
    frozen=root/'frozen';architecture=json.loads((frozen/'designs/initial/construct/architecture.json').read_text())
    alleles=json.loads((root/'case.json').read_text())['hla_class_I']
    assert len(alleles)==len(set(alleles))==6
    assert all(e['patient_allele'] in alleles for e in architecture['elements'].values())
    linkers=json.loads((frozen/'config/linkers.json').read_text())['choices']
    return selected,linkers,architecture['elements']['SEC']['aa'],architecture['elements']['MITD']['aa'],alleles,architecture


def make_request(root,joins,scope):
    frozen=root/'frozen';alleles=data(root)[4]
    expected=json.loads((frozen/'data/checkpoints/provenance/construct-inference.json').read_text())['versions']
    return dict(scope=scope,expected_versions=expected,lengths=LENGTHS,alleles=alleles,joins=joins,
        input_sha256=dict(selected=sha(root/'integrate/selected.tsv'),case=sha(root/'case.json'),
            architecture=sha(frozen/'designs/initial/construct/architecture.json'),linkers=sha(frozen/'config/linkers.json')),
        peptides=sorted({p for join in joins for p in junction_peptides(join['left'],join['linker'],join['right'])}))


def prepare(root):
    done=json.loads((root/'integrate/complete.json').read_text());assert done['status']=='COMPLETED'
    out=root/'construct/prepare';out.mkdir(parents=True,exist_ok=True)
    regions,linkers,sec,mitd,alleles,architecture=data(root);n=len(regions)
    joins=[dict(name=f'{i}:{j}:{linker["name"]}',left=regions[i]['sequence'],right=regions[j]['sequence'],linker=linker['sequence'])
           for i in range(n) for j in range(n) if i!=j for linker in linkers]
    joins += [dict(name=f'SEC:{i}:{linker["name"]}',left=sec,right=regions[i]['sequence'],linker=linker['sequence']) for i in range(n) for linker in linkers]
    joins += [dict(name=f'{i}:MITD:{linker["name"]}',left=regions[i]['sequence'],right=mitd,linker=linker['sequence']) for i in range(n) for linker in linkers]
    request=make_request(root,joins,'Fresh affinity-only inference for every internal and terminal option of the matched-normal selection')
    save(out/'request.json',request);save(root/'construct/architecture.json',architecture)
    receipt(out,[out/'request.json'],regions=n,junction_options=len(joins),peptides=len(request['peptides']),peptide_HLA_pairs=len(request['peptides'])*len(alleles),worker_sha256=sha(__file__))


def predict(root,selected=False):
    name='selected-inference' if selected else 'all-inference';out=root/'construct'/name
    parent=root/'construct/protein' if selected else root/'construct/prepare'
    request=parent/('verification-request.json' if selected else 'request.json')
    frozen=root/'frozen';worker=frozen/'designs/initial/construct/predict.py'
    assert not (out/'predictions.tsv').exists(),'Prior inference outputs require explicit review, not overwrite'
    logs=root/'construct'/(name+'-logs')
    run(logs,'fresh-affinity',[sys.executable,worker,'--inputs',request,
        '--models','/data/results/TN25-112103/openvax-mutations-20261018/references/mhcflurry/2.3.0/models_class1_presentation/models',
        '--model-manifest',frozen/'data/checkpoints/provenance/construct-model-files.json','--output',out])
    meta=json.loads((out/'inference.json').read_text());assert meta['all_queries_freshly_predicted']
    receipt(out,[out/'predictions.tsv',out/'inference.json'],request_sha256=sha(request),
            native_worker_sha256=sha(worker),worker_sha256=sha(__file__),unique_peptides=meta['unique_peptides'],peptide_HLA_pairs=meta['peptide_hla_pairs'])


def predictions(root,selected=False):
    directory=root/'construct'/('selected-inference' if selected else 'all-inference')
    request=root/'construct'/('protein/verification-request.json' if selected else 'prepare/request.json')
    req=json.loads(request.read_text());done=json.loads((directory/'complete.json').read_text());assert done['status']=='COMPLETED'
    meta=json.loads((directory/'inference.json').read_text())
    assert meta['input_sha256']==sha(request) and meta['predictions_sha256']==sha(directory/'predictions.tsv')
    assert meta['versions']==req['expected_versions'] and meta['all_model_files_match_reference'] and meta['native_junction_enumeration_matches']
    result={}
    with (directory/'predictions.tsv').open() as stream:
        for row in csv.DictReader(stream,delimiter='\t'):
            key=row['peptide'],row['allele'];assert key not in result
            a,r=float(row['affinity_nm']),float(row['affinity_percentile']);assert math.isfinite(a) and a>0 and math.isfinite(r) and 0<=r<=100
            result[key]=(a,r)
    assert set(result)=={(p,a) for p in req['peptides'] for a in req['alleles']}
    return result


def solve(root):
    out=root/'construct/protein';out.mkdir(parents=True,exist_ok=True)
    frozen=root/'frozen';sys.path.insert(0,str(frozen/'src'))
    from yash_vaccine.hydrophobicity import summarize,write_protein_screen
    spec=importlib.util.spec_from_file_location('recorded_joint_optimizer',frozen/'designs/initial/construct/optimize.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    regions,linkers,sec,mitd,alleles,architecture=data(root);values=predictions(root);n=len(regions)
    lookup=lambda p,h:values[p,alleles[h]]
    def edge(i,j,linker):
        left=sec if i==n else regions[i]['sequence'];right=mitd if j==n else regions[j]['sequence']
        queries=junction_peptides(left,linker['sequence'],right);scores=[lookup(q,h) for q in queries for h in range(len(alleles))]
        category='SEC' if i==n else 'MITD' if j==n else 'internal';hydro=summarize(left+linker['sequence']+right)
        return dict(left_index=i,right_index=j,left_id='SEC' if i==n else regions[i]['region_id'],right_id='MITD' if j==n else regions[j]['region_id'],
            category=category,left=left,right=right,linker=linker['name'],linker_sequence=linker['sequence'],
            strong_rank=sum(r<=.5 for a,r in scores),mild_rank=sum(r<=2 for a,r in scores),affinity_50=sum(a<=50 for a,r in scores),affinity_500=sum(a<=500 for a,r in scores),
            min_rank=min(r for a,r in scores),min_affinity_nm=min(a for a,r in scores),added_aa=len(linker['sequence']),peptide_occurrences=len(queries),
            max_pair_mean_kd=hydro['max_mean_kd'],hydrophobicity_pass=category!='internal' or hydro['status']=='PASS')
    edges=[edge(i,j,l) for i in range(n) for j in range(n) if i!=j for l in linkers]+[edge(n,i,l) for i in range(n) for l in linkers]+[edge(i,n,l) for i in range(n) for l in linkers]
    eligible=[e for e in edges if e['hydrophobicity_pass']]
    path,solver=module.solve(n,eligible,seconds=600)
    sequence,parts,joins=assemble(path,regions,sec,mitd)
    predicted,totals,categories=assess(sequence,joins,path,alleles,lookup)
    hydro=write_protein_screen(out,sequence,parts)
    (out/'cassette.fasta').write_text('>Yash_matched_normal_SEC_antigens_MITD\n'+sequence+'\n')
    for name,values_ in [('components.tsv',parts),('junctions.tsv',joins),('boundary-predictions.tsv',predicted),('edge-costs.tsv',edges),('by-category.tsv',categories)]:table(out/name,values_)
    order=[regions[e['right_index']] for e in path[:-1]]
    table(out/'antigen-order.tsv',[dict(position=i+1,antigen_id=r['antigen_id'],region_id=r['region_id'],gene=r['gene'],sequence=r['sequence']) for i,r in enumerate(order)])
    save(out/'construct.json',dict(name='matched-normal-20261107',regions=n,protein_length_aa=len(sequence),solver=solver,
        all_options_freshly_predicted=True,hydrophobicity=hydro,totals=totals,categories=categories,
        considered_options=len(edges),hydrophobicity_excluded_options=len(edges)-len(eligible),
        SEC=sec,MITD=mitd,linker_counts=dict(Counter(e['linker'] for e in path)),selected_regions_sha256=sha(root/'integrate/selected.tsv')))
    request=make_request(root,[dict(name=str(j['junction_index']),left=j['left'],linker=j['linker_sequence'],right=j['right']) for j in joins],
        'Independent fresh re-inference of every selected junction after optimization')
    save(out/'verification-request.json',request)
    receipt(out,[out/name for name in ['cassette.fasta','components.tsv','junctions.tsv','boundary-predictions.tsv','edge-costs.tsv','by-category.tsv','antigen-order.tsv','construct.json','verification-request.json','hydrophobicity-windows.tsv','hydrophobicity.json']],
            all_antigens_exactly_once=True,independent_full_chain_scan=True,all_solver_phases_optimal=solver['all_phases_proven_optimal'],worker_sha256=sha(__file__))


def verify(root):
    out=root/'construct/verify';out.mkdir(parents=True,exist_ok=True);protein=root/'construct/protein'
    expected=json.loads((protein/'complete.json').read_text())
    for name,metadata in expected['outputs'].items():assert sha(protein/name)==metadata['sha256']
    regions,_,sec,mitd,alleles,_=data(root);fresh=predictions(root,selected=True)
    sequence=''.join(line for line in (protein/'cassette.fasta').read_text().splitlines() if not line.startswith('>'))
    with (protein/'components.tsv').open() as stream:parts=list(csv.DictReader(stream,delimiter='\t'))
    with (protein/'junctions.tsv').open() as stream:joins=list(csv.DictReader(stream,delimiter='\t'))
    for j in joins:
        for key in ['junction_index','linker_start_aa','linker_end_aa']:j[key]=int(j[key])
        a,b=j['linker_start_aa'],j['linker_end_aa'];assert sequence[a:b]==j['linker_sequence']
        assert sequence[a-10:a]==j['left'][-10:] and sequence[b:b+10]==j['right'][:10]
    assert sequence.startswith(sec) and sequence.endswith(mitd) and len(joins)==len(regions)+1
    assert sorted((p['region_id'],p['aa']) for p in parts if p['kind']=='antigen')==sorted((r['region_id'],r['sequence']) for r in regions)
    sys.path.insert(0,str(root/'frozen/src'))
    from yash_vaccine.hydrophobicity import protein_screen
    _,hydro=protein_screen(sequence,parts);assert hydro['status']=='PASS'
    scan=full_scan(sequence,joins);actual=[dict(**o,allele=a,affinity_nm=fresh[o['peptide'],a][0],affinity_percentile=fresh[o['peptide'],a][1]) for o in scan for a in alleles]
    with (protein/'boundary-predictions.tsv').open() as stream:prior=list(csv.DictReader(stream,delimiter='\t'))
    assert len(prior)==len(actual);aff_diffs=[];rank_diffs=[]
    for before,after in zip(prior,actual):
        for key in ['junction_index','category','start_aa','end_aa','peptide','allele']:assert str(before[key])==str(after[key])
        a,r=float(before['affinity_nm']),float(before['affinity_percentile']);na,nr=after['affinity_nm'],after['affinity_percentile']
        assert [a<=50,a<=500,r<=.5,r<=2]==[na<=50,na<=500,nr<=.5,nr<=2]
        aff_diffs.append(abs(na/a-1));rank_diffs.append(abs(nr-r))
    assert max(aff_diffs)<1e-4 and max(rank_diffs)<1e-4
    result=dict(status='FRESH_FULL_JUNCTION_CHECK_PASSED',protein_sha256=sha(protein/'cassette.fasta'),components_sha256=sha(protein/'components.tsv'),
        protein_length_aa=len(sequence),antigen_regions=len(regions),all_antigens_exactly_once=True,unique_peptide_hla_pairs=len(fresh),
        max_relative_affinity_difference=max(aff_diffs),max_absolute_rank_difference=max(rank_diffs),threshold_classifications_identical=True,
        hydrophobicity=hydro,totals=summarize_predictions(actual),scope='Predicted binding and hydropathy checks; no experimental processing or presentation validation')
    save(root/'construct/final-check.json',result);save(out/'final-check.json',result)
    receipt(out,[out/'final-check.json'],worker_sha256=sha(__file__),all_checks_pass=True)


def encode(root):
    check=json.loads((root/'construct/final-check.json').read_text());assert check['status']=='FRESH_FULL_JUNCTION_CHECK_PASSED'
    frozen=root/'frozen';workspace=root/'rebuild';mrna=workspace/'construct/mrna';mrna.mkdir(parents=True,exist_ok=True)
    original=frozen/'designs/initial/construct/mrna'
    for name in ['run.py','optimizer.py','genbank.py','backbone.json','annotation-sources.json','human-codon-usage.json','requirements.txt']:
        shutil.copyfile(original/name,mrna/name)
    for name in ['run.py','genbank.py']:
        path=mrna/name
        text=path.read_text().replace('Design #1','matched-normal design').replace('YASH_INITIAL','YASH_MATCHED').replace('Yash_Initial','Yash_Matched')
        # Long unbroken paths are split inside GenBank /note qualifiers by the
        # serializer. Use the checksum-verified annotation-input key instead.
        # Routing evidence remains the original measured source; do not relabel
        # it as new matched-normal HLA typing or invent a new evidence path.
        text=text.replace('Generated from designs/initial/construct/mrna/run.py and genbank.py.',
                          'Generated by this design\'s run.py and genbank.py.')
        text=text.replace('records: designs/initial/selection/selected.tsv.',
                          'records: selected_regions in annotation-sources.json.')
        text=text.replace('Source record: designs/initial/selection/selected.tsv',
                          'Source record: selected_regions in annotation-sources.json')
        path.write_text(text)
    (workspace/'selection').mkdir(exist_ok=True)
    shutil.copyfile(root/'integrate/selected.tsv',workspace/'selection/selected.tsv')
    for name in ['cassette.fasta','components.tsv']:
        path=workspace/'construct/protein'/name;path.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(root/'construct/protein'/name,path)
    shutil.copyfile(root/'construct/architecture.json',workspace/'construct/architecture.json')
    shutil.copyfile(root/'construct/final-check.json',workspace/'construct/final-check.json')
    cfg=json.loads((mrna/'backbone.json').read_text());cfg.update(name='matched-normal-poprna-mrna',protein_length_aa=check['protein_length_aa'],
        protein_file_sha256=sha(workspace/'construct/protein/cassette.fasta'),protein_components_sha256=sha(workspace/'construct/protein/components.tsv'))
    protein=''.join(line for line in (workspace/'construct/protein/cassette.fasta').read_text().splitlines() if not line.startswith('>'))
    settings=cfg['optimization'];candidate=mrna/'prepared-coding.txt'
    run(mrna/'execution','record-new-encoding',[sys.executable,mrna/'optimizer.py','--protein','-i',protein,'--output',candidate,
        '--organism',settings['species_table'],'--cutoff',settings['codon_cutoff'],'--seed',settings['seed'],
        '--fold-iterations',settings['structure_iterations_per_pass'],'--sequences',*settings['avoid_sequences_rna']])
    coding=candidate.read_text().strip();from Bio.Seq import Seq
    assert set(coding)<=set('ACGT') and str(Seq(coding).translate())==protein
    cfg['optimization']['recorded_coding_sequence_sha256']=hashlib.sha256(coding.encode()).hexdigest();save(mrna/'backbone.json',cfg)
    annotations=json.loads((mrna/'annotation-sources.json').read_text())
    for key in ['selected_regions','routing_architecture']:
        annotations[key]['sha256']=sha(mrna/annotations[key]['path'])
    annotations['record_date']='07-NOV-2026';save(mrna/'annotation-sources.json',annotations)
    out=root/'construct/mrna'
    run(mrna/'execution','verified-encoding',[sys.executable,mrna/'run.py','--output',out])
    verification=json.loads((out/'verification.json').read_text())
    assert verification['all_constraints_pass'] and verification['exact_parent_protein_translation']
    assert verification['polyA_status']=='unspecified' and not verification['complete_construct_verified']
    receipt(out,[p for p in out.iterdir() if p.is_file() and p.name not in ['complete.json','status.json']],
        worker_sha256=sha(__file__),protein_length_aa=len(protein),new_synonymous_encoding=True,exact_translation_preserved=True,
        seed=settings['seed'],species_table=settings['species_table'],codon_cutoff=settings['codon_cutoff'],
        original_backbone_json_sha256=sha(original/'backbone.json'),new_backbone_json_sha256=sha(mrna/'backbone.json'),
        PolyA_status='unspecified',complete_molecular_sequence=False,source_adaptation='Same pinned encoder and GenBank verifier; new protein, selection hashes and design labels')


if __name__=='__main__':
    root=Path(sys.argv[1]);action=sys.argv[2]
    if action=='prepare':prepare(root)
    elif action=='predict':predict(root)
    elif action=='solve':solve(root)
    elif action=='predict-selected':predict(root,selected=True)
    elif action=='verify':verify(root)
    elif action=='encode':encode(root)
    else:raise ValueError(action)

