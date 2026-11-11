"""Compare complete designs by stable region identity and recorded provenance."""
import json
from pathlib import Path
from .io import rows,sha256,write_json,write_tsv

def read(directory):
    path=Path(directory).resolve()
    if path.is_file():path=path.parent
    if (path/'selection/design.json').is_file():path=path/'selection'
    if any(not (path/name).is_file() for name in ('verification.json','design.json','construct.json','selected.tsv')):
        raise ValueError('This comparison helper expects outputs from our optional recipe runner. For independent designs, compare the relevant tables with your own script; see docs/usage.md.')
    verification=json.loads((path/'verification.json').read_text())
    if sha256(path/'design.json')!=verification['recipe_sha256']:
        raise ValueError('Recipe changed after construction: '+str(path)+'. Run it into a new output before comparing.')
    for name,digest in verification.get('output_sha256',verification['artifact_sha256']).items():
        if sha256(path/name)!=digest:raise ValueError('Recorded design output changed: '+str(path/name))
    selected=list(rows(path/'selected.tsv'));index={row['region_id']:row for row in selected}
    if len(index)!=len(selected):raise ValueError('Duplicate region_id in '+str(path/'selected.tsv'))
    return path,index,json.loads((path/'construct.json').read_text()),json.loads((path/'design.json').read_text())

def _compare_recipe(left,right,out=None):
    lp,li,lc,lr=read(left);rp,ri,rc,rr=read(right)
    common=set(li)&set(ri)
    for identity in common:
        if li[identity]['sequence']!=ri[identity]['sequence']:
            raise ValueError('Same region_id has different sequences: '+identity)
    changes=[]
    for identity in sorted(set(li)|set(ri)):
        a=li.get(identity,{});b=ri.get(identity,{});row=b or a
        changes.append(dict(region_id=identity,status='retained' if a and b else 'removed' if a else 'added',
                            source=row['source'],gene=row['gene'],protein_change=row.get('protein_change',''),
                            left_antigen_id=a.get('antigen_id',''),right_antigen_id=b.get('antigen_id',''),sequence=row['sequence']))
    differences=[]
    if lr['evidence']['sha256']!=rr['evidence']['sha256']:differences.append('evidence bundle')
    if lr['evidence']['context_sha256']!=rr['evidence']['context_sha256']:differences.append('case/HLA context')
    if lr['prediction_cache']['sha256']!=rr['prediction_cache']['sha256']:differences.append('prediction cache')
    if lr['linkers']!=rr['linkers']:differences.append('linker choices')
    if lr['optimization']['strategy']!=rr['optimization']['strategy']:differences.append('optimization strategy')
    if (lp/'run.json').exists() and (rp/'run.json').exists():
        left_code=json.loads((lp/'run.json').read_text())['code']['source_sha256']
        right_code=json.loads((rp/'run.json').read_text())['code']['source_sha256']
        if left_code!=right_code:differences.append('code snapshot')
    def summary(regions,construct):
        return dict(regions=len(regions),protein_length_aa=construct['protein_length_aa'],
                    linker_counts=construct['linker_counts'],objective_tuple=construct['objective_tuple'],
                    boundary_totals=construct['boundary_totals_including_start_M'])
    report=dict(left=str(lp),right=str(rp),matching_key='region_id',retained=len(common),
                added=len(set(ri)-set(li)),removed=len(set(li)-set(ri)),
                differing_inputs=differences,metrics_share_inputs=not differences,
                left_summary=summary(li,lc),right_summary=summary(ri,rc),regions=changes)
    text=['Comparison of recipe-stage cassettes (before SEC/MITD and RNA encoding): '+lr['name']+' → '+rr['name'],
          f'Regions: {len(common)} retained, {report["added"]} added, {report["removed"]} removed.',
          f'Cassette: {lc["protein_length_aa"]} → {rc["protein_length_aa"]} aa.',
          'Regions are matched by region_id; AG labels are local to each selection.']
    for key,label in [('strong_rank','Boundary/HLA hits at rank ≤0.5%'),('mild_rank','Boundary/HLA hits at rank ≤2%'),('added_aa','Added residues')]:
        text.append(f'{label}: {lc["boundary_totals_including_start_M"][key]} → {rc["boundary_totals_including_start_M"][key]}.')
    for name,result in [(lr['name'],lc),(rr['name'],rc)]:
        text.append('Linkers ('+name+'): '+', '.join(f'{key} × {value}' for key,value in sorted(result['linker_counts'].items()))+'.')
    for row in changes:
        if row['status']!='retained':text.append(row['status'].capitalize()+': '+row['gene']+' '+row['protein_change']+' ('+row['region_id']+')')
    if differences:text.append('Interpret metrics with the differing inputs in view: '+', '.join(differences)+'.')
    else:text.append('Both designs use the same evidence, prediction cache, linker choices, optimization strategy and code snapshot.')
    text.append('These are computational design comparisons, not measurements of vaccine performance.')
    report['text']='\n'.join(text)+'\n'
    if out is not None:
        out=Path(out)
        if out.exists():raise ValueError('Choose a new comparison directory: '+str(out))
        out.mkdir(parents=True);(out/'comparison.txt').write_text(report['text'])
        write_json(out/'comparison.json',report);write_tsv(out/'regions.tsv',changes)
    return report


def _complete(path):
    path=Path(path).resolve()
    return (path/'selection/selected.tsv').is_file() and (path/'construct/protein/cassette.fasta').is_file()


def _read_complete(path):
    from .provenance import RecordedSource
    path=Path(path).resolve();repo=path.parents[1]
    if path.name=='design-2':
        from .lead import verify
        source_scope=verify(repo)['source_snapshot_verification']
    elif path.name=='initial':
        manifest=json.loads((path/'candidate.json').read_text());recorded=RecordedSource(repo)
        recorded.verify_files({str((path/name).resolve().relative_to(repo)):dict(sha256=digest)
                               for name,digest in manifest['artifact_sha256'].items()})
        source_scope=recorded.verification_scope()
    else:raise ValueError('Complete-design comparison currently supports designs/initial and designs/design-2')
    selected=list(rows(path/'selection/selected.tsv'))
    index={r['region_id']:r for r in selected}
    if len(index)!=len(selected):raise ValueError('Duplicate region_id in '+str(path))
    construct=json.loads((path/'construct/protein/construct.json').read_text())
    check=json.loads((path/'construct/final-check.json').read_text())
    check['source_snapshot_verification']=source_scope
    if sha256(path/'construct/protein/cassette.fasta')!=check['protein_sha256']:
        raise ValueError('Protein checksum differs: '+str(path))
    return path,index,construct,check


def compare(left,right,out=None):
    """Compare complete recorded designs; preserve explicit recipe-stage scope."""
    if not (_complete(left) and _complete(right)):
        return _compare_recipe(left,right,out)
    lp,li,lc,lv=_read_complete(left);rp,ri,rc,rv=_read_complete(right)
    common=set(li)&set(ri);changes=[]
    for identity in sorted(set(li)|set(ri)):
        a=li.get(identity,{});b=ri.get(identity,{});row=b or a
        if a and b and a['sequence']!=b['sequence']:raise ValueError('Same region_id has different sequences')
        changes.append(dict(region_id=identity,status='retained' if a and b else 'removed' if a else 'added',source=row['source'],gene=row['gene'],protein_change=row.get('protein_change',''),left_antigen_id=a.get('antigen_id',''),right_antigen_id=b.get('antigen_id',''),sequence=row['sequence']))
    def summary(index,c,v):
        return dict(regions=len(index),protein_length_aa=c['protein_length_aa'],linker_counts=c['linker_counts'],junction_totals=v['totals'])
    report=dict(scope='Complete SEC–antigens–MITD proteins',left=str(lp),right=str(rp),matching_key='region_id',retained=len(common),added=len(set(ri)-set(li)),removed=len(set(li)-set(ri)),metrics_share_inputs=lp==rp,differing_inputs=[] if lp==rp else ['DNA evidence and normal/self sequence models','selected regions and junction prediction queries'],left_summary=summary(li,lc,lv),right_summary=summary(ri,rc,rv),regions=changes)
    report['source_snapshot_verification']={'left':lv['source_snapshot_verification'],'right':rv['source_snapshot_verification']}
    lines=[f'Complete recorded designs: {lp.name} → {rp.name}',f'Regions: {report["retained"]} retained, {report["added"]} added, {report["removed"]} removed.',f'Full protein: {lc["protein_length_aa"]} → {rc["protein_length_aa"]} aa.','Regions are matched by region_id; AG labels are local to each selection.']
    for row in changes:
        if row['status']!='retained':lines.append(row['status'].capitalize()+': '+row['gene']+' '+row['protein_change']+' ('+row['region_id']+')')
    if lp!=rp:lines.append('DNA/self evidence, selected antigens and prediction queries differ; junction counts are descriptive, not a controlled comparison of vaccine performance.')
    report['text']='\n'.join(lines)+'\n'
    if out is not None:
        out=Path(out)
        if out.exists():raise ValueError('Choose a new comparison directory: '+str(out))
        out.mkdir(parents=True);(out/'comparison.txt').write_text(report['text']);write_json(out/'comparison.json',report);write_tsv(out/'regions.tsv',changes)
    return report
