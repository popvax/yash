"""Independent design workspaces and optional helpers for our recipe format."""
import collections,copy,importlib.metadata,json,os,re,shutil,subprocess,sys,tempfile
from pathlib import Path
from .io import rows,sha256,write_json
from .settings import SELECTION_KEYS,validate_selection

INITIAL=Path('designs/initial/selection/design.json')
EXPORTS={
    'selection':{'selected.tsv':'selected.tsv','selected.fasta':'selected.fasta',
                 'compact-pool.tsv':'candidate-pool.tsv',
                 'selection-audit.tsv':'selection-audit.tsv',
                 'hydrophobicity-screen.tsv':'hydrophobicity-screen.tsv',
                 'hydrophobicity-summary.json':'hydrophobicity-summary.json',
                 'splice-frame-summary.json':'splice-frame-summary.json',
                 'splice-gene-ranking.tsv':'splice-gene-ranking.tsv',
                 'splice-event-ranking.tsv':'splice-event-ranking.tsv',
                 'splice-rna-verification.json':'splice-rna-verification.json',
                 'selection-verification.json':'selection-verification.json',
                 'mutation-score-verification.json':'mutation-score-verification.json'},
    'construct':{'joint.fasta':'cassette.fasta','joint.json':'construct.json',
                 'joint-components.tsv':'components.tsv','joint-junctions.tsv':'junctions.tsv',
                 'joint-boundary-predictions.tsv':'boundary-predictions.tsv',
                 'hydrophobicity-windows.tsv':'cassette-hydrophobicity-windows.tsv',
                 'hydrophobicity.json':'cassette-hydrophobicity.json','verification.json':'construct-verification.json'}}
STABLE=('mutation-score-verification.json','selected.tsv','selected.fasta','candidate-pool.tsv','selection-audit.tsv',
        'hydrophobicity-screen.tsv','hydrophobicity-summary.json','cassette-hydrophobicity-windows.tsv','cassette-hydrophobicity.json',
        'splice-frame-summary.json','splice-gene-ranking.tsv','splice-event-ranking.tsv','splice-rna-verification.json',
        'cassette.fasta','components.tsv','junctions.tsv','boundary-predictions.tsv')
COUNT_KEYS=('required_peptides','required_peptide_hla_pairs','cache_peptides','cache_peptide_hla_pairs')

def resolve_recipe(path):
    path=Path(path).resolve()
    if path.is_dir() and (path/'selection/design.json').is_file():
        return path/'selection/design.json'
    return path/'design.json' if path.is_dir() else path

def referenced(recipe_path,value):return (Path(recipe_path).parent/value).resolve()

def load(repo,path,*,check_cache=True):
    repo=Path(repo).resolve();path=resolve_recipe(path)
    if not path.is_file():
        raise ValueError('No optional recipe at '+str(path)+'. The lead uses yash-vaccine reproduce designs/design-2. Start an independent design with yash-vaccine init NAME, or copy the Design #1 recipe with --from designs/initial/selection.')
    recipe=json.loads(path.read_text())
    expected={'schema_version','name','description','evidence','selection','linkers','prediction_cache','optimization'}
    if not isinstance(recipe,dict) or set(recipe)!=expected or recipe['schema_version']!=1:
        raise ValueError('Design recipe must use schema_version 1 and the fields in designs/initial/selection/design.json')
    if not isinstance(recipe['name'],str) or not re.fullmatch(r'[a-z0-9][a-z0-9_-]*',recipe['name']):
        raise ValueError('Design name must contain lowercase letters, numbers, hyphens or underscores')
    if not isinstance(recipe['description'],str):raise ValueError('Design description must be text')
    if not isinstance(recipe['selection'],dict) or set(recipe['selection'])!=SELECTION_KEYS:
        raise ValueError('Selection fields must match the documented selection schema')
    validate_selection(recipe['selection'])
    evidence=recipe['evidence']
    if not isinstance(evidence,dict) or set(evidence)!={'manifest','sha256','context','context_sha256'}:
        raise ValueError('Evidence must specify manifest, sha256, context and context_sha256')
    if any(not isinstance(value,str) for value in evidence.values()):raise ValueError('Evidence paths and checksums must be text')
    # The current scientific adapters consume the recorded case bundle. Reject a
    # different bundle instead of appearing to use files the adapters do not read.
    for key,relative,digest in [('manifest','data/checkpoints/manifest.json','sha256'),('context','config/case.json','context_sha256')]:
        target=referenced(path,evidence[key])
        if target!=(repo/relative).resolve():
            raise ValueError('Our optional recipe runner uses the recorded case bundle at '+relative+'. Use an independent design script for different evidence or methods; see docs/usage.md.')
        if sha256(target)!=evidence[digest]:raise ValueError('Evidence checksum mismatch: '+str(target))
    cache=recipe['prediction_cache']
    if not isinstance(cache,dict) or set(cache)!={'path','sha256'}:raise ValueError('Prediction cache needs path and sha256')
    if any(not isinstance(value,str) for value in cache.values()):raise ValueError('Cache path and checksum must be text')
    if check_cache and sha256(referenced(path,cache['path']))!=cache['sha256']:raise ValueError('Prediction-cache checksum mismatch')
    optimization=recipe['optimization']
    if not isinstance(optimization,dict) or set(optimization)!={'strategy','seconds_per_phase'}:
        raise ValueError('Optimization needs strategy and seconds_per_phase')
    if optimization['strategy']!='junction_binding_lexicographic_v1':raise ValueError('Unsupported optimization strategy')
    seconds=optimization['seconds_per_phase']
    if type(seconds) is not int or seconds<=0:raise ValueError('seconds_per_phase must be a positive integer')
    return path,recipe

def rebased(recipe,source,destination):
    result=copy.deepcopy(recipe)
    for section,key in [('evidence','manifest'),('evidence','context'),('prediction_cache','path')]:
        result[section][key]=os.path.relpath(referenced(source,recipe[section][key]),Path(destination).parent)
    return result

def create(repo,name,source=None):
    repo=Path(repo).resolve()
    if not re.fullmatch(r'[a-z0-9][a-z0-9_-]*',name):raise ValueError('Use a lowercase design name with letters, numbers, hyphens or underscores')
    out=repo/'designs'/name
    if out.exists():raise ValueError('Design already exists: '+str(out))
    if source is None:
        out.mkdir(parents=True)
        (out/'README.md').write_text('# '+name+'\n\n'
            'Describe your question, approach and findings here. Use any method or tools you need.\n\n'
            '## Starting point\n\n'
            'The lead is [Design #2](../design-2/README.md), generated with matched-normal DNA.\n'
            'Reuse its [candidate pool](../design-2/exports/integrated-candidate-pool.tsv),\n'
            '[selected regions](../design-2/selection/selected.tsv) and [methods](../design-2/methods.md)\n'
            'where useful. This workspace does not copy either design\'s analysis.\n\n'
            '## Inputs\n\nList the data, software versions and any additional downloads needed.\n'
            'Shared data and access instructions are in ../../data/README.md.\n\n'
            '## Run\n\nFrom the repository root:\n\n'
            '```sh\npython designs/'+name+'/run.py\n```\n\n'
            'Replace the starter script with your analysis, or document another way to run it.\n'
            'The starter does not perform an analysis or create a candidate.\n\n'
            '## Results\n\nKeep useful outputs here: for example selected.tsv and, if assembled, cassette.fasta.\n'
            'Use filenames and formats that suit your approach; no design.json or runner schema is required.\n'
            'Explain changes and evidence so another researcher can reproduce your work.\n')
        (out/'run.py').write_text('"""Starter for '+name+'. Replace main() with your analysis.\n\n'
            'Run from the repository root with: python designs/'+name+'/run.py\n'
            'Use any inputs, methods and output formats appropriate to your question.\n"""\n'
            'from pathlib import Path\n\n'
            'HERE = Path(__file__).resolve().parent\n'
            'REPO = HERE.parents[1]\n'
            'DATA = REPO / "data"\n\n'
            'def main():\n'
            '    print("Workspace: " + str(HERE))\n'
            '    print("Shared evidence: " + str(DATA))\n'
            '    print("Add your analysis to run.py. No analysis has been run and no candidate has been generated.")\n\n'
            'if __name__ == "__main__":\n'
            '    main()\n')
        return out
    source,recipe=load(repo,source,check_cache=False)
    recipe=rebased(recipe,source,out/'design.json');recipe['name']=name
    recipe['description']='Alternative candidate: '+name+'.'
    out.mkdir(parents=True);write_json(out/'design.json',recipe)
    (out/'README.md').write_text('# '+name+'\n\n'
        'This candidate copies the optional Design #1 tumour-only recipe. It varies selection\n'
        'and the internal cassette, before SEC/MITD and RNA encoding. It does not use the\n'
        'matched-normal evidence or full construct workflow of lead Design #2.\n'
        'For an alternative to the lead, use a plain workspace from `yash-vaccine init NAME`\n'
        'and [Design #2\'s methods and evidence](../design-2/README.md) in your own script.\n\n'
        'Describe your question and changes here. You can also replace this approach with your own script.\n\n'
        'Edit the selection or linker choices in design.json, then run from the repository root:\n\n'
        '```sh\nyash-vaccine run designs/'+name+'/design.json\n'
        'yash-vaccine compare designs/initial designs/'+name+'\n```\n\n'
        'Commit this design directory on your branch when it is ready for review.\n')
    return out

def source_record(repo):
    def git(*args):
        try:return subprocess.check_output(['git','-C',str(repo),*args],text=True,stderr=subprocess.DEVNULL).strip()
        except (OSError,subprocess.CalledProcessError):return None
    versions={}
    for name in ('biopython','numpy','ortools'):
        try:versions[name]=importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:versions[name]='unavailable'
    dirty=bool(git('status','--porcelain','--untracked-files=no'))
    return dict(commit=None if dirty else git('rev-parse','HEAD'),tracked_changes=dirty,
                revision_note='Content hashes identify the executed code; a dirty worktree is not attributed to its parent commit.',
                python=sys.version.split()[0],packages=versions,
                source_sha256={str(p.relative_to(repo)):sha256(p) for p in sorted((repo/'src/yash_vaccine').rglob('*.py'))})

def record(repo,out):
    """Attach readable and machine-readable evidence to an already built design."""
    out=Path(out);recipe=json.loads((out/'design.json').read_text())
    regions=list(rows(out/'selected.tsv'));selection=json.loads((out/'selection-verification.json').read_text())
    construct=json.loads((out/'construct.json').read_text());verified=json.loads((out/'construct-verification.json').read_text())
    if selection['status']!='VERIFIED' or verified['status']!='VERIFIED':raise ValueError('Incomplete selection or construct verification')
    if not verified['joint_full_chain_verified'] or not construct['solver']['all_phases_proven_optimal']:
        raise ValueError('Design has not completed full-chain verification and all optimization phases')
    sources=dict(collections.Counter(row['source'] for row in regions))
    report=dict(status='VERIFIED',scope='Local selection and joint construction from recorded sequence/score evidence; no new model inference or raw-read analysis',
                name=recipe['name'],selected_regions=len(regions),protein_length_aa=construct['protein_length_aa'],sources=sources,
                recipe_sha256=sha256(out/'design.json'),evidence_manifest_sha256=recipe['evidence']['sha256'],
                predictions={key:verified[key] for key in COUNT_KEYS},objective_tuple=construct['objective_tuple'],
                all_solver_phases_optimal=True,full_chain_verified=True,
                artifact_sha256={name:sha256(out/name) for name in STABLE})
    report['output_sha256']={name:sha256(out/name) for mapping in EXPORTS.values() for name in mapping.values()}
    write_json(out/'verification.json',report)
    write_json(out/'run.json',dict(code=source_record(Path(repo)),recipe_sha256=report['recipe_sha256'],
                                 evidence_manifest_sha256=recipe['evidence']['sha256'],model_inference_rerun=False))
    return report

def run(repo,path,out=None,prediction_cache=None):
    from .replay import execute
    from .constructs.run import run as construct
    repo=Path(repo).resolve();path,recipe=load(repo,path,check_cache=prediction_cache is None)
    out=Path(out).resolve() if out is not None else path.parent
    if out==repo:raise ValueError('Choose a design directory or a new work/ output directory')
    generated=[name for mapping in EXPORTS.values() for name in mapping.values()]+['verification.json','run.json']
    if out.exists() and (out!=path.parent or any((out/name).exists() for name in generated)):
        raise ValueError('Results already exist at '+str(out)+'. Use reproduce to check them or choose a new --output directory.')
    cache=Path(prediction_cache).resolve() if prediction_cache else referenced(path,recipe['prediction_cache']['path'])
    with tempfile.TemporaryDirectory(prefix='yash-design-') as tmp:
        temp=Path(tmp);selection=temp/'selection';assembly=temp/'construct';linkers=temp/'linkers.json'
        write_json(linkers,recipe['linkers'])
        execute(repo,selection,recipe['selection'])
        construct(repo,assembly,recipe['optimization']['seconds_per_phase'],regions_path=selection/'selected.tsv',cache_path=cache,linkers_path=linkers)
        # Export only after the cache, optimization and full chain have passed.
        out.mkdir(parents=True,exist_ok=True)
        exported=rebased(recipe,path,out/'design.json')
        if prediction_cache:exported['prediction_cache']=dict(path=os.path.relpath(cache,out),sha256=sha256(cache))
        if out!=path.parent or prediction_cache:write_json(out/'design.json',exported)
        for directory,mapping in EXPORTS.items():
            for old,new in mapping.items():shutil.copy2(temp/directory/old,out/new)
    report=record(repo,out)
    if not (out/'README.md').exists():
        (out/'README.md').write_text('# '+recipe['name']+'\n\n'+recipe['description']+'\n\n'
            f'This design contains {report["selected_regions"]} regions and a {report["protein_length_aa"]}-aa protein cassette.\n\n'
            'Read design.json for the recipe, selected.tsv for the antigen selection, cassette.fasta for the protein,\n'
            'and verification.json for checks and file hashes. The tables are tab-separated text.\n')
    return report

def reproduce(repo,path=None,out=None):
    repo=Path(repo).resolve()
    if path is None or Path(path).resolve()==repo/"designs/design-2":
        from .lead import reproduce as reproduce_lead
        return reproduce_lead(repo,out)
    if Path(path).resolve()==repo/"designs/initial":
        from .initial import reproduce as reproduce_complete
        return reproduce_complete(repo,out)
    path,recipe=load(repo,path)
    reference=path.parent;expected=json.loads((reference/'verification.json').read_text())
    if sha256(path)!=expected['recipe_sha256']:raise ValueError('The recorded recipe has changed. Create a named alternative or restore the recipe before reproduction.')
    for name,digest in expected.get('output_sha256',expected['artifact_sha256']).items():
        if sha256(reference/name)!=digest:raise ValueError('Recorded artifact checksum mismatch: '+name)
    def check(destination):
        actual=run(repo,path,destination)
        changed=[name for name,digest in expected['artifact_sha256'].items() if actual['artifact_sha256'].get(name)!=digest]
        if changed:raise ValueError('Reproduced artifacts differ: '+', '.join(changed))
        if actual['objective_tuple']!=expected['objective_tuple']:raise ValueError('Reproduced optimization objective differs')
        return dict(status='VERIFIED',name=recipe['name'],selected_regions=actual['selected_regions'],
                    protein_length_aa=actual['protein_length_aa'],all_recorded_artifacts_identical=True,
                    all_solver_phases_optimal=True,output=str(destination) if out else None)
    if out is not None:return check(Path(out).resolve())
    with tempfile.TemporaryDirectory(prefix='yash-reproduce-') as tmp:return check(Path(tmp)/'design')

def prepare(repo,path,out):
    from .replay import execute
    from .constructs.inputs import request
    repo=Path(repo).resolve();path,recipe=load(repo,path,check_cache=False);out=Path(out).resolve()
    if out.exists():raise ValueError('Choose a new prediction-request directory: '+str(out))
    with tempfile.TemporaryDirectory(prefix='yash-prepare-') as tmp:
        temp=Path(tmp);selection=temp/'selection';linkers=temp/'linkers.json';write_json(linkers,recipe['linkers'])
        execute(repo,selection,recipe['selection'])
        peptides,metadata=request(repo,selection/'selected.tsv',linkers)
        out.mkdir(parents=True)
        shutil.copy2(selection/'selected.tsv',out/'selected.tsv');shutil.copy2(linkers,out/'linkers.json')
        metadata.update(regions=str(out/'selected.tsv'),linkers=str(out/'linkers.json'))
        (out/'peptides.txt').write_text(''.join(p+'\n' for p in peptides))
        write_json(out/'alleles.json',metadata['alleles'])
        metadata.update(peptides_sha256=sha256(out/'peptides.txt'),alleles_sha256=sha256(out/'alleles.json'))
        write_json(out/'prediction-request.json',metadata)
        write_json(out/'design.json',rebased(recipe,path,out/'design.json'))
    return metadata
