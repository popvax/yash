"""Explore shared evidence, start independent designs or reuse our workflow."""
import argparse,json,subprocess,sys
from pathlib import Path

def _main():
    p=argparse.ArgumentParser(prog='yash-vaccine',description=__doc__,epilog='Start with README.md and docs/usage.md. Commands use local evidence unless you explicitly request downloads or remote execution.')
    p.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[2])
    sub=p.add_subparsers(dest='command',required=True)
    reproduce=sub.add_parser('reproduce',help='Rebuild a recorded design and compare its files; local only')
    reproduce.add_argument('recipe',type=Path,nargs='?',help='Design directory or design.json; defaults to lead designs/design-2')
    reproduce.add_argument('--output',type=Path,help='Optional new directory in which to keep reproduced files')
    verify=sub.add_parser('verify',help='Check recorded lead results without rebuilding or starting jobs')
    verify.add_argument('design',type=Path,nargs='?',help='Defaults to designs/design-2')
    init=sub.add_parser('init',help='Start an independent design; --from copies a Design #1 recipe',
        description='Create a plain README/script workspace for a new design, including alternatives to lead Design #2. Optional --from copies a Design #1 tumour-only selection/internal-cassette recipe; it does not copy the matched-normal lead.')
    init.add_argument('name');init.add_argument('--from',dest='source',type=Path,help='Copy a compatible Design #1 tumour-only recipe, not lead Design #2')
    run=sub.add_parser('run',help='Run the optional Design #1 selection/internal-cassette recipe',
        description='Run a Design #1-derived tumour-only selection/internal-cassette recipe. This is separate from matched-normal Design #2; use reproduce to replay the fixed lead or init NAME to develop an alternative.')
    run.add_argument('recipe',type=Path);run.add_argument('--output',type=Path)
    run.add_argument('--prediction-cache',type=Path,help='Use new predictions and record their path/hash in the exported recipe')
    compare=sub.add_parser('compare',help='Compare complete recorded designs or compatible recipe-stage outputs')
    compare.add_argument('left',type=Path);compare.add_argument('right',type=Path)
    compare.add_argument('--output',type=Path,help='Save comparison.txt, regions.tsv and comparison.json in a new directory')
    prepare=sub.add_parser('prepare',help='Prepare prediction inputs for a Design #1-derived recipe; no inference')
    prepare.add_argument('recipe',type=Path);prepare.add_argument('--output',type=Path,required=True)
    predict=sub.add_parser('predict',help='Plan predictions for a prepared request; --execute starts an A10 job')
    predict.add_argument('request',type=Path);predict.add_argument('--output',type=Path,required=True)
    predict.add_argument('--execute',action='store_true')
    data=sub.add_parser('data',help='List source data or plan/execute verified downloads')
    data.add_argument('action',choices=['list','restore'],nargs='?',default='list')
    data.add_argument('--category',choices=['raw','references','analysis','software','workflow','all'])
    data.add_argument('--destination',type=Path);data.add_argument('--key',action='append')
    data.add_argument('--execute',action='store_true');data.add_argument('--anonymous',action='store_true')
    data.add_argument('--credentials-file',type=Path);data.add_argument('--public-base-url')
    upstream=sub.add_parser('upstream',help='Inspect Design #1 DNA/RNA and shared RNA-source adapters; lead controllers are documented separately')
    upstream.add_argument('stage',choices=['dna','rna','lens','regions'])
    upstream.add_argument('--execute',action='store_true')
    args,extra=p.parse_known_args();repo=args.repo.resolve()
    if extra and args.command!='upstream':p.error('unrecognized arguments: '+' '.join(extra))
    if args.command=='verify':
        from .lead import verify
        design=(args.design or repo/'designs/design-2').resolve()
        if design!=repo/'designs/design-2':p.error('verify currently checks the lead designs/design-2; use reproduce for other recorded designs')
        result=verify(repo)
        print(f'Verified lead Design #2: {result["selected_regions"]} regions, {result["protein_aa"]} aa. Recorded evidence, sequences and annotations agree. No rebuild or cloud job started.')
        external=result.get('source_snapshot_verification',{}).get('external_source_hash_records_only',[])
        if external:print(f'{len(external)} upstream source files are referenced by recorded hash, not bundled or byte-verified locally; see docs/external-resources.md.')
    elif args.command=='reproduce':
        from .designs import reproduce
        result=reproduce(repo,args.recipe,args.output)
        print(f'Verified {result["name"]}: {result["selected_regions"]} regions, {result["protein_length_aa"]} aa. All recorded artifacts match; all solver phases are optimal.')
        external=result.get('source_snapshot_verification',{}).get('external_source_hash_records_only',[])
        if external:print(f'{len(external)} upstream source files remain external hash references; the replay does not execute them.')
        if "rna_length_nt" in result:
            if result["rna_length_nt"] is None:
                print('PolyA is unspecified. Known RNA/DNA segments match; complete RNA/template lengths are undefined.')
            else:
                print(f'RNA: {result["rna_length_nt"]:,} nt; DNA template context: {result["template_length_bp"]:,} bp. Junction checks and nucleotide build match.')
        if result['output']:print('Files:',result['output'])
    elif args.command=='init':
        from .designs import create
        out=create(repo,args.name,args.source)
        print('Created',out)
        if args.source:
            print('Design #1 tumour-only recipe: selection/internal cassette only; matched-normal Design #2 is not copied.')
            print('Edit design.json and README.md, then run:')
            print('  yash-vaccine run '+str(out/'design.json'))
        else:
            print('Describe your approach in README.md and implement it in run.py (or your own tools).')
            print('  python '+str(out/'run.py'))
            print('No recipe or prescribed pipeline is required.')
    elif args.command=='run':
        from .designs import run,resolve_recipe
        result=run(repo,args.recipe,args.output,args.prediction_cache)
        print(f'Verified {result["name"]}: {result["selected_regions"]} regions, {result["protein_length_aa"]} aa.')
        print('Files:',(args.output or resolve_recipe(args.recipe).parent).resolve())
        print('Selection: selected.tsv; protein: cassette.fasta; checks: verification.json.')
    elif args.command=='compare':
        from .compare import compare
        print(compare(args.left,args.right,args.output)['text'],end='')
    elif args.command=='prepare':
        from .designs import prepare
        result=prepare(repo,args.recipe,args.output)
        print(f'Prepared {result["unique_peptides"]:,} peptides × {len(result["alleles"])} HLA alleles = {result["peptide_hla_pairs"]:,} predictions. No inference run.')
        print('Request:',args.output.resolve())
    elif args.command=='predict':
        from .io import sha256
        request=args.request.resolve();metadata=json.loads((request/'prediction-request.json').read_text())
        for name,key in [('selected.tsv','regions_sha256'),('linkers.json','linkers_sha256'),('peptides.txt','peptides_sha256'),('alleles.json','alleles_sha256')]:
            if sha256(request/name)!=metadata[key]:raise ValueError('Prepared prediction input changed: '+name+'. Prepare a new request.')
        command=[sys.executable,str(repo/'tools/launch_predictions.py'),'--regions',str(request/'selected.tsv'),
                 '--linkers',str(request/'linkers.json'),'--output',str(args.output.resolve())]
        if args.execute:raise SystemExit(subprocess.call(command+['--execute']))
        result=subprocess.run(command,capture_output=True,text=True)
        if result.returncode:raise ValueError(result.stderr.strip() or result.stdout.strip())
        plan=json.loads(result.stdout)
        print(f'Prediction plan: {plan["unique_peptides"]:,} peptides, {plan["peptide_hla_pairs"]:,} peptide/HLA pairs on {plan["gpu"]}.')
        print('Output:',plan['local_output']);print('No job started. Add --execute to run.')
    elif args.command=='data':
        if args.action=='list':
            from .io import sha256
            access=json.loads((repo/'data/catalog/access.json').read_text());manifest=repo/'data/catalog'/access['bundled_manifest']
            if sha256(manifest)!=access['bundled_manifest_sha256']:raise ValueError('Bundled catalog checksum mismatch')
            category=args.category or 'raw'
            records=[x for x in json.loads(manifest.read_text())['objects'] if category=='all' or x['category']==category]
            print(f'{len(records)} catalog objects ({sum(x["bytes"] for x in records):,} bytes):')
            for item in records:print(f'  {item["bytes"]:>14,} bytes  {item["key"]}')
            if category in ('references','software','workflow'):
                print('Third-party code, models and references are upstream-only; see docs/external-resources.md and data/catalog/external-resources.json.')
            print('External files use the access instructions in data/README.md; listing needs no network.')
        else:
            if args.destination is None:p.error('data restore requires --destination')
            command=[sys.executable,str(repo/'tools/r2/restore.py'),'--destination',str(args.destination)]
            for flag,value in [('--category',args.category),('--credentials-file',args.credentials_file),('--public-base-url',args.public_base_url)]:
                if value is not None:command.extend([flag,str(value)])
            for key in args.key or []:command.extend(['--key',key])
            if args.execute:command.append('--execute')
            if args.anonymous:command.append('--anonymous')
            raise SystemExit(subprocess.call(command))
    elif args.command=='upstream':
        stages=repo/'workflow/stages'
        if args.stage=='dna':
            command=[sys.executable,str(stages/'01_mutation_calling/discovery/calling/run.py'),'calling',*extra]
        elif args.stage=='rna':
            command=[sys.executable,str(stages/'01_mutation_calling/discovery/review/run.py'),'rna',*extra]
        else:
            stage='02_lens' if args.stage=='lens' else '04_region_integration'
            command=[sys.executable,str(stages/stage/'reproduce.py'),*extra]
        if args.execute:
            command.append('--execute')
            raise SystemExit(subprocess.call(command))
        import shlex
        print('Design #1/shared-source adapter; paired Design #2 controllers: designs/design-2/pipeline/README.md')
        print('Plan only. Recorded stage:',args.stage);print('Native command:',shlex.join(command))
        print('Read docs/upstream.md for required data, prepared calling-plan and runtime. Add --execute to launch; no cloud job started.')

def main():
    try:_main()
    except (ValueError,FileNotFoundError,FileExistsError) as error:raise SystemExit(str(error))

if __name__=='__main__':main()
