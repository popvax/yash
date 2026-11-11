"""Plan or run construct predictions for explicit regions/linkers on Modal.

No compute starts without --execute. Outputs go to a new directory, and the
frozen checkpoint cache/manifest are never modified by this launcher.
"""
import argparse,datetime,json,sys
from pathlib import Path


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--regions',type=Path,help='Defaults to designs/initial/selection/selected.tsv')
    p.add_argument('--linkers',type=Path,help='Defaults to config/linkers.json')
    p.add_argument('--output',type=Path,help='New local directory for NPZ and inference provenance')
    p.add_argument('--execute',action='store_true')
    p.add_argument('--image',default='im-DBMghmAQfB3Oz7McdNZeic')
    p.add_argument('--volume',default='research-lens-workspace')
    p.add_argument('--models',default='/data/results/TN25-112103/openvax-mutations-20261018/references/mhcflurry/2.3.0/models_class1_presentation/models')
    args=p.parse_args();root=Path(__file__).resolve().parents[1]
    sys.path.insert(0,str(root/'src'))
    from yash_vaccine.constructs.inputs import request
    from yash_vaccine.io import write_json
    region_file=(args.regions or root/'designs/initial/selection/selected.tsv').resolve()
    linker_file=(args.linkers or root/'config/linkers.json').resolve()
    _,queries=request(root,region_file,linker_file)
    stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    output=(args.output or root/'work'/('predictions-'+stamp)).resolve();remote='/tmp/yash-vaccine-design'
    state=dict(**queries,local_output=str(output),image=args.image,volume=args.volume,models=args.models,
               gpu='A10',cpu=4,memory_mib=16384,execute=args.execute,terminated=False)
    print(json.dumps(state,indent=2),flush=True)
    if not args.execute:return
    if output.exists():p.error('Choose a new --output directory')
    import modal
    output.mkdir(parents=True);state_path=output/'modal-state.json';sandbox=None
    try:
        sandbox=modal.Sandbox.create(app=modal.App.lookup('research-yash-vaccine',create_if_missing=True),
            image=modal.Image.from_id(args.image),gpu='A10',cpu=4,memory=16384,timeout=2400,
            volumes={'/data':modal.Volume.from_name(args.volume)})
        state.update(sandbox_id=sandbox.object_id,started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
        write_json(state_path,state)
        paths=[root/'tools/predict_constructs.py',root/'config/case.json',root/'config/frozen-replay.json',
               root/'data/checkpoints/provenance/construct-model-files.json',
               *root.glob('src/yash_vaccine/*.py'),*root.glob('src/yash_vaccine/constructs/*.py')]
        for local in paths:sandbox.filesystem.write_bytes(local.read_bytes(),remote+'/'+str(local.relative_to(root)))
        sandbox.filesystem.write_bytes(region_file.read_bytes(),remote+'/inputs/regions.tsv')
        sandbox.filesystem.write_bytes(linker_file.read_bytes(),remote+'/inputs/linkers.json')
        proc=sandbox.exec('python',remote+'/tools/predict_constructs.py','--repo',remote,'--models',args.models,
                         '--regions',remote+'/inputs/regions.tsv','--linkers',remote+'/inputs/linkers.json',
                         '--output',remote+'/prediction-output',timeout=1800)
        for line in proc.stdout:print(line,end='',flush=True)
        (output/'stderr.log').write_text(proc.stderr.read());state['exit_code']=proc.wait()
        if state['exit_code']!=0:raise RuntimeError('Prediction failed; inspect '+str(output/'stderr.log'))
        for name in ['construct-predictions.npz','construct-model-files.json','construct-inference.json']:
            (output/name).write_bytes(sandbox.filesystem.read_bytes(remote+'/prediction-output/'+name))
        state['completed']=True
    finally:
        if sandbox:sandbox.terminate(wait=True);state['terminated']=True
        state['finished_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();write_json(state_path,state)


if __name__=='__main__':main()
