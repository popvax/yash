"""Score CTA/self and fusion contexts using staged LENS data and mutation regions.

Execution requires --execute and a fresh output directory. Patient-RNA splice
reconstruction has its own stage under patient_splice/.
"""
import argparse
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
SOURCE='/data/results/TN25-112103/lens-expression-aware-20261030'

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--output-root',required=True)
    p.add_argument('--execute',action='store_true')
    p.add_argument('--mutation-regions',type=Path,required=True,help='Segmented mutation-region TSV, such as data/checkpoints/mutation-regions.tsv')
    p.add_argument('--lens-root',default=SOURCE,help='Completed RNA-only LENS output on the volume')
    p.add_argument('--image',default='im-01M66YBAT7GF00KRAD4BJMX01J',help='Saved runtime with vaxrank 3.46.0; rebuild from environment.freeze.txt if snapshot has expired')
    args=p.parse_args();root=args.output_root;source=args.lens_root
    if not args.mutation_regions.is_file():p.error('--mutation-regions must be an existing TSV')
    assert root.startswith('/data/results/') and root not in [SOURCE,'/data/results/TN25-112103/combined-antigen-regions-20261030']
    print(json.dumps(dict(output_root=root,lens_root=source,mutation_regions=str(args.mutation_regions),sources=['CTA/SELF','FUSION'],execute=args.execute),indent=2))
    if not args.execute:return
    import modal
    volume=modal.Volume.from_name('research-lens-workspace')
    sandbox=modal.Sandbox.create(app=modal.App.lookup('research-lens',create_if_missing=True),image=modal.Image.from_id(args.image),
                                 cpu=4,memory=32768,timeout=14400,runtime='vm',volumes={'/data':volume})
    print('Sandbox:',sandbox.object_id,flush=True)
    def run(*cmd,timeout=3600):
        proc=sandbox.exec(*cmd,timeout=timeout)
        print(proc.stdout.read(),flush=True);print(proc.stderr.read(),flush=True)
        assert proc.wait()==0,cmd
    try:
        run('python3','-m','pip','install','mhcflurry==2.3.14',timeout=300)
        run('python3','-c','from pathlib import Path; p=Path('+repr(root)+'); assert not p.exists(), "Use a new directory"; p.mkdir(parents=True)')
        for name in ['prepare_contexts.py','score_regions.py','verify_regions.py','expand_scoring_targets.py','prepare_scoring_evidence.py','audit_reference_native.py','export_native_reference.py','verify_prediction_runtime.py','consolidate_regions.py']:
            sandbox.filesystem.write_bytes((HERE/name).read_bytes(),root+'/'+name)
        mapping={'final/all-reviewed-candidates.tsv':'lens-all-reviewed.tsv',
                 'prediction/source-metadata.jsonl':'prediction/source-metadata.jsonl',
                 'inputs-checkpoints/original-starfusion-coding-effect.tsv':'inputs-checkpoints/original-starfusion-coding-effect.tsv'}
        for src,dest in mapping.items():
            run('python3','-c','from pathlib import Path; import shutil; p=Path('+repr(root+'/inputs/'+dest)+');p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile('+repr(source+'/'+src)+',p)')
        sandbox.filesystem.write_bytes(args.mutation_regions.read_bytes(),root+'/inputs/mutation-regions.tsv')
        run('python3',root+'/prepare_contexts.py')
        run('python3',root+'/expand_scoring_targets.py',root)
        run('python3',root+'/audit_reference_native.py',root)
        run('python3',root+'/export_native_reference.py',root)
        run('python3',root+'/prepare_scoring_evidence.py',root)
        run('python3',root+'/score_regions.py',root,timeout=7200)
        run('python3',root+'/verify_prediction_runtime.py',root)
        run('python3',root+'/verify_regions.py',root)
        print('Verified output on research-lens-workspace:',root,flush=True)
    finally:
        sandbox.terminate()

if __name__=='__main__':main()
