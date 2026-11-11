"""Wait for selection, then explicitly rebuild and verify the routed construct."""
from pathlib import Path
import argparse
import hashlib
import json
import time

CODE = Path(__file__).resolve().parent
HERE = CODE.parent;REPO=HERE.parents[1]
REMOTE='/data/results/TN25-112103/matched-normal-20261107'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--execute',action='store_true')
    args=parser.parse_args();print(json.dumps(dict(execute=args.execute,stages=['prepare','predict','solve','predict-selected','verify','encode'],
        objective='Same recorded joint junction-affinity and hydropathy recipe; fresh predictions for all new inputs')),flush=True)
    if not args.execute:return
    import modal
    local=REPO/'work/matched-normal-20261107';volume=modal.Volume.from_name('research-lens-workspace')
    while True:
        sf=local/'integrate-state.json'
        if sf.exists():
            state=json.loads(sf.read_text());sandbox=modal.Sandbox.from_id(state['sandbox_id']);rc=sandbox.poll()
            if rc==0:break
            if rc is not None:raise RuntimeError('Integrated selection is terminal: '+str(rc))
        else:
            state=json.loads((local/'alignment-state.json').read_text());sandbox=modal.Sandbox.from_id(state['sandbox_id']);rc=sandbox.poll()
            if rc not in [None,0]:raise RuntimeError('Alignment is terminal: '+str(rc))
        print(json.dumps(dict(waiting_for='integrated matched-normal selection',observed_sandbox=sandbox.object_id,observed_live=rc is None)),flush=True)
        time.sleep(30)
    worker=CODE / 'cassette.py';common=CODE / 'common.py';digest=hashlib.sha256(worker.read_bytes()+common.read_bytes()).hexdigest()
    code=REMOTE+'/code/cassette-'+digest[:12]
    for path in [worker,common]:
        key=(code+'/'+path.name).removeprefix('/data/')
        try:
            with volume.batch_upload() as upload:upload.put_file(str(path),key)
        except FileExistsError:assert b''.join(volume.read_file(key))==path.read_bytes()
    state_file=local/'construct-state.json'
    state=json.loads(state_file.read_text()) if state_file.exists() else dict(bundle_sha256=digest,workers={})
    assert state['bundle_sha256']==digest
    def save_state():
        temporary=state_file.with_suffix('.tmp');temporary.write_text(json.dumps(state,indent=2)+'\n');temporary.replace(state_file)
    solver_image=modal.Image.debian_slim(python_version='3.11').pip_install('numpy==1.26.4','ortools==9.14.6206')
    rna_image=modal.Image.debian_slim(python_version='3.11').pip_install('ViennaRNA==2.7.2','python-codon-tables==0.1.18','biopython==1.86','numpy==1.26.4')
    for phase in ['prepare','predict','solve','predict-selected','verify','encode']:
        folder={'prepare':'prepare','predict':'all-inference','solve':'protein','predict-selected':'selected-inference','verify':'verify','encode':'mrna'}[phase]
        old=state['workers'].get(phase)
        if old:
            sandbox=modal.Sandbox.from_id(old['sandbox_id'])
            if sandbox.poll() not in [None,0]:raise RuntimeError(phase+' is terminal; inspect recorded logs before retry')
            target=local/'construct'/folder
            cached_receipt=target/'complete.json'
            if sandbox.poll()==0 and cached_receipt.exists():
                cached=json.loads(cached_receipt.read_text())
                if cached.get('status')=='COMPLETED' and all(
                        metadata['bytes']>40000000 or ((target/name).exists()
                            and (target/name).stat().st_size==metadata['bytes']
                            and hashlib.sha256((target/name).read_bytes()).hexdigest()==metadata['sha256'])
                        for name,metadata in cached['outputs'].items()):
                    print('DONE_REUSED_VERIFIED',phase,sandbox.object_id,flush=True)
                    continue
        else:
            gpu_phase=phase in ['predict','predict-selected']
            image=modal.Image.from_id('im-DBMghmAQfB3Oz7McdNZeic') if gpu_phase else solver_image if phase=='solve' else rna_image if phase=='encode' else modal.Image.from_id('im-hX1h6HknDMvasm9Z6q530H')
            with modal.enable_output():
                sandbox=modal.Sandbox.create('python','-u',code+'/cassette.py',REMOTE,phase,
                    app=modal.App.lookup('research-yash-matched-normal',create_if_missing=True),image=image,
                    gpu='A10' if gpu_phase else None,cpu=4 if gpu_phase else 2,memory=32768 if gpu_phase else 16384,
                    timeout=21600,volumes={'/data':volume})
            state['workers'][phase]=dict(sandbox_id=sandbox.object_id,image_id=image.object_id,status='running',gpu='A10' if gpu_phase else None);save_state()
            print('START',phase,sandbox.object_id,flush=True)
        stdout=sandbox.stdout.read();stderr=sandbox.stderr.read();sandbox.wait();rc=sandbox.poll()
        (local/('construct-'+phase+'.stdout.log')).write_text(stdout);(local/('construct-'+phase+'.stderr.log')).write_text(stderr)
        state['workers'][phase].update(status='completed' if rc==0 else 'failed',exit_code=rc);save_state()
        if rc:raise RuntimeError(phase+' failed: '+stderr[-4000:])
        remote=REMOTE.removeprefix('/data/')+'/construct/'+folder
        body=b''.join(volume.read_file(remote+'/complete.json'));receipt=json.loads(body);assert receipt['status']=='COMPLETED'
        target=local/'construct'/folder;target.mkdir(parents=True,exist_ok=True);(target/'complete.json').write_bytes(body)
        for name,meta in receipt['outputs'].items():
            if meta['bytes']>40000000:continue
            content=b''.join(volume.read_file(remote+'/'+name));assert hashlib.sha256(content).hexdigest()==meta['sha256']
            (target/name).write_bytes(content)
        print('DONE',phase,flush=True)
    state['status']='COMPLETED';save_state()
    print('MATCHED_NORMAL_CONSTRUCT_AND_KNOWN_RNA_SEGMENTS_VERIFIED; readable final evidence review and comparison remain required',flush=True)


if __name__=='__main__':main()
