"""Run independent paired discovery and native filtering on Modal."""
from pathlib import Path
import argparse
import concurrent.futures
import hashlib
import json
import threading

CODE = Path(__file__).resolve().parent
HERE = CODE.parent;REPO=HERE.parents[1]
REMOTE='/data/results/TN25-112103/matched-normal-20261107'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--execute',action='store_true')
    parser.add_argument('--workers',type=int,default=32,choices=range(1,33))
    args=parser.parse_args();print(json.dumps(dict(stage='paired_discovery',workers=args.workers,execute=args.execute)),flush=True)
    if not args.execute:return
    import modal
    local=REPO/'work/matched-normal-20261107';local.mkdir(parents=True,exist_ok=True)
    volume=modal.Volume.from_name('research-lens-workspace')
    normal=json.loads(b''.join(volume.read_file(REMOTE.removeprefix('/data/')+'/bqsr/complete.json')))
    identity=json.loads(b''.join(volume.read_file(REMOTE.removeprefix('/data/')+'/identity/complete.json')))
    assert normal['status']==identity['status']=='COMPLETED' and identity['identity']=='MATCH'
    worker=CODE / 'paired_calling.py';common=CODE / 'common.py'
    digest=hashlib.sha256(worker.read_bytes()+common.read_bytes()).hexdigest()
    code_root=REMOTE+'/code/paired-'+digest[:12]
    for path in [worker,common]:
        key=(code_root+'/'+path.name).removeprefix('/data/')
        try:
            with volume.batch_upload() as up:up.put_file(str(path),key)
        except FileExistsError:assert b''.join(volume.read_file(key))==path.read_bytes()
    state_path=local/'paired-state.json'
    state=json.loads(state_path.read_text()) if state_path.exists() else dict(bundle_sha256=digest,workers={})
    assert state['bundle_sha256']==digest
    lock=threading.Lock()
    def save_state():
        temporary=state_path.with_suffix('.tmp');temporary.write_text(json.dumps(state,indent=2)+'\n');temporary.replace(state_path)
    def run(action,task=None):
        name=action+('-'+task if task else '')
        with lock:old=state['workers'].get(name)
        if old:
            sb=modal.Sandbox.from_id(old['sandbox_id'])
            if sb.poll() not in [None,0]:raise RuntimeError(name+' is terminal; inspect before a new attempt')
        else:
            command=['python','-u',code_root+'/paired_calling.py',REMOTE,action]+([task] if task else [])
            sb=modal.Sandbox.create(*command,app=modal.App.lookup('research-yash-matched-normal',create_if_missing=True),
                image=modal.Image.from_id('im-sxSffXmGfIo7J1b56AqMFq'),cpu=2 if action=='call' else 4,
                memory=8192 if action=='call' else 16384,timeout=21600,volumes={'/data':volume})
            with lock:state['workers'][name]=dict(sandbox_id=sb.object_id,status='running');save_state()
            print('START',name,sb.object_id,flush=True)
        stdout=sb.stdout.read();stderr=sb.stderr.read();sb.wait();rc=sb.poll()
        (local/(name+'.stdout.log')).write_text(stdout);(local/(name+'.stderr.log')).write_text(stderr)
        with lock:state['workers'][name].update(status='completed' if rc==0 else 'failed',exit_code=rc);save_state()
        if rc:raise RuntimeError(name+' failed: '+stderr[-2500:])
        folder='paired/calling/'+task if action=='call' else 'paired/population/'+task if action=='population' else 'paired/'+action
        body=b''.join(volume.read_file(REMOTE.removeprefix('/data/')+'/'+folder+'/complete.json'))
        result=json.loads(body);assert result['status']=='COMPLETED'
        (local/(name+'.complete.json')).write_bytes(body)
        print('DONE',name,flush=True);return result
    run('prepare')
    plan=json.loads(b''.join(volume.read_file(REMOTE.removeprefix('/data/')+'/calling-plan.json')))
    (local/'paired-calling-plan.json').write_text(json.dumps(plan,indent=2)+'\n')
    # Verify one real paired shard before opening the larger compute scatter.
    pilot=min(plan['shards'],key=lambda s:s['estimated_read_weight'])
    run('call',pilot['id'])
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        list(pool.map(lambda shard:run('call',shard['id']),[s for s in plan['shards'] if s!=pilot]))
    run('gather')
    chromosomes=['chr'+str(c) for c in list(range(1,23))+['X','Y']]
    with concurrent.futures.ThreadPoolExecutor(max_workers=min(args.workers,8)) as pool:
        list(pool.map(lambda chrom:run('population',chrom),chromosomes))
    result=run('filter');state['status']='COMPLETED';save_state()
    print(json.dumps(dict(status='PAIRED_DISCOVERY_AND_FILTERING_COMPLETE',counts=result['counts'])),flush=True)


if __name__=='__main__':main()
