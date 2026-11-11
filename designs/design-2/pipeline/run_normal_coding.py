"""Produce normal coding genotypes after native sample identity is verified."""
from pathlib import Path
import argparse
import concurrent.futures
import hashlib
import json
import threading
import time

CODE = Path(__file__).resolve().parent
HERE = CODE.parent;REPO=HERE.parents[1]
REMOTE='/data/results/TN25-112103/matched-normal-20261107'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--execute',action='store_true')
    parser.add_argument('--wait-identity',action='store_true')
    parser.add_argument('--workers',type=int,default=16,choices=range(1,25))
    args=parser.parse_args();print(json.dumps(dict(stage='normal coding genotypes',execute=args.execute,workers=args.workers)),flush=True)
    if not args.execute:return
    import modal
    local=REPO/'work/matched-normal-20261107';volume=modal.Volume.from_name('research-lens-workspace')
    while True:
        state_path=local/'identity-state.json'
        if state_path.exists():
            state=json.loads(state_path.read_text());sandbox=modal.Sandbox.from_id(state['sandbox_id']);code=sandbox.poll()
            if code==0:break
            if code is not None:raise RuntimeError('Identity stage is terminal: '+str(code))
            observed=[sandbox.object_id]
        else:
            state=json.loads((local/'alignment-state.json').read_text());sandbox=modal.Sandbox.from_id(state['sandbox_id']);code=sandbox.poll()
            if code not in [None,0]:raise RuntimeError('Alignment is terminal: '+str(code))
            observed=[sandbox.object_id] if code is None else []
        if not args.wait_identity:raise RuntimeError('Identity is not complete')
        print(json.dumps(dict(waiting_for='native identity verdict',verified_live_sandboxes=observed)),flush=True)
        time.sleep(30)
    identity=json.loads(b''.join(volume.read_file(REMOTE.removeprefix('/data/')+'/identity/complete.json')))
    assert identity['status']=='COMPLETED' and identity['identity']=='MATCH'
    worker=CODE / 'normal_coding.py';common=CODE / 'common.py';digest=hashlib.sha256(worker.read_bytes()+common.read_bytes()).hexdigest()
    code_root=REMOTE+'/code/normal-coding-'+digest[:12]
    for path in [worker,common]:
        key=(code_root+'/'+path.name).removeprefix('/data/')
        try:
            with volume.batch_upload() as up:up.put_file(str(path),key)
        except FileExistsError:assert b''.join(volume.read_file(key))==path.read_bytes()
    sf=local/'normal-coding-state.json'
    state=json.loads(sf.read_text()) if sf.exists() else dict(bundle_sha256=digest,workers={})
    assert state['bundle_sha256']==digest
    lock=threading.Lock()
    def save_state():
        temporary=sf.with_suffix('.tmp');temporary.write_text(json.dumps(state,indent=2)+'\n');temporary.replace(sf)
    def run(action,task=None):
        name=action+('-'+task if task else '')
        with lock:old=state['workers'].get(name)
        if old:
            sandbox=modal.Sandbox.from_id(old['sandbox_id'])
            if sandbox.poll() not in [None,0]:raise RuntimeError(name+' is terminal; review before retry')
        else:
            command=['python','-u',code_root+'/normal_coding.py',REMOTE,action]+([task] if task else [])
            sandbox=modal.Sandbox.create(*command,app=modal.App.lookup('research-yash-matched-normal',create_if_missing=True),
                image=modal.Image.from_id('im-sxSffXmGfIo7J1b56AqMFq'),cpu=2 if action=='call' else 4,
                memory=8192 if action=='call' else 16384,timeout=21600,volumes={'/data':volume})
            with lock:state['workers'][name]=dict(sandbox_id=sandbox.object_id,status='running');save_state()
            print('START',name,sandbox.object_id,flush=True)
        stdout=sandbox.stdout.read();stderr=sandbox.stderr.read();sandbox.wait();rc=sandbox.poll()
        (local/('normal-coding-'+name+'.stdout.log')).write_text(stdout);(local/('normal-coding-'+name+'.stderr.log')).write_text(stderr)
        with lock:state['workers'][name].update(status='completed' if rc==0 else 'failed',exit_code=rc);save_state()
        if rc:raise RuntimeError('normal coding '+name+' failed: '+stderr[-2000:])
        folder='normal-coding/calling/'+task if action=='call' else 'normal-coding/'+action
        body=b''.join(volume.read_file(REMOTE.removeprefix('/data/')+'/'+folder+'/complete.json'))
        rec=json.loads(body);assert rec['status']=='COMPLETED'
        (local/('normal-coding-'+name+'.complete.json')).write_bytes(body)
        print('DONE',name,flush=True);return rec
    run('prepare')
    plan=json.loads(b''.join(volume.read_file(REMOTE.removeprefix('/data/')+'/normal-coding-plan.json')))
    (local/'normal-coding-plan.json').write_text(json.dumps(plan,indent=2)+'\n')
    pilot=next(s for s in plan['shards'] if s['chrom']=='chr22' and s['bases']>100000)
    run('call',pilot['id'])
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        list(pool.map(lambda shard:run('call',shard['id']),[s for s in plan['shards'] if s!=pilot]))
    rec=run('gather');state['status']='COMPLETED';save_state()
    print(json.dumps(dict(status='NORMAL_CODING_GENOTYPES_COMPLETE',counts=rec['counts'])),flush=True)


if __name__=='__main__':main()
