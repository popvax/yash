"""Continue from completed paired/native-normal jobs through fresh RNA ranking."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import sys
import time

CODE = Path(__file__).resolve().parent
HERE = CODE.parent;REPO=HERE.parents[1]
LOCAL=REPO/'work/matched-normal-20261107'
REMOTE='results/TN25-112103/matched-normal-20261107'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--execute',action='store_true')
    args=parser.parse_args();print(json.dumps(dict(execute=args.execute,through='normal-self models, paired annotation, fresh RNA, native affinity ranking and integrated selection')),flush=True)
    if not args.execute:return
    import modal
    volume=modal.Volume.from_name('research-lens-workspace')
    def wait_worker_group(file,key):
        last=None
        while True:
            path=LOCAL/file
            if path.exists():
                state=json.loads(path.read_text());item=state['workers'].get(key)
                if item:
                    sandbox=modal.Sandbox.from_id(item['sandbox_id']);rc=sandbox.poll()
                    if rc==0:return
                    if rc is not None:raise RuntimeError(key+' terminated: '+str(rc))
                    observation=dict(waiting_for=key,sandbox_id=sandbox.object_id,verified_live=True)
                else:
                    live=[]
                    for name,worker in state['workers'].items():
                        rc=modal.Sandbox.from_id(worker['sandbox_id']).poll()
                        if rc not in [None,0]:raise RuntimeError(name+' terminated: '+str(rc))
                        if rc is None:live.append(worker['sandbox_id'])
                    observation=dict(waiting_for=key,live_prerequisite_sandboxes=live)
            else:
                state=json.loads((LOCAL/'alignment-state.json').read_text());sandbox=modal.Sandbox.from_id(state['sandbox_id']);rc=sandbox.poll()
                if rc not in [None,0]:raise RuntimeError('Alignment terminated: '+str(rc))
                observation=dict(waiting_for=file,alignment_id=sandbox.object_id,alignment_live=rc is None)
            if observation!=last:print(json.dumps(observation),flush=True);last=observation
            time.sleep(30)
    def stage(name):
        subprocess.run([sys.executable,str(CODE / 'run_stage.py'),name,'--execute'],cwd=REPO,check=True)
        state=json.loads((LOCAL/(name+'-state.json')).read_text());sandbox=modal.Sandbox.from_id(state['sandbox_id'])
        while sandbox.poll() is None:
            print(json.dumps(dict(stage=name,sandbox_id=sandbox.object_id,verified_live=True)),flush=True);time.sleep(30)
        folder={'normal_self':'normal-self'}.get(name,name)
        directory=LOCAL/folder
        cached_receipt=directory/'complete.json'
        if sandbox.poll()==0 and cached_receipt.exists():
            cached=json.loads(cached_receipt.read_text())
            if cached.get('status')=='COMPLETED' and all(
                    metadata['bytes']>40000000 or ((directory/filename).exists()
                        and (directory/filename).stat().st_size==metadata['bytes']
                        and hashlib.sha256((directory/filename).read_bytes()).hexdigest()==metadata['sha256'])
                    for filename,metadata in cached['outputs'].items()):
                print(json.dumps(dict(stage=name,status='REUSED_VERIFIED_COMPLETED_STAGE',sandbox_id=sandbox.object_id)),flush=True)
                return
        stdout=sandbox.stdout.read();stderr=sandbox.stderr.read()
        (LOCAL/(name+'.stdout.log')).write_text(stdout);(LOCAL/(name+'.stderr.log')).write_text(stderr)
        if sandbox.poll()!=0:raise RuntimeError(name+' terminated: '+stderr[-3000:])
        body=b''.join(volume.read_file(REMOTE+'/'+folder+'/complete.json'));done=json.loads(body)
        assert done['status']=='COMPLETED';directory=LOCAL/folder;directory.mkdir(parents=True,exist_ok=True)
        (directory/'complete.json').write_bytes(body)
        for filename,metadata in done['outputs'].items():
            if metadata['bytes']>40000000:continue
            path=directory/filename
            if (path.exists() and path.stat().st_size==metadata['bytes']
                    and hashlib.sha256(path.read_bytes()).hexdigest()==metadata['sha256']):
                continue
            content=b''.join(volume.read_file(REMOTE+'/'+folder+'/'+filename))
            assert hashlib.sha256(content).hexdigest()==metadata['sha256']
            path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(content)
        print(json.dumps(dict(stage=name,status='COMPLETED',receipt=done)),flush=True)
    wait_worker_group('normal-coding-state.json','gather')
    stage('normal_self')
    wait_worker_group('paired-state.json','filter')
    for name in ['annotate','rna','rank','integrate']:stage(name)
    print('FRESH_PAIRED_RNA_RANKING_AND_INTEGRATED_SELECTION_COMPLETE; construct verification remains required',flush=True)


if __name__=='__main__':main()
