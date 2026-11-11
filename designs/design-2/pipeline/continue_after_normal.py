"""Continue verified normal processing through paired discovery automatically.

The BQSR controller is started separately. This controller observes its actual
sandbox handles, then runs native identity/QC, the current-target audit and the
new paired-discovery controller. Downstream RNA/construct work remains required.
"""
from pathlib import Path
import argparse
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
    args=parser.parse_args();print(json.dumps(dict(execute=args.execute,through='paired discovery and native filtering')),flush=True)
    if not args.execute:return
    import modal
    volume=modal.Volume.from_name('research-lens-workspace')
    def collect(folder,names):
        completed=json.loads(b''.join(volume.read_file(REMOTE+'/'+folder+'/complete.json')))
        assert completed['status']=='COMPLETED'
        out=LOCAL/folder;out.mkdir(parents=True,exist_ok=True)
        (out/'complete.json').write_text(json.dumps(completed,indent=2)+'\n')
        for name in names:
            body=b''.join(volume.read_file(REMOTE+'/'+folder+'/'+name))
            expected=completed.get('outputs',{}).get(name,{}).get('sha256')
            if expected:
                import hashlib
                assert hashlib.sha256(body).hexdigest()==expected
            (out/name).write_bytes(body)
        return completed
    def wait_stage(stage):
        state=json.loads((LOCAL/(stage+'-state.json')).read_text())
        sandbox=modal.Sandbox.from_id(state['sandbox_id'])
        while True:
            result=sandbox.poll()
            if result==0:break
            if result is not None:
                stdout=sandbox.stdout.read();stderr=sandbox.stderr.read()
                (LOCAL/(stage+'.stdout.log')).write_text(stdout);(LOCAL/(stage+'.stderr.log')).write_text(stderr)
                raise RuntimeError(stage+' is terminal: '+str(result)+' '+stderr[-2000:])
            print(json.dumps(dict(stage=stage,sandbox_id=sandbox.object_id,verified_live=True)),flush=True)
            time.sleep(30)
    def wait_normal():
        while True:
            state_path=LOCAL/'bqsr-state.json'
            if state_path.exists():
                state=json.loads(state_path.read_text())
                apply=state['workers'].get('apply')
                if apply:
                    sb=modal.Sandbox.from_id(apply['sandbox_id']);rc=sb.poll()
                    if rc==0:
                        rec=collect('bqsr',[])
                        assert rec['BQSR_applied']
                        return
                    if rc is not None:raise RuntimeError('BQSR application is terminal: '+str(rc))
                observed=[]
                for name,item in state['workers'].items():
                    rc=modal.Sandbox.from_id(item['sandbox_id']).poll()
                    if rc not in [None,0]:raise RuntimeError('BQSR worker is terminal: '+name+' '+str(rc))
                    if rc is None:observed.append(item['sandbox_id'])
                print(json.dumps(dict(waiting_for='BQSR',verified_live_sandboxes=observed)),flush=True)
            else:
                state=json.loads((LOCAL/'alignment-state.json').read_text())
                sb=modal.Sandbox.from_id(state['sandbox_id']);rc=sb.poll()
                if rc not in [None,0]:raise RuntimeError('Alignment is terminal: '+str(rc))
                print(json.dumps(dict(waiting_for='normal alignment/BQSR controller',alignment_sandbox=sb.object_id,
                                      alignment_live=rc is None,alignment_completed=rc==0)),flush=True)
            time.sleep(30)
    wait_normal()
    subprocess.run([sys.executable,str(CODE / 'run_stage.py'),'identity','--execute'],cwd=REPO,check=True)
    wait_stage('identity')
    identity=collect('identity',['fingerprints.metrics.tsv','fingerprint-verdict.json','normal-contamination.table',
                                 'tumour-contamination.table','allele-fraction-segments.table'])
    assert identity['identity']=='MATCH'
    subprocess.run([sys.executable,str(CODE / 'run_stage.py'),'audit','--execute'],cwd=REPO,check=True)
    # The explicit 17-target audit and unrestricted discovery are independent.
    log=(LOCAL/'paired-controller.log').open('a')
    controller=subprocess.Popen([sys.executable,str(CODE / 'run_paired.py'),'--execute'],cwd=REPO,stdout=log,stderr=log)
    (LOCAL/'paired-controller-process.json').write_text(json.dumps(dict(pid=controller.pid,command=['run_paired.py','--execute']),indent=2)+'\n')
    wait_stage('audit')
    collect('audit',['current-mutation-audit.json','current-mutation-audit.tsv'])
    print('CURRENT_MUTATION_NORMAL_AUDIT_AVAILABLE',flush=True)
    while controller.poll() is None:
        print(json.dumps(dict(paired_controller_pid=controller.pid,verified_live=True)),flush=True)
        time.sleep(30)
    log.close()
    if controller.returncode:raise RuntimeError('Paired discovery controller failed; inspect paired-controller.log')
    collect('paired/filter',['filtering-stats.tsv','scatter-conflicts.json'])
    print('PAIRED_DISCOVERY_FILTERING_AND_CURRENT_TARGET_AUDIT_COMPLETE; downstream RNA and construct work remain',flush=True)


if __name__=='__main__':main()
