"""Read-only status from actual Modal handles; never starts or restarts jobs."""
from pathlib import Path
import argparse
import concurrent.futures
import datetime
import json

CODE = Path(__file__).resolve().parent
HERE = CODE.parent;REPO=HERE.parents[1]
LOCAL=REPO/'work/matched-normal-20261107'
REMOTE='/data/results/TN25-112103/matched-normal-20261107'
NAMES=['resources','alignment','normal_metadata','bqsr','bqsr-scatter','identity_precal','identity','audit','normal-coding',
       'paired','normal_self','annotate','rna','rank','integrate','construct','report']


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    if not any((LOCAL/(name+'-state.json')).exists() for name in NAMES):
        recorded=json.loads((HERE/'execution/completion-audit.json').read_text())
        result=dict(source='recorded_completed_results',design_id='design-2',status=recorded['status'],
                    selected_regions=recorded['selected_regions'],protein_aa=recorded['protein_aa'],
                    live_job_handles_available=False,
                    note='This clone has no local Modal handle files. These are completed recorded results, not live job status. Use yash-vaccine verify to check them.')
        if args.output:
            args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result,indent=2));return
    import modal
    entries=[];not_launched=[]
    for name in NAMES:
        path=LOCAL/(name+'-state.json')
        if not path.exists():not_launched.append(name);continue
        state=json.loads(path.read_text())
        workers=state.get('workers') or {name:state}
        entries.extend((name,key,item['sandbox_id']) for key,item in workers.items())
    def poll(entry):
        stage,key,identifier=entry
        try:
            sandbox=modal.Sandbox.from_id(identifier);rc=sandbox.poll()
            return dict(stage=stage,worker=key,sandbox_id=identifier,exit_code=rc,
                        state='LIVE' if rc is None else 'SUCCEEDED' if rc==0 else 'TERMINAL_FAILURE')
        except Exception as error:
            return dict(stage=stage,worker=key,sandbox_id=identifier,state='OBSERVATION_ERROR',error=repr(error))
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:workers=list(pool.map(poll,entries))
    snapshot=None
    alignment=next((row for row in workers if row['stage']=='alignment'),None)
    if alignment and alignment['state']=='LIVE':
        code='''from pathlib import Path
import json,os
root=Path("REMOTE/alignment")
status=json.loads((root/'status.json').read_text());status.pop('progress_log_tail',None)
inputs={}
for proc in Path('/proc').iterdir():
 if not proc.name.isdigit():continue
 try:
  if not (proc/'cmdline').read_bytes().startswith(b'fastp\\0'):continue
  for fd in (proc/'fd').iterdir():
   try:
    target=fd.resolve()
    if str(target).endswith('.fastq.gz'):
     pos=int(next(line.split()[1] for line in (proc/'fdinfo'/fd.name).read_text().splitlines() if line.startswith('pos:')))
     inputs[str(target)]=[pos,target.stat().st_size]
   except (OSError,StopIteration):pass
 except OSError:pass
v=os.statvfs('/tmp')
print(json.dumps(dict(status=status,input_fraction_consumed=sum(x[0] for x in inputs.values())/sum(x[1] for x in inputs.values()) if inputs else None,
 temporary_available_bytes=v.f_bavail*v.f_frsize,input_progress_is_not_pipeline_completion=True)))
'''.replace('REMOTE',REMOTE)
        try:
            sandbox=modal.Sandbox.from_id(alignment['sandbox_id']);proc=sandbox.exec('python','-c',code)
            body=proc.stdout.read();error=proc.stderr.read()
            snapshot=json.loads(body) if body.strip() else dict(observation_error=error,retry_same_handle=True)
        except Exception as error:snapshot=dict(observation_error=repr(error),retry_same_handle=True)
    counts={stage:{label:sum(r['stage']==stage and r['state']==label for r in workers)
                   for label in ['LIVE','SUCCEEDED','TERMINAL_FAILURE','OBSERVATION_ERROR']} for stage in NAMES if stage not in not_launched}
    result=dict(observed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),stages=counts,
                not_yet_launched=not_launched,alignment=snapshot,workers=workers,
                overall_completion='Not inferred from handles; final receipts and scientific review are still required')
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(result,workers=[r for r in workers if r['state'] in ['LIVE','TERMINAL_FAILURE','OBSERVATION_ERROR']]),indent=2))


if __name__=='__main__':main()
