"""Run independent mutation-evidence stages on Modal; --execute starts compute."""
from pathlib import Path
import argparse,concurrent.futures,datetime,hashlib,json,threading,time
HERE=Path(__file__).resolve().parent
REPO=next(p for p in HERE.parents if (p/'pyproject.toml').is_file())
REMOTE='/data/results/TN25-112103/mutation-evidence-v2-20261102'
IMAGE='im-sxSffXmGfIo7J1b56AqMFq';VOLUME='research-lens-workspace'
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('stage',choices=['population','calling'])
 p.add_argument('--execute',action='store_true');p.add_argument('--workers',type=int)
 p.add_argument('--population-download-parts',type=int,default=0)
 p.add_argument('--output',type=Path,default=REPO/'work/mutation-evidence-v2-20261102')
 a=p.parse_args();a.workers=a.workers or (8 if a.stage=='population' else 32);assert 1<=a.workers<=32
 chroms=['chrX','chr12','chr22']+[f'chr{x}' for x in list(range(1,23))+['Y'] if f'chr{x}' not in ['chrX','chr12','chr22']]
 cp=json.loads((a.output/'calling-plan.json').read_text()) if a.stage=='calling' else None
 tasks=chroms if cp is None else [s['id'] for s in cp['shards']]
 memory=4096 if a.stage=='population' else 8192
 assert a.population_download_parts==0 or 2<=a.population_download_parts<=16
 plan=dict(stage=a.stage,remote_root=REMOTE,volume=VOLUME,image=IMAGE,tasks=tasks,max_workers=a.workers,cpu_per_worker=2,memory_mib_per_worker=memory,population_download_parts=a.population_download_parts,execute=a.execute)
 print(json.dumps({k:v for k,v in plan.items() if k!='tasks'},indent=2),flush=True)
 print('TASKS',len(tasks),flush=True)
 if not a.execute:return
 import modal
 a.output.mkdir(parents=True,exist_ok=True);sf=a.output/(a.stage+'-state.json')
 state=json.loads(sf.read_text()) if sf.exists() else {'plan':plan,'workers':{}}
 assert state['plan']['remote_root']==REMOTE
 state.setdefault('controller_history',[]).append(dict(started_utc=now(),max_workers=a.workers))
 state['plan']=plan
 lock=threading.Lock()
 def save():
  t=sf.with_suffix('.tmp');t.write_text(json.dumps(state,indent=2)+'\n');t.replace(sf)
 volume=modal.Volume.from_name(VOLUME);image=modal.Image.from_id(IMAGE);app=modal.App.lookup('yash-mutation-evidence-v2',create_if_missing=True)
 worker_name='population.py' if a.stage=='population' else 'call.py'
 with volume.batch_upload(force=True) as upload:
  upload.put_file(str(HERE/worker_name),REMOTE.removeprefix('/data/')+'/code/'+worker_name)
  if cp:upload.put_file(str(a.output/'calling-plan.json'),REMOTE.removeprefix('/data/')+'/calling-plan.json')
 def worker(task):
  old=state['workers'].get(task);s=None
  if old:
   try:
    candidate=modal.Sandbox.from_id(old['sandbox_id']);rc=candidate.poll()
    if rc is None or rc==0:s=candidate
   except modal.exception.NotFoundError:pass
  if s is None:
   extra=[str(a.population_download_parts)] if a.stage=='population' else []
   s=modal.Sandbox.create('python','-u',REMOTE+'/code/'+worker_name,REMOTE,task,*extra,app=app,image=image,cpu=2,memory=memory,timeout=14400,volumes={'/data':volume})
   with lock:
    state['workers'][task]=dict(sandbox_id=s.object_id,status='running',started_utc=now());save()
   print('START',task,s.object_id,flush=True)
  stdout=s.stdout.read();stderr=s.stderr.read();s.wait();rc=s.poll()
  (a.output/(task+'.'+a.stage+'.stdout.log')).write_text(stdout);(a.output/(task+'.'+a.stage+'.stderr.log')).write_text(stderr)
  with lock:state['workers'][task].update(status='completed' if rc==0 else 'failed',exit_code=rc,completed_utc=now());save()
  if rc!=0:raise RuntimeError(task+' '+a.stage+' failed; inspect logs before retrying')
  folder='population' if a.stage=='population' else 'calling'
  data=b''.join(volume.read_file(REMOTE.removeprefix('/data/')+'/'+folder+'/'+task+'/complete.json'))
  (a.output/(task+'.'+a.stage+'.json')).write_bytes(data)
  print('DONE',task,flush=True);return json.loads(data)
 errors=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
  if a.stage=='population':
   remaining=[]
   for c in chroms:
    receipt=a.output/(c+'.population.json')
    if receipt.exists():
     saved=json.loads(receipt.read_text());assert saved['status']=='VERIFIED' and saved['chrom']==c
    else:remaining.append(c)
   remaining.sort(key=lambda c:state['workers'].get(c,{}).get('status')=='running',reverse=True)
   futures={pool.submit(worker,c):c for c in remaining}
   for f in concurrent.futures.as_completed(futures):
    try:f.result()
    except Exception as e:errors.append((futures[f],repr(e)));print('FAILED',errors[-1],flush=True)
  else:
   control=False;pending={s['id']:s for s in cp['shards']};active={};done=0;last_print=0
   for task in list(pending):
    receipt=a.output/(task+'.calling.json')
    if receipt.exists():
     saved=json.loads(receipt.read_text());assert saved['status']=='COMPLETED' and saved['shard']==pending[task]
     del pending[task];done+=1
   while pending or active:
    if not control:
     paths=[a.output/(c+'.population.json') for c in ['chrX','chr12','chr22']]
     if all(p.is_file() for p in paths):
      reports=[json.loads(p.read_text()) for p in paths];checks={k:v for r in reports for k,v in r['spot_checks'].items()}
      assert set(checks)=={'ATRX','CYP2D6','KMT2D'} and all(r['status']=='VERIFIED' for r in reports)
      control=True;(a.output/'germline-resource-validation.json').write_text(json.dumps(dict(status='VERIFIED',checks=checks),indent=2)+'\n')
      print('GERMLINE_RESOURCE_CONTROLS_VERIFIED',flush=True)
    slots=(min(4,a.workers) if done<2 else a.workers)-len(active)
    ready=[s for s in pending.values() if (a.output/(s['chrom']+'.population.json')).is_file()] if control else []
    ready.sort(key=lambda s:(state['workers'].get(s['id'],{}).get('status')=='running',s['estimated_read_weight']),reverse=True)
    for s in ready[:max(0,slots)]:active[pool.submit(worker,s['id'])]=s['id'];del pending[s['id']]
    finished=[f for f in active if f.done()]
    for f in finished:
     task=active.pop(f)
     try:f.result();done+=1
     except Exception as e:errors.append((task,repr(e)));print('FAILED',errors[-1],flush=True)
    if errors:raise RuntimeError(errors)
    if time.time()-last_print>30:
     print('PROGRESS',json.dumps(dict(completed=done,running=len(active),pending=len(pending),control_sites_verified=control)),flush=True);last_print=time.time()
    if pending and not active:
     ps=a.output/'population-state.json'
     if ps.exists():
      failed=[k for k,v in json.loads(ps.read_text())['workers'].items() if v['status']=='failed']
      if failed:raise RuntimeError('Population preparation failed for '+str(failed))
    if pending or active:time.sleep(3)
 if errors:raise RuntimeError(errors)
 state['status']='COMPLETED';state['completed_utc']=now();save();print(a.stage.upper()+'_COMPLETE',flush=True)
if __name__=='__main__':main()
