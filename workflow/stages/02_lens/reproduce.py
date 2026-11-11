"""Reproduce all three sources in a NEW Modal output directory.
Requires modal==1.6.3, an authenticated account, and the saved research-lens-workspace
volume. The saved reference/input/container artifacts are intentionally reused.
Usage: python reproduce.py --output-root /data/results/TN25-112103/NAME
The original execution involved recorded repairs; this clean entrypoint has not
been separately run end to end. It uses the same scripts, inputs and patches.
"""
import argparse,datetime,hashlib,json,re,time
from pathlib import Path
import modal
R=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--output-root',required=True);a=p.parse_args()
assert re.fullmatch(r'/data/results/TN25-112103/[A-Za-z0-9_-]+',a.output_root)
assert a.output_root!='/data/results/TN25-112103/lens-expression-aware-20261030'
state={'root':a.output_root,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'compute_terminated':False};cpu=gpu=None
volume=modal.Volume.from_name('research-lens-workspace')
statefile=R/'reproduction-state.json'
def save():statefile.write_text(json.dumps(state,indent=2)+'\n')
def committed_bytes(path, expected_hash):
 for _ in range(60):
  try:
   data=b''.join(volume.read_file(path.removeprefix('/data/')))
   if hashlib.sha256(data).hexdigest()==expected_hash:return data
  except (FileNotFoundError,modal.exception.NotFoundError):pass
  time.sleep(5)
 raise TimeoutError('Committed volume output has not matched its writer: '+path)
def execute(s,*args,timeout=43200):
 q=s.exec(*args,timeout=timeout)
 for line in q.stdout:print(line,end='',flush=True)
 err=q.stderr.read()
 if err:print(err,flush=True)
 if q.wait():raise RuntimeError('Command failed: '+repr(args))
try:
 cpu=modal.Sandbox.create(app=modal.App.lookup('research-lens'),image=modal.Image.from_id('im-01M51KXRB2T0C1YGNP5JTEPDVR'),cpu=8,memory=65536,runtime='vm',timeout=43200,volumes={'/data':modal.Volume.from_name('research-lens-workspace')})
 state['sandbox_id']=cpu.object_id;save()
 execute(cpu,'test','!','-e',a.output_root)
 execute(cpu,'mkdir','-p','/work/lens-expression-aware/patches')
 for f in R.rglob('*'):
  if f.is_file() and f.suffix in ('.py','.sh','.patch','.json'):
   cpu.filesystem.write_bytes(f.read_bytes(),'/work/lens-expression-aware/'+str(f.relative_to(R)))
 execute(cpu,'env','LENS_RUN_ROOT='+a.output_root,'bash','/work/lens-expression-aware/canonical_run.sh','pre')
 input_path=a.output_root+'/prediction/unique-peptides.txt'
 committed_bytes(input_path,hashlib.sha256(cpu.filesystem.read_bytes(input_path)).hexdigest())
 gpu=modal.Sandbox.create(app=modal.App.lookup('research-lens'),image=modal.Image.from_id('im-DBMghmAQfB3Oz7McdNZeic'),gpu='A10',cpu=2,memory=8192,timeout=3600,volumes={'/data':modal.Volume.from_name('research-lens-workspace')})
 state['gpu_sandbox_id']=gpu.object_id;save()
 gpu.filesystem.write_text((R/'predict_affinity.py').read_text(),'/tmp/predict_affinity.py')
 execute(gpu,'env','LENS_RUN_ROOT='+a.output_root,'python3','/tmp/predict_affinity.py',timeout=3300)
 model_bytes=gpu.filesystem.read_bytes(a.output_root+'/prediction/model-run.json')
 model_run=json.loads(model_bytes)
 data=committed_bytes(a.output_root+'/prediction/affinity-all-unique.tsv',model_run['affinity_output_sha256'])
 gpu.terminate();gpu=None;state['gpu_terminated']=True;save()
 cpu.filesystem.write_bytes(data,'/work/lens-expression-aware/affinity-all-unique.tsv')
 state['affinity_output_sha256']=hashlib.sha256(data).hexdigest();save()
 cpu.filesystem.write_bytes(data,a.output_root+'/prediction/affinity-all-unique.tsv')
 cpu.filesystem.write_bytes(model_bytes,a.output_root+'/prediction/model-run.json')
 cpu.filesystem.write_text(json.dumps({'sha256':state['affinity_output_sha256'],'bytes':len(data)},indent=2)+'\n',a.output_root+'/provenance/prediction-transfer.json')
 execute(cpu,'env','LENS_RUN_ROOT='+a.output_root,'bash','/work/lens-expression-aware/canonical_run.sh','post')
 state['completed']=True
finally:
 if gpu:gpu.terminate();state['gpu_terminated']=True
 if cpu:cpu.terminate();state['compute_terminated']=True
 state['finished_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();save()
