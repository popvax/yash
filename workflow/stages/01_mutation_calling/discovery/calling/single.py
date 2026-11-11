"""Launch or reattach a single stage. Requires explicit --execute."""
from pathlib import Path
import argparse,datetime,hashlib,json
HERE=Path(__file__).resolve().parent
REPO=next(p for p in HERE.parents if (p/'pyproject.toml').is_file())
REMOTE='/data/results/TN25-112103/mutation-evidence-v2-20261102'
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['preflight','support','gather']);p.add_argument('--execute',action='store_true');p.add_argument('--output',type=Path,default=REPO/'work/mutation-evidence-v2-20261102');a=p.parse_args()
 script=HERE/(a.stage+'.py');assert script.is_file()
 image='im-sxSffXmGfIo7J1b56AqMFq';cpu=4;memory=16384
 plan=dict(stage=a.stage,root=REMOTE,volume='research-lens-workspace',image=image,cpu=cpu,memory_mib=memory,execute=a.execute,worker_sha256=hashlib.sha256(script.read_bytes()).hexdigest())
 print(json.dumps(plan,indent=2),flush=True)
 if not a.execute:return
 import modal
 a.output.mkdir(parents=True,exist_ok=True);sf=a.output/(a.stage+'-state.json');v=modal.Volume.from_name(plan['volume']);app=modal.App.lookup('yash-mutation-evidence-v2',create_if_missing=True)
 s=None
 if sf.exists():
  old=json.loads(sf.read_text());existing=modal.Sandbox.from_id(old['sandbox_id']);rc=existing.poll()
  if rc is None or rc==0:
   assert old['worker_sha256']==plan['worker_sha256'],'Do not overwrite code for a live/completed stage';s=existing
 if s is None:
  with v.batch_upload(force=True) as upload:
   upload.put_file(str(script),REMOTE.removeprefix('/data/')+'/code/'+script.name)
   if a.stage!='preflight':upload.put_file(str(a.output/'calling-plan.json'),REMOTE.removeprefix('/data/')+'/calling-plan.json')
   if (a.output/'calling-plan-original.json').exists():upload.put_file(str(a.output/'calling-plan-original.json'),REMOTE.removeprefix('/data/')+'/calling-plan-original.json')
   upload.put_file(str(HERE/'policy.json'),REMOTE.removeprefix('/data/')+'/policy.json')
   upload.put_file(str(REPO/'config/case.json'),REMOTE.removeprefix('/data/')+'/case.json')
   upload.put_file(str(REPO/'data/checkpoints/provenance/construct-model-files.json'),REMOTE.removeprefix('/data/')+'/affinity-model-manifest.json')
  s=modal.Sandbox.create('python','-u',REMOTE+'/code/'+script.name,REMOTE,app=app,image=modal.Image.from_id(image),cpu=cpu,memory=memory,timeout=14400,volumes={'/data':v})
  plan.update(sandbox_id=s.object_id,started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat());sf.write_text(json.dumps(plan,indent=2)+'\n')
 print('SANDBOX',s.object_id,flush=True)
 out=s.stdout.read();err=s.stderr.read();s.wait();rc=s.poll()
 (a.output/(a.stage+'.stdout.log')).write_text(out);(a.output/(a.stage+'.stderr.log')).write_text(err)
 if rc!=0:
  print(out[-3000:],err[-3000:],flush=True)
  raise RuntimeError((a.stage,rc))
 data=b''.join(v.read_file(REMOTE.removeprefix('/data/')+'/'+a.stage+'/complete.json'));(a.output/(a.stage+'.complete.json')).write_bytes(data)
 report=json.loads(data)
 fields=['status','outcomes','tiers','input_variants','usable_primary_proteins','primary_DNA_and_RNA_eligible','primary_ranked_regions','primary_genes','new_peptide_HLA_predictions','unique_peptides','checks','all_verified']
 summary={k:report[k] for k in fields if k in report}
 if 'counts' in report:summary['counts']={k:v for k,v in report['counts'].items() if not k.startswith('native_') or k=='native_PASS'}
 print(json.dumps(summary,indent=2),flush=True)
 if a.stage=='preflight':
  for name in ['preflight.json','population-source-probe.json']:(a.output/name).write_bytes(b''.join(v.read_file(REMOTE.removeprefix('/data/')+'/'+name)))
 print('COMPLETED',a.stage,flush=True)
if __name__=='__main__':main()
