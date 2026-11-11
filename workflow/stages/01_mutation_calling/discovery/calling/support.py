"""Fresh sample-wide QC and contamination inputs; no Caris nomination data."""
from pathlib import Path
import csv,datetime,hashlib,json,subprocess,sys,time
R=Path(sys.argv[1]);D=R/'support';D.mkdir(parents=True,exist_ok=True)
plan=json.loads((R/'calling-plan.json').read_text())
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def save(name,value):
 p=D/name;t=p.with_suffix('.tmp');t.write_text(json.dumps(value,indent=2)+'\n');t.replace(p)
def run(name,args,outputs):
 marker=D/(name+'.done.json')
 if marker.exists():
  v=json.loads(marker.read_text());assert v['command']==list(map(str,args)) and all(sha(p)==h for p,h in v['outputs'].items());return
 save('status.json',dict(phase=name,updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
 with (D/(name+'.log')).open('w') as log:subprocess.run(list(map(str,args)),stdout=log,stderr=log,check=True)
 save(marker.name,dict(command=list(map(str,args)),outputs={str(p):sha(p) for p in outputs}))
start=time.time();bam=Path(plan['bam']);ref=Path(plan['reference'])
assert sha(bam)==plan['bam_sha256'] and sha(ref)==plan['reference_sha256']
run('bam-quickcheck',['samtools','quickcheck','-v',bam],[])
run('bam-stats',['samtools','stats','-@','2',bam],[])
common=Path('/data/raft/references/homo_sapiens/vcfs/small_exac_common_3.hg38.vcf.gz')
def gatk(tool,*args):return ['/gatk/gatk','--java-options','-Xmx10g -XX:ParallelGCThreads=2',tool,*args]
pileups=D/'tumour-pileups.table';cont=D/'contamination.table';seg=D/'allele-fraction-segments.table'
run('pileup-summaries',gatk('GetPileupSummaries','-R',ref,'-I',bam,'-V',common,'-L',common,'-O',pileups),[pileups])
run('contamination',gatk('CalculateContamination','-I',pileups,'-O',cont,'--tumor-segmentation',seg),[cont,seg])
report=dict(status='COMPLETED',elapsed_seconds=time.time()-start,input_sha256=dict(bam=plan['bam_sha256'],reference=plan['reference_sha256'],common_sites=sha(common)),
 common_sites_scope='Standard Small ExAC common-site panel, used only for contamination estimation; not the new mutation germline resource.',
 commands=[json.loads(p.read_text()) for p in sorted(D.glob('*.done.json'))],
 outputs={p.name:dict(bytes=p.stat().st_size,sha256=sha(p)) for p in [pileups,cont,seg,D/'bam-stats.log']},
 completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
save('complete.json',report);save('status.json',report);print(json.dumps(report),flush=True)
