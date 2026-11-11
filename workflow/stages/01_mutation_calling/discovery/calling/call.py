"""One independent Mutect2 discovery shard; no given alleles/Caris input."""
from pathlib import Path
import datetime,hashlib,json,subprocess,sys,time
R=Path(sys.argv[1]);sid=sys.argv[2];plan=json.loads((R/'calling-plan.json').read_text())
s=next(x for x in plan['shards'] if x['id']==sid);D=R/'calling'/sid;D.mkdir(parents=True,exist_ok=True)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def save(name,v):
 p=D/name;t=p.with_suffix('.tmp');t.write_text(json.dumps(v,indent=2)+'\n');t.replace(p)
if (D/'complete.json').exists():
 c=json.loads((D/'complete.json').read_text());assert all(sha(D/k)==v['sha256'] for k,v in c['outputs'].items());print(json.dumps(c));sys.exit(0)
source=R/'population'/s['chrom']/'complete.json';assert source.is_file(),'Population resource must finish before discovery'
pop=json.loads(source.read_text());resource=source.parent/'caller-pass.vcf.gz'
assert pop['status']=='VERIFIED' and resource.stat().st_size==pop['outputs'][resource.name]['bytes']
assert plan['given_alleles'] is None and plan['matched_normal'] is None and plan['panel_of_normals'] is None
command=['/gatk/gatk','--java-options','-Xmx6g -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1','Mutect2',
 '-R',plan['reference'],'-I',plan['bam'],'--germline-resource',str(resource),'-L',s['interval'],
 '--native-pair-hmm-threads','1','--f1r2-tar-gz',str(D/'f1r2.tar.gz'),'-O',str(D/'unfiltered.vcf.gz')]
assert '--alleles' not in command and '--genotype-germline-sites' not in command
save('command.json',dict(command=command,discovery=True,Caris_input=False,shard=s,population_sha256=pop['outputs']['caller-pass.vcf.gz']['sha256']))
start=time.time();save('status.json',dict(phase='calling',started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),command=command))
with (D/'gatk.log').open('w') as log:subprocess.run(command,stdout=log,stderr=log,check=True)
outputs={p.name:dict(bytes=p.stat().st_size,sha256=sha(p)) for p in [D/'unfiltered.vcf.gz',D/'unfiltered.vcf.gz.tbi',D/'unfiltered.vcf.gz.stats',D/'f1r2.tar.gz']}
report=dict(status='COMPLETED',shard=s,elapsed_seconds=time.time()-start,completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),outputs=outputs)
save('complete.json',report);save('status.json',report);print(json.dumps(report),flush=True)
