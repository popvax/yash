"""Gather independent discovery shards and filter jointly across the sample."""
from pathlib import Path
import collections,datetime,hashlib,json,pysam,subprocess,sys,time
R=Path(sys.argv[1]);D=R/'gather';D.mkdir(parents=True,exist_ok=True)
plan=json.loads((R/'calling-plan.json').read_text());support=json.loads((R/'support/complete.json').read_text())
assert support['status']=='COMPLETED' and len(plan['shards'])>=24 and plan['given_alleles'] is None
chromosomes={s['chrom'] for s in plan['shards']}
assert chromosomes=={'chr'+str(c) for c in list(range(1,23))+['X','Y']}
reference=pysam.FastaFile(plan['reference'])
for chrom in chromosomes:
 pieces=sorted((s for s in plan['shards'] if s['chrom']==chrom),key=lambda s:s['start'])
 assert pieces[0]['start']==1 and pieces[-1]['end']==reference.get_reference_length(chrom) and all(a['end']+1==b['start'] for a,b in zip(pieces,pieces[1:])),chrom
reference.close()
records=[json.loads((R/'calling'/s['id']/'complete.json').read_text()) for s in plan['shards']]
assert all(r['status']=='COMPLETED' and r['shard']==s for r,s in zip(records,plan['shards']))
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def save(name,value):
 p=D/name;t=p.with_suffix('.tmp');t.write_text(json.dumps(value,indent=2)+'\n');t.replace(p)
def run(name,args,outputs):
 args=list(map(str,args));mark=D/(name+'.done.json')
 if mark.exists():
  m=json.loads(mark.read_text());assert m['command']==args and all(sha(p)==h for p,h in m['outputs'].items());return
 save('status.json',dict(phase=name,updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
 with (D/(name+'.log')).open('w') as log:subprocess.run(args,stdout=log,stderr=log,check=True)
 save(mark.name,dict(command=args,outputs={str(p):sha(p) for p in outputs}))
def gatk(tool,*args):return ['/gatk/gatk','--java-options','-Xmx11g -XX:ParallelGCThreads=2',tool,*args]
started=time.time();vcfs=[R/'calling'/s['id']/'unfiltered.vcf.gz' for s in plan['shards']]
# Before merging, verify every shard's content digest, not just its marker.
for s,receipt in zip(plan['shards'],records):
 for name,v in receipt['outputs'].items():assert sha(R/'calling'/s['id']/name)==v['sha256']
merged=D/'unfiltered.vcf.gz';stats=D/'unfiltered.vcf.gz.stats';orientation=D/'orientation-model.tar.gz'
run('merge-vcfs',gatk('MergeVcfs',*[a for p in vcfs for a in ['-I',p]],'-O',merged),[merged])
run('merge-stats',gatk('MergeMutectStats',*[a for p in vcfs for a in ['--stats',str(p)+'.stats']],'-O',stats),[stats])
run('orientation-model',gatk('LearnReadOrientationModel',*[a for s in plan['shards'] for a in ['-I',R/'calling'/s['id']/'f1r2.tar.gz']],'-O',orientation),[orientation])
# Caris is normalized for post-discovery annotation; it does not nominate calls.
caris_source=Path('/data/inputs/TN25-112103/DNA_TN25-112103.vcf')
caris_tagged=D/'caris-comparison-input.vcf.gz';caris_meta={}
with pysam.VariantFile(str(caris_source)) as inp:
 with pysam.VariantFile(str(caris_tagged),'wz',header=inp.header) as out:
  for i,r in enumerate(inp,1):
   cid='caris_'+str(i);r.id=cid
   caris_meta[cid]=dict(id=cid,chrom=r.contig,position=r.pos,ref=r.ref,alts=r.alts,filter=list(r.filter),info=dict(r.info),samples={name:dict(v.items()) for name,v in r.samples.items()})
   out.write(r)
save('caris-comparison-metadata.json',caris_meta)
caris_normalized=D/'caris-comparison-normalized.vcf.gz'
run('normalize-Caris-for-comparison-only',['bcftools','norm','-f',plan['reference'],'-m','-any','-Oz','-o',caris_normalized,caris_tagged],[caris_normalized])
run('index-Caris-comparison',['bcftools','index','-t','-f',caris_normalized],[Path(str(caris_normalized)+'.tbi')])

with pysam.VariantFile(str(merged)) as calls:raw_records=sum(1 for _ in calls)
report=dict(status='COMPLETED',Caris_nominations=False,shards=len(records),counts=dict(raw_records=raw_records),
 elapsed_seconds=time.time()-started,
 input_sha256={'calling-plan.json':sha(R/'calling-plan.json'),'support-receipt':sha(R/'support/complete.json')},
 outputs={p.name:dict(bytes=p.stat().st_size,sha256=sha(p)) for p in [merged,stats,orientation,caris_normalized,Path(str(caris_normalized)+'.tbi'),D/'caris-comparison-metadata.json']},completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
save('complete.json',report);save('status.json',report);print(json.dumps(report),flush=True)
