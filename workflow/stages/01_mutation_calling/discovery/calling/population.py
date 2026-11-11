"""Prepare pinned gnomAD v4.1 genome references on the research volume.

Run inside the GATK image. Source VCF fields/QC are preserved in a compact
annotation file; the separate calling resource contains PASS records only.
"""
from pathlib import Path
import concurrent.futures,datetime,hashlib,json,subprocess,sys,tempfile,time,urllib.request

ROOT=Path(sys.argv[1]);CHROM=sys.argv[2];PARTS=int(sys.argv[3]) if len(sys.argv)>3 else 0
assert PARTS==0 or 2<=PARTS<=16
D=ROOT/'population'/CHROM;D.mkdir(parents=True,exist_ok=True)
BASE='https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/gnomad/v4.1/genomes/'
URL=BASE+f'gnomad.genomes.v4.1.sites.{CHROM}_trimmed.vcf.bgz'
KEEP={'AF','AC','AN','grpmax','AF_grpmax','AN_grpmax','AC_grpmax','AF_sas','AN_sas','AC_sas'}
probe=json.loads((ROOT/'population-source-probe.json').read_text())
import re
fields={re.search(r'ID=([^,]+)',s).group(1) for s in probe['genomes']['INFO']}
assert KEEP<=fields
remove=['ID','QUAL']+['INFO/'+k for k in sorted(fields-KEEP)]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def write(name,obj):
 p=D/name;t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(obj,indent=2)+'\n');t.replace(p)
def run(args):
 with (D/'commands.log').open('a') as log:
  log.write(json.dumps(args)+'\n');log.flush();subprocess.run(args,stdout=log,stderr=log,check=True)
start=time.time();started=datetime.datetime.now(datetime.timezone.utc).isoformat()
if (D/'complete.json').exists():
 x=json.loads((D/'complete.json').read_text())
 assert all(sha(D/f)==v['sha256'] for f,v in x['outputs'].items())
 print(json.dumps(x));sys.exit(0)
urlmeta={}
with urllib.request.urlopen(urllib.request.Request(URL,method='HEAD'),timeout=60) as r:
 urlmeta={k:r.headers.get(k) for k in ['Content-Length','ETag','Last-Modified']}
write('status.json',dict(phase='downloading_and_extracting',url=URL,started_utc=started,source=urlmeta))
out=D/'annotation.vcf.gz';tmp=D/'annotation.partial.vcf.gz'
parts=[];part_dir=None
if PARTS:
 size=int(urlmeta['Content-Length']);part_dir=tempfile.TemporaryDirectory(prefix='gnomad-ranges-')
 def download(i):
  lo=size*i//PARTS;hi=size*(i+1)//PARTS-1;path=Path(part_dir.name)/str(i)
  for attempt in range(3):
   try:
    request=urllib.request.Request(URL,headers={'Range':f'bytes={lo}-{hi}','If-Range':urlmeta['ETag']})
    with urllib.request.urlopen(request,timeout=90) as source,path.open('wb') as dest:
     assert source.status==206 and source.headers.get('Content-Range')==f'bytes {lo}-{hi}/{size}'
     assert source.headers.get('ETag')==urlmeta['ETag']
     while True:
      block=source.read(1024*1024)
      if not block:break
      dest.write(block)
    assert path.stat().st_size==hi-lo+1
    return path
   except Exception:
    if attempt==2:raise
    time.sleep(2)
 with concurrent.futures.ThreadPoolExecutor(max_workers=min(8,PARTS)) as pool:parts=list(pool.map(download,range(PARTS)))
 write('status.json',dict(phase='extracting_verified_ranges',parts=PARTS,source_bytes=size,source=urlmeta))
for attempt in range(3):
 total=0;h=hashlib.sha256()
 try:
  with (D/'commands.log').open('a') as log:
   args=['bcftools','annotate','-x',','.join(remove),'-Oz','-o',str(tmp),'-'];log.write(json.dumps(args)+'\n');log.flush()
   p=subprocess.Popen(args,stdin=subprocess.PIPE,stdout=log,stderr=log)
   try:
    if parts:
     for path in parts:
      with path.open('rb') as source:
       while True:
        b=source.read(1024*1024)
        if not b:break
        h.update(b);total+=len(b);p.stdin.write(b)
    else:
     with urllib.request.urlopen(URL,timeout=180) as source:
      while True:
       b=source.read(1024*1024)
       if not b:break
       h.update(b);total+=len(b);p.stdin.write(b)
    p.stdin.close();rc=p.wait()
    assert rc==0,f'bcftools annotate exited {rc}'
   except BaseException:
    p.kill();p.wait();raise
  if urlmeta['Content-Length']:assert total==int(urlmeta['Content-Length'])
  tmp.replace(out);break
 except Exception as e:
  write('status.json',dict(phase='source_retry',attempt=attempt+1,error=repr(e),bytes=total))
  if attempt==2:raise
  time.sleep(2)
run(['bcftools','index','-t','-f',str(out)])
if part_dir:part_dir.cleanup()
caller=D/'caller-pass.vcf.gz'
run(['bcftools','view','-f','PASS','-Oz','-o',str(caller),str(out)])
run(['bcftools','index','-t','-f',str(caller)])
import pysam
v=pysam.VariantFile(str(caller));spot={}
for gene,c,pos,ref,alt,expected in [('ATRX','chrX',77682471,'C','G',.5158322228567802),('CYP2D6','chr22',42130692,'G','A',.19041485180070117),('KMT2D','chr12',49052257,'C','T',.0891850591779953)]:
 if c!=CHROM:continue
 rows=[r for r in v.fetch(c,pos-1,pos) if r.pos==pos and r.ref==ref and alt in r.alts]
 assert len(rows)==1,(gene,'missing or duplicate frequency record')
 r=rows[0];af=float(r.info['AF'][r.alts.index(alt)]);assert abs(af-expected)<1e-4,(gene,af,expected)
 spot[gene]=dict(chrom=c,position=pos,ref=ref,alt=alt,AF=af,QC=list(r.filter))
outputs={p.name:dict(bytes=p.stat().st_size,sha256=sha(p)) for p in [out,Path(str(out)+'.tbi'),caller,Path(str(caller)+'.tbi')]}
report=dict(status='VERIFIED',chrom=CHROM,source_url=URL,source_http_metadata=urlmeta,source_compressed_bytes=total,
 source_compressed_sha256=h.hexdigest(),retained_INFO=sorted(KEEP),annotation_policy='All source FILTER labels retained',
 calling_policy='Only source PASS records',spot_checks=spot,outputs=outputs,elapsed_seconds=time.time()-start,
 started_utc=started,completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
write('complete.json',report);write('status.json',report);print(json.dumps(report),flush=True)
