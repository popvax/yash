"""Read-only raw-input and population-header inspection inside the GATK image."""
from pathlib import Path
import datetime,gzip,json,pysam,subprocess,sys,urllib.request
R=Path(sys.argv[1]);R.mkdir(parents=True,exist_ok=True);D=R/'preflight';D.mkdir(exist_ok=True)
bam=Path('/data/inputs/TN25-112103/DNA_TN25-112103.bam')
ref=Path('/data/results/TN25-112103/mutect2-tumour-only-20261029/resources/reference.fa')
old=Path('/data/raft/references/homo_sapiens/vcfs/af-only-gnomad.hg38.vcf.gz')
b=pysam.AlignmentFile(str(bam),'rb');b.check_index();fa=pysam.FastaFile(str(ref));v=pysam.VariantFile(str(old));checks={}
for gene,c,pos in [('ATRX','chrX',77682471),('CYP2D6','chr22',42130692),('KMT2D','chr12',49052257)]:
 checks[gene]=[str(x).strip() for x in v.fetch(c,pos-1,pos) if x.pos==pos]
report=dict(started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),bam_bytes=bam.stat().st_size,
 bam_header_programs=b.header.to_dict().get('PG'),bam_samples=b.header.to_dict().get('RG'),
 bam_reference_agreement=all(fa.get_reference_length(c)==n for c,n in zip(b.references,b.lengths)),
 bam_index=[dict(contig=s.contig,mapped=s.mapped,unmapped=s.unmapped,length=b.get_reference_length(s.contig)) for s in b.get_index_statistics()],
 old_germline_header=str(v.header),old_resource_locus_records=checks,
 tools={name:subprocess.run(args,text=True,errors='replace',capture_output=True,check=True).stdout[:1000] for name,args in [('gatk',['/gatk/gatk','--version']),('bcftools',['bcftools','--version']),('samtools',['samtools','--version'])]})
assert report['bam_reference_agreement']
(R/'preflight.json').write_text(json.dumps(report,indent=2)+'\n')
base='https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/gnomad/v4.1/';headers={}
for cohort in ['genomes','exomes']:
 url=base+cohort+'/gnomad.'+cohort+'.v4.1.sites.chrX_trimmed.vcf.bgz'
 with urllib.request.urlopen(url,timeout=90) as stream:
  z=gzip.GzipFile(fileobj=stream);lines=[]
  for line in z:
   s=line.decode();lines.append(s)
   if s.startswith('#CHROM'):break
  z.close()
 headers[cohort]=dict(url=url,header=''.join(lines),INFO=[s.strip() for s in lines if s.startswith('##INFO')],FILTER=[s.strip() for s in lines if s.startswith('##FILTER')])
(R/'population-source-probe.json').write_text(json.dumps(headers,indent=2)+'\n')
(D/'complete.json').write_text(json.dumps(dict(status='COMPLETED',bam_reference_agreement=True,source_headers=list(headers)),indent=2)+'\n')
print(json.dumps(dict(status='COMPLETED',bam_bytes=report['bam_bytes'],bam_reference_agreement=True)),flush=True)
