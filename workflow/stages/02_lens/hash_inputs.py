"""Hash the actual RNA inputs, discovery checkpoints and required references.
Streams files; shared reference targets are hashed once but every path is recorded.
"""
import argparse,hashlib,json,os,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--root',required=True);a=p.parse_args();r=Path(a.root)
ref=Path('/data/raft/references/homo_sapiens');old=Path('/data/results/TN25-112103/original-fastqs/outputs');sample=old/'samples/TN25_research/TN25-112103/ar-TN25-112103';inp=Path('/data/inputs/TN25-112103')
paths=[inp/'vendor-transformed-expression.tsv',inp/'RNA_TN25-112103.bam',inp/'RNA_TN25-112103.bam.bai',inp/'original/RNA_TN25-112103_S45.R1.fastq.gz',inp/'original/RNA_TN25-112103_S45.R2.fastq.gz',sample/'salmon_aln_quant/TN25_research-TN25-112103-ar-TN25-112103.quant.sf',sample/'starfusion/TN25_research-TN25-112103-ar-TN25-112103.star-fusion.fusion_predictions.abridged.coding_effect.tsv',old/'samples/multisample/msnaf/counts.msnaf.input.txt',ref/'fasta/Homo_sapiens.assembly38.fa',ref/'fasta/Homo_sapiens.assembly38.fa.fai',ref/'annot/gencode.v37.annotation.with.hervs.gtf',ref/'cta_self/cta.pirlygenes.06APR2026.gene_list',ref/'cta_self/canonical.annotated.mtec.norm_tiss.subcell.21APR2026.tsv',ref/'protein/gencode.v37.pc_translations.fa',ref/'pyensembl/GRCh38/ensembl103/Homo_sapiens.GRCh38.pep.all.fa.gz',Path('/data/raft/references/qc/review_scoring_defaults.yaml')]
snaf=ref/'splice_variants/snaf/snaf-data'
paths += [snaf/'Alt91_db'/x for x in ['Hs_Ensembl_exon_add_col.txt','mRNA-ExonIDs.txt','Hs_gene-seq-2000_flank.fa','Homo_sapiens.GRCh38.91.gtf','df_start_codon.txt']]
paths += [snaf/'controls'/x for x in ['GTEx_junction_counts.h5ad','tcga_matched_control_junction_count.h5ad']]
ctat=ref/'fusions/starfusion/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play'
for directory in [ctat,ctat/'ref_genome.fa.star.idx']:
 for base,dirs,files in os.walk(directory,followlinks=False):
  paths += [Path(base)/x for x in files]
paths += list((Path('/data/results/TN25-112103/openvax-mutations-20261018/references/mhcflurry/2.3.0/models_class1_presentation/models')).rglob('*'))
paths += [r/'review_scoring_defaults.original.yaml',r/'review_scoring_defaults.expression-aware.yaml']
cache={};rows=[];missing=[];start=time.time()
for path in sorted(set(paths)):
 if not path.exists():missing.append(str(path));continue
 if not path.is_file():continue
 resolved=str(path.resolve())
 if resolved not in cache:
  h=hashlib.sha256()
  with path.open('rb') as f:
   while chunk:=f.read(8*1024*1024):h.update(chunk)
  cache[resolved]=h.hexdigest()
 row={'path':str(path),'resolved_path':resolved,'bytes':path.stat().st_size,'sha256':cache[resolved]};rows.append(row)
 print(json.dumps({'hashed':len(rows),'bytes':row['bytes'],'path':str(path)}),flush=True)
out=r/'provenance';out.mkdir(exist_ok=True)
summary={'files':rows,'missing_declared_paths':missing,'unique_resolved_files':len(cache),'seconds':time.time()-start,'scope':'RNA inputs; original full discovery checkpoints; CTA/GTF/proteome references; native MSNAF inputs and normal controls; CTAT library including STAR index; MHCflurry weights; review configuration'}
(out/'input-reference-hashes.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({'hashed_paths':len(rows),'missing_paths':missing,'seconds':time.time()-start}),flush=True)
assert not missing,missing
