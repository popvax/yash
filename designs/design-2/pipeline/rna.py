"""Fresh Isovar reconstruction from paired DNA evidence, with explicit normal-evidence tiers."""
import os
os.environ['OMP_NUM_THREADS']='4';os.environ['MKL_NUM_THREADS']='4'
from pathlib import Path
import collections,csv,datetime,hashlib,importlib.metadata,inspect,json,sys,time,traceback
import pysam,varcode,isovar
from pyensembl import Genome
from isovar import ReadCollector,ProteinSequenceCreator
from isovar.isovar_result import fragment_ids
R=Path(sys.argv[1]);D=R/'rna';D.mkdir(parents=True,exist_ok=True)
inputs=json.loads((R/'annotate/rna-input.json').read_text());cfg=json.loads((R/'analysis-genome-config.json').read_text());g=Genome(**cfg)
assert json.loads((R/'annotate/complete.json').read_text())['status']=='COMPLETED'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def save(name,v):
 p=D/name;t=p.with_suffix('.tmp');t.write_text(json.dumps(v,indent=2,default=str)+'\n');t.replace(p)
def table(name,rows,empty=('variant_id','status')):
 fields=list(dict.fromkeys(k for r in rows for k in r)) or list(empty)
 with (D/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
  for r in rows:w.writerow({k:json.dumps(v) if isinstance(v,(tuple,list,dict)) else v for k,v in r.items()})
def status(phase,**kw):save('status.json',dict(phase=phase,updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),**kw));print(phase,kw,flush=True)
start_time=time.time();fa=pysam.FastaFile('/data/results/TN25-112103/mutect2-tumour-only-20261029/resources/reference.fa')
bam_path='/data/inputs/TN25-112103/RNA_TN25-112103.bam';bam=pysam.AlignmentFile(bam_path,'rb');bam.check_index()
collector=ReadCollector(min_mapping_quality=20,use_secondary_alignments=False,use_duplicate_reads=False,
 read_filter=lambda r:not r.is_supplementary and (not r.has_tag('NH') or r.get_tag('NH')==1))
creator=ProteinSequenceCreator(max_protein_sequences_per_variant=0)
save('settings.json',dict(versions={n:importlib.metadata.version(n) for n in ['isovar','varcode','pyensembl','pysam']},
 read_policy=dict(min_mapping_quality=20,secondary=False,supplementary=False,duplicates=False,NH_when_present=1,fragment_unit='Native Isovar fragment-ID union; no UMI molecule claim'),
 native_fragment_ID_implementation=inspect.getsource(fragment_ids),
 protein_creator_parameters=vars(creator),genome=cfg,scope='Noncommon paired-native-PASS coding calls; normal-present alleles excluded, unresolved-normal cases retained separately from primary evidence',
 RNA_BAM=bam_path,normal_DNA='/data/results/TN25-112103/matched-normal-20261107/bqsr/apply/normal.recal.bam',gene_expression_from_read_counts=False))
variants=[];audit=[]
for x in inputs:
 c,pos,rf,al=x['normalized_allele'];assert fa.fetch(c,pos-1,pos-1+len(rf)).upper()==rf
 v=varcode.Variant(c,pos,rf,al,genome=g,convert_ucsc_contig_names=True);start,ref,alt=v.start,v.ref,v.alt
 if not ref and alt:
  while fa.fetch(c,start,start+1).upper()==alt[0]:
   alt=alt[1:]+alt[0];start+=1;assert start-v.start<10000
 elif ref and not alt:
  while fa.fetch(c,start+len(ref)-1,start+len(ref)).upper()==ref[0]:
   ref=ref[1:]+ref[0];start+=1;assert start-v.start<10000
 if not ref:rp=start;rr=fa.fetch(c,rp-1,rp).upper();ra=rr+alt
 elif not alt:rp=start-1;anchor=fa.fetch(c,rp-1,rp).upper();rr=anchor+ref;ra=anchor
 else:rp=start;rr=ref;ra=alt
 lo=max(0,min(pos,rp)-26);hi=max(pos+len(rf),rp+len(rr))+25;window=fa.fetch(c,lo,hi).upper()
 def apply(p,r,a):
  off=p-1-lo;assert window[off:off+len(r)]==r;return window[:off]+a+window[off+len(r):]
 assert apply(pos,rf,al)==apply(rp,rr,ra)
 rv=varcode.Variant(c,rp,rr,ra,genome=g,convert_ucsc_contig_names=True);variants.append((x,rv))
 audit.append(dict(variant_id=x['variant_id'],normalized_DNA_allele=x['normalized_allele'],RNA_representation=[c,rp,rr,ra],repeat_shift=start-v.start,exact_local_haplotype_equivalence=True))
save('indel-representation-audit.json',audit)
evidence=[];proteins=[];errors=[];accounting=[]
save('rna-evidence.json',evidence);save('protein-hypotheses.json',proteins);save('errors.json',errors)
for i,(x,v) in enumerate(variants,1):
 vid=x['variant_id'];status('reconstructing',completed=i-1,total=len(variants),variant_id=vid,gene=x['gene'])
 try:
  result=list(isovar.run_isovar([v],bam,read_collector=collector,protein_sequence_creator=creator))[0]
  raw=collector.read_evidence_for_variant(v,bam);a=set(fragment_ids(raw.alt_reads));f=set(fragment_ids(raw.ref_reads));o=set(fragment_ids(raw.other_reads));union=a|f|o
  assert (len(a),len(f),len(o),len(union))==(result.num_alt_fragments,result.num_ref_fragments,result.num_other_fragments,result.num_total_fragments)
  frac=len(a)/len(union) if union else None
  if union:assert abs(frac-result.fraction_alt_fragments)<1e-10
  rec=dict(variant_id=vid,gene=x['gene'],protein_change=x['protein_change'],DNA_normalized_allele=x['normalized_allele'],
   DNA_evidence_tier=x['evidence_tier'],primary_DNA_eligible=x['primary_DNA_eligible'],somatic_status=x['somatic_status'],
   rna_ref_fragments=len(f),rna_alt_fragments=len(a),rna_other_fragments=len(o),rna_total_fragments=len(union),rna_alt_fraction=frac,
   rna_alt_only_fragments=len(a-f-o),passes_isovar_filters=bool(result.passes_all_filters),filter_values=dict(result.filter_values),protein_hypotheses=len(result.sorted_protein_sequences))
  if not union:rec['outcome']='NO_ELIGIBLE_LOCAL_RNA_COVERAGE'
  elif not a:rec['outcome']='NO_ALT_RNA_WITH_LOCAL_COVERAGE'
  elif not result.sorted_protein_sequences:rec['outcome']='ALT_RNA_WITHOUT_RECONSTRUCTED_PROTEIN'
  elif not result.passes_all_filters:rec['outcome']='RECONSTRUCTION_REQUIRES_ISOVAR_FILTER_REVIEW'
  else:rec['outcome']='RNA_SUPPORTED_PROTEIN'
  evidence.append(rec);accounting.append(dict(variant_id=vid,alt=len(a),ref=len(f),other=len(o),union=len(union),ref_alt_overlap=len(a&f),alt_other_overlap=len(a&o)))
  for j,p in enumerate(result.sorted_protein_sequences,1):
   ts=sorted(p.transcript_ids);genes=sorted({g.transcript_by_id(t).gene_name for t in ts})
   gene=genes[0] if len(genes)==1 else x['gene']
   usable=bool(result.passes_all_filters and p.contains_mutation and p.amino_acids and 0<=p.mutation_start_idx<p.mutation_end_idx<=len(p.amino_acids))
   proteins.append(dict(protein_id=vid+f'_protein_{j:04d}',variant_id=vid,gene=gene,RNA_transcript_gene_names=genes,gene_assignment_unambiguous=len(genes)==1,
    DNA_annotation_gene=x['gene'],DNA_annotation_protein_change=x['protein_change'],DNA_annotation_transcript_id=x['transcript_id'],hypothesis_rank=j,sequence=p.amino_acids,
    mutation_start_0based=p.mutation_start_idx,mutation_end_0based_exclusive=p.mutation_end_idx,frameshift=p.frameshift,
    contains_mutation=p.contains_mutation,ends_with_stop_codon=p.ends_with_stop_codon,supporting_reads=p.num_supporting_reads,supporting_fragments=p.num_supporting_fragments,
    transcript_ids=ts,variant_passes_isovar_filters=bool(result.passes_all_filters),usable_mutant_context=usable,
    primary_DNA_eligible=x['primary_DNA_eligible'],DNA_evidence_tier=x['evidence_tier'],somatic_status=x['somatic_status']))
 except Exception as e:errors.append(dict(variant_id=vid,error=traceback.format_exc()));print(errors[-1]['error'],flush=True)
 save('rna-evidence.json',evidence);save('protein-hypotheses.json',proteins);save('errors.json',errors)
assert not errors and len(evidence)==len(inputs),(len(errors),len(evidence),len(inputs))
primary=[p for p in proteins if p['hypothesis_rank']==1 and p['usable_mutant_context']]
save('primary-proteins.json',primary);save('fragment-accounting.json',accounting)
table('rna-evidence.tsv',evidence);table('protein-hypotheses.tsv',proteins);table('primary-proteins.tsv',primary)
(D/'primary-proteins.fasta').write_text(''.join('>'+p['protein_id']+' '+p['gene']+'\n'+p['sequence']+'\n' for p in primary))
summary=dict(status='COMPLETED',input_variants=len(inputs),all_inputs_accounted_for=True,outcomes=dict(collections.Counter(r['outcome'] for r in evidence)),
 usable_primary_proteins=len(primary),primary_DNA_and_RNA_eligible=sum(p['primary_DNA_eligible'] and p['gene_assignment_unambiguous'] for p in primary),
 protein_hypotheses=len(proteins),all_repeat_shifts_haplotype_equivalent=True,errors=len(errors),elapsed_seconds=time.time()-start_time,
 input_sha256=dict(variants=sha(R/'annotate/rna-input.json'),genome_config=sha(R/'analysis-genome-config.json')),
 outputs={p.name:dict(bytes=p.stat().st_size,sha256=sha(p)) for p in [D/'rna-evidence.json',D/'rna-evidence.tsv',D/'primary-proteins.json',D/'primary-proteins.tsv',D/'primary-proteins.fasta',D/'indel-representation-audit.json']},completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
save('complete.json',summary);status('completed',**summary)
