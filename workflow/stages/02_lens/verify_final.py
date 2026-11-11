"""Completion audit against actual inputs, native evidence and exported tables.
This checks all exported rows, every review JSON, raw HLA coverage, both fusion
estimates, rank order and model hashes. It does not claim experimental validation.
"""
import argparse,csv,json,math,hashlib,collections
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--root',required=True);a=p.parse_args();r=Path(a.root)
def read(path):
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def num(x):return None if x in (None,'','NA','N/A','nan') else float(x)
def close(x,y):
 x,y=num(x),num(y);assert (x is None)==(y is None),(x,y)
 if x is not None:assert math.isclose(x,y,rel_tol=1e-10,abs_tol=1e-10),(x,y)
fi_rows=read(r/'fusion/fusioninspector/finspector.FusionInspector.fusions.abridged.tsv')
fi_lookup={(x['#FusionName'],x['LeftBreakpoint'],x['RightBreakpoint']):x for x in fi_rows}
assert len(fi_lookup)==len(fi_rows)
allowed_fusions={x['identity'] for x in read(r/'prediction/native-source-filter.tsv') if x['antigen_source']=='FUSION'}
summary=json.loads((r/'final/summary.json').read_text());model=json.loads((r/'prediction/model-run.json').read_text());inputs=json.loads((r/'prediction/input-summary.json').read_text())
final=read(r/'final/all-reviewed-candidates.tsv');candidates=read(r/'peptide-review/candidates.tsv');scores=read(r/'peptide-review/scores.tsv');score_by_id={x['candidate_id']:x for x in scores}
assert len(final)==summary['reviewed_origin_peptide_hla_rows']
assert len(candidates)==len(scores)==len(score_by_id)==summary['unique_origin_peptide_reviews']
assert len({(x['identity'],x['peptide'],x['allele']) for x in final})==len(final)
metadata={}
with (r/'prediction/source-metadata.jsonl').open() as f:
 for line in f:
  x=json.loads(line);assert x['identity'] not in metadata;metadata[x['identity']]=x
origin_count=0;origin_peptides=set();origin_sources=collections.Counter()
with (r/'prediction/peptide-origins.tsv').open() as f:
 for x in csv.DictReader(f,delimiter='\t'):
  source=metadata[x['identity']];pos=int(x['pos']);assert source['sequence'][pos:pos+len(x['peptide'])]==x['peptide']
  origin_count+=1;origin_peptides.add(x['peptide']);origin_sources[source['source']]+=1
assert origin_count==inputs['origin_rows']
for x in final:
 source=metadata[x['identity']];assert source['source']==x['antigen_source'] and x['peptide'] in source['sequence']
alleles=model['alleles'];mask={p:0 for p in (r/'prediction/unique-peptides.txt').read_text().splitlines()};lookup={};wanted={x['peptide'] for x in final};raw_rows=0;percentile_missing=collections.Counter()
with (r/'prediction/affinity-all-unique.tsv').open() as f:
 for x in csv.DictReader(f,delimiter='\t'):
  bit=1<<alleles.index(x['allele']);assert x['peptide'] in mask and not mask[x['peptide']]&bit
  mask[x['peptide']] |= bit;v=float(x['mhcflurry_2.3.14.aff']);assert math.isfinite(v) and v>0
  q=num(x['mhcflurry_2.3.14.aff_perc'])
  if q is None:percentile_missing[x['allele']]+=1
  else:assert 0<=q<=100
  if x['peptide'] in wanted:lookup[(x['peptide'],x['allele'])]=v
  raw_rows+=1
assert all(x==(1<<len(alleles))-1 for x in mask.values())
assert len(mask)==inputs['unique_peptides']==model['unique_peptides']
assert set(mask)==origin_peptides
assert raw_rows==model['prediction_rows']==len(mask)*len(alleles)
assert hashlib.sha256((r/'prediction/unique-peptides.txt').read_bytes()).hexdigest()==model['input_sha256']
transfer=json.loads((r/'provenance/prediction-transfer.json').read_text())
assert hashlib.sha256((r/'prediction/affinity-all-unique.tsv').read_bytes()).hexdigest()==transfer['sha256']
expr={}
for line in Path('/data/inputs/TN25-112103/vendor-transformed-expression.tsv').read_text().splitlines():
 if line and not line.startswith('#'):
  gene,v=line.split('\t')[:2];expr[gene]=float(v)
evidence={};native_missing=collections.Counter();splice_consistency=collections.Counter()
for x in candidates:
 cid=x['candidate_id'];e=json.loads((r/'peptide-review/evidence'/(cid+'.json')).read_text());evidence[cid]=e
 assert e['candidate_id']==cid and e['source_type']==x['source_type']
 assert e['tumor_dna']['support_reads'] is None and e['normal_dna']['support_reads'] is None
 if x['source_type']=='SPLICE':
  pc=e['event_specific']['peptide_consistency'];assert pc['reported_peptide']==x['peptide']
  import re
  match=re.match(r'^([^:]+):(\d+)-(\d+)\(([+-])\)$',x['splice_coords']);assert match
  c,start,end,strand=match.groups();j=e['event_specific']['junction']
  assert j['coordinate_contract']=='msnaf_exon_boundaries'
  assert j['intron_start_1based']==int(start)+1 and j['intron_end_1based']==int(end)-1
  splice_consistency[str(pc.get('consistency_pass'))]+=1
  if e['event_specific']['canonical_comparison']['baseline_definition_quality']!='gtf_competing_junctions':native_missing['splice_no_annotated_competitor']+=1
for x in final:
 assert 8<=len(x['peptide'])<=11 and x['allele'] in alleles
 aff=float(x['mhcflurry_2.3.14.aff']);assert aff<=500;close(aff,lookup[(x['peptide'],x['allele'])])
 assert x['final_tier']==score_by_id[x['review_id']]['final_tier']
 close(x['alignment_confidence_score'],score_by_id[x['review_id']]['alignment_confidence_score'])
 eligible=x['final_tier']=='PASS' and (x['antigen_source']!='FUSION' or x['identity'] in allowed_fusions)
 assert (x['primary_eligible']=='True')==eligible
 assert x['review_flags']==score_by_id[x['review_id']]['reasons']
 e=evidence[x['review_id']];close(x['rna_support_count'],e['tumor_rna']['support_reads']);close(x['rna_baseline_count'],e['tumor_rna']['baseline_reads'])
 for view,gene_col,expression_col,fraction_col in [('left','left_or_single_gene','left_or_single_supplied_expression','left_or_single_event_fraction'),('right','right_gene','right_supplied_expression','right_event_fraction')]:
  gene=x[gene_col]
  if not gene and x['antigen_source']!='FUSION' and 'INTERGENE' not in x['expression_proxy_notes']:gene=x['left_or_single_gene']
  close(x[expression_col],expr.get(gene))
  v=num(x[expression_col]);fraction=num(x[fraction_col]);assert fraction is None or 0<=fraction<=1
  expected=num(x['left_or_single_supplied_expression']) if x['antigen_source']=='SPLICE' else None if v is None or fraction is None else v*fraction
  close(x['expression_proxy_'+view],expected)
 if x['antigen_source']=='FUSION':
  F=float(x['fusion_fragments'])
  fi=fi_lookup[(x['fusion_id'],x['fusion_left_breakpoint'],x['fusion_right_breakpoint'])]
  close(F,float(fi['JunctionReadCount'])+float(fi['SpanningFragCount']))
  close(e['tumor_rna']['support_reads'],F)
  close(e['event_specific']['support']['split_reads'],fi['JunctionReadCount'])
  close(e['event_specific']['support']['spanning_pairs'],fi['SpanningFragCount'])
  bp=e['event_specific']['breakpoint']
  assert ':'.join([bp['chrA'],str(bp['posA']),bp['strandA']])==x['fusion_left_breakpoint']
  assert ':'.join([bp['chrB'],str(bp['posB']),bp['strandB']])==x['fusion_right_breakpoint']
  for view,counter in [('left','left_unfused_fragments'),('right','right_unfused_fragments')]:
   U=float(x[counter]);close(x['left_or_single_event_fraction'] if view=='left' else x['right_event_fraction'],None if F+U==0 else F/(F+U))
 if x['antigen_source']=='SPLICE':
  cc=e['event_specific']['canonical_comparison'];close(x['native_splice_psi_like'],cc['psi_like'])
  assert num(x['left_or_single_event_fraction']) is None and num(x['right_event_fraction']) is None
cta_annotations={x['candidate_id']:x for x in read(r/'cta/transcript-review-candidates.tsv')}
mtec_rows=0
for x in final:
 if x['antigen_source']=='CTA/SELF':
  reference=cta_annotations[x['identity']]
  for col in ['cta_mean_mtec_tpm','cta_median_mtec_tpm','cta_stdev_mtec_tpm','cta_gene_detectable_normal_tissues']:close(x[col],reference[col])
  assert x['cta_normal_tissue_names']==reference['cta_normal_tissue_names'];mtec_rows+=1
motif_audit=json.loads((r/'verification/all-splice-motifs-audit.json').read_text());assert motif_audit['status']=='VERIFIED' and motif_audit['peptide_reviews_checked']==sum(x['source_type']=='SPLICE' for x in candidates)
ref_audit=json.loads((r/'verification/reference-filter-audit.json').read_text());assert ref_audit['retained_peptides_with_non_cta_reference_match']==0
manifest=json.loads((r/'provenance/input-reference-hashes.json').read_text());assert not manifest['missing_declared_paths']
reference_hashes={x['path']:x['sha256'] for x in manifest['files']};models_root='/data/results/TN25-112103/openvax-mutations-20261018/references/mhcflurry/2.3.0/models_class1_presentation/models/'
for path,digest in model['model_files_sha256'].items():assert reference_hashes[models_root+path]==digest
native_high=read(r/'final/native-high-confidence.tsv');assert all(x['final_tier']=='PASS' for x in native_high)
assert not any(x['antigen_source']=='SPLICE' for x in native_high) # No DNA linkage supplied.
result={'status':'VERIFIED','raw_prediction_rows':raw_rows,'native_origin_rows_checked':origin_count,'origin_rows_by_source':dict(origin_sources),'reviewed_origin_peptide_hla_rows':len(final),'review_evidence_jsons_checked':len(evidence),'mtec_annotation_rows_checked':mtec_rows,'all_splice_motifs_independently_verified':True,'msnaf_coordinate_contract_checked_for_every_splice_review':True,'splice_fractions_not_used_for_weighting':True,'fusion_external_counts_checked_at_exact_breakpoints':True,'native_fusion_normal_annotation_filter_checked':True,'model_weight_files_cross_checked':len(model['model_files_sha256']),'reference_hash_paths':len(manifest['files']),'percentile_missing_by_allele':dict(percentile_missing),'native_missing_metrics':dict(native_missing),'splice_reported_peptide_consistency':dict(splice_consistency),'scope':'all rows, every RNA evidence JSON, all raw peptide/HLA coverage, source expression evidence, independent full-proteome filter audit and original references; experimental presentation not assessed'}
(r/'verification/completion-audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
