"""Export reviewed LENS eligibility and expression evidence for region scoring."""
import argparse,csv,json,math,re,sys,collections,hashlib
from functools import lru_cache
from pathlib import Path
sys.path.insert(0,'/opt/lenstools')
from expression_weights import event_expression_proxies

def read(path):
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(path,rows,columns=None):
 columns=columns or list(dict.fromkeys(k for row in rows for k in row))
 with path.open('w') as f:
  w=csv.DictWriter(f,fieldnames=columns,delimiter='\t');w.writeheader();w.writerows(rows)
def number(x):return None if x in (None,'','NA','N/A','nan') else float(x)
def main():
 p=argparse.ArgumentParser();p.add_argument('--root',required=True);a=p.parse_args();r=Path(a.root);out=r/'final';out.mkdir(exist_ok=True)
 candidates={(x['identity'],x['peptide']):x for x in read(r/'peptide-review/candidates.tsv')}
 scores={x['candidate_id']:x for x in read(r/'peptide-review/scores.tsv')}
 assert len(scores)==len(candidates)
 fusion={x['candidate_id']:x for x in read(r/'fusion/event-review-candidates.tsv')}
 allowed_fusions={x['identity'] for x in read(r/'prediction/native-source-filter.tsv') if x['antigen_source']=='FUSION'}
 expression={}
 for line in Path('/data/inputs/TN25-112103/vendor-transformed-expression.tsv').read_text().splitlines():
  if line and not line.startswith('#'):
   gene,v=line.split('\t')[:2];expression[gene]=float(v)
 gene_names={}
 with (r/'splice/review.exons.complete.gtf').open() as f:
  for line in f:
   at=dict(re.findall(r'(\w+) "([^"]*)"',line));gene_names[at['gene_id']]=at.get('gene_name','')
 @lru_cache(None)
 def evidence(cid):
  e=json.loads((r/'peptide-review/evidence'/(cid+'.json')).read_text());assert e['candidate_id']==cid;return e
 rows=[];fractions={};stage_counts=collections.Counter()
 for row in read(r/'prediction/lens-self-filtered.tsv'):
  affinity=float(row['mhcflurry_2.3.14.aff'])
  if affinity>500:continue
  candidate=candidates[(row['identity'],row['peptide'])];cid=candidate['candidate_id'];score=scores[cid];e=evidence(cid)
  source=row['antigen_source'];notes=[];left_gene=candidate['gene'];right_gene=''
  left_expression=expression.get(left_gene);right_expression=left_expression;left_fraction=right_fraction=1.0
  base_quality='not_applicable';event_id=left_gene
  row.update({k:v for k,v in candidate.items() if k.startswith('cta_') or k=='annotation_match'})
  if source=='SPLICE':
   event_id=row['splice_description'];cc=e['event_specific']['canonical_comparison'];raw_psi=cc['psi_like'];left_fraction=right_fraction=None;base_quality=cc['baseline_definition_quality'];row['native_splice_psi_like']=raw_psi
   if base_quality!='gtf_competing_junctions':
    notes.append('NO_ANNOTATED_COMPETING_JUNCTION_FOUND')
   previous=fractions.setdefault(event_id,raw_psi);assert previous==raw_psi
   gene_ids=re.findall(r'ENSG[0-9]+',event_id)
   if len(set(gene_ids))>1:
    left_gene=gene_names.get(gene_ids[0],'');right_gene=gene_names.get(gene_ids[-1],'')
    left_expression=expression.get(left_gene);right_expression=expression.get(right_gene)
    notes.append('INTERGENE_USES_ASSIGNED_PRIMARY_GENE_ONLY')
   else:right_expression=left_expression
  elif source=='FUSION':
   f=fusion[row['identity']];left_gene=row['fusion_left_gene'].split('^')[0];right_gene=row['fusion_right_gene'].split('^')[0]
   left_expression=number(f.get('left_supplied_expression'));right_expression=number(f.get('right_supplied_expression'))
   left_fraction=number(f.get('left_fraction'));right_fraction=number(f.get('right_fraction'))
   row.update({k:v for k,v in f.items() if k in ['fusion_fragments','left_unfused_fragments','right_unfused_fragments','FAR_left_reported','FAR_right_reported','fusion_inspector_match']})
   event_id='|'.join([row['fusion_id'],row['fusion_left_breakpoint'],row['fusion_right_breakpoint']]);base_quality=f['fusion_inspector_match']
  elif source!='CTA/SELF':raise ValueError(source)
  lp,rp=event_expression_proxies(source,left_expression,right_expression,left_fraction,right_fraction)
  if left_expression is None:notes.append('LEFT_OR_SINGLE_GENE_EXPRESSION_MISSING')
  if right_expression is None:notes.append('RIGHT_GENE_EXPRESSION_MISSING')
  if left_expression==0 or right_expression==0:notes.append('SUPPLIED_GENE_EXPRESSION_ZERO')
  if lp is None or rp is None:notes.append('EXPRESSION_WEIGHT_INCOMPLETE')
  row.update(review_id=cid,source_type=source,gene=candidate['gene'],event_id=event_id,alignment_confidence_score=score['alignment_confidence_score'],final_tier=score['final_tier'],review_flags=score['reasons'],evidence_file='peptide-review/evidence/'+cid+'.json',left_or_single_gene=left_gene,right_gene=right_gene,left_or_single_supplied_expression=left_expression,right_supplied_expression=right_expression,left_or_single_event_fraction=left_fraction,right_event_fraction=right_fraction,fraction_baseline_quality=base_quality,expression_proxy_left=lp,expression_proxy_right=rp,expression_proxy_notes=';'.join(notes),rna_support_count=e['tumor_rna'].get('support_reads'),rna_baseline_count=e['tumor_rna'].get('baseline_reads'),rna_count_scope={'CTA/SELF':'transcript_vs_library','SPLICE':'competing_junctions','FUSION':'local_breakpoint_windows'}[source])
  source_ok=source!='FUSION' or row['identity'] in allowed_fusions
  row.update(native_source_filter_status=('FUSION_ANNOTATION_PASS' if source_ok else 'FUSION_NORMAL_ANNOTATION_EXCLUDED') if source=='FUSION' else 'RNA_ONLY_SPLICE_NO_DNA_LINKAGE_GATE' if source=='SPLICE' else 'NATIVE_CTA_SELF_FILTER',primary_eligible=score['final_tier']=='PASS' and source_ok)
  rows.append(row);stage_counts[(source,score['final_tier'])]+=1
 assert len(rows)==len({(x['identity'],x['peptide'],x['allele']) for x in rows})
 columns=list(dict.fromkeys(k for row in rows for k in row));write(out/'all-reviewed-candidates.tsv',rows,columns)
 write(out/'pass-candidates.tsv',[x for x in rows if x['final_tier']=='PASS'],columns)
 events=collections.defaultdict(list)
 for x in rows:events[(x['antigen_source'],x['event_id'])].append(x)
 event_rows=[]
 for (source,event),xs in sorted(events.items()):
  passed=[x for x in xs if x['primary_eligible']]
  event_rows.append(dict(source=source,event_id=event,gene=xs[0]['gene'],origin_sequences=len({x['identity'] for x in xs}),unique_peptides=len({x['peptide'] for x in xs}),peptide_hla_rows=len(xs),rna_pass_peptides=len({x['peptide'] for x in xs if x['final_tier']=='PASS'}),primary_eligible_peptides=len({x['peptide'] for x in passed}),review_peptides=len({x['peptide'] for x in xs if x['final_tier']=='REVIEW'}),fail_peptides=len({x['peptide'] for x in xs if x['final_tier']=='FAIL'}),primary_eligible_hla_alleles=';'.join(sorted({x['allele'] for x in passed})),best_affinity_nm=min(float(x['mhcflurry_2.3.14.aff']) for x in xs),best_primary_affinity_nm=min((float(x['mhcflurry_2.3.14.aff']) for x in passed),default=None)))
 write(out/'event-summary.tsv',event_rows)
 summary={'reviewed_origin_peptide_hla_rows':len(rows),'unique_origin_peptide_reviews':len(candidates),
  'unique_peptides':len({x['peptide'] for x in rows}),
  'source_tiers':{source:{tier:stage_counts[(source,tier)] for tier in ['PASS','REVIEW','FAIL']} for source in ['CTA/SELF','SPLICE','FUSION']},
  'events':dict(collections.Counter(x['source'] for x in event_rows)),
  'primary_eligible_rows':dict(collections.Counter(x['antigen_source'] for x in rows if x['primary_eligible'])),
  'splice_expression':'assigned primary gene only; native fraction remains QC evidence',
  'fusion_expression':'left partner proxy used for region scoring; both raw partner measurements retained',
  'scoring':'Native Vaxrank region scores are computed downstream; no per-peptide expression ranking here'}
 (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
