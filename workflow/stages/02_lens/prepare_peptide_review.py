"""Prepare native peptide-level review, sharing only peptide-independent evidence.
CTA/fusion native extractors do not use peptide sequence; reuse their event
metrics with a new candidate ID. Splice evidence is recomputed with each actual
peptide/CDS, including native coding-read and sequence-consistency checks.
HLA does not affect read QC: one review per origin/peptide is shared by its alleles.
"""
import csv,json,argparse,hashlib,re
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--root',required=True);a=p.parse_args();r=Path(a.root)
out=r/'peptide-review';out.mkdir(exist_ok=True);evdir=out/'evidence';evdir.mkdir(exist_ok=True)
genes={};txgenes={}
for name in ['splice/review.exons.complete.gtf','cta/cta.selected.gtf']:
 with (r/name).open() as f:
  for l in f:
   d=dict(re.findall(r'(\w+) "([^"]*)"',l));g=d.get('gene_id','').split('.')[0]
   genes[g]=d.get('gene_name','');txgenes[d.get('transcript_id','')]=d.get('gene_name','')
with (r/'cta/transcript-review-candidates.tsv').open() as f:cta={x['candidate_id']:x for x in csv.DictReader(f,delimiter='\t')}
rows={}
with (r/'prediction/lens-self-filtered.tsv').open() as f:
 for row in csv.DictReader(f,delimiter='\t'):
  if float(row['mhcflurry_2.3.14.aff'])>500:continue
  key=row['identity'],row['peptide']
  if key not in rows or float(row['mhcflurry_2.3.14.aff'])<float(rows[key]['mhcflurry_2.3.14.aff']):rows[key]=row
allrows=[];splice=[];reused=0
for (identity,peptide),row in rows.items():
 source=row['antigen_source'];cid='PEP_'+hashlib.sha256((identity+'|'+peptide).encode()).hexdigest()[:20]
 row.update(candidate_id=cid,source_type=source,peptide_length=len(peptide),hla_allele=row['allele'],lens_row_id=identity)
 if source=='CTA/SELF':
  row.update({k:v for k,v in cta[identity].items() if k.startswith('cta_') or k=='annotation_match'})
  row['gene']=txgenes.get(row.get('transcript_id',''),'')
 elif source=='SPLICE':row['gene']=genes.get(row.get('splice_gene_id','').split('.')[0],'')
 elif source=='FUSION':row['gene']=row.get('fusion_id','')
 else:raise ValueError(source)
 if row.get('transcript_id') in (None,'NA','N/A'):row['transcript_id']=''
 if source=='SPLICE':splice.append(row)
 else:
  original=r/('cta/review-evidence' if source=='CTA/SELF' else 'fusion/event-review-evidence')/(identity+'.json')
  e=json.loads(original.read_text());assert e['source_type']==source;e['candidate_id']=cid
  (evdir/(cid+'.json')).write_text(json.dumps(e,indent=2)+'\n');reused+=1
 allrows.append(row)
assert len({x['candidate_id'] for x in allrows})==len(allrows)
columns=list(dict.fromkeys(k for x in allrows for k in x))
def write(path,rows):
 with path.open('w') as f:
  w=csv.DictWriter(f,fieldnames=columns,delimiter='\t');w.writeheader();w.writerows(rows)
write(out/'candidates.tsv',allrows)
for i in range(4):write(out/f'splice-candidates.{i}.tsv',splice[i::4])
summary={'unique_origin_peptide_reviews':len(allrows),'actual_splice_peptide_reviews':len(splice),'native_peptide_independent_evidence_reused':reused,'affinity_threshold_nm_inclusive':500,'raw_predictions_retained':True,'HLA_QC_sharing':'same origin/peptide RNA evidence shared across all HLA predictions'}
(out/'input-summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))
