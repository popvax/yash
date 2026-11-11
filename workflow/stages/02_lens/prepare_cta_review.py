"""Prepare one native transcript-level CTA review row per translated isoform.
Native CTA read evidence is transcript-level, independent of peptide/HLA.
Preserve mTEC means as means; do not relabel them as 95th percentiles.
"""
import csv,json,re,argparse
from pathlib import Path
from prepare_pmhc_inputs import fasta
p=argparse.ArgumentParser();p.add_argument('--root',required=True);a=p.parse_args();r=Path(a.root)
txgene={}
for line in (r/'cta/cta.selected.gtf').read_text().splitlines():
 d=dict(re.findall(r'(\w+) "([^"]*)"',line));txgene[d['transcript_id']]=d.get('gene_name','')
with (r/'annotations/cta_annotation_ref.collapsed.tsv').open() as f:
 annotations={x['transcript']:x for x in csv.DictReader(f,delimiter='\t')}
cols=['candidate_id','source_type','peptide','hla_allele','peptide_length','gene','transcript_id','lens_row_id','cta_gene_detectable_normal_tissues','cta_normal_tissue_names','cta_mean_mtec_tpm','cta_median_mtec_tpm','cta_stdev_mtec_tpm','annotation_match']
rows=[]
for identity,seq,m in fasta(r/'cta/cta.pep.fa'):
 tx=m['transcript_id'];ann=annotations.get(tx,{})
 names=ann.get('gene_detectable_normal_tissues','')
 missing=names in ('','NA','N/A','None','nan')
 rows.append(dict(candidate_id=identity,source_type='CTA/SELF',peptide='',hla_allele='',peptide_length=0,gene=txgene[tx],transcript_id=tx,lens_row_id=identity,cta_gene_detectable_normal_tissues='' if missing else len({x.strip() for x in names.split(',') if x.strip()}),cta_normal_tissue_names=names,cta_mean_mtec_tpm=ann.get('mean_mtec_tpm',''),cta_median_mtec_tpm=ann.get('median_mtec_tpm',''),cta_stdev_mtec_tpm=ann.get('stdev_mtec_tpm',''),annotation_match='exact_transcript' if ann else 'missing'))
with (r/'cta/transcript-review-candidates.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows(rows)
print(json.dumps({'translated_isoforms':len(rows),'annotated_isoforms':sum(x['annotation_match']!='missing' for x in rows),'review_scope':'transcript-level native CTA metrics; empty peptide denotes transcript-only review; peptide-specific flags must be derived later'}))
