"""Prepare native junction-level evidence and a compatible exon reference.
MSNAF UIDs use unversioned ENSG IDs. Native LENSTools does not strip GTF gene
versions during competing-junction lookup, so normalize ONLY that identifier
attribute in the copied exon reference. Coordinates and transcript IDs are kept.
Peptide/CDS are empty: this preliminary pass measures junctions, not coding frames.
"""
import argparse,csv,re,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--root',required=True);a=p.parse_args();r=Path(a.root)
with (r/'splice/msnaf_peptides.csv').open() as f:
 events={}
 for x in csv.DictReader(f):
  if x['uid'] in events:assert events[x['uid']]['coord']==x['coord']
  else:events[x['uid']]=x
genes={g for uid in events for g in re.findall(r'ENSG[0-9]+',uid)};names={};n=0
pattern=re.compile(r'gene_id "(ENSG[0-9]+)(?:\.[0-9]+)?"')
with open('/data/raft/references/homo_sapiens/annot/gencode.v37.annotation.with.hervs.gtf') as fi,(r/'splice/review.exons.complete.gtf').open('w') as fo:
 for line in fi:
  if line.startswith('#'):continue
  c=line.split('\t')
  if len(c)!=9 or c[2]!='exon':continue
  m=pattern.search(c[8])
  if m and m[1] in genes:
   gn=re.search(r'gene_name "([^"]+)"',c[8]);names[m[1]]=gn[1] if gn else ''
   fo.write(pattern.sub(lambda m:'gene_id "'+m[1]+'"',line));n+=1
rows=[]
for i,(uid,x) in enumerate(sorted(events.items())):
 rows.append(dict(candidate_id='JUNCTION_'+str(i),source_type='SPLICE',peptide='',hla_allele='',peptide_length=0,gene=names.get(uid.split(':')[0],''),transcript_id='',lens_row_id=i,splice_description=uid,splice_coords=x['coord'],splice_coding_sequence=''))
with (r/'splice/junction-review-candidates.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
summary={'junctions':len(rows),'gene_ids':len(genes),'gene_ids_with_gtf_exons':len(names),'exon_rows':n,'identifier_change':'strip version suffix of gene_id only; preserve original GTF separately','review_scope':'junction evidence/fractions only; no peptide/frame validation in this preliminary table'}
(r/'splice/junction-review-input-summary.complete.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))
