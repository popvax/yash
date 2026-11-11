"""Attach supplied gene expression to native transcript quant rows.
Does not overwrite Salmon TPM or infer isoform-specific abundance. Outputs only
positive, measured CTA-gene rows to the native selection input; the audit keeps
zero/missing genes. Native LENS >=0 threshold would otherwise admit zero values.
"""
import argparse,csv,json,re,hashlib
from pathlib import Path
p=argparse.ArgumentParser()
for x in ('expression','quants','gtf','genes','out'):p.add_argument('--'+x,required=True)
a=p.parse_args(); out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
expr={}
with open(a.expression) as f:
 for line in f:
  if line.startswith('#') or not line.strip():continue
  gene,value=line.rstrip().split('\t')[:2]
  if gene in expr:raise ValueError('Duplicate expression gene '+gene)
  expr[gene]=float(value)
genes={x.split('\t')[0].strip() for x in Path(a.genes).read_text().splitlines() if x.strip()}
txgene={}; attr=re.compile(r'(\w+) "([^"]*)"')
with open(a.gtf) as f:
 for line in f:
  if line.startswith('#'):continue
  c=line.rstrip().split('\t')
  if len(c)!=9 or c[2]!='transcript':continue
  d=dict(attr.findall(c[8]));g=d.get('gene_name');t=d.get('transcript_id')
  if g in genes and d.get('transcript_type',d.get('transcript_biotype'))=='protein_coding':txgene[t]=g
selected=set()
with open(a.quants) as fi,(out/'expression-adapted.quant.sf').open('w') as fo:
 r=csv.DictReader(fi,delimiter='\t');w=csv.DictWriter(fo,fieldnames=r.fieldnames+['SuppliedGeneExpression'],delimiter='\t');w.writeheader()
 for row in r:
  t=row['Name'];g=txgene.get(t);v=expr.get(g)
  if v is not None and v>0:
   row['SuppliedGeneExpression']=v;w.writerow(row);selected.add(t)
with open(a.gtf) as fi,(out/'cta.selected.gtf').open('w') as fo:
 for line in fi:
  if line.startswith('#'):continue
  c=line.rstrip().split('\t')
  if len(c)==9 and dict(attr.findall(c[8])).get('transcript_id') in selected:fo.write(line)
with (out/'cta-gene-selection-audit.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['gene','supplied_expression','status','selected_transcripts'])
 for g in sorted(genes):
  v=expr.get(g);w.writerow([g,'' if v is None else v,'missing' if v is None else 'zero' if v==0 else 'positive',sum(txgene[t]==g for t in selected)])
manifest={'selection':'measured supplied gene expression >0; all matching protein-coding transcript rows eligible; no isoform abundance claim','cta_genes':len(genes),'positive_genes_with_transcripts':len({txgene[t] for t in selected}),'selected_transcripts':len(selected),'input_sha256':{k:hashlib.sha256(Path(getattr(a,k)).read_bytes()).hexdigest() for k in ('expression','quants','gtf','genes')}}
(out/'expression-adapter-summary.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(manifest))
