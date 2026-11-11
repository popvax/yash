"""Format native pre-HLA sequence outputs for one common affinity screen.
Native splice peptides are already enumerated at 8-11 aa. Native CTA/fusion
FASTA sequences are scanned at those lengths, retaining every source identity.
This adapter changes neither translated sequences nor frame assignments.
"""
import argparse,csv,hashlib,json
from pathlib import Path

def fasta(path):
    identity=None;seq=[];meta={}
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            if identity is not None:yield identity,''.join(seq),meta
            identity=line[1:].split()[0];seq=[];meta={}
        elif line.startswith(';'):
            meta.update(x.split(':',1) for x in line[1:].split() if ':' in x)
        elif line.strip():seq.append(line.strip())
    if identity is not None:yield identity,''.join(seq),meta

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',required=True);a=p.parse_args()
    root=Path(a.root);out=root/'prediction';out.mkdir(exist_ok=True)
    origins=[];peptides=set();counts={};metadata=[]
    allowed=set('ACDEFGHIKLMNPQRSTVWY')
    for source,path in [('CTA/SELF',root/'cta/cta.pep.fa'),('FUSION',root/'fusion/fusion.pep.fa')]:
        n=0
        for identity,seq,meta in fasta(path):
            metadata.append({'identity':identity,'source':source,'sequence':seq,**meta});n+=1
            for k in range(8,12):
                for i in range(len(seq)-k+1):
                    pep=seq[i:i+k]
                    if set(pep)<=allowed:origins.append((identity,pep,i));peptides.add(pep)
        counts[source]=n
    n=0
    with (root/'splice/msnaf_peptides.csv').open() as f,(out/'splice.pep.fa').open('w') as fa:
        for row in csv.DictReader(f):
            pep=row['peptide'];assert len(pep) in (8,9,10,11)
            identity=hashlib.md5(json.dumps(row,sort_keys=True).encode()).hexdigest()[:15]
            meta={'antigen_source':'SPLICE','splice_description':row['uid'],'splice_gene_id':row['uid'].split(':')[0],'splice_coords':row['coord'],'splice_coding_sequence':row['coding_sequence'],'pep_context':row['peptide_context']}
            metadata.append({'identity':identity,'source':'SPLICE','sequence':pep,**meta})
            fa.write('>'+identity+'\n;'+ ' '.join(k+':'+v for k,v in meta.items())+'\n'+pep+'\n')
            if set(pep)<=allowed:origins.append((identity,pep,0));peptides.add(pep)
            n+=1
    counts['SPLICE']=n
    with (out/'peptide-origins.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['identity','peptide','pos']);w.writerows(origins)
    (out/'unique-peptides.txt').write_text('\n'.join(sorted(peptides))+'\n')
    with (out/'source-metadata.jsonl').open('w') as f:
        for row in metadata:f.write(json.dumps(row)+'\n')
    counts.update(unique_peptides=len(peptides),origin_rows=len(origins))
    (out/'input-summary.json').write_text(json.dumps(counts,indent=2)+'\n');print(json.dumps(counts))
if __name__=='__main__':main()
