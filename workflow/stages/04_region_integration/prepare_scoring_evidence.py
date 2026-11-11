"""Reuse cached affinities and make an exact reference index for all queries.

Multi-pattern matching is an execution optimization of exact substring matching;
membership for the original event queries is compared to native Vaxrank results.
"""
import collections,csv,dataclasses,gzip,hashlib,json,re,sys
from pathlib import Path
import ahocorasick
from Bio import SeqIO

R=Path(sys.argv[1]);L=Path('/data/results/TN25-112103/lens-expression-aware-20261030')
queries=set((R/'expanded-target-peptides.txt').read_text().splitlines());found=set();rows=0
with (L/'prediction/affinity-all-unique.tsv').open() as src,(R/'scoring-affinities.tsv').open('w') as out:
    reader=csv.DictReader(src,delimiter='\t');writer=csv.DictWriter(out,fieldnames=reader.fieldnames,delimiter='\t');writer.writeheader()
    for x in reader:
        if x['peptide'] in queries:writer.writerow(x);found.add(x['peptide']);rows+=1
missing=sorted(queries-found)
print('CACHE_COVERAGE',len(queries),'queries;',len(missing),'require additional affinity predictions',flush=True)
if missing:
    import torch,importlib.metadata
    from mhcflurry import Class1PresentationPredictor
    assert importlib.metadata.version('mhcflurry')=='2.3.14'
    torch.set_num_threads(4);torch.set_float32_matmul_precision('highest')
    model=Class1PresentationPredictor.load('/data/results/TN25-112103/openvax-mutations-20261018/references/mhcflurry/2.3.0/models_class1_presentation/models').affinity_predictor
    alleles=['HLA-A*01:01','HLA-A*33:03','HLA-B*08:01','HLA-B*15:08','HLA-C*01:02','HLA-C*07:02']
    extra=[]
    for allele in alleles:
        affinities=model.predict(peptides=missing,allele=allele,model_kwargs={'batch_size':4096})
        ranks=model.percentile_ranks(affinities,allele=allele)
        extra.extend(dict(zip(['peptide','allele','mhcflurry_2.3.14.aff','mhcflurry_2.3.14.aff_perc'],values)) for values in zip(missing,[allele]*len(missing),affinities,ranks))
    with (R/'scoring-affinities.tsv').open('a') as out:
        writer=csv.DictWriter(out,fieldnames=['peptide','allele','mhcflurry_2.3.14.aff','mhcflurry_2.3.14.aff_perc'],delimiter='\t');writer.writerows(extra)
    with (R/'additional-affinities.tsv').open('w') as out:
        writer=csv.DictWriter(out,fieldnames=list(extra[0]),delimiter='\t');writer.writeheader();writer.writerows(extra)
    rows+=len(extra)
assert rows==len(queries)*6
print('PREDICTIONS',len(queries),'missing',len(missing),flush=True)
ctas={s.strip() for s in Path('/data/raft/references/homo_sapiens/cta_self/cta.pirlygenes.06APR2026.gene_list').read_text().splitlines() if s.strip()}
cta_gene_ids=set()
annotation=Path('/data/raft/references/homo_sapiens/cta_self/canonical.annotated.mtec.norm_tiss.subcell.21APR2026.tsv')
with annotation.open() as f:
    for row in csv.DictReader(f,delimiter='\t'):
        if row['gene_name'] in ctas:cta_gene_ids.add(row['gene_id'].split('.')[0])
automaton=ahocorasick.Automaton()
for p in queries:automaton.add_word(p,p)
automaton.make_automaton();all_hits={};noncta_hits={};proteins=0
reference=R/'native-reference-snapshot.jsonl'
sources_count=0
with reference.open() as f:
    for line in f:
        record=json.loads(line);proteins+=1;sources_count+=len(record['sources'])
        hits={p for _,p in automaton.iter(record['sequence'])}
        for source in record['sources']:
            gene=source['gene_id'];name=source['gene_name']
            if name in ctas:cta_gene_ids.add(gene)
            for p in hits:
                all_hits.setdefault(p,source)
                if gene not in cta_gene_ids and name not in ctas:noncta_hits.setdefault(p,source)
# This is membership evidence with a documented witness, not an exhaustive list
# of matching reference sources. Native provenance-complete status is not reused.
result={p:dict(occurs=p in all_hits,occurs_in_non_CTA=p in noncta_hits,
               witness=all_hits.get(p),non_CTA_witness=noncta_hits.get(p)) for p in queries}
native=json.loads((R/'event-peptide-self-reference.json').read_text())
for p,m in native.items():assert result[p]['occurs']==m['occurs'],p
(R/'scoring-reference-membership.json').write_text(json.dumps(result)+'\n')
summary=dict(unique_scoring_peptides=len(queries),prediction_rows=rows,new_prediction_peptides=len(missing),new_prediction_rows=len(missing)*6,
             reference_proteins=proteins,reference_sources=sources_count,exact_reference_peptides=len(all_hits),non_CTA_reference_peptides=len(noncta_hits),
             native_Vaxrank_membership_parity_queries=len(native),native_membership_parity=True,
             reference_path=str(reference),reference_sha256=hashlib.sha256(reference.read_bytes()).hexdigest(),
             reference_matching='Exact substring, standard amino-acid alphabet, no I/L equivalence; one witness per hit; all proteins scanned',
             cta_exception_genes=sorted(ctas),cta_excluded_gene_ids=sorted(cta_gene_ids))
(R/'scoring-evidence-summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({k:v for k,v in summary.items() if k!='cta_exception_genes'},indent=2),flush=True)
