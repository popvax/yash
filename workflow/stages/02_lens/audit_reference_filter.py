"""Independent full-output check of the native reference-proteome filter.
Uses multi-pattern exact matching, rather than repeating native per-peptide scans.
The curated CTA exception is the same stated policy, not an implied novelty claim.
"""
import sys,csv,json,collections
from pathlib import Path
sys.path.insert(0,'/work/lens-expression-aware/verification-deps')
import ahocorasick
r=Path(sys.argv[1]);out=r/'verification';out.mkdir(exist_ok=True)
def rows(path):
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
original=rows(r/'prediction/lens-annotated.tsv');retained=rows(r/'prediction/lens-self-filtered.tsv')
keys=lambda xs:{(x['identity'],x['peptide'],x['allele']) for x in xs}
assert len(keys(original))==len(original) and len(keys(retained))==len(retained)
assert keys(retained)<=keys(original)
ctas={x.strip() for x in Path('/data/raft/references/homo_sapiens/cta_self/cta.pirlygenes.06APR2026.gene_list').read_text().splitlines() if x.strip()}
auto=ahocorasick.Automaton()
for pep in {x['peptide'] for x in retained}:auto.add_word(pep,pep)
auto.make_automaton();noncta_hits={};ref_count=0;noncta_count=0
header=None;seq=[]
def scan(header,seq):
 global ref_count,noncta_count
 ref_count+=1;gene=header.split('|')[-2]
 if gene in ctas:return
 noncta_count+=1
 for _,pep in auto.iter(''.join(seq)):noncta_hits.setdefault(pep,gene)
with open('/data/raft/references/homo_sapiens/protein/gencode.v37.pc_translations.fa') as f:
 for line in f:
  if line.startswith('>'):
   if header is not None:scan(header,seq)
   header=line[1:].strip();seq=[]
  else:seq.append(line.strip())
 if header is not None:scan(header,seq)
summary={'working_input_rows':len(original),'native_retained_rows':len(retained),'native_removed_rows':len(original)-len(retained),'retained_unique_peptides':len(auto),'reference_proteins':ref_count,'non_cta_reference_proteins':noncta_count,'retained_peptides_with_non_cta_reference_match':len(noncta_hits),'independent_method':'pyahocorasick 2.2.0 exact multi-pattern scan of every non-CTA reference protein','cta_exception':'matches to curated CTA proteins remain allowed; reference matching is not universal novelty evidence'}
(out/'reference-filter-audit.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
assert not noncta_hits,list(noncta_hits.items())[:10]
