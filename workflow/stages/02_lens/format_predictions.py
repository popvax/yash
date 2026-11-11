"""Join native origins to raw MHCflurry results for native LENS downstream tools.
The 1000 nM working table covers the positive native LENS binding-contribution
range. Every raw prediction remains archived; this is not a biological cutoff.
"""
import csv,json,collections,argparse
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--root',required=True);p.add_argument('--max-affinity',type=float,default=1000);p.add_argument('--affinity-file');a=p.parse_args();r=Path(a.root);d=r/'prediction'
by_pep=collections.defaultdict(list);raw=0
with Path(a.affinity_file or d/'affinity-all-unique.tsv').open() as f:
 reader=csv.DictReader(f,delimiter='\t');hdr=reader.fieldnames
 for x in reader:
  raw+=1
  if float(x['mhcflurry_2.3.14.aff'])<a.max_affinity:by_pep[x['peptide']].append(x)
seen=set();n=0
with (d/'peptide-origins.tsv').open() as f,(d/'mhcflurry.working.tsv').open('w') as out:
 w=csv.DictWriter(out,fieldnames=['identity']+hdr+['pos'],delimiter='\t');w.writeheader()
 for x in csv.DictReader(f,delimiter='\t'):
  key=x['identity'],x['peptide']
  if key in seen:continue
  seen.add(key)
  for pred in by_pep.get(x['peptide'],[]):w.writerow(dict(identity=x['identity'],pos=x['pos'],**pred));n+=1
with (d/'all-sources.pep.fa').open('w') as out:
 for name in ['cta/cta.pep.fa','fusion/fusion.pep.fa','prediction/splice.pep.fa']:out.write((r/name).read_text())
summary={'raw_prediction_rows':raw,'working_max_affinity_exclusive_nm':a.max_affinity,'working_unique_peptides':len(by_pep),'working_origin_hla_rows':n,'repeated_positions':'first position retained in native working table; all occurrences in peptide-origins.tsv'}
(d/'working-table-summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))
