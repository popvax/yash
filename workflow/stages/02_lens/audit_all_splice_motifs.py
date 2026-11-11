"""Independently verify every reviewed splice motif against indexed genome bases."""
import csv,json,sys,re,collections
from pathlib import Path
import pysam
from Bio.Seq import Seq
r=Path(sys.argv[1]);fa=pysam.FastaFile('/data/raft/references/homo_sapiens/fasta/Homo_sapiens.assembly38.fa');counts=collections.Counter();uids={}
accepted={'GT-AG','GC-AG','AT-AC','CT-AC','CT-GC','GT-AT'}
with (r/'peptide-review/candidates.tsv').open() as f:
 for x in csv.DictReader(f,delimiter='\t'):
  if x['source_type']!='SPLICE':continue
  e=json.loads((r/'peptide-review/evidence'/(x['candidate_id']+'.json')).read_text());j=e['event_specific']['junction'];c,a,b,st=re.match(r'^([^:]+):(\d+)-(\d+)\(([+-])\)$',x['splice_coords']).groups();a,b=int(a),int(b)
  assert (j['intron_start_1based'],j['intron_end_1based'])==(a+1,b-1)
  ends=fa.fetch(c,a,a+2).upper()+fa.fetch(c,b-3,b-1).upper()
  oriented=ends if st=='+' else str(Seq(ends).reverse_complement());expected=oriented[:2]+'-'+oriented[2:]
  assert j['motif']==expected and j['motif_valid']==(expected in accepted),(x['candidate_id'],j,expected)
  category='forward_consensus' if expected in {'GT-AG','GC-AG','AT-AC'} else 'reverse_consensus_allowed_by_native_policy' if expected in accepted else 'non_consensus'
  counts[category]+=1;uids[x['splice_description']]=category
result={'status':'VERIFIED','peptide_reviews_checked':sum(counts.values()),'unique_junctions_checked':len(uids),'peptide_review_categories':dict(counts),'unique_junction_categories':dict(collections.Counter(uids.values())),'method':'independent genome fetch of intron ends; reverse-complement the entire four-base end sequence for minus strand; compare every stored motif and validity flag','native_policy':'accepted six motifs retained, including both orientations; motif matching alone does not establish transcript strand or tumour specificity'}
(r/'verification/all-splice-motifs-audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
