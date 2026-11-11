"""Check geometry with synthetic plus/minus introns and real RNA CIGAR support."""
import sys,json,csv,re,collections
from pathlib import Path
sys.path.insert(0,'/opt/lenstools')
from src import metadata as m
from Bio.Seq import Seq
import pysam
r=Path(sys.argv[1]);work=Path('/work/lens-expression-aware/contract-test');work.mkdir(exist_ok=True)
# Expected minus motif comes from reverse-complementing the ENTIRE intron,
# independently of the helper's individual-end handling.
seqs={'plus':'A'*10+'GT'+'A'*20+'AG'+'A'*10,'minus':'A'*10+'CT'+'A'*20+'AC'+'A'*10,'bad':'A'*10+'CC'+'A'*20+'CC'+'A'*10}
fa_path=work/'motifs.fa';fa_path.write_text(''.join('>'+k+'\n'+v+'\n' for k,v in seqs.items()));pysam.faidx(str(fa_path));fa=pysam.FastaFile(str(fa_path))
checks=[]
for name,strand in [('plus','+'),('minus','-'),('bad','+')]:
 motif,valid=m._junction_motif(fa,name,11,34,strand)
 intron=seqs[name][10:34];transcribed=intron if strand=='+' else str(Seq(intron).reverse_complement())
 expected=transcribed[:2]+'-'+transcribed[-2:];assert motif==expected
 assert valid==(name!='bad');checks.append({'name':name,'motif':motif,'valid':valid})
rows=list(csv.DictReader((r/'peptide-review/candidates.tsv').open(),delimiter='\t'));by_uid={}
for x in rows:
 if x['source_type']=='SPLICE':by_uid.setdefault(x['splice_description'],x)
# Use independently recorded actual-CIGAR examples, including negative strands.
audit=json.loads((r/'verification/splice-coordinate-contract-audit.json').read_text())
bam=pysam.AlignmentFile('/data/inputs/TN25-112103/RNA_TN25-112103.bam','rb');ref=pysam.FastaFile('/data/raft/references/homo_sapiens/fasta/Homo_sapiens.assembly38.fa');gtf=m._load_gtf_junctions(str(r/'splice/review.exons.complete.gtf'))
real=[]
for example in audit:
 x=by_uid[example['uid']];c,a,b,st=re.match(r'^([^:]+):(\d+)-(\d+)\(([+-])\)$',x['splice_coords']).groups();a,b=int(a),int(b)
 sample,event=m._extract_splice_evidence(bam,c,a,b,st,gtf,x['splice_description'],20,200,'fragment',x['splice_coding_sequence'])
 norm=event['normalized_junction'];assert (norm['start'],norm['end'])==(a+1,b-1)
 fragments=set()
 for read in bam.fetch(c,max(0,a-201),b+200):
  if read.is_unmapped or read.is_secondary or read.mapping_quality<20:continue
  if any((j['start'],j['end'])==(a+1,b-1) for j in m._extract_junctions_from_read(read)):fragments.add(read.query_name)
 assert sample['support_reads']==len(fragments),(sample['support_reads'],len(fragments))
 intron=ref.fetch(c,a,b-1).upper();transcribed=intron if st=='+' else str(Seq(intron).reverse_complement());expected=transcribed[:2]+'-'+transcribed[-2:]
 motif,valid=m._junction_motif(ref,c,a+1,b-1,st);assert motif==expected
 real.append({'gene':x['gene'],'support_fragments':sample['support_reads'],'intron_start':a+1,'intron_end':b-1,'motif':motif,'motif_valid':valid})
result={'synthetic_motifs':checks,'actual_RNA_CIGAR_cases':real,'status':'PASSED'};(r/'verification/splice-contract-test.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
