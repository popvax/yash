"""Supplementary read-only checks of the recorded 21 selected antigens.

Counts describe alignments, not a replacement for Mutect2 or Isovar. No design
eligibility or sequence is changed. Run against already mounted patient BAMs.
"""
import collections
import datetime
import hashlib
import json
from pathlib import Path
import re
import statistics
import sys

import pysam
from Bio.Seq import Seq

root = Path(sys.argv[1])
request = json.loads((root / 'request.json').read_text())
dna = pysam.AlignmentFile('/data/inputs/TN25-112103/DNA_TN25-112103.bam', 'rb')
rna = pysam.AlignmentFile('/data/inputs/TN25-112103/RNA_TN25-112103.bam', 'rb')

def eligible(r):
    return not (r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail) and r.mapping_quality >= 20 and (not r.has_tag('NH') or r.get_tag('NH') == 1)

def med(values):
    return statistics.median(values) if values else None

def snv(bam, allele):
    chrom, pos, ref, alt = allele
    counts = collections.Counter(); names = collections.defaultdict(set)
    qualities = collections.defaultdict(list); strands = collections.defaultdict(collections.Counter)
    for r in bam.fetch(chrom, pos-1, pos):
        if not eligible(r): continue
        for q, p in r.get_aligned_pairs(matches_only=True):
            if p != pos-1: continue
            bq = r.query_qualities[q] if r.query_qualities is not None else 0
            if bq < 20: continue
            base = r.query_sequence[q]
            counts[base] += 1; names[base].add(r.query_name)
            qualities[base].append(bq); strands[base]['reverse' if r.is_reverse else 'forward'] += 1
    return dict(read_counts=dict(counts), fragment_names_per_base={k:len(v) for k,v in names.items()},
                alt_BQ_median=med(qualities[alt]), alt_strand=dict(strands[alt]))

def gap(bam, chrom, low, high):
    # low and high are 1-based exon boundary coordinates. A CIGAR N/D
    # spanning (low, high-1) in zero-based coordinates removes high-low-1 bases.
    names = collections.defaultdict(set); anchors=[]; tags=collections.defaultdict(collections.Counter)
    ends=0; clipped=0; nreads=0
    for r in bam.fetch(chrom, low-1, high):
        if not eligible(r): continue
        p=r.reference_start
        for i,(op,n) in enumerate(r.cigartuples or []):
            if op in (2,3) and p == low and p+n == high-1:
                key='N' if op==3 else 'D'; names[key].add(r.query_name)
                if op==3:
                    nreads+=1
                    a=r.cigartuples[i-1][1] if i else 0
                    b=r.cigartuples[i+1][1] if i+1<len(r.cigartuples) else 0
                    anchors.append(min(a,b))
                    if min(a,b)<10: ends+=1
                    if any(o==4 for o,l in r.cigartuples):clipped+=1
                    for t in ['XS','ts','NH']:
                        tags[t][str(r.get_tag(t)) if r.has_tag(t) else 'ABSENT']+=1
                    tags['read1_reverse' if r.is_read1 else 'read2_reverse'][str(r.is_reverse)]+=1
            if op in (0,2,3,7,8):p+=n
    return dict(exact_N_reads=nreads, exact_N_fragments=len(names['N']),exact_D_fragments=len(names['D']),
                shorter_anchor_median=med(anchors),N_reads_with_anchor_lt10=ends,
                N_reads_softclipped=clipped,tags={k:dict(v) for k,v in tags.items()})

result={'reviewed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'request_sha256':hashlib.sha256((root/'request.json').read_bytes()).hexdigest(),
        'worker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'pysam_version':pysam.__version__, 'policy':'primary, nonduplicate, QC-pass alignments; MQ >=20; NH=1 when present; SNV/codon BQ >=20',
        'mutation':{},'splice':{},'cta':{}}
result['RNA_header_provenance']={k:rna.header.to_dict().get(k) for k in ['RG','PG','CO']}
for x in request['mutations']:
    result['mutation'][x['antigen_id']]={'DNA':snv(dna,x['allele']),'RNA':snv(rna,x['allele'])}
    print('mutation',x['antigen_id'],flush=True)
for x in request['splice']:
    result['splice'][x['antigen_id']]={'RNA':gap(rna,x['chrom'],x['low'],x['high']),
                                      'DNA':gap(dna,x['chrom'],x['low'],x['high'])}
    print('splice',x['antigen_id'],flush=True)

# A nearby patient sequence difference in the PPIC reconstruction.
h=collections.defaultdict(set)
for r in rna.fetch('chr5',123023934,123023966):
    if not eligible(r):continue
    positions={p:q for q,p in r.get_aligned_pairs(matches_only=True)}
    if not all(p-1 in positions for p in [123023945,123023955]):continue
    if min(r.query_qualities[positions[p-1]] for p in [123023945,123023955])<20:continue
    key=''.join(r.query_sequence[positions[p-1]] for p in [123023945,123023955])
    h[key].add(r.query_name)
result['PPIC_two_site_RNA']={'positions_1based':[123023945,123023955],
                            'fragment_names_by_genomic_haplotype':{k:len(v) for k,v in h.items()}}

want={x['transcript_id']:x for x in request['cta']}; cds=collections.defaultdict(list)
gtf=Path('/data/raft/references/homo_sapiens/annot/gencode.v37.annotation.with.hervs.gtf')
genome=pysam.FastaFile('/data/results/TN25-112103/mutect2-tumour-only-20261029/resources/reference.fa')
for line in gtf.open():
    if line.startswith('#'):continue
    c=line.rstrip().split('\t')
    if len(c)!=9 or c[2]!='CDS':continue
    m=re.search(r'transcript_id "([^"]+)"',c[8])
    if m and m[1] in want:cds[m[1]].append((c[0],int(c[3]),int(c[4]),c[6],c[7]))
for tx,x in want.items():
    parts=cds[tx];assert parts,(tx,'missing CDS')
    chrom=parts[0][0];strand=parts[0][3];parts.sort(key=lambda c:c[1],reverse=strand=='-')
    # All selected transcripts have a complete annotated start codon.
    assert parts[0][4]=='0',(tx,parts[0])
    coords=[]
    for _,a,b,_,_ in parts:coords.extend(range(a-1,b) if strand=='+' else range(b-1,a-2,-1))
    coords=coords[x['window_start']*3:x['window_end']*3];assert len(coords)==75
    reference_nt=''.join(genome.fetch(chrom,a,a+1) for a in coords).upper()
    if strand=='-':reference_nt=str(Seq(reference_nt).complement())
    reference_aa=str(Seq(reference_nt).translate())
    assert reference_aa==x['sequence'],(tx,reference_aa,x['sequence'])
    covered=collections.defaultdict(set);full=collections.defaultdict(set); mismatches=collections.Counter()
    fragment_bases=collections.defaultdict(lambda:collections.defaultdict(set))
    read1_orientation=collections.Counter()
    for r in rna.fetch(chrom,min(coords),max(coords)+1):
        if not eligible(r):continue
        p={p:q for q,p in r.get_aligned_pairs(matches_only=True)}
        if r.is_read1:read1_orientation['reverse' if r.is_reverse else 'forward']+=1
        for a in coords:
            if a in p and r.query_qualities[p[a]]>=20:fragment_bases[r.query_name][a].add(r.query_sequence[p[a]])
        for i in range(25):
            ps=coords[3*i:3*i+3]
            if all(a in p and r.query_qualities[p[a]]>=20 for a in ps):
                nt=''.join(r.query_sequence[p[a]] for a in ps)
                if strand=='-':nt=str(Seq(nt).complement())
                aa=str(Seq(nt).translate())
                if aa==x['sequence'][i]:covered[i].add(r.query_name)
                else:mismatches[(i,aa)]+=1
        if all(a in p and r.query_qualities[p[a]]>=20 for a in coords):
            nt=''.join(r.query_sequence[p[a]] for a in coords)
            if strand=='-':nt=str(Seq(nt).complement())
            full[str(Seq(nt).translate())].add(r.query_name)
    merged_full=collections.Counter();best_match=set()
    best_start=x['sequence'].index(x['best_peptide']);best_positions=coords[3*best_start:3*(best_start+len(x['best_peptide']))]
    for name,bases in fragment_bases.items():
        if all(len(bases[a])==1 for a in coords):
            nt=''.join(next(iter(bases[a])) for a in coords)
            if strand=='-':nt=str(Seq(nt).complement())
            merged_full[str(Seq(nt).translate())]+=1
        if all(len(bases[a])==1 for a in best_positions):
            nt=''.join(next(iter(bases[a])) for a in best_positions)
            if strand=='-':nt=str(Seq(nt).complement())
            if str(Seq(nt).translate())==x['best_peptide']:best_match.add(name)
    result['cta'][x['antigen_id']]={'transcript_id':tx,'chrom':chrom,'strand':strand,
        'genomic_reference_translation_matches':True,
        'protein_window_start_0based':x['window_start'],'selected_sequence':x['sequence'],
        'genomic_codon_positions_0based':coords,
        'reference_matching_fragments_per_codon':[len(covered[i]) for i in range(25)],
        'full_25aa_matching_fragments':len(full[x['sequence']]),
        'merged_pair_25aa_matching_fragments':merged_full[x['sequence']],
        'merged_pair_window_observed_sequences':dict(merged_full),
        'best_predicted_peptide':x['best_peptide'],'best_peptide_matching_fragments':len(best_match),
        'read1_orientation_over_window':dict(read1_orientation),
        'full_window_observed_sequences':{k:len(v) for k,v in full.items()},
        'nonreference_codon_observations':[dict(codon_index_0based=i,aa=a,reads=n) for (i,a),n in mismatches.items()]}
    print('cta',x['antigen_id'],'full-window fragments',len(full[x['sequence']]),flush=True)

for x in request['splice']:
    chrom,low,high,strand=x['chrom'],x['low'],x['high'],x['strand']
    left=x['left_nt'];right=75-left
    if strand=='+':coords=list(range(low-left,low))+list(range(high-1,high-1+right))
    else:coords=list(range(high+left-2,high-2,-1))+list(range(low-1,low-1-right,-1))
    nt=''.join(genome.fetch(chrom,a,a+1) for a in coords).upper()
    if strand=='-':nt=str(Seq(nt).complement())
    aa=str(Seq(nt).translate());assert aa==x['sequence'],(x['antigen_id'],aa,x['sequence'])
    candidates=list(rna.fetch(chrom,min(coords),max(coords)+1));names=set()
    for r in candidates:
        if not eligible(r):continue
        p=r.reference_start
        for op,n in r.cigartuples or []:
            if op==3 and p==low and p+n==high-1:names.add(r.query_name)
            if op in (0,2,3,7,8):p+=n
    bases=collections.defaultdict(lambda:collections.defaultdict(set))
    for r in candidates:
        if not eligible(r) or r.query_name not in names:continue
        p={p:q for q,p in r.get_aligned_pairs(matches_only=True)}
        for a in coords:
            if a in p and r.query_qualities[p[a]]>=20:bases[r.query_name][a].add(r.query_sequence[p[a]])
    covered=collections.Counter();full=collections.Counter();best=set()
    start=x['sequence'].index(x['best_peptide']);best_coords=coords[start*3:(start+len(x['best_peptide']))*3]
    def translation(b,positions):
        if not all(len(b[a])==1 for a in positions):return None
        nt=''.join(next(iter(b[a])) for a in positions)
        if strand=='-':nt=str(Seq(nt).complement())
        return str(Seq(nt).translate())
    for name,b in bases.items():
        full_aa=translation(b,coords)
        if full_aa:full[full_aa]+=1
        if translation(b,best_coords)==x['best_peptide']:best.add(name)
        for i in range(25):
            if translation(b,coords[i*3:i*3+3])==x['sequence'][i]:covered[i]+=1
    result['splice'][x['antigen_id']]['sequence_check']={
        'reference_join_translation_matches':True,'junction_fragment_names':len(names),
        'matching_fragments_per_codon':[covered[i] for i in range(25)],
        'full_25aa_matching_fragments':full[x['sequence']],
        'best_peptide':x['best_peptide'],'best_peptide_matching_fragments':len(best),
        'full_window_observed_sequences':dict(full)}
    print('splice sequence',x['antigen_id'],'best peptide fragments',len(best),flush=True)

result['status']='COMPLETED'
(root/'raw-read-review.json').write_text(json.dumps(result,indent=2)+'\n')
print('COMPLETED',flush=True)
