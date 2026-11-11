"""Build a bounded, genotype-compatible normal coding-sequence search library.

Unphased local combinations are possibilities, not measured haplotypes. A hit
is conservatively removed from altered-target scoring. A negative search is
limited to the modelled coding sources; it is not proof of normal absence.
"""
from pathlib import Path
import bisect
import collections
import itertools
import json
import re
import sys

from Bio.Seq import Seq
from pyensembl import Genome
from common import save,receipt,sha

MAX_COMBINATIONS=4096
AA=set('ACDEFGHIKLMNPQRSTVWY')


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGTN','TGCAN'))[::-1]


def project(site,transcript,cds):
    chrom,pos,ref=site['chrom'],site['pos'],site['ref']
    offsets=[transcript.spliced_offset(p)-transcript.first_start_codon_spliced_offset
             for p in range(pos,pos+len(ref))]
    lo,hi=min(offsets),max(offsets)+1
    if lo<0 or hi>len(cds) or hi-lo!=len(offsets):
        raise ValueError('Reference allele does not lie wholly in contiguous coding sequence')
    reverse=transcript.strand=='-'
    expected=reverse_complement(ref) if reverse else ref
    if cds[lo:hi]!=expected:raise ValueError('CDS/reference allele mismatch')
    alts={index:reverse_complement(alt) if reverse else alt for index,alt in site['choices'].items()}
    return dict(site,lo=lo,hi=hi,reference_cds_allele=expected,choices=alts)


def phase_compatible(sites,choices):
    # A chosen local chromosome must have a consistent index within every
    # measured phase block. Different phase blocks remain unlinked.
    possibilities={}
    for site,choice in zip(sites,choices):
        if site['phase_block'] is None:continue
        block=site['phase_block'];gt=site['phased_GT']
        allowed={i for i,value in enumerate(gt) if value==choice}
        possibilities[block]=possibilities.get(block,set(range(len(gt))))&allowed
        if not possibilities[block]:return False
    return True


def main():
    root=Path(sys.argv[1]);out=root/'normal-self';out.mkdir(parents=True,exist_ok=True)
    source=root/'normal-coding/gather/called-normal-alleles.json'
    complete=json.loads((source.parent/'complete.json').read_text());assert complete['status']=='COMPLETED'
    plan=json.loads((root/'normal-coding-plan.json').read_text());genome=Genome(**plan['genome'])
    alleles=json.loads(source.read_text());sites={};issues=[];counts=collections.Counter()
    quality_by_site={}
    for row in alleles:
        key=(*row['allele'][:3],tuple(row['source_alts']))
        quality_by_site[key]=quality_by_site.get(key,True) and row['quality_gate']
    for row in alleles:
        chrom,pos,ref,alt=row['allele'];key=(chrom,pos,ref,tuple(row['source_alts']))
        if key in sites:continue
        gt=row['native_genotype']['GT'];options={}
        for index in sorted(set(gt)):
            value=ref if index==0 else row['source_alts'][index-1]
            if set(value)-set('ACGTN'):continue
            options[index]=value
        phased=bool(row.get('native_GT_phased'));phased_gt=gt
        # HaplotypeCaller often records PGT/PID while GT remains unphased.
        pgt=row['native_genotype'].get('PGT');pid=row['native_genotype'].get('PS') or row['native_genotype'].get('PID')
        if pgt and '|' in pgt:
            parsed=[int(x) for x in pgt.split('|')]
            if sorted(parsed)==sorted(gt):phased=True;phased_gt=parsed
        quality=quality_by_site[key]
        sites[key]=dict(chrom=chrom,pos=pos,ref=ref,choices=options,GT=gt,
                        phase_block=(chrom+':'+str(pid)) if phased and pid is not None else None,
                        phased_GT=phased_gt,quality_gate=quality and len(options)==len(set(gt)),
                        native_genotype=row['native_genotype'],native_GT_phased=row.get('native_GT_phased',False))
    # Sweep coding intervals rather than inventing a transcript for an allele.
    coding=collections.defaultdict(list)
    for line in Path(plan['annotation']).open():
        if line.startswith('#'):continue
        fields=line.split('\t')
        if len(fields)!=9 or fields[2] not in ['CDS','start_codon','stop_codon']:continue
        match=re.search(r'transcript_id "([^"]+)"',fields[8])
        if not match:continue
        chrom=fields[0] if fields[0].startswith('chr') else 'chr'+fields[0]
        coding[chrom].append((int(fields[3]),int(fields[4]),match.group(1).split('.')[0]))
    indexes={}
    for chrom,ranges in coding.items():
        ranges.sort();ends=[];maximum=0
        for lo,hi,tid in ranges:maximum=max(maximum,hi);ends.append(maximum)
        indexes[chrom]=(ranges,[r[0] for r in ranges],ends)
    by_transcript=collections.defaultdict(list)
    for site in sites.values():
        ranges,starts,maxends=indexes.get(site['chrom'],([],[],[]));end=site['pos']+len(site['ref'])-1
        i=bisect.bisect_right(starts,end)-1;tids=set()
        while i>=0 and maxends[i]>=site['pos']:
            lo,hi,tid=ranges[i]
            if hi>=site['pos']:tids.add(tid)
            i-=1
        counts['sites_overlapping_annotated_CDS']+=bool(tids)
        for tid in tids:by_transcript[tid].append(site)
    records=[]
    for number,(tid,raw_sites) in enumerate(sorted(by_transcript.items()),1):
        transcript=genome.transcript_by_id(tid)
        try:cds=transcript.coding_sequence
        except Exception as error:
            issues.append(dict(transcript_id=tid,gene=transcript.gene_name,reason='ANNOTATED_CDS_UNAVAILABLE: '+repr(error),sites=raw_sites));continue
        if not cds or len(cds)%3:
            issues.append(dict(transcript_id=tid,gene=transcript.gene_name,reason='INCOMPLETE_ANNOTATED_CDS',sites=raw_sites));continue
        projected=[]
        for site in raw_sites:
            try:
                value=project(site,transcript,cds)
                if not site['quality_gate']:raise ValueError('Normal genotype failed the recorded quality gate')
                projected.append(value)
            except Exception as error:
                issues.append(dict(transcript_id=tid,gene=transcript.gene_name,reason=str(error),site=site))
        clusters=[]
        for site in sorted(projected,key=lambda r:(r['lo'],r['hi'])):
            frameshift=any((len(a)-(site['hi']-site['lo']))%3 for a in site['choices'].values())
            if clusters and (clusters[-1]['frameshift'] or site['lo']-clusters[-1]['end']<=33):
                clusters[-1]['sites'].append(site);clusters[-1]['end']=max(clusters[-1]['end'],site['hi']);clusters[-1]['frameshift']|=frameshift
            else:clusters.append(dict(sites=[site],end=site['hi'],frameshift=frameshift))
        for cluster in clusters:
            ss=cluster['sites'];combinations=1
            for site in ss:combinations*=len(site['choices'])
            reason=None
            if combinations>MAX_COMBINATIONS:reason='LOCAL_HAPLOTYPE_ENUMERATION_CAP'
            if any(a['hi']>b['lo'] for a,b in zip(ss,ss[1:])):reason='OVERLAPPING_NORMAL_ALLELE_REPRESENTATIONS'
            if reason:
                issues.append(dict(transcript_id=tid,gene=transcript.gene_name,reason=reason,sites=ss,possible_combinations=combinations));continue
            left=max(0,(ss[0]['lo']-33)//3*3)
            right=len(cds) if cluster['frameshift'] else min(len(cds),((ss[-1]['hi']+35)//3)*3)
            for choices in itertools.product(*(sorted(site['choices']) for site in ss)):
                if not phase_compatible(ss,choices):counts['phase_incompatible_combinations']+=1;continue
                parts=[];position=left
                for site,index in zip(ss,choices):
                    parts.extend([cds[position:site['lo']],site['choices'][index]]);position=site['hi']
                parts.append(cds[position:right]);sequence=''.join(parts)
                protein=str(Seq(sequence[:len(sequence)//3*3]).translate(to_stop=True))
                if not protein or set(protein)-AA:
                    issues.append(dict(transcript_id=tid,gene=transcript.gene_name,reason='AMBIGUOUS_OR_EMPTY_LOCAL_TRANSLATION',sites=ss));continue
                identifier='NSELF_'+__import__('hashlib').sha256((tid+':'+str(left)+':'+str(choices)+':'+protein).encode()).hexdigest()[:24]
                records.append(dict(context_id=identifier,gene=transcript.gene_name,transcript_id=tid,sequence=protein,
                    reference_CDS_start_0based=left,reference_CDS_end_exclusive=right,selected_normal_alleles=list(choices),
                    normal_sites=ss,interpretation='GENOTYPE_COMPATIBLE_LOCAL_SEQUENCE_POSSIBILITY',
                    measured_complete_haplotype=False,healthy_translation_or_presentation_proven=False))
        if number%100==0:save(out/'status.json',dict(phase='normal_coding_sequence_models',transcripts_completed=number,total=len(by_transcript),contexts=len(records),unresolved=len(issues)))
    typed_hla=json.loads((root/'inputs/patient-hla-proteins.json').read_text())
    records.extend(typed_hla['contexts'])
    save(out/'typed-HLA-source.json',typed_hla)
    save(out/'normal-sequence-contexts.json',records);save(out/'unresolved-normal-sequence-models.json',issues)
    (out/'normal-sequence-contexts.fasta').write_text(''.join('>'+r['context_id']+' '+r['gene']+' '+r.get('transcript_id','typed_HLA_reference')+'\n'+r['sequence']+'\n' for r in records))
    save(out/'scope.json',dict(maximum_local_combinations=MAX_COMBINATIONS,maximum_target_peptide_length=11,
        local_cluster_separation_coding_nt=33,frameshifts_include_downstream_coding_variants=True,
        reference_proteome_searched_separately=True,unphased_combinations='Potential self; conservatively exclude hits, do not claim measured phase',
        negative_result='No hit in the searched genotyped coding models; not proof of normal-tissue absence',
        unresolved_cases_retained=True,normal_expression_or_HLA_presentation_measured=False))
    receipt(out,[out/name for name in ['normal-sequence-contexts.json','normal-sequence-contexts.fasta','unresolved-normal-sequence-models.json','scope.json','typed-HLA-source.json']],
            normal_genotype_source_sha256=sha(source),worker_sha256=sha(__file__),contexts=len(records),
            transcripts_considered=len(by_transcript),unresolved_models=len(issues),counts=dict(counts))


if __name__=='__main__':main()
