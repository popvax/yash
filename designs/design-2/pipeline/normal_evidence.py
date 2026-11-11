"""Read-only native normal-genotype and focal-allele evidence helpers."""
from collections import Counter,defaultdict
from pathlib import Path
import json

import pysam
from isovar import ReadCollector
from isovar.base_quality import BaseQualityPolicy


class NormalEvidence:
    def __init__(self,root):
        root=Path(root)
        done=json.loads((root/'normal-coding/gather/complete.json').read_text())
        assert done['status']=='COMPLETED'
        plan=json.loads((root/'normal-coding-plan.json').read_text())
        self.bam=pysam.AlignmentFile(plan['normal_bam'],'rb');self.bam.check_index()
        self.fasta=pysam.FastaFile(plan['reference'])
        self.gvcf=pysam.VariantFile(root/'normal-coding/gather/normal.coding.g.vcf.gz')
        assert list(self.gvcf.header.samples)==[plan['normal_sample']]
        self.sample=plan['normal_sample']
        self.collector=ReadCollector(use_secondary_alignments=False,use_duplicate_reads=False,
            min_mapping_quality=20,merge_overlapping_fragments=False,use_reads_without_base_qualities=False,
            read_filter=lambda read:not read.is_supplementary and not read.is_qcfail
                and (not read.has_tag('NH') or read.get_tag('NH')==1))
        self.quality=BaseQualityPolicy(20)
        self.scatter_conflicts=json.loads((root/'normal-coding/gather/scatter-conflicts.json').read_text())

    def equivalent(self,allele,position,ref,alt):
        chrom,pos,rf,al=allele
        if (pos,rf,al)==(position,ref,alt):return True
        if set(rf+al+ref+alt)-set('ACGTN'):return False
        lo=min(pos,position)-1;hi=max(pos+len(rf),position+len(ref))+1
        sequence=self.fasta.fetch(chrom,lo,hi).upper()
        def apply(p,r,a):
            offset=p-1-lo
            return sequence[:offset]+a+sequence[offset+len(r):] if sequence[offset:offset+len(r)]==r else None
        left=apply(pos,rf,al)
        return left is not None and left==apply(position,ref,alt)

    def read_counts(self,variant):
        evidence=self.collector.read_evidence_for_variant(variant,self.bam)
        fragments=defaultdict(set);raw_counts=Counter();quality_status=Counter()
        for category,reads in [('REF',evidence.ref_reads),('ALT',evidence.alt_reads),('OTHER',evidence.other_reads)]:
            for read in reads:
                assessments=self.quality.assess(read)
                passing=any(item.status=='passed' for item in assessments)
                quality_status['passed' if passing else 'failed' if any(item.status=='failed' for item in assessments) else 'unassessed']+=1
                if passing:
                    fragments[read.name].add(category);raw_counts[category]+=1
        counts=Counter(next(iter(labels)) if len(labels)==1 else 'DISAGREEING_MATES' for labels in fragments.values())
        return dict(qualifying_fragment_names=len(fragments),usable_fragment_names=sum(v for k,v in counts.items() if k!='DISAGREEING_MATES'),
                    fragment_counts=dict(counts),quality_passing_observations=dict(raw_counts),quality_assessment_counts=dict(quality_status),
                    base_quality_policy=self.quality.description(),fragment_unit='Original read names; no UMI molecule claim')

    def genotype(self,allele):
        chrom,pos,rf,al=allele;start=pos-1;end=start+len(rf)
        conflicting=[r for r in self.scatter_conflicts if r['allele'][0]==chrom and r['allele'][1]-1<end
                     and r['allele'][1]-1+len(r['allele'][2])>start]
        observations=[];positive=[];reference_ranges=[];complex_records=[]
        for record in self.gvcf.fetch(chrom,max(0,start-101),end+100):
            if record.start>=end or record.stop<=start:continue
            sample=record.samples[self.sample];gt=sample.get('GT');gq=sample.get('GQ')
            complete=bool(gt and all(i is not None for i in gt))
            explicit_filters=set(record.filter)-{'PASS','.'}
            high_quality=bool(complete and gq is not None and gq>=20 and not explicit_filters)
            row=dict(position=record.pos,ref=record.ref,alts=list(record.alts or []),stop=record.stop,
                     genotype=dict(sample.items()),native_FILTER=list(record.filter),high_quality=high_quality)
            observations.append(row)
            if not high_quality:continue
            called=[record.ref if i==0 else record.alts[i-1] for i in gt]
            target_indexes=[i for i,a in enumerate(record.alts or [],1) if a not in ['<NON_REF>','*']
                            and self.equivalent(allele,record.pos,record.ref,a)]
            if any(i in target_indexes for i in gt):
                positive.append(row);continue
            # Equal-length native haplotypes can contain a nominated SNV along
            # with another difference. Examine the called base, not only exact VCF keys.
            if len(rf)==len(al)==1 and record.start<=start<record.start+len(record.ref) \
                    and all(len(a)==len(record.ref) and set(a)<=set('ACGTN') for a in called):
                bases={a[start-record.start] for a in called}
                if al in bases:positive.append(row)
                else:reference_ranges.append((start,end,gq))
                continue
            if all(i==0 for i in gt):
                reference_ranges.append((max(start,record.start),min(end,record.stop),gq))
            else:
                complex_records.append(row)
        ranges=[]
        for lo,hi,gq in sorted(reference_ranges):
            if ranges and lo<=ranges[-1][1]:ranges[-1][1]=max(ranges[-1][1],hi)
            else:ranges.append([lo,hi])
        covered=bool(ranges and ranges[0][0]<=start and ranges[-1][1]>=end and len(ranges)==1)
        if conflicting or (positive and covered):status='NORMAL_GENOTYPE_CONFLICT'
        elif positive:status='NORMAL_TARGET_ALLELE_PRESENT'
        elif complex_records:status='NORMAL_COMPLEX_OVERLAP_REQUIRES_REVIEW'
        elif covered:status='NORMAL_TARGET_ALLELE_NOT_CALLED'
        else:status='NORMAL_GENOTYPE_UNRESOLVED'
        return dict(status=status,high_quality_reference_span_complete=covered,
                    positive_native_genotypes=positive,observations=observations,scatter_conflicts=conflicting,
                    note='Reference genotype evidence is conditional on the native model; missing records are not absence')

    def assess(self,allele,variant):
        reads=self.read_counts(variant);genotypes=self.genotype(allele)
        alt=reads['fragment_counts'].get('ALT',0);adequate=reads['usable_fragment_names']>=20
        status=genotypes['status']
        if status=='NORMAL_TARGET_ALLELE_PRESENT' and adequate and alt>0:
            classification='PRESENT_IN_NORMAL'
        elif status=='NORMAL_TARGET_ALLELE_NOT_CALLED' and adequate and alt<3:
            classification='NORMAL_REFERENCE_SUPPORTED'
        else:
            classification='UNRESOLVED_NORMAL_EVIDENCE'
        return dict(classification=classification,adequate_qualifying_normal_fragments=adequate,
                    target_alternate_normal_fragments=alt,read_evidence=reads,genotype_evidence=genotypes,
                    no_normal_ALT_call_is_not_a_zero_coverage_rule=True,
                    healthy_tissue_expression_or_presentation_proven=False)
