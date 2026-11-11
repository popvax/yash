"""Normal-only coding-region genotypes for callability and self-sequence review.

These are genotype measurements, not a clinical inherited-disease report.
"""
from pathlib import Path
import collections
import bisect
import json
import sys

import pysam
from common import gatk,run,receipt,save,sha

CHROMS=['chr'+str(c) for c in list(range(1,23))+['X','Y']]
OLD=Path('/data/results/TN25-112103/mutation-evidence-v2-20261102')


def prepare(root):
    out=root/'normal-coding/prepare';out.mkdir(parents=True,exist_ok=True)
    normal=json.loads((root/'bqsr/complete.json').read_text())
    identity=json.loads((root/'identity/complete.json').read_text())
    assert normal['status']==identity['status']=='COMPLETED' and identity['identity']=='MATCH'
    genome=json.loads((OLD/'analysis-genome-config.json').read_text())
    annotation=Path(genome['gtf_path_or_url']);fasta=pysam.FastaFile(normal['reference'])
    ranges=collections.defaultdict(list)
    with annotation.open() as stream:
        for line in stream:
            if line.startswith('#'):continue
            fields=line.rstrip().split('\t')
            if len(fields)!=9 or fields[2] not in ['CDS','start_codon','stop_codon']:continue
            chrom=fields[0] if fields[0].startswith('chr') else 'chr'+fields[0]
            if chrom not in CHROMS:continue
            lo=max(0,int(fields[3])-1-100);hi=min(fasta.get_reference_length(chrom),int(fields[4])+100)
            ranges[chrom].append((lo,hi))
    shards=[];total_bases=0
    for chrom in CHROMS:
        merged=[]
        for lo,hi in sorted(ranges[chrom]):
            if merged and lo<=merged[-1][1]:merged[-1][1]=max(merged[-1][1],hi)
            else:merged.append([lo,hi])
        # Partition by disjoint genomic bins. Each base belongs to one BED shard.
        buckets=collections.defaultdict(list)
        for lo,hi in merged:
            while lo<hi:
                bucket=lo//20000000;end=min(hi,(bucket+1)*20000000)
                buckets[bucket].append([lo,end]);lo=end
        for bucket,intervals in sorted(buckets.items()):
            sid=f'{chrom}-{bucket:03d}';bed=out/(sid+'.bed')
            bed.write_text(''.join(f'{chrom}\t{lo}\t{hi}\n' for lo,hi in intervals))
            bases=sum(hi-lo for lo,hi in intervals);total_bases+=bases
            shards.append(dict(id=sid,chrom=chrom,bed=str(bed),bed_sha256=sha(bed),bases=bases,
                               sample_ploidy=2))
    assert shards and total_bases>20000000
    plan=dict(normal_bam=normal['normal_bam'],normal_sample='TN25-112103-N',reference=normal['reference'],
              normal_bam_sha256=normal['outputs'][Path(normal['normal_bam']).name]['sha256'],
              annotation=str(annotation),annotation_sha256=sha(annotation),genome=genome,
              scope='GENCODE/Ensembl 103 primary nuclear CDS and start/stop codons, with 100-bp flanks',
              excluded='Mitochondrial and unplaced/alternate contigs; these are unassessed, not reference-negative',
              ploidy_note='Diploid genotype model; sex-chromosome dosage is not clinically interpreted',
              total_target_bases=total_bases,shards=shards)
    save(root/'normal-coding-plan.json',plan)
    receipt(out,[out/(s['id']+'.bed') for s in shards],shards=len(shards),target_bases=total_bases,
            normal_bam_sha256=plan['normal_bam_sha256'],annotation_sha256=plan['annotation_sha256'],worker_sha256=sha(__file__))


def call(root,sid):
    plan=json.loads((root/'normal-coding-plan.json').read_text())
    shard=next(s for s in plan['shards'] if s['id']==sid)
    bed=Path(shard['bed']);assert sha(bed)==shard['bed_sha256']
    out=root/'normal-coding/calling'/sid;out.mkdir(parents=True,exist_ok=True)
    gvcf=out/'normal.g.vcf.gz'
    run(out,'haplotype-caller',gatk('HaplotypeCaller','-R',plan['reference'],'-I',plan['normal_bam'],
        '-L',bed,'-ERC','GVCF','--native-pair-hmm-threads','1','--sample-ploidy',shard['sample_ploidy'],'-O',gvcf))
    with pysam.VariantFile(gvcf) as records:
        assert list(records.header.samples)==[plan['normal_sample']]
    receipt(out,[gvcf,Path(str(gvcf)+'.tbi')],shard=shard,normal_sample=plan['normal_sample'],
            plan_sha256=sha(root/'normal-coding-plan.json'),worker_sha256=sha(__file__))


def gather(root):
    plan=json.loads((root/'normal-coding-plan.json').read_text())
    out=root/'normal-coding/gather';out.mkdir(parents=True,exist_ok=True)
    inputs=[]
    target_intervals=collections.defaultdict(list)
    for shard in plan['shards']:
        directory=root/'normal-coding/calling'/shard['id'];rec=json.loads((directory/'complete.json').read_text())
        assert rec['status']=='COMPLETED' and rec['shard']==shard
        path=directory/'normal.g.vcf.gz';assert sha(path)==rec['outputs'][path.name]['sha256'];inputs.append(path)
        for line in Path(shard['bed']).read_text().splitlines():
            chrom,lo,hi=line.split('\t');target_intervals[chrom].append((int(lo),int(hi)))
    # Merge rather than assume file order; a source interval can emit an event
    # extending over a BED boundary. Duplicates are reconciled explicitly below.
    merged=out/'normal.coding.g.vcf.gz'
    run(out,'merge-normal-gvcfs',gatk('MergeVcfs',*[a for p in inputs for a in ['-I',p]],'-O',merged,heap='11g'))
    target_starts={}
    for chrom,intervals in target_intervals.items():
        intervals.sort();target_starts[chrom]=[lo for lo,hi in intervals]
    coverage={label:collections.defaultdict(list) for label in ['genotyped','GQ20_DP10','GQ20_DP20','GQ20_DP30']}
    def add_coverage(label,chrom,start,end):
        intervals=target_intervals.get(chrom,[]);starts=target_starts.get(chrom,[])
        index=max(0,bisect.bisect_right(starts,start)-1)
        while index<len(intervals) and intervals[index][0]<end:
            lo=max(start,intervals[index][0]);hi=min(end,intervals[index][1]);index+=1
            if hi<=lo:continue
            ranges=coverage[label][chrom]
            if ranges and ranges[-1][0]<=lo<=ranges[-1][1]:ranges[-1][1]=max(ranges[-1][1],hi)
            else:ranges.append([lo,hi])
    variants=[];signatures={};variant_indexes={};conflicts=[];counts=collections.Counter()
    with pysam.VariantFile(merged) as records:
        for record in records:
            counts['records']+=1
            sample=record.samples[plan['normal_sample']];gt=sample.get('GT')
            if not gt or any(i is None for i in gt):continue
            add_coverage('genotyped',record.contig,record.start,record.stop)
            gq=sample.get('GQ');depth=sample.get('MIN_DP') if sample.get('MIN_DP') is not None else sample.get('DP')
            if gq is not None and gq>=20 and not (set(record.filter)-{'PASS','.'}):
                for cutoff in [10,20,30]:
                    if depth is not None and depth>=cutoff:add_coverage(f'GQ20_DP{cutoff}',record.contig,record.start,record.stop)
            for index in sorted(set(gt)-{0}):
                alt=record.alts[index-1]
                if alt in ['<NON_REF>','*'] or set(alt+record.ref)-set('ACGTN'):continue
                counts['called_non_reference_alleles']+=1
                key=(record.contig,record.pos,record.ref,alt)
                evidence=dict(allele=key,source_alts=list(record.alts),native_genotype=dict(sample.items()),
                              native_GT_phased=sample.phased,
                              native_FILTER=list(record.filter),native_QUAL=record.qual,
                              target_alt_index=index,target_present=True,source_stop=record.stop)
                signature=json.dumps(evidence,sort_keys=True)
                if key in signatures:
                    counts['duplicate_called_alleles']+=1
                    if signatures[key]!=signature:
                        conflicts.append(dict(allele=key,evidence=evidence))
                        variants[variant_indexes[key]]['quality_gate']=False
                        variants[variant_indexes[key]]['scatter_conflict']=True
                    continue
                signatures[key]=signature
                depth=sample.get('DP');quality=sample.get('GQ')
                explicit_filters=set(record.filter)-{'PASS','.'}
                evidence['quality_gate']=bool(depth is not None and depth>=20 and quality is not None and quality>=20 and not explicit_filters)
                evidence['scatter_conflict']=False
                variant_indexes[key]=len(variants)
                variants.append(evidence)
    # Long events or duplicate scatter records can make clipped ranges arrive
    # out of order. Canonicalize their union before counting covered bases.
    for chroms in coverage.values():
        for chrom,ranges in chroms.items():
            merged_ranges=[]
            for lo,hi in sorted(ranges):
                if merged_ranges and lo<=merged_ranges[-1][1]:merged_ranges[-1][1]=max(merged_ranges[-1][1],hi)
                else:merged_ranges.append([lo,hi])
            chroms[chrom]=merged_ranges
    coverage_bases={label:sum(hi-lo for ranges in chroms.values() for lo,hi in ranges) for label,chroms in coverage.items()}
    assert all(value<=plan['total_target_bases'] for value in coverage_bases.values())
    save(out/'callability-summary.json',dict(target_bases=plan['total_target_bases'],genotype_covered_bases=coverage_bases,
         target_scope=plan['scope'],reference_blocks_use_MIN_DP=True,source='Native GATK GVCF genotype evidence, clipped to requested intervals',
         individual_variant_read_evidence_still_required=True))
    save(out/'called-normal-alleles.json',variants);save(out/'scatter-conflicts.json',conflicts)
    receipt(out,[merged,Path(str(merged)+'.tbi'),out/'called-normal-alleles.json',out/'scatter-conflicts.json',out/'callability-summary.json'],
            counts=dict(counts),unique_called_alleles=len(variants),quality_passing_alleles=sum(r['quality_gate'] for r in variants),
            scatter_conflicts=len(conflicts),normal_bam_sha256=plan['normal_bam_sha256'],
            scope=plan['scope'],not_a_complete_normal_transcriptome=True,worker_sha256=sha(__file__))


def main():
    root=Path(sys.argv[1]);action=sys.argv[2]
    if action=='prepare':prepare(root)
    elif action=='call':call(root,sys.argv[3])
    elif action=='gather':gather(root)
    else:raise ValueError(action)


if __name__=='__main__':main()
