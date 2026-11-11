"""Independent paired Mutect2 discovery, complete population review and filtering.

This worker reuses verified reference assets, not old tumour-only calls. Every
newly emitted ALT is normalized and reviewed before native sample-wide filtering.
"""
from pathlib import Path
import bisect
import collections
import gzip
import json
import math
import shutil
import sys

import pysam
from common import gatk,run,receipt,save,sha

OLD_DISCOVERY=Path('/data/results/TN25-112103/mutation-evidence-v2-20261102')
OLD_REVIEW=Path('/data/results/TN25-112103/mutation-evidence-v3-20261102')
CHROMS=['chr'+str(c) for c in list(range(1,23))+['X','Y']]
KEEP={'AF','AC','AN','grpmax','AF_grpmax','AN_grpmax','AC_grpmax','AF_sas','AN_sas','AC_sas','lcr','segdup','non_par'}


def scalar(value,index=0):
    if isinstance(value,(tuple,list)):
        return value[index] if len(value)>index else None
    return value


def prepare(root):
    out=root/'paired/prepare';out.mkdir(parents=True,exist_ok=True)
    normal=json.loads((root/'bqsr/complete.json').read_text())
    identity=json.loads((root/'identity/complete.json').read_text())
    assert normal['status']==identity['status']=='COMPLETED' and identity['identity']=='MATCH'
    old=json.loads((OLD_DISCOVERY/'calling-plan.json').read_text())
    assert old['given_alleles'] is None and old['Caris_used_to_define_intervals'] is False
    ref=pysam.FastaFile(old['reference'])
    for chrom in CHROMS:
        parts=sorted((s for s in old['shards'] if s['chrom']==chrom),key=lambda s:s['start'])
        assert parts[0]['start']==1 and parts[-1]['end']==ref.get_reference_length(chrom)
        assert all(a['end']+1==b['start'] for a,b in zip(parts,parts[1:]))
    normal_path=Path(normal['normal_bam'])
    plan=dict(old,normal_bam=str(normal_path),normal_sample='TN25-112103-N',tumour_sample='TN25-112103',
              normal_bam_sha256=normal['outputs'][normal_path.name]['sha256'],normal_bam_bytes=normal_path.stat().st_size,
              matched_normal=str(normal_path),given_alleles=None,Caris_nominations=False,
              analysis_kind='New independent paired tumour-normal discovery',
              normal_processing='BWA-MEM2, paired fixmate, duplicate marking and native GATK BQSR',
              population_review='Every emitted ALT against cached full-chromosome gnomAD v4.1 genomes and exomes',
              source_interval_plan_sha256=sha(OLD_DISCOVERY/'calling-plan.json'))
    for key in ['samples','sample']:
        plan.pop(key,None)
    for path,sample in [(Path(plan['bam']),plan['tumour_sample']),(normal_path,plan['normal_sample'])]:
        with pysam.AlignmentFile(path,'rb') as bam:
            bam.check_index()
            assert {r['SM'] for r in bam.header.to_dict()['RG']}=={sample}
            assert all(ref.get_reference_length(c)==n for c,n in zip(bam.references,bam.lengths))
    assert sha(Path(plan['bam']))==plan['bam_sha256']
    assert sha(Path(plan['reference']))==plan['reference_sha256']
    populations={}
    controls={}
    for chrom in CHROMS:
        directory=OLD_DISCOVERY/'population'/chrom
        previous=json.loads((directory/'complete.json').read_text())
        assert previous['status']=='VERIFIED'
        path=directory/'caller-pass.vcf.gz'
        assert path.stat().st_size==previous['outputs'][path.name]['bytes']
        populations[chrom]=dict(path=str(path),**previous['outputs'][path.name],source_receipt_sha256=sha(directory/'complete.json'))
        controls.update(previous['spot_checks'])
    # Independently repeat the three controls that caught the earlier resource gap.
    for gene,chrom,pos,rf,alt,expected in [
        ('ATRX','chrX',77682471,'C','G',.5158322228567802),
        ('CYP2D6','chr22',42130692,'G','A',.19041485180070117),
        ('KMT2D','chr12',49052257,'C','T',.0891850591779953)]:
        with pysam.VariantFile(populations[chrom]['path']) as variants:
            rows=[r for r in variants.fetch(chrom,pos-1,pos) if r.pos==pos and r.ref==rf and alt in r.alts]
            assert len(rows)==1
            value=float(scalar(rows[0].info['AF'],rows[0].alts.index(alt)))
            assert abs(value-expected)<1e-4,(gene,value)
            controls[gene]=dict(chromosome=chrom,position=pos,ref=rf,alt=alt,AF=value,verified_again=True)
    plan['population_resources']=populations
    save(root/'calling-plan.json',plan)
    save(out/'population-controls.json',controls)
    case=json.loads((OLD_DISCOVERY/'case.json').read_text())
    case['matched_normal_DNA']=True;case['kind']='paired_analysis_context'
    save(root/'case.json',case)
    for name in ['analysis-genome-config.json','analysis-environment-verified.json','affinity-model-manifest.json']:
        shutil.copyfile(OLD_DISCOVERY/name,root/name)
    receipt(out,[out/'population-controls.json'],normal_bam_sha256=plan['normal_bam_sha256'],
            calling_plan_sha256=sha(root/'calling-plan.json'),shards=len(plan['shards']),
            independent_discovery=True,Caris_nominations=False,worker_sha256=sha(__file__))


def call(root,sid):
    plan=json.loads((root/'calling-plan.json').read_text())
    shard=next(s for s in plan['shards'] if s['id']==sid)
    out=root/'paired/calling'/sid;out.mkdir(parents=True,exist_ok=True)
    assert plan['given_alleles'] is None and plan['matched_normal']==plan['normal_bam']
    population=plan['population_resources'][shard['chrom']]
    assert sha(population['path'])==population['sha256']
    normal=Path(plan['normal_bam']);assert normal.stat().st_size==plan['normal_bam_bytes']
    command=gatk('Mutect2','-R',plan['reference'],'-I',plan['bam'],'-I',normal,
        '-normal',plan['normal_sample'],'--germline-resource',population['path'],
        '-L',shard['interval'],'--native-pair-hmm-threads','1',
        '--f1r2-tar-gz',out/'f1r2.tar.gz','-O',out/'unfiltered.vcf.gz')
    assert '--alleles' not in command
    run(out,'mutect2-paired',command)
    receipt(out,[out/n for n in ['unfiltered.vcf.gz','unfiltered.vcf.gz.tbi','unfiltered.vcf.gz.stats','f1r2.tar.gz']],
            shard=shard,normal_sample=plan['normal_sample'],tumour_sample=plan['tumour_sample'],
            independent_discovery=True,Caris_nominations=False,worker_sha256=sha(__file__))


def gather(root):
    plan=json.loads((root/'calling-plan.json').read_text());out=root/'paired/gather';out.mkdir(parents=True,exist_ok=True)
    paths=[]
    for shard in plan['shards']:
        directory=root/'paired/calling'/shard['id'];completed=json.loads((directory/'complete.json').read_text())
        assert completed['status']=='COMPLETED' and completed['shard']==shard
        for name,data in completed['outputs'].items():assert sha(directory/name)==data['sha256']
        paths.append(directory/'unfiltered.vcf.gz')
    merged=out/'merged.vcf.gz';stats=out/'unfiltered.vcf.gz.stats';orientation=out/'orientation-model.tar.gz'
    run(out,'merge-vcfs',gatk('MergeVcfs',*[a for p in paths for a in ['-I',p]],'-O',merged,heap='11g'))
    run(out,'merge-statistics',gatk('MergeMutectStats',*[a for p in paths for a in ['--stats',str(p)+'.stats']],'-O',stats,heap='11g'))
    run(out,'read-orientation-model',gatk('LearnReadOrientationModel',*[a for p in paths for a in ['-I',p.parent/'f1r2.tar.gz']],'-O',orientation,heap='11g'))
    tagged=out/'tagged.vcf.gz';writers={};record_count=0;alt_count=0
    try:
        with pysam.VariantFile(merged) as source:
            header=source.header.copy()
            assert set(header.samples)=={plan['tumour_sample'],plan['normal_sample']}
            header.info.add('MN_ALT_INDEX','A','Integer','Original ALT index before population normalization')
            for chrom in CHROMS:
                directory=root/'paired/population'/chrom;directory.mkdir(parents=True,exist_ok=True)
                writers[chrom]=pysam.VariantFile(directory/'raw.vcf.gz','wz',header=header)
            with pysam.VariantFile(tagged,'wz',header=header) as destination:
                for i,record in enumerate(source,1):
                    record.translate(header);record.id='mnraw_'+str(i)
                    assert record.contig in writers and len(record.info['POPAF'])==len(record.alts)
                    record.info['MN_ALT_INDEX']=tuple(range(len(record.alts)))
                    destination.write(record);writers[record.contig].write(record)
                    record_count+=1;alt_count+=len(record.alts)
    finally:
        for writer in writers.values():writer.close()
    normalized=out/'normalized-all.vcf.gz'
    run(out,'normalize-all-emitted-alleles',['bcftools','norm','-f',plan['reference'],'-m','-any','-Oz','-o',normalized,tagged])
    handles={c:(root/'paired/population'/c/'alleles.jsonl').open('w') for c in CHROMS}
    counts=collections.Counter()
    try:
        with pysam.VariantFile(normalized) as records:
            for r in records:
                assert len(r.alts)==1
                handles[r.contig].write(json.dumps(dict(raw_id=r.id,raw_alt_index=r.info['MN_ALT_INDEX'][0],
                    normalized_allele=[r.contig,r.pos,r.ref,r.alts[0]],original_POPAF=r.info['POPAF'][0]))+'\n')
                counts[r.contig]+=1
    finally:
        for handle in handles.values():handle.close()
    assert sum(counts.values())==alt_count
    receipt(out,[merged,tagged,normalized,stats,orientation],records=record_count,all_emitted_alleles=alt_count,
            counts_by_chromosome=dict(counts),normal_and_tumour_explicit=True,Caris_nominations=False,worker_sha256=sha(__file__))


def equivalent(fasta,key,record):
    chrom,pos,rf,al=key
    if (pos,rf,al)==(record.pos,record.ref,record.alts[0]):return True
    if len(rf)==len(al) or len(record.ref)==len(record.alts[0]):return False
    if set(rf+al+record.ref+record.alts[0])-set('ACGTN'):return False
    lo=min(pos,record.pos)-1;hi=max(pos+len(rf),record.pos+len(record.ref))+1
    sequence=fasta.fetch(chrom,lo,hi).upper()
    def apply(p,r,a):
        off=p-1-lo
        return sequence[:off]+a+sequence[off+len(r):] if sequence[off:off+len(r)]==r else None
    left=apply(pos,rf,al)
    return left is not None and left==apply(record.pos,record.ref,record.alts[0])


def population(root,chrom):
    out=root/'paired/population'/chrom;out.mkdir(parents=True,exist_ok=True)
    plan=json.loads((root/'calling-plan.json').read_text());fasta=pysam.FastaFile(plan['reference'])
    genome=OLD_DISCOVERY/'population'/chrom/'annotation.vcf.gz'
    exome=OLD_REVIEW/'population'/chrom/'exomes-annotation.vcf.gz'
    genome_receipt=json.loads((genome.parent/'complete.json').read_text())
    exome_receipt=json.loads((exome.parent/'complete.json').read_text())
    assert genome_receipt['status']=='VERIFIED' and exome_receipt['status']=='COMPLETED'
    assert sha(genome)==genome_receipt['outputs'][genome.name]['sha256']
    assert sha(exome)==exome_receipt['exome_compact_sha256']
    bykey=collections.defaultdict(list)
    for line in (out/'alleles.jsonl').open():
        r=json.loads(line);bykey[tuple(r['normalized_allele'])].append(r)
    positions=sorted({key[1] for key in bykey});keys_by_position=collections.defaultdict(list)
    for key in bykey:keys_by_position[key[1]].append(key)
    intervals=[]
    for pos in positions:
        lo=max(0,pos-101);hi=max(pos+len(key[2])+100 for key in keys_by_position[pos])
        if intervals and lo<=intervals[-1][1]:intervals[-1][1]=max(intervals[-1][1],hi)
        else:intervals.append([lo,hi])
    observations={key:{'genomes':[],'exomes':[]} for key in bykey}
    for cohort,path in [('genomes',genome),('exomes',exome)]:
        seen=set()
        with pysam.VariantFile(path) as source:
            fields=set(source.header.info)
            for index,(lo,hi) in enumerate(intervals):
                for record in source.fetch(chrom,lo,hi):
                    assert len(record.alts)==1
                    signature=(record.pos,record.ref,record.alts[0])
                    if signature in seen:continue
                    seen.add(signature)
                    left=bisect.bisect_left(positions,max(1,record.pos-100))
                    right=bisect.bisect_right(positions,record.pos+len(record.ref)+100)
                    for pos in positions[left:right]:
                        for key in keys_by_position[pos]:
                            if equivalent(fasta,key,record):
                                observations[key][cohort].append(dict(position=record.pos,ref=record.ref,alt=record.alts[0],
                                    FILTER=list(record.filter),**{name:scalar(record.info.get(name)) if name in fields else None for name in sorted(KEEP)}))
                if index%1000==0:
                    save(out/'status.json',dict(phase='matching_'+cohort,intervals=index,total_intervals=len(intervals)))
    priors={};counts=collections.Counter()
    evidence=out/'population-observations.jsonl.gz'
    with gzip.open(evidence,'wt') as stream:
        for key,cohorts in observations.items():
            good=[r for records in cohorts.values() for r in records if r['FILTER']==['PASS'] and r.get('AF') is not None and r.get('AN')]
            global_af=max([r['AF'] for r in good],default=None)
            for original in bykey[key]:
                native=10**-original['original_POPAF'];updated=max(native,global_af if global_af is not None else native)
                location=original['raw_id'],original['raw_alt_index']
                assert location not in priors and 0<updated<=1
                priors[location]=-math.log10(updated);counts['prior_increased']+=updated>native*1.01
            stream.write(json.dumps(dict(normalized_allele=key,population=[dict(cohort=c,status='MATCHED' if records else 'NOT_FOUND',records=records) for c,records in cohorts.items()]))+'\n')
    updated=out/'updated.vcf.gz';seen=set()
    with pysam.VariantFile(out/'raw.vcf.gz') as source:
        with pysam.VariantFile(updated,'wz',header=source.header) as writer:
            for record in source:
                values=list(record.info['POPAF'])
                for i in range(len(record.alts)):
                    key=record.id,i;assert key in priors and key not in seen
                    seen.add(key);values[i]=priors[key]
                record.info['POPAF']=tuple(values);writer.write(record)
    assert seen==set(priors)
    run(out,'index-updated',['bcftools','index','-t','-f',updated])
    receipt(out,[updated,Path(str(updated)+'.tbi'),evidence],chromosome=chrom,all_emitted_alleles=len(priors),
            unique_alleles=len(bykey),counts=dict(counts),genome_sha256=sha(genome),exome_sha256=sha(exome),
            policy='max(native, QC-PASS genome global AF, QC-PASS exome global AF); missing is unknown, not measured zero',
            genome_context_note='Compact genome resource omits lcr/segdup; restore original source flags for coding candidates during annotation',
            worker_sha256=sha(__file__))


def filter_calls(root):
    plan=json.loads((root/'calling-plan.json').read_text());out=root/'paired/filter';out.mkdir(parents=True,exist_ok=True)
    paths=[];expected=0
    for chrom in CHROMS:
        directory=root/'paired/population'/chrom;completed=json.loads((directory/'complete.json').read_text())
        assert completed['status']=='COMPLETED';expected+=completed['all_emitted_alleles']
        path=directory/'updated.vcf.gz';assert sha(path)==completed['outputs'][path.name]['sha256'];paths.append(path)
    merged=out/'unfiltered.vcf.gz';filtered=out/'filtered.vcf.gz'
    run(out,'merge-population-reviewed-calls',gatk('MergeVcfs',*[a for p in paths for a in ['-I',p]],'-O',merged,heap='11g'))
    run(out,'filter-paired-calls',gatk('FilterMutectCalls','-R',plan['reference'],'-V',merged,
        '--stats',root/'paired/gather/unfiltered.vcf.gz.stats',
        '--contamination-table',root/'identity/tumour-contamination.table',
        '--tumor-segmentation',root/'identity/allele-fraction-segments.table',
        '--ob-priors',root/'paired/gather/orientation-model.tar.gz',
        '--threshold-strategy','FALSE_DISCOVERY_RATE','--false-discovery-rate','0.05',
        '--filtering-stats',out/'filtering-stats.tsv','-O',filtered,heap='11g'))
    tagged=out/'tagged.vcf.gz'
    with pysam.VariantFile(filtered) as source:
        header=source.header.copy()
        header.info.add('MN_EFFECTIVE_FILTER','A','String','Native site plus per-ALT filters, preserved before normalization')
        with pysam.VariantFile(tagged,'wz',header=header) as writer:
            for r in source:
                r.translate(header);raw=r.info.get('AS_FilterStatus','SITE')
                if isinstance(raw,tuple):raw=','.join(raw)
                states=raw.split('|');assert len(states)==len(r.alts)
                site=set(r.filter)-{'PASS','.'}
                r.info['MN_EFFECTIVE_FILTER']=tuple('&'.join(sorted(site|(set(x.split(','))-{'SITE','PASS','.',''}))) or 'PASS' for x in states)
                writer.write(r)
    normalized=out/'normalized.vcf.gz'
    run(out,'normalize-filtered',['bcftools','norm','-f',plan['reference'],'-m','-any','-Oz','-o',normalized,tagged])
    run(out,'index-filtered',['bcftools','index','-t','-f',normalized])
    counts=collections.Counter();records={};conflicts=[]
    with pysam.VariantFile(normalized) as source:
        header=source.header.copy()
        header.info.add('MN_SCATTER_CONFLICT',0,'Flag','Conflicting native evidence among duplicate normalized scatter calls; held for review')
        for r in source:
            assert len(r.alts)==1;counts['all_emitted_alleles']+=1
            effective=r.info['MN_EFFECTIVE_FILTER'][0];counts['native_'+effective]+=1
            key=r.contig,r.pos,r.ref,r.alts[0]
            r.translate(header)
            signature=lambda x:json.dumps(dict(samples={name:dict(x.samples[name].items()) for name in [plan['tumour_sample'],plan['normal_sample']]},GERMQ=x.info.get('GERMQ'),effective=x.info['MN_EFFECTIVE_FILTER'][0]),sort_keys=True)
            if key in records:
                counts['duplicate_normalized_alleles']+=1
                if signature(records[key])!=signature(r):
                    records[key].info['MN_SCATTER_CONFLICT']=True
                    conflicts.append(dict(allele=key,first=signature(records[key]),other=signature(r)))
            else:records[key]=r.copy()
    assert counts['all_emitted_alleles']==expected
    order={c:i for i,c in enumerate(header.contigs)}
    all_unique=out/'all-unique.vcf.gz';passed=out/'pass.vcf.gz'
    with pysam.VariantFile(all_unique,'wz',header=header) as all_writer,pysam.VariantFile(passed,'wz',header=header) as pass_writer:
        for key in sorted(records,key=lambda k:(order[k[0]],k[1],k[2],k[3])):
            r=records[key];r.id='mn_'+'_'.join(map(str,key));all_writer.write(r)
            if r.info['MN_EFFECTIVE_FILTER'][0]=='PASS' and not r.info.get('MN_SCATTER_CONFLICT'):
                pass_writer.write(r);counts['unique_native_PASS_without_scatter_conflict']+=1
    for path in [all_unique,passed]:run(out,'index-'+path.stem,['bcftools','index','-t','-f',path])
    save(out/'scatter-conflicts.json',conflicts)
    receipt(out,[merged,filtered,normalized,Path(str(normalized)+'.tbi'),all_unique,Path(str(all_unique)+'.tbi'),passed,Path(str(passed)+'.tbi'),out/'filtering-stats.tsv',out/'scatter-conflicts.json'],
            counts=dict(counts),unique_alleles=len(records),scatter_conflicts=len(conflicts),
            scope='All newly emitted paired calls; no old variant whitelist',native_strategy='FALSE_DISCOVERY_RATE',model_FDR_target=.05,
            empirical_FDR_established=False,worker_sha256=sha(__file__))


def main():
    root=Path(sys.argv[1]);action=sys.argv[2]
    if action=='prepare':prepare(root)
    elif action=='call':call(root,sys.argv[3])
    elif action=='gather':gather(root)
    elif action=='population':population(root,sys.argv[3])
    elif action=='filter':filter_calls(root)
    else:raise ValueError(action)


if __name__=='__main__':main()
