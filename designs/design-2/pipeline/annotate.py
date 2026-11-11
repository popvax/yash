"""Annotate new paired PASS alleles and attach measured normal evidence."""
from pathlib import Path
import bisect
import collections
import csv
import gzip
import json
import math
import sys
import urllib.request

import pysam
import varcode
from pyensembl import Genome
from common import sha,save,receipt
from normal_evidence import NormalEvidence

OLD=Path('/data/results/TN25-112103/mutation-evidence-v2-20261102')
OLD_REVIEW=Path('/data/results/TN25-112103/mutation-evidence-v3-20261102')


def scalar(value,index=0):
    return value[index] if isinstance(value,(tuple,list)) and len(value)>index else None if isinstance(value,(tuple,list)) else value


def table(path,rows):
    fields=list(dict.fromkeys(key for row in rows for key in row)) or ['variant_id','status']
    with path.open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,delimiter='\t');writer.writeheader()
        for row in rows:writer.writerow({key:json.dumps(value) if isinstance(value,(tuple,list,dict)) else value for key,value in row.items()})


def main():
    root=Path(sys.argv[1]);out=root/'annotate';out.mkdir(parents=True,exist_ok=True)
    native=json.loads((root/'paired/filter/complete.json').read_text());assert native['status']=='COMPLETED'
    plan=json.loads((root/'calling-plan.json').read_text())
    policy=json.loads((root/'inputs/paired-policy.json').read_text());save(root/'policy.json',policy)
    assert json.loads((root/'audit/complete.json').read_text())['status']=='COMPLETED'
    focal_audit={tuple(row['allele']):row for row in json.loads((root/'audit/current-mutation-audit.json').read_text())}
    genome_config=json.loads((root/'analysis-genome-config.json').read_text())
    genome=Genome(**genome_config)
    normal=NormalEvidence(root)
    chromosomes={shard['chrom'] for shard in plan['shards']};raw=collections.defaultdict(list)
    with Path(genome_config['gtf_path_or_url']).open() as stream:
        for line in stream:
            if line.startswith('#'):continue
            fields=line.rstrip().split('\t')
            if len(fields)!=9 or fields[2] not in ['CDS','start_codon','stop_codon']:continue
            chrom=fields[0] if fields[0].startswith('chr') else 'chr'+fields[0]
            if chrom in chromosomes:raw[chrom].append((max(1,int(fields[3])-2),int(fields[4])+2))
    intervals={};starts={}
    for chrom,ranges in raw.items():
        merged=[]
        for lo,hi in sorted(ranges):
            if merged and lo<=merged[-1][1]+1:merged[-1][1]=max(merged[-1][1],hi)
            else:merged.append([lo,hi])
        intervals[chrom]=merged;starts[chrom]=[lo for lo,hi in merged]
    caris_meta=json.loads((OLD/'gather/caris-comparison-metadata.json').read_text());caris=collections.defaultdict(list)
    with pysam.VariantFile(OLD/'gather/caris-comparison-normalized.vcf.gz') as records:
        for r in records:caris[r.contig,r.pos,r.ref,r.alts[0]].append(r.id)
    candidates=[];errors=[];counts=collections.Counter()
    with pysam.VariantFile(root/'paired/filter/pass.vcf.gz') as records:
        assert set(records.header.samples)=={plan['tumour_sample'],plan['normal_sample']}
        for r in records:
            counts['all_unique_native_PASS']+=1
            key=(r.contig,r.pos,r.ref,r.alts[0]);chrom,pos,rf,alt=key
            index=bisect.bisect_right(starts.get(chrom,[]),pos+len(rf)-1)-1
            if index<0 or intervals[chrom][index][1]<pos:continue
            counts['overlaps_coding_or_2bp_boundary']+=1
            try:
                variant=varcode.Variant(chrom,pos,rf,alt,genome=genome,convert_ucsc_contig_names=True)
                effects=variant.effects();altered=[e for e in effects if getattr(e,'modifies_protein_sequence',False)]
                if not altered:
                    counts['no_predicted_protein_change']+=1;continue
                effect=type(effects)(altered).top_priority_effect()
                tumour=r.samples[plan['tumour_sample']];blood=r.samples[plan['normal_sample']]
                assert r.info['MN_EFFECTIVE_FILTER'][0]=='PASS'
                matches=caris.get(key,[])
                evidence=normal.assess(list(key),variant)
                row=dict(variant_id=r.id,normalized_allele=list(key),gene=effect.gene_name,
                    gene_id=getattr(effect,'gene_id',None),protein_change=getattr(effect,'short_description',None),
                    consequence=type(effect).__name__,transcript_id=getattr(effect,'transcript_id',None) or getattr(getattr(effect,'transcript',None),'id',None),
                    other_effects=[dict(gene=e.gene_name,consequence=type(e).__name__,protein_change=getattr(e,'short_description',None),transcript_id=getattr(e,'transcript_id',None)) for e in altered],
                    native_FILTER=list(r.filter),effective_allele_filter=r.info['MN_EFFECTIVE_FILTER'][0],
                    DNA_DP=tumour.get('DP'),DNA_AD=list(tumour.get('AD') or []),DNA_model_AF=scalar(tumour.get('AF')),
                    GERMQ=scalar(r.info.get('GERMQ')),native_population_AF_prior=10**-float(scalar(r.info['POPAF'])),
                    native_INFO=dict(r.info),native_tumour_sample=dict(tumour.items()),native_normal_sample=dict(blood.items()),
                    Caris_exact_matches=matches,Caris_flags=sorted({flag for ident in matches for flag in caris_meta[ident]['filter']}),
                    normal_evidence=evidence,independent_current_target_normal_verdict=focal_audit.get(key,{}).get('normal_verdict'))
                candidates.append(row)
            except Exception as error:errors.append(dict(variant_id=r.id,allele=key,error=repr(error)))
            if len(candidates)%25==0:save(out/'status.json',dict(phase='annotating_and_assessing_normal',native_PASS_seen=counts['all_unique_native_PASS'],coding=len(candidates),errors=len(errors)))
    save(out/'annotation-errors.json',errors)
    assert not errors,errors[:5]
    by_key={tuple(row['normalized_allele']):row for row in candidates}
    for chrom in chromosomes:
        with gzip.open(root/'paired/population'/chrom/'population-observations.jsonl.gz','rt') as stream:
            for line in stream:
                row=json.loads(line);key=tuple(row['normalized_allele'])
                if key in by_key:
                    assert 'population' not in by_key[key]
                    by_key[key]['population']=row['population']
    assert all('population' in row for row in candidates)
    # Compact genomes omit regional flags. Restore those from their exact pinned
    # original source for the new coding candidates, retaining the old policy.
    by_chrom=collections.defaultdict(list)
    for row in candidates:by_chrom[row['normalized_allele'][0]].append(row)
    source_metadata=[]
    for chrom,rows in sorted(by_chrom.items()):
        prior=json.loads((OLD/'population'/chrom/'complete.json').read_text());url=prior['source_url']
        with urllib.request.urlopen(urllib.request.Request(url,method='HEAD'),timeout=90) as response:
            assert response.headers.get('ETag')==prior['source_http_metadata']['ETag']
        original_index=OLD_REVIEW/'annotate/population-indices'/(chrom+'.genomes.tbi')
        index=out/(chrom+'.genomes.tbi')
        if original_index.exists():index.write_bytes(original_index.read_bytes())
        else:
            with urllib.request.urlopen(url+'.tbi',timeout=90) as response:index.write_bytes(response.read())
        source_metadata.append(dict(chromosome=chrom,url=url,etag=prior['source_http_metadata']['ETag'],index_sha256=sha(index)))
        with pysam.VariantFile(url,index_filename=str(index)) as genomes,pysam.VariantFile(OLD_REVIEW/'population'/chrom/'exomes-annotation.vcf.gz') as exomes:
            for row in rows:
                c,p,rf,al=row['normalized_allele']
                for cohort,source in [('genomes',genomes),('exomes',exomes)]:
                    observation=next(item for item in row['population'] if item['cohort']==cohort)
                    flags=set()
                    for record in source.fetch(c,max(0,p-101),p+len(rf)+100):
                        if record.pos<=p+len(rf)-1 and record.pos+len(record.ref)-1>=p:
                            flags.update(name for name in ['lcr','segdup'] if name in source.header.info and record.info.get(name))
                    observation['overlapping_population_site_regional_flags']=sorted(flags)
    for row in candidates:
        dna=[];population=[];context=[];rule=policy['DNA_screen']
        dp,ad,af=row['DNA_DP'],row['DNA_AD'],row['DNA_model_AF']
        for name,failed in [('LOW_DEPTH',dp is None or dp<rule['minimum_sample_depth']),
                            ('LOW_ALT_AD',len(ad)<2 or ad[1]<rule['minimum_alternate_AD']),
                            ('LOW_TUMOUR_AF',af is None or af<rule['minimum_model_tumour_AF']),
                            ('LOW_GERMQ',row['GERMQ'] is None or row['GERMQ']<rule['minimum_GERMQ'])]:
            if failed:dna.append(name)
        good=[];any_match=False
        for cohort in row['population']:
            if cohort['overlapping_population_site_regional_flags']:context.append('LOW_COMPLEXITY_OR_SEGMENTAL_DUPLICATION')
            for entry in cohort['records']:
                any_match=True
                if entry['FILTER']!=['PASS'] or entry.get('AF') is None or not entry.get('AN'):continue
                good.append(dict(cohort=cohort['cohort'],**entry))
                if entry['AF']>=policy['population_screen']['global_AF_exclusion_at_or_above']:population.append('COMMON_GLOBAL_'+cohort['cohort'].upper())
                if entry.get('AF_grpmax') is not None and entry['AF_grpmax']>=policy['population_screen']['ancestry_group_AF_exclusion_at_or_above']:population.append('COMMON_ANCESTRY_GROUP_'+cohort['cohort'].upper())
        chrom,pos,rf,al=row['normalized_allele']
        if chrom in policy['context_review']['sex_chromosomes']:context.append('SEX_CHROMOSOME_CONTEXT')
        if set(row['Caris_flags'])&set(policy['context_review']['exact_Caris_quality_flags']):context.append('CARIS_FP_OR_LOW_QUALITY_FLAG')
        if any_match and not good:context.append('INSUFFICIENT_OR_FILTERED_POPULATION_RECORD')
        if len(rf)>1 and len(al)>1 and not good:context.append('MNV_WITHOUT_EXACT_POPULATION_RECORD')
        if af is not None and af>=policy['context_review']['high_tumour_model_AF_at_or_above']:context.append('HIGH_AF_OR_PLOIDY_AMBIGUITY')
        maximum=max([x['AF'] for x in good],default=None)
        if maximum is not None:assert maximum<=row['native_population_AF_prior']*1.03
        normal_status=row['normal_evidence']['classification']
        if row['independent_current_target_normal_verdict']=='TARGET_ALLELE_PRESENT_IN_NORMAL' and normal_status!='PRESENT_IN_NORMAL':
            context.append('INDEPENDENT_NORMAL_GENOTYPE_DISAGREEMENT')
        primary=not(dna or population or context) and normal_status=='NORMAL_REFERENCE_SUPPORTED'
        row.update(DNA_screen_failures=sorted(set(dna)),population_exclusions=sorted(set(population)),manual_review_flags=sorted(set(context)),
                   population_status='QC_PASS_OBSERVATION' if good else 'QC_FILTERED_ONLY' if any_match else 'UNOBSERVED_NOT_MEASURED_ZERO',
                   QC_PASS_combined_global_AF=maximum,primary_DNA_eligible=primary,
                   RNA_review_eligible=not(dna or population or context) and normal_status!='PRESENT_IN_NORMAL',
                   somatic_status='SUPPORTED_PAIRED_SOMATIC' if primary else 'PRESENT_IN_NORMAL' if normal_status=='PRESENT_IN_NORMAL' else 'RESEARCH_REVIEW_REQUIRED',
                   evidence_tier='PAIRED_NORMAL_SUPPORTED' if primary else 'EXCLUDED_NORMAL_PRESENT' if normal_status=='PRESENT_IN_NORMAL' else 'REVIEW_NORMAL_EVIDENCE' if not(dna or population or context) else 'REVIEW_DNA_POPULATION_OR_CONTEXT')
    candidates.sort(key=lambda row:row['variant_id'])
    save(out/'variant-evidence.json',candidates);table(out/'variant-evidence.tsv',candidates)
    inputs=[row for row in candidates if row['RNA_review_eligible']]
    save(out/'rna-input.json',inputs);save(out/'population-source-metadata.json',source_metadata)
    receipt(out,[out/name for name in ['variant-evidence.json','variant-evidence.tsv','rna-input.json','annotation-errors.json','population-source-metadata.json']],
        counts=dict(counts),protein_changing_native_PASS=len(candidates),primary_DNA_eligible=sum(row['primary_DNA_eligible'] for row in candidates),
        normal_classifications=dict(collections.Counter(row['normal_evidence']['classification'] for row in candidates)),
        rna_review_variants=len(inputs),all_samples_named_explicitly=True,Caris_nominations=False,
        policy_sha256=sha(root/'policy.json'),worker_sha256=sha(__file__))


if __name__=='__main__':main()
