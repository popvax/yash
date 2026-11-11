"""Declare independent primary-nuclear calling scope and read-weighted interval groups.

No Caris VCF is read. Interval allocations use BAM-index mapped-read counts by
chromosome; boundaries within each chromosome split its full span evenly.
"""
import argparse,hashlib,json,math
from pathlib import Path

def main():
 a=argparse.ArgumentParser(description=__doc__)
 a.add_argument('--preflight',type=Path,required=True);a.add_argument('--output',type=Path,required=True)
 a.add_argument('--shards',type=int,default=128);x=a.parse_args();p=json.loads(x.preflight.read_text())
 assert p['bam_reference_agreement'] and p['bam_bytes']==13260535753
 order=['chr'+str(c) for c in list(range(1,23))+['X','Y']]
 stats={r['contig']:r for r in p['bam_index']};weights={c:stats[c]['mapped'] for c in order}
 assert all(n>0 for n in weights.values()) and x.shards>=len(order)
 extra=x.shards-len(order);total=sum(weights.values());raw={c:extra*weights[c]/total for c in order}
 sizes={c:1+math.floor(raw[c]) for c in order}
 for c in sorted(order,key=lambda c:raw[c]-math.floor(raw[c]),reverse=True)[:x.shards-sum(sizes.values())]:sizes[c]+=1
 shards=[]
 for c in order:
  n=sizes[c];length=stats[c]['length'];last=0
  for i in range(n):
   start=length*i//n+1;end=length*(i+1)//n
   assert start==last+1 and end>=start;last=end
   shards.append(dict(id=f's{len(shards):03d}_{c}',chrom=c,start=start,end=end,interval=f'{c}:{start}-{end}',estimated_read_weight=weights[c]/n))
  assert last==length
 assert len(shards)==x.shards
 plan=dict(schema_version=1,sample='TN25-112103',gatk_version='4.6.2.0',
 bam='/data/inputs/TN25-112103/DNA_TN25-112103.bam',bam_sha256='4f12c8dbca50e26eda82bbc05166286797231227323e27b03f113906e19076e5',
 reference='/data/results/TN25-112103/mutect2-tumour-only-20261029/resources/reference.fa',
 reference_sha256='c3e61ad6d770db22a6f27fdbcd3dba93d0d80f649018c8887ec803f12ea90450',
 scope='Entire primary nuclear chromosomes, independently discovered. chrM and nonprimary/unplaced contigs excluded.',
 capture_BED=None,reference_caveat='Original hg38.pU2AF1_Y_PAR_masked.noHap.fa alignment reference is unavailable. Compatible GRCh38 is used; original alignment-mask effects persist.',
 reprocessing='Uses supplied coordinate-sorted, Sentieon deduplicated/realigned BAM. No unverified second deduplication or BQSR.',
 matched_normal=None,panel_of_normals=None,given_alleles=None,Caris_used_to_define_intervals=False,
 germline_resource='gnomAD v4.1 genomes, source PASS only; chr-specific files',
 filtering=dict(threshold_strategy='FALSE_DISCOVERY_RATE',model_false_discovery_rate=.05,scope='Joint sample after gathering all shard calls/statistics',empirical_FDR_established=False),
 population_review=dict(global_AF_exclusion_at_or_above=.001,ancestry_group_AF_exclusion_at_or_above=.001,QC_required='PASS',absence_policy='Unobserved/unknown, never a measured AF of zero'),
 shards=shards,excluded_contigs=[r for r in p['bam_index'] if r['contig'] not in order],
 preflight_sha256=hashlib.sha256(x.preflight.read_bytes()).hexdigest())
 x.output.write_text(json.dumps(plan,indent=2)+'\n')
 print(json.dumps(dict(shards=len(shards),concurrent_containers=32,cpu_per_container=2,shards_per_chromosome=sizes),indent=2))
if __name__=='__main__':main()
