#!/usr/bin/env python3
"""Compare local RNA BAM and normal gVCF at common SNP positions; no network or upload.

The --normal-pileups input is a GATK GetPileupSummaries TSV with contig,
position, ref_count, alt_count, other_alt_count and allele_frequency columns.
It is not a Picard haplotype-map file. Output contains aggregate counts only;
--marker-output explicitly opts into a separate detailed genotype file.
"""
from pathlib import Path
import argparse,csv,hashlib,json,collections,time
import pysam


def file_sha256(path):
    digest=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(4*1024*1024),b''):
            digest.update(block)
    return digest.hexdigest()


def summarize_markers(markers):
    """Compute aggregate counts and prose; no dataset-specific result literals."""
    total=len(markers)
    matches=sum(c['threshold_genotypes_match'] for c in markers)
    uncertain=sum(c['RNA_called_genotype_by_threshold']=='uncertain' for c in markers)
    het_hom=sum(0 in c['normal_genotype'] and c['RNA_called_genotype_by_threshold'] in ('0/0','1/1') for c in markers)
    homo_alt=[c for c in markers if c['normal_genotype']==[1,1]]
    homo_alt_matches=sum(c['RNA_called_genotype_by_threshold']=='1/1' for c in homo_alt)
    opposite=sum(c['RNA_called_genotype_by_threshold']=='0/0' for c in homo_alt)
    other=total-matches-uncertain-het_hom
    return {
        'comparison_status':'BOUNDED_COMPARISON_COMPLETE' if total else 'NO_INFORMATIVE_MARKERS',
        'markers_retained':total,
        'exact_threshold_genotype_matches':matches,
        'RNA_uncertain_genotypes':uncertain,
        'normal_heterozygote_RNA_homozygote_count':het_hom,
        'other_nonmatching_genotype_calls':other,
        'normal_homozygous_ALT_marker_count':len(homo_alt),
        'normal_homozygous_ALT_RNA_match_count':homo_alt_matches,
        'normal_homozygous_ALT_RNA_homozygous_REF_conflicts':opposite,
        'opposite_homozygotes':opposite,
        'opposite_homozygote_eligible_marker_count':len(homo_alt),
        'opposite_homozygote_scope':'Only normal 1/1 versus RNA 0/0 among selected normal homozygous-ALT markers. Normal 0/0 markers are not selected, and normal 0/1 with homozygous RNA is reported separately.',
        'interpretation':(
            f'{matches} of {total} retained markers have exact threshold-genotype matches; '
            f'{uncertain} RNA calls are uncertain, {het_hom} normal heterozygotes have homozygous RNA calls, '
            f'and {other} other calls do not match. '
            f'{homo_alt_matches} of {len(homo_alt)} normal homozygous-ALT markers have matching RNA calls; '
            f'{opposite} have homozygous-REF RNA calls. '
            'The causes of nonmatching or uncertain calls are unresolved. '
            'No claim that all markers match, or that samples are authenticated to a named patient.'
        ),
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rna-bam',type=Path,required=True)
    parser.add_argument('--rna-index',type=Path,required=True)
    parser.add_argument('--normal-gvcf',type=Path,required=True)
    parser.add_argument('--normal-pileups','--fingerprint-map',dest='normal_pileups',type=Path,required=True,help='GATK GetPileupSummaries TSV; --fingerprint-map is a compatibility alias, not a Picard map')
    parser.add_argument('--output',type=Path,required=True,help='New aggregate JSON output')
    parser.add_argument('--marker-output',type=Path,help='Optional new private detailed genotype output')
    args=parser.parse_args()
    for path in (args.rna_bam,args.rna_index,args.normal_gvcf,args.normal_pileups):
        if not path.is_file():parser.error('A local input file is required: '+str(path))
    for path in (args.output,args.marker_output):
        if path and path.exists():parser.error('Choose a new output path: '+str(path))
    if args.marker_output and args.output.resolve()==args.marker_output.resolve():
        parser.error('Aggregate and detailed outputs must have different paths')
    RNA=args.rna_bam
    RNA_INDEX=args.rna_index
    NORMAL=args.normal_gvcf
    PILEUPS=args.normal_pileups
    common={}
    with PILEUPS.open() as f:
     for r in csv.DictReader((s for s in f if not s.startswith('#')),delimiter='\t'):
      if r['contig'] not in {f'chr{i}' for i in range(1,23)}:continue
      depth=sum(int(r[k]) for k in ('ref_count','alt_count','other_alt_count'))
      if depth>=20:common[r['contig'],int(r['position'])]=r
    candidates=[];start=time.time();records=0
    with pysam.VariantFile(str(NORMAL)) as v:
     for r in v:
      records+=1
      key=r.contig,r.pos
      if key not in common or len(r.ref)!=1 or r.ref not in 'ACGT':continue
      s=r.samples[0];gt=s.get('GT')
      if gt not in [(0,1),(1,0),(1,1)] or (s.get('DP') or 0)<20 or (s.get('GQ') or 0)<30:continue
      if not r.alts or len(r.alts[0])!=1 or r.alts[0] not in 'ACGT':continue
      if any(x not in ('PASS','.') for x in r.filter):continue
      candidates.append({'chromosome':r.contig,'position_1based':r.pos,'reference':r.ref,'alternate':r.alts[0],'normal_genotype':list(gt),'normal_DP':s['DP'],'normal_GQ':s['GQ'],'population_AF':float(common[key]['allele_frequency'])})
    print(json.dumps({'gvcf_records_scanned':records,'common_covered_positions':len(common),'eligible_candidates':len(candidates),'seconds':round(time.time()-start,1)}),flush=True)
    candidates.sort(key=lambda r:hashlib.sha256(f"{r['chromosome']}:{r['position_1based']}".encode()).hexdigest())
    chosen=[];perchr=collections.Counter();positions=collections.defaultdict(list);attempted=0
    with pysam.AlignmentFile(str(RNA),'rb',index_filename=str(RNA_INDEX)) as b:
     for c in candidates:
      chrom=c['chromosome'];pos=c['position_1based']-1
      if perchr[chrom]>=5 or any(abs(pos-x)<1000000 for x in positions[chrom]):continue
      attempted+=1
      calls=collections.defaultdict(set)
      for r in b.fetch(chrom,pos,pos+1):
       if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail or r.mapping_quality<30:continue
       if r.has_tag('NH') and r.get_tag('NH')!=1:continue
       for q,ref in r.get_aligned_pairs(matches_only=True):
        if ref==pos:
         if r.query_qualities is not None and r.query_qualities[q]>=30:calls[r.query_name].add(r.query_sequence[q])
         break
      counts=collections.Counter(next(iter(v)) for v in calls.values() if len(v)==1)
      focal=counts[c['reference']]+counts[c['alternate']]
      if focal<10:continue
      af=counts[c['alternate']]/focal
      called='1/1' if af>=0.95 else '0/0' if af<=0.05 else '0/1' if 0.2<=af<=0.8 and min(counts[c['reference']],counts[c['alternate']])>=2 else 'uncertain'
      normal='0/1' if 0 in c['normal_genotype'] else '1/1'
      c.update({'RNA_fragment_counts':dict(counts),'RNA_focal_fragments':focal,'RNA_alternate_fraction':af,'RNA_called_genotype_by_threshold':called,'threshold_genotypes_match':called==normal,'opposite_homozygote':normal=='1/1' and called=='0/0','conflicting_mate_names_removed':sum(len(x)>1 for x in calls.values())})
      chosen.append(c);perchr[chrom]+=1;positions[chrom].append(pos)
    summary={
        'normal_gVCF_sha256':hashlib.sha256(NORMAL.read_bytes()).hexdigest(),
        'normal_pileups_sha256':hashlib.sha256(PILEUPS.read_bytes()).hexdigest(),
        'RNA_BAM_sha256':file_sha256(RNA),
        'RNA_index_sha256':file_sha256(RNA_INDEX),
        'gvcf_records_scanned':records,
        'eligible_candidate_count':len(candidates),
        'RNA_candidate_positions_attempted':attempted,
        'markers_retained':len(chosen),
        'chromosome_counts':dict(perchr),
        'exact_threshold_genotype_matches':sum(c['threshold_genotypes_match'] for c in chosen),
        'RNA_uncertain_genotypes':sum(c['RNA_called_genotype_by_threshold']=='uncertain' for c in chosen),
        'opposite_homozygotes':sum(c['opposite_homozygote'] for c in chosen),
        'normal_genotype_counts':dict(collections.Counter('0/1' if 0 in c['normal_genotype'] else '1/1' for c in chosen)),
        'scope':'Bounded comparison of the supplied RNA BAM against the supplied normal coding gVCF at common-marker positions. Not Picard CrosscheckFingerprints, not a calibrated identity probability, not named-patient authentication, and not exhaustive swap/contamination testing.',
        'filters':'Normal DP>=20 GQ>=30, single-base REF and first ALT, diploid GT 0/1, 1/0 or 1/1, with native PASS/unfiltered; RNA primary nonduplicate non-QCfail MAPQ>=30 NH=1 if provided BQ>=30, consistent fragment/read-name calls; >=10 REF+ALT fragments. Deterministic SHA-256 ordering; max five per autosome, minimum 1 Mb between retained loci.',
        'elapsed_seconds':round(time.time()-start,2),
    }
    summary.update(summarize_markers(chosen))
    summary['software']={'python':__import__('sys').version.split()[0],'pysam':pysam.__version__}
    summary['script_sha256']=file_sha256(Path(__file__))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(summary,indent=2)+'\n')
    if args.marker_output:
        args.marker_output.parent.mkdir(parents=True,exist_ok=True)
        args.marker_output.write_text(json.dumps({'scope':'Private detailed marker genotypes; no upload performed.','markers':chosen},indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    main()
