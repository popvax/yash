"""Annotate segmented native-PASS coding alleles and apply the research policy."""
from pathlib import Path
import bisect
import collections
import concurrent.futures
import csv
import gzip
import json
import math
import sys
import time
import urllib.request

import pysam
import varcode
from pyensembl import Genome
from common import DISCOVERY, sha, save, receipt, scalar

ROOT = Path(sys.argv[1])
OUT = ROOT / 'annotate'
OUT.mkdir(parents=True, exist_ok=True)
policy = json.loads((ROOT / 'policy.json').read_text())
cfg = json.loads((ROOT / 'analysis-genome-config.json').read_text())
genome = Genome(**cfg)
plan = json.loads((ROOT / 'calling-plan.json').read_text())
assert json.loads((ROOT / 'gather/complete.json').read_text())['status'] == 'COMPLETED'
started = time.time()


def table(name, values):
    fields = list(dict.fromkeys(k for r in values for k in r)) or ['variant_id', 'status']
    with (OUT / name).open('w') as f:
        writer = csv.DictWriter(f, delimiter='\t', fieldnames=fields)
        writer.writeheader()
        for r in values:
            writer.writerow({k: json.dumps(v) if isinstance(v, (dict, list, tuple)) else v for k, v in r.items()})


raw = collections.defaultdict(list)
chromosomes = {s['chrom'] for s in plan['shards']}
with open(cfg['gtf_path_or_url']) as f:
    for line in f:
        if line.startswith('#'):
            continue
        x = line.rstrip().split('\t')
        if len(x) != 9 or x[2] not in ['CDS', 'start_codon', 'stop_codon']:
            continue
        chrom = x[0] if x[0].startswith('chr') else 'chr' + x[0]
        if chrom in chromosomes:
            raw[chrom].append((max(1, int(x[3]) - 2), int(x[4]) + 2))
intervals, starts = {}, {}
for chrom, rr in raw.items():
    merged = []
    for a, b in sorted(rr):
        if merged and a <= merged[-1][1] + 1:
            merged[-1] = merged[-1][0], max(b, merged[-1][1])
        else:
            merged.append((a, b))
    intervals[chrom] = merged
    starts[chrom] = [a for a, b in merged]
del raw
caris_meta = json.loads((DISCOVERY / 'gather/caris-comparison-metadata.json').read_text())
caris_keys = collections.defaultdict(list)
with pysam.VariantFile(DISCOVERY / 'gather/caris-comparison-normalized.vcf.gz') as src:
    for r in src:
        caris_keys[r.contig, r.pos, r.ref, r.alts[0]].append(r.id)
model_evidence = {}
with gzip.open(ROOT / 'gather/model-evidence.jsonl.gz', 'rt') as f:
    for line in f:
        r = json.loads(line)
        model_evidence[r['variant_id']] = r
counts = collections.Counter()
candidates, errors = [], []
with pysam.VariantFile(ROOT / 'gather/pass.vcf.gz') as src:
    for r in src:
        counts['all_native_PASS'] += 1
        c, p, ref, alt = r.contig, r.pos, r.ref, r.alts[0]
        index = bisect.bisect_right(starts.get(c, []), p + len(ref) - 1) - 1
        if index < 0 or intervals[c][index][1] < p:
            continue
        counts['overlaps_CDS_or_2bp_boundary'] += 1
        try:
            variant = varcode.Variant(c, p, ref, alt, genome=genome, convert_ucsc_contig_names=True)
            effects = variant.effects()
            altered = [e for e in effects if getattr(e, 'modifies_protein_sequence', False)]
            if not altered:
                counts['no_predicted_protein_change'] += 1
                continue
            e = type(effects)(altered).top_priority_effect()
            sample = next(iter(r.samples.values()))
            ad = sample.get('AD')
            matches = caris_keys.get((c, p, ref, alt), [])
            metadata = model_evidence[r.id]
            assert metadata['normalized_allele'] == [c, p, ref, alt]
            candidates.append(dict(
                variant_id=r.id, normalized_allele=[c, p, ref, alt], gene=e.gene_name,
                gene_id=getattr(e, 'gene_id', None), protein_change=getattr(e, 'short_description', None),
                consequence=type(e).__name__, transcript_id=getattr(e, 'transcript_id', None) or getattr(getattr(e, 'transcript', None), 'id', None),
                other_effects=[dict(gene=x.gene_name, consequence=type(x).__name__, protein_change=getattr(x, 'short_description', None), transcript_id=getattr(x, 'transcript_id', None)) for x in altered],
                native_FILTER=list(r.filter), effective_allele_filter=r.info['V3_EFFECTIVE_FILTER'][0],
                DNA_DP=sample.get('DP'), DNA_AD=list(ad) if ad else None, DNA_model_AF=scalar(sample.get('AF')),
                GERMQ=r.info.get('GERMQ'), MBQ=list(r.info.get('MBQ', [])), MMQ=list(r.info.get('MMQ', [])),
                native_population_AF_prior=10 ** -float(scalar(r.info['POPAF'])), native_INFO=dict(r.info), native_sample=dict(sample.items()),
                Caris_exact_matches=matches, Caris_flags=sorted({v for cid in matches for v in caris_meta[cid]['filter']}),
                model_evidence=metadata['models']))
        except Exception as error:
            errors.append(dict(variant_id=r.id, allele=[c, p, ref, alt], error=repr(error)))
        if len(candidates) % 100 == 0:
            save(OUT / 'status.json', dict(phase='annotating_native_PASS', native_PASS_seen=counts['all_native_PASS'], coding=len(candidates)))
save(OUT / 'annotation-errors.json', errors)
assert not errors, errors[:5]
by_key = {tuple(r['normalized_allele']): r for r in candidates}
for chrom in chromosomes:
    with gzip.open(ROOT / 'population' / chrom / 'population-observations.jsonl.gz', 'rt') as src:
        for line in src:
            r = json.loads(line)
            k = tuple(r['normalized_allele'])
            if k in by_key:
                assert 'population' not in by_key[k], ('Duplicate normalized population record', k)
                by_key[k]['population'] = r['population']
assert all('population' in r for r in candidates)
by_chrom = collections.defaultdict(list)
for r in candidates:
    by_chrom[r['normalized_allele'][0]].append(r)


def equivalent(fasta, c, p, ref, alt, r):
    if (p, ref, alt) == (r.pos, r.ref, r.alts[0]):
        return True
    if len(ref) == len(alt) or len(r.ref) == len(r.alts[0]):
        return False
    lo, hi = min(p, r.pos) - 1, max(p + len(ref), r.pos + len(r.ref)) + 1
    sequence = fasta.fetch(c, lo, hi).upper()
    def apply(pos, rf, al):
        off = pos - 1 - lo
        return sequence[:off] + al + sequence[off + len(rf):] if sequence[off:off + len(rf)] == rf else None
    left = apply(p, ref, alt)
    return left is not None and left == apply(r.pos, r.ref, r.alts[0])


def enrich(item):
    # The reused genome compact files omit lcr/segdup. Recover those fields
    # from the exact original source, after confirming the pinned source ETag.
    chrom, records = item
    source = json.loads((DISCOVERY / 'population' / chrom / 'complete.json').read_text())
    url = source['source_url']
    with urllib.request.urlopen(urllib.request.Request(url, method='HEAD'), timeout=90) as response:
        assert response.headers.get('ETag') == source['source_http_metadata']['ETag']
    indices = OUT / 'population-indices'
    indices.mkdir(exist_ok=True)
    index = indices / (chrom + '.genomes.tbi')
    with urllib.request.urlopen(url + '.tbi', timeout=90) as response:
        index.write_bytes(response.read())
    meta = dict(chromosome=chrom, genome_url=url, genome_ETag=source['source_http_metadata']['ETag'], index_sha256=sha(index))
    fasta = pysam.FastaFile(plan['reference'])
    try:
        with pysam.VariantFile(url, index_filename=str(index)) as genomes, pysam.VariantFile(ROOT / 'population' / chrom / 'exomes-annotation.vcf.gz') as exomes:
            for record in records:
                c, p, ref, alt = record['normalized_allele']
                for cohort, handle in [('genomes', genomes), ('exomes', exomes)]:
                    original = next(x for x in record['population'] if x['cohort'] == cohort)
                    flags, matches = set(), []
                    for r in handle.fetch(c, max(0, p - 101), p + max(1, len(ref)) + 100):
                        assert len(r.alts) == 1
                        if r.pos <= p + len(ref) - 1 and r.pos + len(r.ref) - 1 >= p:
                            flags.update(k for k in ['lcr', 'segdup'] if k in handle.header.info and r.info.get(k))
                        if equivalent(fasta, c, p, ref, alt, r):
                            matches.append(dict(position=r.pos, ref=r.ref, alt=r.alts[0], FILTER=list(r.filter),
                                **{k: scalar(r.info.get(k)) if k in handle.header.info else None for k in ['AF', 'AC', 'AN', 'grpmax', 'AF_grpmax', 'AN_grpmax', 'AF_sas', 'AN_sas', 'AC_sas', 'lcr', 'segdup', 'non_par']}))
                    sig = lambda x: (x['position'], x['ref'], x['alt'], tuple(x['FILTER']))
                    assert {sig(x) for x in original['records']} == {sig(x) for x in matches}, ('Population matching changed', record['variant_id'], cohort)
                    original['records'] = matches
                    original['overlapping_population_site_regional_flags'] = sorted(flags)
        return meta
    finally:
        fasta.close()


save(OUT / 'status.json', dict(phase='enriching_population_context', protein_changing=len(candidates), chromosomes=len(by_chrom)))
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    source_metadata = list(pool.map(enrich, sorted(by_chrom.items())))
for r in candidates:
    dna, context, exclusions = [], [], []
    rule = policy['DNA_screen']
    dp, ad, af = r['DNA_DP'], r['DNA_AD'], r['DNA_model_AF']
    for name, fails in [
        ('LOW_DEPTH', dp is None or dp < rule['minimum_sample_depth']),
        ('LOW_ALT_AD', not ad or ad[1] < rule['minimum_alternate_AD']),
        ('LOW_TUMOUR_AF', af is None or af < rule['minimum_model_tumour_AF']),
        ('LOW_GERMQ', r['GERMQ'] is None or r['GERMQ'] < rule['minimum_GERMQ'])]:
        if fails:
            dna.append(name)
    good, any_match = [], False
    for cohort in r['population']:
        if cohort['overlapping_population_site_regional_flags']:
            context.append('LOW_COMPLEXITY_OR_SEGMENTAL_DUPLICATION')
        for x in cohort['records']:
            any_match = True
            if x['FILTER'] != ['PASS'] or x['AF'] is None or not x['AN']:
                continue
            good.append(dict(cohort=cohort['cohort'], **x))
            if x['AF'] >= policy['population_screen']['global_AF_exclusion_at_or_above']:
                exclusions.append('COMMON_GLOBAL_' + cohort['cohort'].upper())
            if (x['AF_grpmax'] is not None and x['AF_grpmax']
                    >= policy['population_screen']['ancestry_group_AF_exclusion_at_or_above']):
                exclusions.append('COMMON_ANCESTRY_GROUP_' + cohort['cohort'].upper())
            if x.get('lcr') or x.get('segdup'):
                context.append('LOW_COMPLEXITY_OR_SEGMENTAL_DUPLICATION')
    c, p, ref, alt = r['normalized_allele']
    if c in policy['context_review']['sex_chromosomes']:
        context.append('SEX_CHROMOSOME_CONTEXT')
    if set(r['Caris_flags']) & set(policy['context_review']['exact_Caris_quality_flags']):
        context.append('CARIS_FP_OR_LOW_QUALITY_FLAG')
    if any_match and not good:
        context.append('INSUFFICIENT_OR_FILTERED_POPULATION_RECORD')
    if len(ref) > 1 and len(alt) > 1 and not good:
        context.append('MNV_WITHOUT_EXACT_POPULATION_RECORD')
    if af is not None and af >= policy['context_review']['high_tumour_model_AF_at_or_above']:
        context.append('HIGH_AF_OR_PLOIDY_AMBIGUITY')
    global_af = max([x['AF'] for x in good], default=None)
    if global_af is not None:
        # GATK rounds INFO/POPAF on export. Allow that small numeric difference,
        # but stop if the updated resource still materially understates AF.
        assert global_af <= r['native_population_AF_prior'] * 1.03, ('Updated prior missing a reviewed observation', r['variant_id'])
    model_eligible = {}
    common_failures = set(dna) - {'LOW_GERMQ'}
    for model, data in r['model_evidence'].items():
        model_eligible[model] = bool(data['effective_filter'] == 'PASS' and (data['GERMQ'] or 0) >= rule['minimum_GERMQ'] and not (common_failures or context or exclusions))
    r.update(DNA_screen_failures=sorted(set(dna)), population_exclusions=sorted(set(exclusions)),
             manual_review_flags=sorted(set(context)), population_status='QC_PASS_OBSERVATION' if good else 'QC_FILTERED_ONLY' if any_match else 'UNOBSERVED_NOT_MEASURED_ZERO',
             somatic_status='UNCONFIRMED_NO_MATCHED_NORMAL', primary_DNA_eligible=any(model_eligible.values()),
             primary_DNA_eligible_by_model=model_eligible, QC_PASS_combined_global_AF=global_af,
             RNA_review_eligible=model_eligible['segmented'])
    r['evidence_tier'] = 'EXCLUDED_POPULATION_COMMON' if exclusions else 'REVIEW_DNA_QUALITY' if dna else 'REVIEW_CONTEXT' if context else 'PUTATIVE_TUMOUR_ONLY'
    assert r['primary_DNA_eligible'] == (not (dna or exclusions or context))
candidates.sort(key=lambda x: x['variant_id'])
save(OUT / 'variant-evidence.json', candidates)
table('variant-evidence.tsv', candidates)
save(OUT / 'rna-input.json', [r for r in candidates if r['RNA_review_eligible']])
save(OUT / 'population-source-metadata.json', source_metadata)
counts.update(protein_changing=len(candidates), primary_DNA_eligible=sum(r['primary_DNA_eligible'] for r in candidates),
              RNA_review_eligible=sum(r['RNA_review_eligible'] for r in candidates))
receipt(OUT, [OUT / n for n in ['variant-evidence.json', 'variant-evidence.tsv', 'rna-input.json', 'population-source-metadata.json', 'annotation-errors.json']],
        counts=dict(counts), tiers=dict(collections.Counter(r['evidence_tier'] for r in candidates)),
        primary_DNA_by_model={m: sum(r['primary_DNA_eligible_by_model'][m] for r in candidates) for m in ['segmented']},
        annotation='Varcode / Ensembl 103 on segmented native-PASS coding calls',
        elapsed_seconds=time.time() - started, worker_sha256=sha(__file__), errors=0,
        input_sha256=dict(union_PASS=sha(ROOT / 'gather/pass.vcf.gz'), policy=sha(ROOT / 'policy.json'), GTF=sha(cfg['gtf_path_or_url'])))
