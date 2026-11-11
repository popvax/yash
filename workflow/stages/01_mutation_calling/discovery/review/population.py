"""Prepare local exomes and update priors for every emitted allele on one chromosome."""
from pathlib import Path
import bisect
import collections
import concurrent.futures
import gzip
import hashlib
import json
import math
import re
import shutil
import sys
import tempfile
import time
import urllib.request

import pysam
from common import DISCOVERY, sha, save, run, receipt, scalar

ROOT = Path(sys.argv[1])
CHROM = sys.argv[2]
OUT = ROOT / 'population' / CHROM
OUT.mkdir(parents=True, exist_ok=True)
assert json.loads((ROOT / 'setup/complete.json').read_text())['status'] == 'COMPLETED'
probe = {x['chrom']: x for x in json.loads((ROOT / 'exome-source-probe.json').read_text())}[CHROM]
URL = probe['url']
KEEP = {'AF', 'AC', 'AN', 'grpmax', 'AF_grpmax', 'AN_grpmax', 'AC_grpmax', 'AF_sas', 'AN_sas', 'AC_sas', 'lcr', 'segdup', 'non_par'}
genomes = DISCOVERY / 'population' / CHROM / 'annotation.vcf.gz'
genome_receipt = json.loads((genomes.parent / 'complete.json').read_text())
assert sha(genomes) == genome_receipt['outputs'][genomes.name]['sha256']
exomes = OUT / 'exomes-annotation.vcf.gz'
source_receipt = OUT / 'exomes-source.json'


def download_source():
    with urllib.request.urlopen(urllib.request.Request(URL, method='HEAD'), timeout=90) as response:
        size = int(response.headers['Content-Length'])
        etag = response.headers.get('ETag')
    assert size == probe['bytes'] and etag == probe['etag'], 'Source changed since preflight'
    full = OUT / 'exomes-source.vcf.bgz'
    if source_receipt.exists():
        prior = json.loads(source_receipt.read_text())
        assert prior['source_bytes'] == size and prior['etag'] == etag
        assert full.stat().st_size == size and sha(full) == prior['source_sha256']
        return prior
    save(OUT / 'status.json', dict(phase='downloading_exome_ranges', source_bytes=size, url=URL))
    parts = 16
    with tempfile.TemporaryDirectory(prefix='all-exomes-') as temp:
        def part(i):
            lo, hi = size * i // parts, size * (i + 1) // parts - 1
            path = Path(temp) / str(i)
            for attempt in range(4):
                try:
                    req = urllib.request.Request(URL, headers={'Range': f'bytes={lo}-{hi}', 'If-Range': etag})
                    with urllib.request.urlopen(req, timeout=90) as src, path.open('wb') as dst:
                        assert src.status == 206 and src.headers.get('Content-Range') == f'bytes {lo}-{hi}/{size}'
                        assert src.headers.get('ETag') == etag
                        shutil.copyfileobj(src, dst, 1024 * 1024)
                    assert path.stat().st_size == hi - lo + 1
                    return path
                except Exception:
                    if attempt == 3:
                        raise
                    time.sleep(2)
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
            paths = list(pool.map(part, range(parts)))
        digest = hashlib.sha256()
        partial = OUT / 'exomes-source.partial.vcf.bgz'
        with partial.open('wb') as dst:
            for p in paths:
                with p.open('rb') as src:
                    for chunk in iter(lambda: src.read(1024 * 1024), b''):
                        digest.update(chunk)
                        dst.write(chunk)
        assert partial.stat().st_size == size
        partial.replace(full)
    info = dict(url=URL, source_bytes=size, etag=etag, source_sha256=digest.hexdigest(),
                transport='16 verified byte ranges, four connections; concatenated in source order')
    save(source_receipt, info)
    return info


source = download_source()
if not exomes.exists():
    with pysam.VariantFile(OUT / 'exomes-source.vcf.bgz') as f:
        fields = set(f.header.info)
    assert {'AF', 'AC', 'AN', 'AF_grpmax', 'AN_grpmax'} <= fields
    remove = ['ID', 'QUAL'] + ['INFO/' + k for k in sorted(fields - KEEP)]
    partial = OUT / 'exomes-annotation.partial.vcf.gz'
    run(OUT, 'extract-exome-fields', ['bcftools', 'annotate', '-x', ','.join(remove), '-Oz', '-o', partial, OUT / 'exomes-source.vcf.bgz'])
    partial.replace(exomes)
run(OUT, 'index-exomes', ['bcftools', 'index', '-t', '-f', exomes])
plan = json.loads((ROOT / 'calling-plan.json').read_text())
fasta = pysam.FastaFile(plan['reference'])
alleles = [json.loads(line) for line in (OUT / 'alleles.jsonl').open()]
bykey = collections.defaultdict(list)
for a in alleles:
    bykey[tuple(a['normalized_allele'])].append(a)
observations = {k: {'genomes': [], 'exomes': []} for k in bykey}
positions = sorted({k[1] for k in bykey})
keys_by_position = collections.defaultdict(list)
intervals = []
for k in bykey:
    keys_by_position[k[1]].append(k)
for p in positions:
    end = max(p + len(k[2]) + 100 for k in keys_by_position[p])
    start = max(0, p - 101)
    if intervals and start <= intervals[-1][1]:
        intervals[-1][1] = max(intervals[-1][1], end)
    else:
        intervals.append([start, end])


def equivalent(key, rec):
    chrom, pos, ref, alt = key
    if (pos, ref, alt) == (rec.pos, rec.ref, rec.alts[0]):
        return True
    if len(ref) == len(alt) or len(rec.ref) == len(rec.alts[0]):
        return False
    if set(ref + alt + rec.ref + rec.alts[0]) - set('ACGTN'):
        return False
    lo = min(pos, rec.pos) - 1
    hi = max(pos + len(ref), rec.pos + len(rec.ref)) + 1
    sequence = fasta.fetch(chrom, lo, hi).upper()
    def apply(p, r, a):
        off = p - 1 - lo
        return sequence[:off] + a + sequence[off + len(r):] if sequence[off:off + len(r)] == r else None
    left = apply(pos, ref, alt)
    return left is not None and left == apply(rec.pos, rec.ref, rec.alts[0])


save(OUT / 'status.json', dict(phase='matching_all_emitted_alleles', alleles=len(alleles), unique_alleles=len(bykey), intervals=len(intervals)))
for cohort, path in [('genomes', genomes), ('exomes', exomes)]:
    with pysam.VariantFile(path) as src:
        available_fields = set(src.header.info)
        seen = set()
        for index, (lo, hi) in enumerate(intervals):
            for r in src.fetch(CHROM, lo, hi):
                assert len(r.alts) == 1
                signature = (r.pos, r.ref, r.alts[0])
                if signature in seen:
                    continue
                seen.add(signature)
                left = bisect.bisect_left(positions, max(1, r.pos - 100))
                right = bisect.bisect_right(positions, r.pos + len(r.ref) + 100)
                for p in positions[left:right]:
                    for key in keys_by_position[p]:
                        if equivalent(key, r):
                            record = dict(position=r.pos, ref=r.ref, alt=r.alts[0], FILTER=list(r.filter),
                                          **{k: scalar(r.info.get(k)) if k in available_fields else None for k in sorted(KEEP)})
                            observations[key][cohort].append(record)
            if index % 1000 == 0:
                save(OUT / 'status.json', dict(phase='matching_' + cohort, intervals_completed=index, intervals=len(intervals), alleles=len(alleles)))
    print('MATCHED', cohort, CHROM, sum(bool(x[cohort]) for x in observations.values()), flush=True)

priors = {}
counts = collections.Counter()
with gzip.open(OUT / 'population-observations.jsonl.gz', 'wt') as f:
    for key, maps in observations.items():
        good = [r for cohort in maps.values() for r in cohort if r['FILTER'] == ['PASS'] and r.get('AF') is not None and r.get('AN')]
        global_af = max([r['AF'] for r in good], default=None)
        counts['matched_QC_PASS'] += bool(good)
        for original in bykey[key]:
            native = 10 ** -original['original_POPAF']
            updated = max(native, global_af if global_af is not None else native)
            assert 0 < updated <= 1
            loc = (original['raw_id'], original['raw_alt_index'])
            assert loc not in priors
            priors[loc] = -math.log10(updated)
            counts['prior_increased'] += updated > native * 1.01
        f.write(json.dumps(dict(normalized_allele=key, population=[dict(cohort=c, status='MATCHED' if r else 'NOT_FOUND', records=r) for c, r in maps.items()])) + '\n')
updated = OUT / 'updated.vcf.gz'
seen = set()
with pysam.VariantFile(OUT / 'raw.vcf.gz') as src:
    with pysam.VariantFile(updated, 'wz', header=src.header) as dst:
        for r in src:
            old = list(r.info['POPAF'])
            for ai in range(len(r.alts)):
                location = (r.id, ai)
                assert location in priors and location not in seen
                seen.add(location)
                old[ai] = priors[location]
            r.info['POPAF'] = tuple(old)
            dst.write(r)
assert seen == set(priors) and len(seen) == len(alleles)
run(OUT, 'index-updated', ['bcftools', 'index', '-t', '-f', updated])
receipt(OUT, [updated, Path(str(updated) + '.tbi'), OUT / 'population-observations.jsonl.gz', source_receipt],
        chromosome=CHROM, all_emitted_alleles=len(alleles), unique_normalized_alleles=len(bykey), counts=dict(counts),
        genome_source_receipt_sha256=sha(genomes.parent / 'complete.json'), genome_compact_sha256=sha(genomes),
        exome_source=source, exome_compact_sha256=sha(exomes), worker_sha256=sha(__file__),
        matching='Exact alleles or verified equal local haplotypes for indels within a 100-base flank',
        population_policy='max(native, QC-PASS genome global AF, QC-PASS exome global AF); absent observations are not measured zero')
