"""Verify and map EVERY original ALT to a normalized allele, before filtering."""
from pathlib import Path
import collections
import json
import shutil
import sys
import pysam
from common import DISCOVERY, sha, save, run, receipt

ROOT = Path(sys.argv[1])
OUT = ROOT / 'setup'
OUT.mkdir(parents=True, exist_ok=True)
original = DISCOVERY / 'gather/unfiltered.vcf.gz'
previous = json.loads((DISCOVERY / 'gather/complete.json').read_text())
assert sha(original) == previous['outputs']['unfiltered.vcf.gz']['sha256']
assert previous['Caris_nominations'] is False
plan = json.loads((DISCOVERY / 'calling-plan.json').read_text())
chromosomes = ['chr' + str(c) for c in list(range(1, 23)) + ['X', 'Y']]
assert {s['chrom'] for s in plan['shards']} == set(chromosomes)
for name in ['case.json', 'analysis-genome-config.json', 'analysis-environment-verified.json', 'calling-plan.json', 'affinity-model-manifest.json']:
    shutil.copyfile(DISCOVERY / name, ROOT / name)
for name in ['unfiltered.vcf.gz.stats', 'orientation-model.tar.gz']:
    assert sha(DISCOVERY / 'gather' / name) == previous['outputs'][name]['sha256']
support = ROOT / 'support'
support.mkdir(exist_ok=True)
for name in ['contamination.table', 'allele-fraction-segments.table']:
    shutil.copyfile(DISCOVERY / 'support' / name, support / name)
    assert sha(support / name) == sha(DISCOVERY / 'support' / name)
policy = json.loads((DISCOVERY / 'policy.json').read_text())
for key in ['minimum_alternate_median_mapping_quality', 'minimum_alternate_median_base_quality']:
    policy['DNA_screen'].pop(key, None)
policy['DNA_screen']['mapping_and_base_quality'] = 'Native GATK FILTER only; no additional median MQ/BQ cutoff'
policy['population_screen'].pop('minimum_group_AN', None)
policy['population_screen']['ancestry_group_sampling'] = 'No minimum group allele count; apply the 0.1% cutoff to recorded QC-PASS grpmax AF regardless of group size'
policy['manual_review_flags'] = [flag for flag in policy['manual_review_flags'] if flag != 'SMALL_GROUP_FREQUENCY_REVIEW']
policy['RNA_screen'] = 'Native Isovar fragment-ID unions (read names), MQ>=20, no secondary/supplementary/duplicates, NH=1 when present; no additional three-fragment cutoff'
save(ROOT / 'policy.json', policy)
tagged = OUT / 'tagged-original.vcf.gz'
writers = {}
count = 0
try:
    with pysam.VariantFile(original) as src:
        h = src.header.copy()
        h.info.add('V3_ALT_INDEX', 'A', 'Integer', 'Original ALT index before normalized population lookup')
        for c in chromosomes:
            d = ROOT / 'population' / c
            d.mkdir(parents=True, exist_ok=True)
            writers[c] = pysam.VariantFile(d / 'raw.vcf.gz', 'wz', header=h)
        with pysam.VariantFile(tagged, 'wz', header=h) as dst:
            for index, r in enumerate(src):
                assert r.contig in writers and len(r.info['POPAF']) == len(r.alts)
                r.translate(h)
                r.id = 'all_' + str(index)
                r.info['V3_ALT_INDEX'] = tuple(range(len(r.alts)))
                dst.write(r)
                writers[r.contig].write(r)
                count += 1
finally:
    for w in writers.values():
        w.close()
normalized = OUT / 'normalized-all.vcf.gz'
command = run(OUT, 'normalize-all', ['bcftools', 'norm', '-f', plan['reference'], '-m', '-any', '-Oz', '-o', normalized, tagged])
maps = {}
counts = collections.Counter()
try:
    for c in chromosomes:
        maps[c] = (ROOT / 'population' / c / 'alleles.jsonl').open('w')
    with pysam.VariantFile(normalized) as src:
        for r in src:
            assert len(r.alts) == 1
            entry = dict(raw_id=r.id, raw_alt_index=r.info['V3_ALT_INDEX'][0],
                         normalized_allele=[r.contig, r.pos, r.ref, r.alts[0]],
                         original_POPAF=r.info['POPAF'][0])
            maps[r.contig].write(json.dumps(entry) + '\n')
            counts[r.contig] += 1
finally:
    for f in maps.values():
        f.close()
assert sum(counts.values()) == 888490
assert count == previous['counts'].get('raw_records', count)
save(OUT / 'commands.json', [command])
files = [tagged, normalized, OUT / 'commands.json']
files += [ROOT / 'population' / c / name for c in chromosomes for name in ['raw.vcf.gz', 'alleles.jsonl']]
# Keep the setup output paths relative to ROOT so all shards are covered.
receipt(ROOT, files + [ROOT / 'policy.json'], stage='setup', original_records=count,
        normalized_alleles=sum(counts.values()), counts_by_chromosome=dict(counts),
        original_unfiltered_sha256=sha(original), DNA_discovery_repeated=False, Caris_nominations=False,
        source_root=str(DISCOVERY), worker_sha256=sha(__file__))
shutil.copyfile(ROOT / 'complete.json', OUT / 'complete.json')
