"""Filter the complete discovery set using population priors and tumour segmentation."""
from pathlib import Path
import collections
import json
import sys
import pysam
from common import DISCOVERY, sha, save, run, receipt

ROOT = Path(sys.argv[1])
MODEL = 'segmented'
OUT = ROOT / 'filter' / MODEL
OUT.mkdir(parents=True, exist_ok=True)
chromosomes = ['chr' + str(c) for c in list(range(1, 23)) + ['X', 'Y']]
inputs = []
expected = 0
for chrom in chromosomes:
    d = ROOT / 'population' / chrom
    r = json.loads((d / 'complete.json').read_text())
    assert r['status'] == 'COMPLETED'
    path = d / 'updated.vcf.gz'
    assert sha(path) == r['outputs'][path.name]['sha256']
    inputs.append(path)
    expected += r['all_emitted_alleles']
assert expected == 888490
merged = OUT / 'unfiltered.vcf.gz'
def gatk(tool, *args):
    return ['/gatk/gatk', '--java-options', '-Xmx11g -XX:ParallelGCThreads=2', tool, *args]
commands = [run(OUT, 'merge-updated-calls', gatk('MergeVcfs', *[a for p in inputs for a in ['-I', p]], '-O', merged))]
plan = json.loads((ROOT / 'calling-plan.json').read_text())
filtered = OUT / 'filtered.vcf.gz'
args = gatk('FilterMutectCalls', '-R', plan['reference'], '-V', merged,
            '--stats', DISCOVERY / 'gather/unfiltered.vcf.gz.stats',
            '--contamination-table', ROOT / 'support/contamination.table',
            '--ob-priors', DISCOVERY / 'gather/orientation-model.tar.gz',
            '--threshold-strategy', 'FALSE_DISCOVERY_RATE', '--false-discovery-rate', '0.05',
            '--filtering-stats', OUT / 'filtering-stats.tsv', '-O', filtered)
if MODEL == 'segmented':
    args += ['--tumor-segmentation', ROOT / 'support/allele-fraction-segments.table']
commands.append(run(OUT, 'filter-mutect-calls', args))
tagged = OUT / 'tagged.vcf.gz'
with pysam.VariantFile(filtered) as src:
    h = src.header.copy()
    h.info.add('V3_EFFECTIVE_FILTER', 'A', 'String', 'Native site and per-ALT filters before biallelic normalization')
    with pysam.VariantFile(tagged, 'wz', header=h) as dst:
        for r in src:
            r.translate(h)
            raw = r.info.get('AS_FilterStatus', 'SITE')
            if isinstance(raw, tuple):
                raw = ','.join(raw)
            statuses = raw.split('|')
            assert len(statuses) == len(r.alts)
            site = set(r.filter) - {'PASS', '.'}
            r.info['V3_EFFECTIVE_FILTER'] = tuple('&'.join(sorted(site | (set(x.split(',')) - {'SITE', 'PASS', '.', ''}))) or 'PASS' for x in statuses)
            dst.write(r)
normalized = OUT / 'normalized.vcf.gz'
commands.append(run(OUT, 'normalize-filtered', ['bcftools', 'norm', '-f', plan['reference'], '-m', '-any', '-Oz', '-o', normalized, tagged]))
commands.append(run(OUT, 'index-normalized', ['bcftools', 'index', '-t', '-f', normalized]))
counts = collections.Counter()
with pysam.VariantFile(normalized) as src:
    for r in src:
        assert len(r.alts) == 1
        counts['normalized_alleles'] += 1
        counts['native_' + r.info['V3_EFFECTIVE_FILTER'][0]] += 1
assert counts['normalized_alleles'] == expected
save(OUT / 'commands.json', commands)
receipt(OUT, [merged, filtered, normalized, Path(str(normalized) + '.tbi'), OUT / 'filtering-stats.tsv', OUT / 'commands.json'],
        model=MODEL, counts=dict(counts), every_emitted_allele_population_reviewed=True,
        population_prior='max(native, genome global AF, exome global AF), source QC PASS only',
        filtering_strategy='FALSE_DISCOVERY_RATE', target=0.05, empirical_FDR_established=False,
        worker_sha256=sha(__file__), DNA_discovery_repeated=False)
