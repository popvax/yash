"""Audit actual full-pool inputs, priors, downstream accounting and score arithmetic."""
from pathlib import Path
import collections
import csv
import gzip
import itertools
import json
import math
import struct
import sys

import pysam
from common import DISCOVERY, sha, save, receipt

ROOT = Path(sys.argv[1])
OUT = ROOT / 'audit'
OUT.mkdir(parents=True, exist_ok=True)
checks = []


def check(name, condition, **evidence):
    assert condition, (name, evidence)
    checks.append(dict(check=name, verified=True, **evidence))


def load(path):
    return json.loads(Path(path).read_text())


def rows(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def equal(a, b):
    if isinstance(a, (tuple, list)) and isinstance(b, (tuple, list)):
        return len(a) == len(b) and all(equal(x, y) for x, y in zip(a, b))
    if isinstance(a, dict) and isinstance(b, dict):
        return set(a) == set(b) and all(equal(a[k], b[k]) for k in a)
    if isinstance(a, float) or isinstance(b, float):
        return a is not None and b is not None and math.isclose(a, b, rel_tol=1e-6, abs_tol=1e-8)
    return a == b


setup = load(ROOT / 'setup/complete.json')
check('RNA_BAM_matches_shared_data_catalog', sha('/data/inputs/TN25-112103/RNA_TN25-112103.bam') == '4e48531afb8b1ecd9ae64493847def981ab0535c53444e1b456f48519287a819',
      catalog_source='data/catalog/object-manifest.json: raw/TN25-112103/RNA_TN25-112103.bam')
check('RNA_index_matches_shared_data_catalog', sha('/data/inputs/TN25-112103/RNA_TN25-112103.bam.bai') == '27e9aa7a225f6fffd4c918f0c7fcd9d0caefe1ad40d05b7382e9767535f3c943')
check('complete_original_discovery_scope', setup['normalized_alleles'] == 888490 and setup['original_records'] == 879912 and not setup['DNA_discovery_repeated'] and not setup['Caris_nominations'])
for name, value in setup['outputs'].items():
    check('setup_hash:' + name, sha(ROOT / name) == value['sha256'])
original_count = 0
tlod_roundings = 0
maximum_tlod_rounding_absolute = 0.0
maximum_tlod_rounding_relative = 0.0
with pysam.VariantFile(DISCOVERY / 'gather/unfiltered.vcf.gz') as old, pysam.VariantFile(ROOT / 'filter/segmented/unfiltered.vcf.gz') as new:
    for a, b in itertools.zip_longest(old, new):
        assert a is not None and b is not None, 'Original raw records were lost or added'
        assert (a.contig, a.pos, a.ref, a.alts) == (b.contig, b.pos, b.ref, b.alts)
        assert equal({k: dict(v.items()) for k, v in a.samples.items()}, {k: dict(v.items()) for k, v in b.samples.items()})
        assert {k: v.phased for k, v in a.samples.items()} == {k: v.phased for k, v in b.samples.items()}
        assert equal(a.qual, b.qual)
        old_info = {k: v for k, v in a.info.items() if k != 'POPAF'}
        new_info = {k: v for k, v in b.info.items() if k not in ['POPAF', 'V3_ALT_INDEX']}
        old_tlod, new_tlod = old_info.pop('TLOD'), new_info.pop('TLOD')
        assert len(old_tlod) == len(new_tlod)
        for before, after in zip(old_tlod, new_tlod):
            if before != after:
                # HTSlib writes INFO floats with six significant digits. Match
                # that exact serialization, rather than accepting an arbitrary
                # tolerance or claiming every original TLOD bit is unchanged.
                expected = struct.unpack('f', struct.pack('f', float(format(before, '.6g'))))[0]
                assert after == expected, ('Unexpected TLOD change', a.contig, a.pos, before, after, expected)
                tlod_roundings += 1
                maximum_tlod_rounding_absolute = max(maximum_tlod_rounding_absolute, abs(before - after))
                maximum_tlod_rounding_relative = max(maximum_tlod_rounding_relative, abs(before - after) / abs(before))
        assert equal(old_info, new_info)
        assert tuple(a.filter) == tuple(b.filter)
        original_count += 1
check('all_raw_alleles_genotypes_and_counts_preserved_with_documented_float_serialization', original_count == 879912,
      records=original_count, TLOD_values_with_six_significant_digit_rounding=tlod_roundings,
      maximum_TLOD_absolute_rounding=maximum_tlod_rounding_absolute, maximum_TLOD_relative_rounding=maximum_tlod_rounding_relative,
      original_unrounded_values_preserved_in=str(DISCOVERY / 'gather/unfiltered.vcf.gz'))
total, increased = 0, 0
for chrom in ['chr' + str(c) for c in list(range(1, 23)) + ['X', 'Y']]:
    directory = ROOT / 'population' / chrom
    completed = load(directory / 'complete.json')
    for name, meta in completed['outputs'].items():
        check('population_hash:' + chrom + ':' + name, sha(directory / name) == meta['sha256'])
    original = [json.loads(line) for line in (directory / 'alleles.jsonl').open()]
    observed = {}
    with gzip.open(directory / 'population-observations.jsonl.gz', 'rt') as f:
        for line in f:
            x = json.loads(line)
            key = tuple(x['normalized_allele'])
            assert key not in observed
            good = [r for cohort in x['population'] for r in cohort['records'] if r['FILTER'] == ['PASS'] and r.get('AF') is not None and r.get('AN')]
            observed[key] = max([r['AF'] for r in good], default=None)
    expected = {}
    for x in original:
        k = tuple(x['normalized_allele'])
        assert k in observed
        native = 10 ** -x['original_POPAF']
        af = max(native, observed[k] if observed[k] is not None else native)
        expected[x['raw_id'], x['raw_alt_index']] = af
        increased += af > native * 1.01
    seen = set()
    with pysam.VariantFile(directory / 'updated.vcf.gz') as src:
        for r in src:
            for i, popaf in enumerate(r.info['POPAF']):
                loc = r.id, i
                assert loc in expected and loc not in seen
                # HTSlib's VCF float serialization need not preserve every bit.
                assert math.isclose(10 ** -popaf, expected[loc], rel_tol=1e-4, abs_tol=1e-12), (chrom, loc, popaf, expected[loc])
                seen.add(loc)
    check('every_emitted_ALT_has_expected_dual_cohort_prior:' + chrom, len(seen) == len(original) == completed['all_emitted_alleles'], alleles=len(seen))
    total += len(seen)
check('every_emitted_allele_population_reviewed', total == 888490, alleles=total, priors_increased_more_than_1_percent=increased)
for stage in ['filter/segmented', 'gather', 'annotate', 'rna', 'rank', 'summarize']:
    directory = ROOT / stage
    completed = load(directory / 'complete.json')
    check('completed:' + stage, completed['status'] == 'COMPLETED')
    for name, meta in completed['outputs'].items():
        check('stage_hash:' + stage + ':' + name, sha(directory / name) == meta['sha256'])
    if stage.startswith('filter/'):
        check('full_filter_scope:' + stage, completed['counts']['normalized_alleles'] == 888490 and completed['every_emitted_allele_population_reviewed'])
        commands = load(directory / 'commands.json')
        args = commands[1]
        check('unchanged_FDR_target:' + stage, args[args.index('--threshold-strategy') + 1] == 'FALSE_DISCOVERY_RATE' and args[args.index('--false-discovery-rate') + 1] == '0.05')
        check('balance_model:' + stage, ('--tumor-segmentation' in args) == stage.endswith('/segmented'))

evidence = {r['variant_id']: r for r in load(ROOT / 'annotate/variant-evidence.json')}
input_rows = load(ROOT / 'annotate/rna-input.json')
input_ids = {r['variant_id'] for r in input_rows}
expected_ids = {r['variant_id'] for r in evidence.values() if r['RNA_review_eligible']}
check('RNA_pool_rebuilt_from_all_new_survivors', input_ids == expected_ids and len(input_rows) == len(input_ids), variants=len(input_ids))
rna = {r['variant_id']: r for r in load(ROOT / 'rna/rna-evidence.json')}
check('all_new_RNA_inputs_accounted_for', set(rna) == input_ids)
for r in rna.values():
    if r['rna_total_fragments']:
        assert math.isclose(r['rna_alt_fraction'], r['rna_alt_fragments'] / r['rna_total_fragments'], rel_tol=1e-10)
    else:
        assert r['rna_alt_fraction'] is None
proteins = {r['variant_id']: r for r in load(ROOT / 'rna/primary-proteins.json')}
check('all_primary_proteins_native_RNA_usable', all(p['usable_mutant_context'] and p['variant_passes_isovar_filters'] and rna[p['variant_id']]['passes_isovar_filters'] for p in proteins.values()))
predictions = rows(ROOT / 'rank/peptide-HLA-affinity.tsv')
affinity = {(r['peptide'], r['allele']): float(r['affinity_nm']) for r in predictions}
alleles = load(ROOT / 'case.json')['hla_class_I']
peptides = {p for p, a in affinity}
check('fresh_prediction_coverage', len(affinity) == len(predictions) == len(peptides) * 6 and set(affinity) == {(p, a) for p in peptides for a in alleles}, predictions=len(predictions), peptides=len(peptides))
reference = load(ROOT / 'rank/reference-matches.json')
check('reference_check_coverage', set(reference) == peptides and all(x['source_provenance_complete'] and not x['excluded_gene_ids'] for x in reference.values()))
def contribution(nm):
    score = 0 if nm >= 5000 else (1 + math.exp(-350 / 150)) / (1 + math.exp((nm - 350) / 150))
    return score if score >= 1e-5 else 0
peptide_scores = {p: sum(contribution(affinity[p, a]) for a in alleles) for p in peptides}
occurrences = collections.defaultdict(list)
for r in rows(ROOT / 'rank/peptide-occurrences.tsv'):
    a, b = int(r['start']), int(r['end'])
    p = proteins[r['variant_id']]
    assert p['sequence'][a:b] == r['peptide'] and a < p['mutation_end_0based_exclusive'] and b > p['mutation_start_0based']
    if not reference[r['peptide']]['occurs']:
        occurrences[r['variant_id']].append((a, b, r['peptide']))
windows = rows(ROOT / 'rank/all-scored-windows.tsv')
maxima = collections.defaultdict(lambda: -1.0)
for r in windows:
    vid = r['variant_id']
    a, b = int(r['window_start_0based']), int(r['window_end_exclusive'])
    assert b - a == 25 and proteins[vid]['sequence'][a:b] == r['sequence']
    inside = [p for x, y, p in occurrences[vid] if a <= x and y <= b]
    score = sum(peptide_scores[p] for p in inside)
    assert math.isclose(score, float(r['affinity_contribution_sum']), rel_tol=1e-8, abs_tol=1e-8)
    best = min([affinity[p, h] for p in inside for h in alleles], default=None)
    assert (best is None and not r['best_target_affinity_nm']) or math.isclose(best, float(r['best_target_affinity_nm']), rel_tol=1e-10)
    if r['weighted_score']:
        assert math.isclose(float(r['weighted_score']), score * math.sqrt(float(r['supplied_gene_value']) * rna[vid]['rna_alt_fraction']), rel_tol=1e-8, abs_tol=1e-8)
    maxima[vid] = max(maxima[vid], score)
check('all_window_sequences_scores_and_weights_recounted', True, windows=len(windows))
all_ranked = rows(ROOT / 'rank/all-ranked-regions.tsv')
check('all_maximum_windows_recounted', all(math.isclose(float(r['affinity_contribution_sum']), maxima[r['variant_id']], rel_tol=1e-8, abs_tol=1e-8) for r in all_ranked), regions=len(all_ranked))
primary = load(ROOT / 'rank/primary-ranked-regions.json')
for r in primary:
    d = evidence[r['variant_id']]
    assert d['primary_DNA_eligible'] and r['gene_assignment_unambiguous'] and r['passes_500nm_gate'] and r['weighted_score'] is not None
for model in ['segmented']:
    values = load(ROOT / 'summarize' / (model + '-primary-ranked-regions.json'))
    expected = [r['variant_id'] for r in primary if evidence[r['variant_id']]['primary_DNA_eligible_by_model'][model]]
    check('model_specific_ranking:' + model, [r['variant_id'] for r in values] == expected and [r['model_rank'] for r in values] == list(range(1, len(values) + 1)), regions=len(values))
summary = load(ROOT / 'summarize/summary.json')
check('segmented_primary_selection', summary['primary_model'] == 'segmented')
save(OUT / 'checks.json', checks)
receipt(OUT, [OUT / 'checks.json'], checks=len(checks), all_verified=True,
        full_emitted_allele_priors_verified=total, raw_record_preservation_verified=original_count,
        newly_rebuilt_RNA_inputs=len(input_ids), fresh_peptide_HLA_pairs=len(predictions), all_window_scores_verified=len(windows),
        ranked_regions_by_model=summary['ranked_regions_by_model'], worker_sha256=sha(__file__),
        scientific_scope='Accounting and computational consistency; no experimental or matched-normal somatic validation')
