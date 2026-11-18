#!/usr/bin/env python3
"""Validate complete appendix records against five authenticated Design 4 sources.

Ordinary verification derives all core records and inherited Tier 1 fields, checks
paired affinity nM against stored exact peptide/HLA values, and checks finite
ranges and derived flags of percentiles. Only --recompute checks percentile
VALUES and repeats affinity inference with the recorded assets. No calculation
here establishes biological or clinical meaning.
"""
import argparse
import csv
import math
from pathlib import Path
import re
import sys

from validation import (INPUTS, PINNED_INPUT_DIGESTS, SHORT, GUARDANT_COVERED, SOURCE_LABEL, DIGEST_NOTE,
                        LOOKUP_METHOD, FRESH_NOTE, INTERPRETATION, CROSS_MODEL_NOTE,
                        SHARED_PREDICATE, ENCODED, REGION, RESIDUAL, TIER, Invalid,
                        require, sha256, keys, equal, number, boolean, read_json,
                        read_tsv, unique, compare_records, aggregate, hash_assets,
                        validate_asset_receipt, residual_low_dai)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--design4', required=True, type=Path)
    ap.add_argument('--recompute', action='store_true')
    ap.add_argument('--models-dir', type=Path, help='local affinity directory; recorded asset identity remains enforced')
    args = ap.parse_args()
    base = args.design4
    data = Path(__file__).resolve().parent.parent / 'data'
    receipt = read_json(data / 'design4-input-digests.json')
    keys(receipt, ('design4_inputs_read_by_this_appendix', 'note'), 'input digest receipt')
    equal(receipt['note'], DIGEST_NOTE, 'input digest receipt.note')
    recorded = receipt['design4_inputs_read_by_this_appendix']
    keys(recorded, INPUTS, 'input digest receipt.files')
    equal(recorded, PINNED_INPUT_DIGESTS, 'pinned Design 4 source identities')
    for rel, digest in recorded.items():
        require(type(digest) is str and re.fullmatch(r'[0-9a-f]{64}', digest), f'{rel}: invalid recorded SHA-256')
        require(sha256(base / rel) == digest, f'{rel}: source SHA-256 differs from receipt')
    print('PASS: all five declared Design 4 source identities match')

    selected = read_tsv(base / INPUTS[0], required=('gene', 'sequence', 'best_target_peptide', 'best_target_allele', 'best_target_affinity_nm'))
    selected_map = unique(selected, ('gene',), 'selected source')
    require(bool(selected), 'selected source is empty')
    source = {k[0]: row for k, row in selected_map.items()}
    for gene, row in source.items():
        require(bool(gene) and re.fullmatch(r'[ACDEFGHIKLMNPQRSTVWY]+', row['sequence']), f'{gene}: invalid encoded sequence')
        require((row['best_target_peptide'] == '' and row['best_target_allele'] == '')
                or (row['best_target_peptide'] != '' and row['best_target_peptide'] in row['sequence']
                    and row['best_target_allele'] in SHORT), f'{gene}: invalid source lead peptide/allele')
        number(row['best_target_affinity_nm'], f'{gene} source lead affinity', positive=True)
    contrasts = read_tsv(base / INPUTS[1], required=('scenario', 'gene', 'mutant_peptide', 'counterfactual_peptide', 'allele', 'mutant_affinity_nm', 'counterfactual_affinity_nm'))
    comparisons = [r for r in contrasts if r['scenario'] == 'single_aa_WT_reversion'
                   and r['gene'] in source and r['mutant_peptide'] in source[r['gene']]['sequence']]
    unique(comparisons, ('gene', 'mutant_peptide', 'allele'), 'counterfactual source')
    encoded, by_gene = [], {g: [] for g in source}
    for row in comparisons:
        gene, peptide, allele = row['gene'], row['mutant_peptide'], row['allele']
        wt = row['counterfactual_peptide']
        require(allele in SHORT and 8 <= len(peptide) <= 11 and len(wt) == len(peptide)
                and sum(a != b for a, b in zip(peptide, wt)) == 1,
                f'{gene} {peptide}: invalid single-AA WT comparison')
        m = number(row['mutant_affinity_nm'], f'{gene} {peptide} source mutant nM', positive=True)
        w = number(row['counterfactual_affinity_nm'], f'{gene} {peptide} source WT nM', positive=True)
        d = w / m
        require(math.isfinite(d), f'{gene} {peptide}: non-finite DAI')
        lead = source[gene]
        encoded.append(dict(zip(ENCODED, (gene, peptide, wt, allele, round(m, 4), round(w, 4), round(d, 4),
            (peptide, allele) == (lead['best_target_peptide'], lead['best_target_allele']), m <= 500, True))))
        by_gene[gene].append((d, m, peptide, allele))
    compare_records(data / 'encoded-peptide-dai.tsv', ENCODED, encoded, ('gene', 'mutant_peptide', 'allele'),
                    numeric=('mutant_affinity_nm', 'wildtype_affinity_nm', 'dai_wt_over_mut'),
                    booleans=('is_design4_lead', 'binds_le_500nM', 'encoded_in_design4_window'))
    residual = []
    for gene, comparisons in by_gene.items():
        for d, m, peptide, allele in comparisons:
            if residual_low_dai(d, m):
                n = len(source[gene]['sequence'])
                residual.append(dict(zip(RESIDUAL, (gene, peptide, allele, round(m, 4), round(d, 4), n,
                    f"inside the encoded {n}-aa window; relabelling a region's lead does NOT remove this peptide from the construct"))))
    compare_records(data / 'residual-low-dai-peptides.tsv', RESIDUAL, residual, ('gene', 'peptide', 'allele'),
                    numeric=('mutant_affinity_nm', 'dai'), integers=('encoded_window_aa',))
    regions = []
    for gene, lead in source.items():
        c = by_gene[gene]
        lp, la = lead['best_target_peptide'], lead['best_target_allele']
        ld = next((d for d, m, p, al in c if (p, al) == (lp, la)), None)
        eligible = [x for x in c if x[1] <= 500]
        best = max(eligible, key=lambda x: x[0]) if eligible else None
        better = bool(best and ld is not None and best[0] > ld + 1e-9)
        regions.append(dict(zip(REGION, (gene, len(lead['sequence']), lp, la,
            float(lead['best_target_affinity_nm']), round(ld, 4) if ld is not None else None,
            round(best[0], 4) if best else None, best[2] if best else '', best[3] if best else '',
            better, ('allele_only' if best[2] == lp else 'peptide_change') if better else '',
            sum(residual_low_dai(d, m) for d, m, p, al in c), gene in GUARDANT_COVERED, bool(c)))))
    compare_records(data / 'region-dai-summary.tsv', REGION, regions, ('gene',),
                    numeric=('design4_lead_nm', 'design4_lead_dai', 'best_dai_in_window_le_500nM'),
                    integers=('window_aa', 'n_mutation_spanning_pep_allele_dai_le_1_and_le500nM'),
                    booleans=('higher_dai_choice_exists', 'guardant360_covered', 'dai_defined'))
    equal(read_json(data / 'receipt.json'), {
        'regions_examined': len(source), 'peptide_allele_rows': len(encoded),
        'residual_low_dai_pep_allele_pairs': len(residual), 'source': SOURCE_LABEL,
        'model_inference_run': False}, 'core receipt')
    print('PASS: all core table fields, blanks, types, unique keys and receipt derive from source')

    mix = {}
    for row in read_tsv(base / INPUTS[2], required=('Peptide', *('%Rank_' + s for s in SHORT.values())), comments=True):
        for al, short in SHORT.items():
            value = row['%Rank_' + short]
            if value != '':
                key = row['Peptide'], al
                require(key not in mix, f'MixMHCpred source: duplicate key {key}')
                mix[key] = number(value, f'MixMHCpred source {key}', percentile=True)
    hom = [r for r in read_tsv(base / INPUTS[3], required=('model_priority_category', 'query_peptide', 'homolog_peptide', 'shared_binding_allele', 'homolog_genes', 'shared_allele_homolog_nm'))
           if r['model_priority_category'] == 'TIER1_1MM_PARALOG_SAME_HLA' and r['shared_binding_allele']]
    hom_map = unique(hom, ('query_peptide', 'homolog_peptide', 'shared_binding_allele'), 'TIER1 source')
    cache = {}
    for row in read_tsv(base / INPUTS[4], required=('peptide', 'allele', 'affinity_nm')):
        key = row['peptide'], row['allele']
        require(key not in cache, f'junction affinity source: duplicate key {key}')
        cache[key] = number(row['affinity_nm'], f'junction affinity source {key}', positive=True)
    homolog_values = {}
    for row in hom:
        key = row['homolog_peptide'], row['shared_binding_allele']
        val = number(row['shared_allele_homolog_nm'], f'homolog affinity source {key}', positive=True)
        require(key not in homolog_values or homolog_values[key] == val, f'homolog source: conflicting duplicate {key}')
        homolog_values[key] = val
    tier = read_tsv(data / 'tier1-direction.tsv', TIER)
    tier_map = unique(tier, ('query_peptide', 'homolog_peptide', 'shared_allele'), 'TIER1 table')
    require(bool(tier) and set(tier_map) == set(hom_map), 'TIER1 table: missing or unexpected source keys')
    per_source = {'all-junction-affinity.tsv.gz': 0, 'offtarget-homologs.tsv.gz': 0}
    maxdiff, worst = 0.0, None
    for key, row in tier_map.items():
        src = hom_map[key]
        require(row['shared_allele'] in SHORT, f'{key}: unknown shared HLA')
        equal(row['homolog_genes'], src['homolog_genes'], f'{key} homolog_genes')
        for side in ('query', 'homolog'):
            nm = number(row[f'mhcflurry_{side}_nm'], f'{key} {side} nM', positive=True)
            percentile = number(row[f'mhcflurry_{side}_percentile'], f'{key} {side} percentile', percentile=True)
            mk = row[f'{side}_peptide'], row['shared_allele']
            want = mix.get(mk)
            rank = number(row[f'mixmhcpred_{side}_rank'], f'{key} {side} Mix rank', blank=want is None, percentile=True)
            require(rank == want, f'{key} {side}: inherited Mix rank differs from source')
            require(nm <= 500 or percentile <= 2 or (rank is not None and rank <= 2),
                    f'{key} {side}: original shared-HLA disjunction is not satisfied')
            if mk in cache:
                value, name = cache[mk], 'all-junction-affinity.tsv.gz'
            else:
                require(mk in homolog_values, f'{key} {side}: no exact peptide/HLA affinity source')
                value, name = homolog_values[mk], 'offtarget-homologs.tsv.gz'
            delta = abs(value - nm)
            require(delta <= 0.05, f'{key} {side}: affinity differs from stored source by {delta} nM')
            per_source[name] += 1
            if delta > maxdiff:
                maxdiff, worst = delta, [mk[0], mk[1], value, nm]
        q, h = float(row['mhcflurry_query_nm']), float(row['mhcflurry_homolog_nm'])
        qp, hp = float(row['mhcflurry_query_percentile']), float(row['mhcflurry_homolog_percentile'])
        qm, hm = mix.get((key[0], key[2])), mix.get((key[1], key[2]))
        equal(row['mhcflurry_relative'], 'homolog_stronger_or_equal' if h <= q else 'homolog_weaker', f'{key} direction')
        for field, expected in (('both_exceed_500nM_mhcflurry', q > 500 and h > 500),
                                ('both_le_2pct_mixmhcpred', qm is not None and hm is not None and qm <= 2 and hm <= 2),
                                ('both_le_2pct_mhcflurry_percentile', qp <= 2 and hp <= 2)):
            equal(boolean(row[field], f'{key} {field}'), expected, f'{key} {field}')
    corroboration = {'cells_checked': len(tier) * 2, 'cells_corroborated': len(tier) * 2,
        'cells_uncovered_or_disagreeing': 0, 'uncovered_detail': [], 'max_abs_difference_nm': round(maxdiff, 8),
        'largest_difference_cell': worst, 'tolerance_nm': 0.05, 'lookup_method': LOOKUP_METHOD,
        'per_source_corroborated': per_source, 'model_inference_run': False}
    equal(read_json(data / 'cell-corroboration.json'), corroboration, 'corroboration receipt')
    flag = lambda row, field: row[field] == 'True'
    stronger = lambda row: row['mhcflurry_relative'] == 'homolog_stronger_or_equal'
    over = [row for row in tier if flag(row, 'both_exceed_500nM_mhcflurry')]
    comparable = [row for row in tier if row['mixmhcpred_query_rank'] != '' and row['mixmhcpred_homolog_rank'] != '']
    derived = {'tier1_rows_with_shared_allele': len(tier),
        'mhcflurry_homolog_stronger_or_equal': sum(map(stronger, tier)),
        'mhcflurry_homolog_weaker': len(tier) - sum(map(stronger, tier)),
        'both_exceed_500nM_mhcflurry': len(over),
        'of_which_homolog_stronger_or_equal': sum(map(stronger, over)),
        'of_which_homolog_weaker': len(over) - sum(map(stronger, over)),
        'both_le_2pct_mixmhcpred': sum(flag(row, 'both_le_2pct_mixmhcpred') for row in tier),
        'both_le_2pct_mixmhcpred_among_the_over500nM': sum(flag(row, 'both_le_2pct_mixmhcpred') for row in over),
        'both_le_2pct_mhcflurry_percentile_among_the_over500nM': sum(flag(row, 'both_le_2pct_mhcflurry_percentile') for row in over),
        'pairs_with_both_mixmhcpred_ranks': len(comparable),
        'cross_model_direction_reversed': sum(stronger(row) != (float(row['mixmhcpred_homolog_rank']) <= float(row['mixmhcpred_query_rank'])) for row in comparable),
        'tier1_table_sha256': sha256(data / 'tier1-direction.tsv'),
        'fresh_inference_cells_note': FRESH_NOTE, 'interpretation': INTERPRETATION,
        'cross_model_note': CROSS_MODEL_NOTE, 'shared_hla_predicate': SHARED_PREDICATE}
    summary = read_json(data / 'tier1-direction-summary.json')
    keys(summary, (*derived, 'models'), 'Tier 1 summary')
    for name, expected in derived.items():
        equal(summary[name], expected, f'Tier 1 summary.{name}')
    models = summary['models']
    keys(models, ('fresh', 'model_bundle', 'read_from_design4_not_rerun', 'not_consulted'), 'models')
    for name, expected in {'fresh': 'mhcflurry 2.3.14', 'read_from_design4_not_rerun': 'MixMHCpred 3.0 (mixmhcpred3.tsv)',
                           'not_consulted': 'PRIME (not part of the shared-HLA predicate)'}.items():
        equal(models[name], expected, 'models.' + name)
    bundle = models['model_bundle']
    keys(bundle, ('package_version', 'num_supported_alleles', 'assets'), 'model bundle')
    equal(bundle['package_version'], '2.3.14', 'model bundle.package_version')
    equal(bundle['num_supported_alleles'], 20205, 'model bundle.num_supported_alleles')
    assets = bundle['assets']
    validate_asset_receipt(assets)
    print('PASS: complete Tier 1 source membership, inherited fields, numeric bounds, shared-HLA predicate, labels, summaries and receipts')
    print('PASS: model receipt schema and internal hash/count/calibration consistency')

    if args.recompute:
        directory = (args.models_dir or Path(assets['resolved_affinity_dir'])).expanduser().resolve()
        require(directory.is_dir(), f'model directory is absent: {directory}')
        print(f"MODEL DIRECTORY: {directory}" + (' (--models-dir override)' if args.models_dir else ''))
        print(f"RECORDED GENERATION PATH: {assets['resolved_affinity_dir']}")
        actual = hash_assets(directory)
        equal(actual, assets['per_file_sha256'], 'local model assets')
        equal(aggregate(actual), assets['aggregate_sha256'], 'local model aggregate')
        import mhcflurry
        from mhcflurry import Class1AffinityPredictor
        equal(mhcflurry.__version__, bundle['package_version'], 'loaded mhcflurry version')
        predictor = Class1AffinityPredictor.load(models_dir=str(directory))
        equal(len(predictor.supported_alleles), bundle['num_supported_alleles'], 'loaded supported allele count')
        peptides = sorted({row[s + '_peptide'] for row in tier for s in ('query', 'homolog')})
        predictions = {}
        for allele in sorted({row['shared_allele'] for row in tier}):
            frame = predictor.predict_to_dataframe(peptides=peptides, allele=allele, include_percentile_ranks=True, throw=False)
            for _, row in frame.iterrows():
                predictions[(row['peptide'], allele)] = float(row['prediction']), float(row['prediction_percentile'])
        for key, row in tier_map.items():
            for side in ('query', 'homolog'):
                nm, pc = predictions[(row[side + '_peptide'], key[2])]
                require(math.isfinite(nm) and nm > 0 and math.isfinite(pc) and 0 <= pc <= 100,
                        f'{key} {side}: invalid fresh model output')
                require(abs(nm - float(row[f'mhcflurry_{side}_nm'])) <= 0.05
                        and abs(pc - float(row[f'mhcflurry_{side}_percentile'])) <= 0.001,
                        f'{key} {side}: fresh affinity/percentile differs')
        print('PASS: fresh affinity and percentile replay from that same authenticated local model directory')
    else:
        print('NOT REPLAYED: model inference and percentile values (use --recompute)')
    print('VERIFIED: declared arithmetic and source correspondence; no biological or clinical conclusion')


if __name__ == '__main__':
    try:
        main()
    except (Invalid, OSError, ValueError, KeyError, TypeError, ZeroDivisionError, ImportError,
            csv.Error, OverflowError, EOFError) as error:
        print(f'FAIL: {error}', file=sys.stderr)
        sys.exit(1)
