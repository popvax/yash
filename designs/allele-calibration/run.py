"""Allele calibration and HLA coverage audit of lead Design #3.

Offline. Reads only recorded repository evidence, runs no inference, starts no
cloud job and changes no shared file. With no flag it recomputes every result
and compares it with the committed files; `--write` refreshes them.

Three questions:

1. Is the recorded selection objective allele-neutral? It scores a peptide from
   its predicted affinity in nM. MHCflurry's nM scale differs between alleles,
   so the same nM is a different percentile for each of Yash's six class-I
   alleles. The committed prediction tables are enough to measure that offset.
2. Does that offset explain the design's concentration on HLA-C, or is the
   concentration a property of the antigens? Answered by comparing the loci at
   matched stringency on both scales.
3. What does the selected set depend on: which regions rest on one allele, what
   a haplotype loss would cost, and how far the order moves when the nM offset
   is removed or HLA-C is down-weighted for its lower cell-surface density.
"""
from pathlib import Path
import argparse
import collections
import csv
import io
import itertools
import json
import math
import random
import re
import sys

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO / 'src'))

from calibration import ALLELES, CALIBRATION_SOURCES, Calibration, open_text
from yash_vaccine.io import require_assertions, sha256, write_json, write_tsv

RESULTS = HERE / 'audit'
REGION_ANNOTATIONS = 'designs/design-3/review/candidate-annotations/region-annotations.tsv'
PEPTIDE_ANNOTATIONS = 'designs/design-3/review/candidate-annotations/peptide-annotations.tsv.gz'
SELECTED = 'designs/design-3/selection/selected.tsv'
HLA_CLASS_I = 'data/hla/class-I.tsv'
SCORE_CONFIG = 'designs/initial/inputs/mutations/rank/score-config.json'
CONTRIBUTIONS = 'designs/initial/inputs/mutations/rank/affinity-contributions.tsv'

# The recorded objective. These constants are asserted against SCORE_CONFIG at
# run time, and the score they produce against every row of CONTRIBUTIONS.
MIDPOINT, WIDTH, MIN_SCORE, CUTOFF = 350.0, 150.0, 1e-5, 5000.0
GATE_NM = 500.0
# Gates reported side by side. The nM gate is the design's own. The rank gates
# are NetMHCpan's conventional strong/weak thresholds, applied here to
# MHCflurry's percentile, which has its own reference distribution.
RANK_GATES = (0.5, 2.0)
# HLA-C surface density relative to HLA-A/B. 1.00 is the recorded assumption;
# 0.25 and 0.20 are the total-cellular ratios and 0.06 sits inside the surface
# range reported by Apps et al. 2015 (J Immunol 194:3594-3600).
HLA_C_WEIGHTS = (1.00, 0.25, 0.20, 0.06)


def epitope_score(affinity_nm):
    """The recorded per-peptide contribution for a predicted affinity."""
    if affinity_nm >= CUTOFF:
        return 0.0
    logistic = lambda value: 1.0 / (1.0 + math.exp((value - MIDPOINT) / WIDTH))
    score = logistic(affinity_nm) / logistic(0.0)
    return score if score >= MIN_SCORE else 0.0


def read_rows(relative):
    with open_text(REPO / relative) as stream:
        yield from csv.DictReader(stream, delimiter='\t')


def locus(allele):
    return allele[4]


def check_score_config():
    """Confirm the objective constants used here are the recorded ones."""
    text = json.loads((REPO / SCORE_CONFIG).read_text())['native_epitope_config']
    found = {name: float(value) for name, value in re.findall(
        r'(logistic_epitope_score_midpoint|logistic_epitope_score_width'
        r'|min_epitope_score|binding_affinity_cutoff)=([0-9.e+-]+)', text)}
    expected = {'logistic_epitope_score_midpoint': MIDPOINT,
                'logistic_epitope_score_width': WIDTH,
                'min_epitope_score': MIN_SCORE,
                'binding_affinity_cutoff': CUTOFF}
    if found != expected:
        raise ValueError(f'{SCORE_CONFIG} records {found}, not {expected}')
    if "scoring_mode='affinity'" not in text:
        raise ValueError(SCORE_CONFIG + ' is not in affinity scoring mode')
    return dict(source=SCORE_CONFIG, **dict(sorted(found.items())))


def check_peptide_scores():
    """Confirm epitope_score reproduces every recorded per-peptide contribution."""
    worst, checked, scored, unscored = 0.0, 0, [], []
    for row in read_rows(CONTRIBUTIONS):
        affinity = float(row['affinity_nm'])
        recorded = float(row['contribution'])
        worst = max(worst, abs(epitope_score(affinity) - recorded))
        checked += 1
        (scored if recorded > 0 else unscored).append(affinity)
    if not (checked > 4000 and worst < 1e-12):
        raise ValueError(f'per-peptide mismatch: {checked} rows, worst {worst:.3e}')
    return dict(source=CONTRIBUTIONS, rows_checked=checked,
                max_absolute_difference=float(f'{worst:.3e}'),
                largest_scored_affinity_nm=round(max(scored), 4),
                smallest_unscored_affinity_nm=round(min(unscored), 4))


def formula_weight(region, selected):
    """sqrt of the recorded expression proxy, where that proxy is recorded.

    Splice, CTA/self and fusion regions use sqrt(supplied gene value); mutation
    regions multiply it by the observed RNA alternate fragment fraction, as
    score-config.json records. None when supplied_gene_value was left blank.
    """
    supplied = region['supplied_gene_value']
    if supplied == '':
        return None
    value = float(supplied)
    if region['source'] == 'MUTATION':
        fraction = selected.get(region['region_id'], {}).get('rna_alternate_fraction')
        if fraction:
            value *= float(fraction)
    return math.sqrt(value)


def collect(calibration):
    """Per-region, per-allele target-peptide evidence under every nM scale."""
    regions = {row['annotation_id']: row for row in read_rows(REGION_ANNOTATIONS)}
    selected = {row['region_id']: row for row in read_rows(SELECTED)}
    recorded = collections.Counter()
    neutral = collections.Counter()
    by_reference = {allele: collections.Counter() for allele in ALLELES}
    weighted_by_c = {weight: collections.Counter() for weight in HLA_C_WEIGHTS}
    best_nm = collections.defaultdict(lambda: math.inf)
    best_pct = collections.defaultdict(lambda: math.inf)
    best_neutral = collections.defaultdict(lambda: math.inf)
    peptides = collections.Counter()
    interpolated = 0
    for row in read_rows(PEPTIDE_ANNOTATIONS):
        identity, allele = row['annotation_id'], row['allele']
        if identity not in regions or allele not in ALLELES:
            continue
        affinity = float(row['affinity_nm'])
        peptide = row['peptide']
        if (peptide, allele) not in calibration.observed:
            interpolated += 1
        percentile = calibration.percentile(allele, affinity, peptide=peptide,
                                            length=len(peptide))
        equivalent = calibration.neutral_affinity(percentile)
        recorded[identity] += epitope_score(affinity)
        neutral[identity] += epitope_score(equivalent)
        for reference in ALLELES:
            by_reference[reference][identity] += epitope_score(
                calibration.affinity_at(reference, percentile))
        for weight in HLA_C_WEIGHTS:
            factor = weight if locus(allele) == 'C' else 1.0
            weighted_by_c[weight][identity] += factor * epitope_score(equivalent)
        peptides[identity] += 1
        key = (identity, allele)
        best_nm[key] = min(best_nm[key], affinity)
        best_pct[key] = min(best_pct[key], percentile)
        best_neutral[key] = min(best_neutral[key], equivalent)
    return dict(regions=regions, selected=selected, recorded=recorded, neutral=neutral,
                by_reference=by_reference, weighted_by_c=weighted_by_c, best_nm=best_nm,
                best_pct=best_pct, best_neutral=best_neutral, peptides=peptides,
                interpolated=interpolated)


def recorded_weight(data, identity):
    """The expression weight implied by the recorded score and contribution sum.

    Every region carrying a recorded `original_weighted_score` gets its exact
    recorded weight, so the recorded column of the re-scoring table is the
    design's own number rather than a reconstruction. Six plan regions have a
    recorded score but a blank `supplied_gene_value`; this recovers their
    weight without joining to another design's candidate pool, where one of
    them is listed under a different gene symbol.
    """
    region = data['regions'][identity]
    if region['original_weighted_score'] == '':
        return None
    total = data['recorded'][identity]
    if total <= 0:
        return None
    return float(region['original_weighted_score']) / total


def check_recorded_objective(data):
    """Check the implied weights against the recorded expression formula."""
    worst, checked, implied_only = 0.0, 0, []
    for identity, region in data['regions'].items():
        weight = recorded_weight(data, identity)
        if weight is None:
            continue
        checked += 1
        expected = formula_weight(region, data['selected'])
        if expected is None:
            implied_only.append(dict(gene=region['gene'], source=region['source'],
                                     region_id=region['region_id'],
                                     implied_expression_proxy=round(weight ** 2, 6)))
            continue
        worst = max(worst, abs(weight - expected) / max(abs(expected), 1e-12))
    if not (checked >= 185 and worst < 1e-5):
        raise ValueError(f'objective mismatch: {checked} regions, worst {worst:.3e}')
    return dict(regions_with_recorded_score=checked,
                regions_checked_against_formula=checked - len(implied_only),
                max_relative_difference=float(f'{worst:.3e}'),
                residual_note=('Consistent with single-precision arithmetic in the '
                               'recorded pipeline: rows sharing one expression proxy '
                               'carry different residuals, which rounding of the proxy '
                               'alone cannot produce. The per-peptide score itself '
                               'agrees to the tolerance in recorded_objective.per_peptide.'),
                weight_implied_from_recorded_score_only=sorted(
                    implied_only, key=lambda row: row['gene']))


def check_case_alleles():
    """Confirm the six audited alleles are the recorded class-I calls."""
    recorded = tuple(sorted(row['allele'] for row in read_rows(HLA_CLASS_I)
                            if row['mhc_class'] == 'I'))
    if recorded != tuple(sorted(ALLELES)):
        raise ValueError('audited alleles differ from ' + HLA_CLASS_I
                         + ': ' + ';'.join(recorded))
    return list(recorded)


def calibration_table(calibration):
    for allele in calibration.alleles:
        for length in sorted(calibration.by_length[allele]):
            curve = calibration.by_length[allele][length]
            yield dict(
                allele=allele, peptide_length=length, recorded_points=len(curve),
                percentile_at_50nm=round(calibration.percentile(allele, 50.0, length=length), 4),
                percentile_at_500nm=round(calibration.percentile(allele, 500.0, length=length), 4),
                nm_at_0p5_percentile=round(calibration.affinity_at(allele, 0.5, length=length), 1),
                nm_at_2_percentile=round(calibration.affinity_at(allele, 2.0, length=length), 1),
                allele_level_percentile_at_500nm=round(calibration.percentile(allele, 500.0), 4),
                allele_level_nm_at_0p5_percentile=round(calibration.affinity_at(allele, 0.5), 1),
                allele_level_nm_at_2_percentile=round(calibration.affinity_at(allele, 2.0), 1))


def plan_identities(data):
    return sorted(identity for identity, region in data['regions'].items()
                  if region['current_plan'] == 'True')


def coverage_matrix(data, plan):
    for identity in plan:
        region = data['regions'][identity]
        for allele in ALLELES:
            key = (identity, allele)
            if data['best_nm'][key] == math.inf:
                continue
            yield dict(
                antigen_gene=region['gene'], source=region['source'],
                protein_change=region['protein_change'], region_id=region['region_id'],
                annotation_id=identity, allele=allele, locus=locus(allele),
                best_target_affinity_nm=round(data['best_nm'][key], 3),
                best_target_percentile=round(data['best_pct'][key], 4),
                allele_neutral_affinity_nm=round(data['best_neutral'][key], 3),
                passes_recorded_500nm_gate=data['best_nm'][key] <= GATE_NM,
                passes_allele_neutral_500nm_gate=data['best_neutral'][key] <= GATE_NM,
                strong_binder_0p5_percentile=data['best_pct'][key] <= 0.5,
                weak_binder_2_percentile=data['best_pct'][key] <= 2.0)


def covered(data, identity, alleles, gate):
    field, threshold = gate
    return any(data[field][(identity, allele)] <= threshold for allele in alleles)


def gates():
    """The gates every structural claim is reported against."""
    yield 'recorded_500nm', ('best_nm', GATE_NM)
    for rank in RANK_GATES:
        yield f'rank_{str(rank).replace(".", "p")}_percent', ('best_pct', rank)


HAPLOTYPE_PAIRS = {'A': ('HLA-A*01:01', 'HLA-A*33:03'),
                   'B': ('HLA-B*08:01', 'HLA-B*15:08'),
                   'C': ('HLA-C*01:02', 'HLA-C*07:02')}


def escape_rows(data, plan):
    for name, gate in gates():
        baseline = sum(1 for identity in plan if covered(data, identity, ALLELES, gate))
        scenarios = [('single_allele_loss', (allele,)) for allele in ALLELES]
        scenarios += [('haplotype_loss', combination) for combination in itertools.product(
            HAPLOTYPE_PAIRS['A'], HAPLOTYPE_PAIRS['B'], HAPLOTYPE_PAIRS['C'])]
        for scenario, combination in scenarios:
            kept = [allele for allele in ALLELES if allele not in combination]
            lost = [data['regions'][i]['gene'] for i in plan
                    if covered(data, i, ALLELES, gate)
                    and not covered(data, i, kept, gate)]
            yield dict(gate=name, scenario=scenario,
                       lost_alleles=';'.join(combination),
                       retained_alleles=';'.join(kept),
                       baseline_regions_covered=baseline,
                       regions_covered=baseline - len(lost), regions_lost=len(lost),
                       lost_region_genes=';'.join(sorted(lost)))


def zero_coverage_regions(data, plan):
    """Regions inside no allele's gate, which a coverage count only implies."""
    result = {}
    for name, gate in gates():
        result[name] = sorted(data['regions'][i]['gene'] for i in plan
                              if not covered(data, i, ALLELES, gate))
    return result


def single_allele_rows(data, plan):
    for name, gate in gates():
        for identity in plan:
            region = data['regions'][identity]
            inside = [a for a in ALLELES if data[gate[0]][(identity, a)] <= gate[1]]
            if len(inside) != 1:
                continue
            others = [a for a in ALLELES if a != inside[0]]
            # The next-best allele by affinity and by percentile need not be the
            # same allele, so each is reported with the allele it belongs to.
            by_nm = min(others, key=lambda a: data['best_nm'][(identity, a)])
            by_pct = min(others, key=lambda a: data['best_pct'][(identity, a)])
            yield dict(gate=name, antigen_gene=region['gene'], source=region['source'],
                       protein_change=region['protein_change'],
                       region_id=region['region_id'], only_allele=inside[0],
                       only_allele_affinity_nm=round(data['best_nm'][(identity, inside[0])], 3),
                       only_allele_percentile=round(data['best_pct'][(identity, inside[0])], 4),
                       next_best_allele_by_affinity=by_nm,
                       next_best_affinity_nm=round(data['best_nm'][(identity, by_nm)], 3),
                       next_best_affinity_percentile=round(data['best_pct'][(identity, by_nm)], 4),
                       next_best_allele_by_percentile=by_pct,
                       next_best_percentile=round(data['best_pct'][(identity, by_pct)], 4),
                       next_best_percentile_affinity_nm=round(data['best_nm'][(identity, by_pct)], 3),
                       next_best_alleles_agree=by_nm == by_pct)


def ranks(values):
    order = sorted(values, key=lambda item: (-item[1], item[0]))
    return {identity: position for position, (identity, _) in enumerate(order, start=1)}


def rescoring_rows(data, plan):
    weights = {i: recorded_weight(data, i) for i in plan}
    scored = [i for i in plan if weights[i] is not None]
    recorded_rank = ranks([(i, data['recorded'][i] * weights[i]) for i in scored])
    neutral_rank = ranks([(i, data['neutral'][i] * weights[i]) for i in scored])
    stressed_rank = ranks([(i, data['weighted_by_c'][0.20][i] * weights[i]) for i in scored])
    reference_ranks = {a: ranks([(i, data['by_reference'][a][i] * weights[i])
                                 for i in scored]) for a in ALLELES}
    for identity in plan:
        region = data['regions'][identity]
        weight = weights[identity]
        ranked = weight is not None
        best_allele = min(ALLELES, key=lambda a: data['best_nm'][(identity, a)])
        best_pct_allele = min(ALLELES, key=lambda a: data['best_pct'][(identity, a)])
        row = dict(
            antigen_gene=region['gene'], source=region['source'],
            protein_change=region['protein_change'], region_id=region['region_id'],
            annotation_id=identity, target_peptide_allele_pairs=data['peptides'][identity],
            in_recorded_score_population=ranked,
            excluded_reason='' if ranked else 'no recorded original_weighted_score',
            best_target_affinity_nm=round(min(data['best_nm'][(identity, a)]
                                              for a in ALLELES), 3),
            best_target_allele=best_allele, best_target_locus=locus(best_allele),
            best_target_percentile=round(min(data['best_pct'][(identity, a)]
                                             for a in ALLELES), 4),
            best_percentile_allele=best_pct_allele,
            best_percentile_locus=locus(best_pct_allele),
            allele_neutral_best_affinity_nm=round(min(data['best_neutral'][(identity, a)]
                                                      for a in ALLELES), 3),
            recorded_contribution_sum=round(data['recorded'][identity], 6),
            allele_neutral_contribution_sum=round(data['neutral'][identity], 6),
            expression_weight=round(weight, 6) if ranked else '',
            recorded_weighted_score=round(data['recorded'][identity] * weight, 6) if ranked else '',
            allele_neutral_weighted_score=round(data['neutral'][identity] * weight, 6) if ranked else '',
            hla_c_0p20_weighted_score=round(data['weighted_by_c'][0.20][identity] * weight, 6) if ranked else '',
            recorded_rank=recorded_rank.get(identity, ''),
            allele_neutral_rank=neutral_rank.get(identity, ''),
            hla_c_0p20_rank=stressed_rank.get(identity, ''),
            rank_change_allele_neutral=(recorded_rank[identity] - neutral_rank[identity]
                                        if ranked else ''),
            rank_change_hla_c_0p20=(recorded_rank[identity] - stressed_rank[identity]
                                    if ranked else ''),
            worst_rank_change_over_reference_scales=(
                max(abs(recorded_rank[identity] - reference_ranks[a][identity])
                    for a in ALLELES) if ranked else ''),
            alleles_within_500nm=sum(1 for a in ALLELES
                                     if data['best_nm'][(identity, a)] <= GATE_NM),
            alleles_within_0p5_percentile=sum(1 for a in ALLELES
                                              if data['best_pct'][(identity, a)] <= 0.5),
            alleles_within_2_percentile=sum(1 for a in ALLELES
                                            if data['best_pct'][(identity, a)] <= 2.0))
        yield row


def reference_scale_rows(data, plan):
    """How far the order depends on which scale is called allele-neutral."""
    weights = {i: recorded_weight(data, i) for i in plan}
    scored = [i for i in plan if weights[i] is not None]
    neutral_rank = ranks([(i, data['neutral'][i] * weights[i]) for i in scored])
    recorded_rank = ranks([(i, data['recorded'][i] * weights[i]) for i in scored])
    for reference in ALLELES:
        order = ranks([(i, data['by_reference'][reference][i] * weights[i]) for i in scored])
        moved = [i for i in scored if order[i] != neutral_rank[i]]
        worst = max(scored, key=lambda i: (abs(order[i] - neutral_rank[i]), i))
        yield dict(reference_scale=reference, regions_ranked=len(scored),
                   regions_moved_vs_geometric_mean=len(moved),
                   max_rank_change_vs_geometric_mean=abs(order[worst] - neutral_rank[worst]),
                   max_rank_change_vs_recorded=max(abs(order[i] - recorded_rank[i])
                                                   for i in scored),
                   largest_mover_vs_geometric_mean=data['regions'][worst]['gene'],
                   largest_mover_geometric_mean_rank=neutral_rank[worst],
                   largest_mover_reference_rank=order[worst])


def splice_rows(data):
    """Re-rank the eligible splice windows and genes under the three scales."""
    rows = []
    for identity, region in data['regions'].items():
        if region['source'] != 'SPLICE':
            continue
        weight = recorded_weight(data, identity)
        if weight is None:
            continue
        rows.append(dict(
            gene=region['gene'], window_id=region['window_id'], region_id=region['region_id'],
            annotation_id=identity, sequence=region['sequence'],
            in_design_3=region['current_plan'] == 'True',
            recorded_weighted_score=round(data['recorded'][identity] * weight, 6),
            allele_neutral_weighted_score=round(data['neutral'][identity] * weight, 6),
            hla_c_0p20_weighted_score=round(data['weighted_by_c'][0.20][identity] * weight, 6)))
    best = {}
    for key, label in (('recorded_weighted_score', 'recorded'),
                       ('allele_neutral_weighted_score', 'allele_neutral'),
                       ('hla_c_0p20_weighted_score', 'hla_c_0p20')):
        per_gene = {}
        for row in rows:
            if row['gene'] not in per_gene or row[key] > per_gene[row['gene']][key]:
                per_gene[row['gene']] = row
        order = sorted(per_gene.values(), key=lambda row: (-row[key], row['gene']))
        best[label] = [(row['gene'], row['window_id'], round(row[key], 4)) for row in order]
    for row in rows:
        for label, key in (('recorded', 'recorded_weighted_score'),
                           ('allele_neutral', 'allele_neutral_weighted_score'),
                           ('hla_c_0p20', 'hla_c_0p20_weighted_score')):
            gene_order = [gene for gene, _, _ in best[label]]
            row[label + '_gene_rank'] = gene_order.index(row['gene']) + 1
            row[label + '_is_best_window_for_gene'] = any(
                window == row['window_id'] for gene, window, _ in best[label]
                if gene == row['gene'])
    rows.sort(key=lambda row: (-row['recorded_weighted_score'], row['gene']))
    return rows, best


def sensitivity_rows(data, plan):
    weights = {i: recorded_weight(data, i) for i in plan}
    scored = [i for i in plan if weights[i] is not None]
    orders = {w: ranks([(i, data['weighted_by_c'][w][i] * weights[i]) for i in scored])
              for w in HLA_C_WEIGHTS}
    for identity in scored:
        region = data['regions'][identity]
        row = dict(antigen_gene=region['gene'], source=region['source'],
                   region_id=region['region_id'], annotation_id=identity)
        for weight in HLA_C_WEIGHTS:
            label = str(weight).replace('.', 'p')
            row['score_hla_c_' + label] = round(
                data['weighted_by_c'][weight][identity] * weights[identity], 6)
            row['rank_hla_c_' + label] = orders[weight][identity]
        row['rank_change_1p0_to_0p06'] = orders[1.00][identity] - orders[0.06][identity]
        yield row


def neutral_equivalent_percentile(calibration, affinity_nm):
    """The percentile at which the allele-neutral scale reads this affinity."""
    low, high = 1e-4, 100.0
    for _ in range(200):
        middle = (low + high) / 2
        if calibration.neutral_affinity(middle) < affinity_nm:
            low = middle
        else:
            high = middle
    return (low + high) / 2


def locus_rows(data, plan, calibration):
    """Coverage per allele at matched stringency on both scales."""
    equivalents = {nm: neutral_equivalent_percentile(calibration, nm)
                   for nm in (50.0, GATE_NM)}
    for allele in ALLELES:
        row = dict(allele=allele, locus=locus(allele),
                   percentile_at_50nm=round(calibration.percentile(allele, 50.0), 4),
                   percentile_at_500nm=round(calibration.percentile(allele, 500.0), 4),
                   nm_at_0p5_percentile=round(calibration.affinity_at(allele, 0.5), 1),
                   nm_at_2_percentile=round(calibration.affinity_at(allele, 2.0), 1),
                   regions_within_50nm=sum(1 for i in plan
                                           if data['best_nm'][(i, allele)] <= 50.0),
                   regions_within_500nm=sum(1 for i in plan
                                            if data['best_nm'][(i, allele)] <= GATE_NM))
        for nm, percentile in sorted(equivalents.items()):
            row[f'neutral_{int(nm)}nm_equivalent_percentile'] = round(percentile, 4)
            row[f'regions_within_neutral_{int(nm)}nm_equivalent_rank'] = sum(
                1 for i in plan if data['best_pct'][(i, allele)] <= percentile)
        for rank in RANK_GATES:
            row[f'regions_within_{str(rank).replace(".", "p")}_percentile'] = sum(
                1 for i in plan if data['best_pct'][(i, allele)] <= rank)
        row['regions_with_best_affinity_here'] = sum(
            1 for i in plan if min(ALLELES, key=lambda a: data['best_nm'][(i, a)]) == allele)
        row['regions_with_best_percentile_here'] = sum(
            1 for i in plan if min(ALLELES, key=lambda a: data['best_pct'][(i, a)]) == allele)
        yield row


BACKGROUND_SEED, BACKGROUND_DRAWS = 20261117, 4000


def locus_background_rows(data, plan, calibration):
    """What the locus split looks like before this selection is applied.

    The 23-of-30 HLA-C share only means something against a baseline. Three
    are available from the recorded tables: single arbitrary peptides, random
    groups of arbitrary peptides matched to the size of a region's target set,
    and the candidate regions this selection did not take.
    """
    scored = collections.defaultdict(dict)
    for (peptide, allele), percentile in calibration.observed.items():
        scored[peptide][allele] = (calibration.observed_nm[(peptide, allele)], percentile)
    complete = sorted(p for p, d in scored.items() if len(d) == len(ALLELES))

    def best_locus(peptides, index):
        return locus(min(((scored[p][a][index], a) for p in peptides
                          for a in ALLELES))[1])

    def shares(groups, index):
        counted = collections.Counter(best_locus(group, index) for group in groups)
        total = sum(counted.values())
        return {key: 100.0 * counted[key] / total for key in 'ABC'}, total

    peptide_groups = [(p,) for p in complete]
    target_peptides = collections.defaultdict(set)
    for row in read_rows(PEPTIDE_ANNOTATIONS):
        if row['annotation_id'] in data['regions'] and row['allele'] in ALLELES:
            target_peptides[row['annotation_id']].add(row['peptide'])
    group_size = round(sum(len(target_peptides[i]) for i in plan) / len(plan))
    sampler = random.Random(BACKGROUND_SEED)
    sampled = [sampler.sample(complete, group_size) for _ in range(BACKGROUND_DRAWS)]
    non_plan = sorted(i for i in data['regions'] if i not in set(plan))

    populations = [
        ('single_arbitrary_peptides', peptide_groups, ''),
        (f'random_groups_of_{group_size}_arbitrary_peptides', sampled,
         f'seed {BACKGROUND_SEED}, {BACKGROUND_DRAWS} draws, size matched to the '
         f'mean distinct target peptides per selected region'),
        ('non_selected_candidate_regions',
         [target_peptides[i] for i in non_plan if target_peptides[i]], ''),
        ('selected_regions', [target_peptides[i] for i in plan], ''),
    ]
    for name, groups, note in populations:
        by_nm, total = shares(groups, 0)
        by_pct, _ = shares(groups, 1)
        yield dict(population=name, groups=total, note=note,
                   **{f'best_affinity_locus_{k}_percent': round(by_nm[k], 1) for k in 'ABC'},
                   **{f'best_percentile_locus_{k}_percent': round(by_pct[k], 1) for k in 'ABC'})


def background_allele_shares(calibration):
    """Per-allele share of best-affinity calls on arbitrary peptides."""
    scored = collections.defaultdict(dict)
    for (peptide, allele), percentile in calibration.observed.items():
        scored[peptide][allele] = calibration.observed_nm[(peptide, allele)]
    complete = [d for d in scored.values() if len(d) == len(ALLELES)]
    counted = collections.Counter(min(ALLELES, key=lambda a: d[a]) for d in complete)
    return {a: round(100.0 * counted[a] / len(complete), 1) for a in ALLELES}, len(complete)


def build():
    require_assertions('Design #3 allele-calibration audit')
    config = check_score_config()
    peptide_check = check_peptide_scores()
    calibration = Calibration(REPO)
    data = collect(calibration)
    plan = plan_identities(data)
    if len(plan) != 30:
        raise ValueError(f'expected the 30-region plan, found {len(plan)}')
    objective = check_recorded_objective(data)
    case_alleles = check_case_alleles()

    coverage = sorted(coverage_matrix(data, plan),
                      key=lambda row: (row['antigen_gene'], row['allele']))
    escape = list(escape_rows(data, plan))
    single = list(single_allele_rows(data, plan))
    rescoring = sorted(rescoring_rows(data, plan),
                       key=lambda row: (row['recorded_rank'] == '', row['recorded_rank'],
                                        row['antigen_gene']))
    splice, splice_best = splice_rows(data)
    sensitivity = sorted(sensitivity_rows(data, plan), key=lambda row: row['rank_hla_c_1p0'])
    calib_rows = list(calibration_table(calibration))
    loci = list(locus_rows(data, plan, calibration))
    references = list(reference_scale_rows(data, plan))
    background = list(locus_background_rows(data, plan, calibration))
    background_alleles, background_peptides = background_allele_shares(calibration)
    zero_coverage = zero_coverage_regions(data, plan)

    ranked = [row for row in rescoring if row['in_recorded_score_population']]
    unranked = [dict(antigen_gene=row['antigen_gene'], source=row['source'],
                     region_id=row['region_id'],
                     recorded_contribution_sum=row['recorded_contribution_sum'])
                for row in rescoring if not row['in_recorded_score_population']]
    locus_by_affinity = collections.Counter()
    locus_by_percentile = collections.Counter()
    for identity in plan:
        locus_by_affinity[locus(min(ALLELES, key=lambda a: data['best_nm'][(identity, a)]))] += 1
        locus_by_percentile[locus(min(ALLELES, key=lambda a: data['best_pct'][(identity, a)]))] += 1
    gate_flips_out = [row['antigen_gene'] for row in rescoring
                      if row['best_target_affinity_nm'] <= GATE_NM
                      < row['allele_neutral_best_affinity_nm']]
    gate_flips_in = [row['antigen_gene'] for row in rescoring
                     if row['allele_neutral_best_affinity_nm'] <= GATE_NM
                     < row['best_target_affinity_nm']]
    by_gate = collections.defaultdict(list)
    for row in escape:
        by_gate[row['gate']].append(row)
    stringency = {row['allele']: row['percentile_at_500nm'] for row in loci}

    summary = dict(
        status='AUDIT_COMPLETE_NO_SELECTION_CHANGE_PROPOSED',
        scope=('Offline re-scoring of recorded Design #3 annotations. No new inference, '
               'no construct, no change to shared evidence or to the selected regions.'),
        design='designs/design-3', regions_audited=len(plan),
        case_alleles=case_alleles, case_allele_source=HLA_CLASS_I,
        recorded_objective=dict(
            formula='logistic(affinity_nm; midpoint=350, width=150) / logistic(0), '
                    'zeroed below 1e-5 or at/above 5000 nM, summed over target '
                    'peptide/allele pairs, times the recorded expression weight',
            config=config, per_peptide=peptide_check, per_region=objective),
        calibration=dict(
            sources=len(CALIBRATION_SOURCES),
            recorded_peptide_allele_percentiles=len(calibration.observed),
            affinity_points_per_allele={a: len(calibration.by_affinity[a]) for a in ALLELES},
            raw_rows_per_allele=calibration.raw_row_count(),
            monotonicity_inversions_on_raw_rows=calibration.inversions(),
            tied_affinity_percentile_disagreements=calibration.tied_affinity_disagreements(),
            worst_relative_percentile_spread_within_one_affinity=float(
                f'{calibration.worst_tied_affinity_spread():.3e}'),
            length_dependence=('Close to, but not exactly, length-independent; the '
                               'per-length columns of allele-calibration.tsv differ, '
                               'most visibly for HLA-B*08:01 at 50 nM. Every target '
                               'percentile used here is a recorded value, so the '
                               'per-length curves are never relied on.'),
            percentile_at_500nm=stringency,
            nm_at_2_percentile={row['allele']: row['nm_at_2_percentile'] for row in loci},
            nm_gate_stringency_spread=round(max(stringency.values())
                                            / min(stringency.values()), 2),
            target_rows_needing_interpolated_percentile=data['interpolated']),
        locus_concentration=dict(
            question='Does the per-allele nM offset explain the concentration on HLA-C?',
            answer=('No. The locus holding a region\'s best target peptide is identical on '
                    'the nM scale and on the percentile scale, and HLA-C still leads at '
                    'every threshold tested. The offset does flatten coverage counts at '
                    'the gate, mostly by redistributing within HLA-A and HLA-B. HLA-C\'s '
                    'lead is widest at the design\'s own 500 nM gate and narrowest at '
                    'the allele-neutral and 2 percent rank gates.'),
            background=dict(
                note=('The 23-of-30 share is above every available baseline, but most of '
                      'it is already present before this selection is applied, so the '
                      'excess attributable to the selection is smaller than the raw '
                      'share suggests.'),
                arbitrary_peptides_scored_for_all_six_alleles=background_peptides,
                best_affinity_share_by_allele_on_arbitrary_peptides=background_alleles,
                populations=background),
            regions_by_best_affinity_locus=dict(sorted(locus_by_affinity.items())),
            regions_by_best_percentile_locus=dict(sorted(locus_by_percentile.items())),
            most_lenient_allele_at_500nm=max(stringency, key=stringency.get),
            per_allele=loci),
        structural_dependence=dict(
            gates_reported=[name for name, _ in gates()],
            note=('Which regions rest on a single allele is gate-specific, so every gate '
                  'is reported rather than one being preferred.'),
            single_allele_regions_by_gate={
                gate: sorted(row['antigen_gene'] for row in single if row['gate'] == gate)
                for gate in sorted({row['gate'] for row in single})},
            baseline_regions_covered={
                gate: rows[0]['baseline_regions_covered'] for gate, rows in by_gate.items()},
            regions_inside_no_allele_gate=zero_coverage,
            zero_coverage_note=('A region can leave the single-allele set either by '
                                'gaining a second allele or by failing the gate '
                                'entirely; PLEC and SYT12 leave it between the 2 and '
                                '0.5 percent gates by failing entirely.'),
            worst_single_allele_loss={
                gate: min(row['regions_covered'] for row in rows
                          if row['scenario'] == 'single_allele_loss')
                for gate, rows in by_gate.items()},
            haplotype_loss_range={
                gate: [min(row['regions_covered'] for row in rows
                           if row['scenario'] == 'haplotype_loss'),
                       max(row['regions_covered'] for row in rows
                           if row['scenario'] == 'haplotype_loss')]
                for gate, rows in by_gate.items()},
            haplotype_configurations_tested=sum(
                1 for row in escape if row['scenario'] == 'haplotype_loss'
                and row['gate'] == 'recorded_500nm')),
        allele_neutral_rescoring=dict(
            population=('the regions carrying a recorded original_weighted_score; the rest '
                        'are reported with contribution sums only, as the design left them'),
            regions_ranked=len(ranked), regions_not_ranked=unranked,
            rank_change_convention='recorded_rank minus alternative rank; '
                                   'positive means the region moves up the order',
            largest_rank_rise=max(row['rank_change_allele_neutral'] for row in ranked),
            largest_rank_fall=min(row['rank_change_allele_neutral'] for row in ranked),
            regions_unmoved=sum(1 for row in ranked
                                if row['rank_change_allele_neutral'] == 0),
            total_contribution_recorded=round(sum(data['recorded'][i] for i in plan), 4),
            total_contribution_allele_neutral=round(sum(data['neutral'][i] for i in plan), 4),
            passes_500nm_gate_only_via_allele_offset=gate_flips_out,
            newly_passes_500nm_gate_when_allele_neutral=gate_flips_in,
            allele_neutral_500nm_equivalent_percentile=round(
                neutral_equivalent_percentile(calibration, GATE_NM), 4),
            gate_flip_caveat=('500 nM on the allele-neutral scale is neither conventional '
                              'rank threshold. Under a 2 percent rank gate no region flips; '
                              'see escape-robustness.tsv and single-allele-dependence.tsv.'),
            scale_dependence=('The geometric-mean scale is one choice of neutral scale and '
                              'the order is not invariant to it; see '
                              'reference-scale-sensitivity.tsv.'),
            reference_scale_max_rank_change_vs_geometric_mean=max(
                row['max_rank_change_vs_geometric_mean'] for row in references),
            reference_scale_max_rank_change_vs_recorded=max(
                row['max_rank_change_vs_recorded'] for row in references)),
        splice_gene_order=dict(
            recorded=[gene for gene, _, _ in splice_best['recorded']],
            allele_neutral=[gene for gene, _, _ in splice_best['allele_neutral']],
            hla_c_0p20=[gene for gene, _, _ in splice_best['hla_c_0p20']],
            recorded_top5=[gene for gene, _, _ in splice_best['recorded'][:5]],
            allele_neutral_top5=[gene for gene, _, _ in splice_best['allele_neutral'][:5]],
            hla_c_0p20_top5=[gene for gene, _, _ in splice_best['hla_c_0p20'][:5]],
            windows_compared=len(splice), genes_compared=len(splice_best['recorded']),
            fifth_to_sixth_score_ratio=round(
                splice_best['recorded'][4][2] / splice_best['recorded'][5][2], 3),
            strength_of_test=('Weak by construction: ten genes, and the fifth-to-sixth gap '
                              'is wide, so the top five survive most monotone rescalings. '
                              'These scales are also not independent of each other.')),
        hla_c_sensitivity=dict(
            weights_scanned=list(HLA_C_WEIGHTS),
            applied_on='the allele-neutral scale',
            basis='Apps et al. 2015 J Immunol 194:3594-3600 (doi:10.4049/jimmunol.1403234): '
                  'in normal peripheral blood lymphocytes of one haplotype, HLA-A was 5x and '
                  'HLA-B 4x the level of HLA-C by mass spectrometry of total cellular protein, '
                  'and HLA-A 13-16x and HLA-B 17-18x by flow cytometry of surface protein. '
                  'Not measured for this tumour.',
            rank_change_convention='rank at weight 1.00 minus rank at weight 0.06; '
                                   'positive means the region moves up when HLA-C is '
                                   'down-weighted',
            largest_rank_rise=max(row['rank_change_1p0_to_0p06'] for row in sensitivity),
            largest_rank_fall=min(row['rank_change_1p0_to_0p06'] for row in sensitivity)),
        limitations=[
            'Predicted binding only. No processing, presentation, immunogenicity or '
            'T-cell measurement, and no experimental validation of any peptide.',
            'MHCflurry percentile ranks are internal to the pinned model bundle. They '
            'normalise its nM scale across alleles; they are not measured ligand yields, '
            'and the 0.5 and 2 percent thresholds are NetMHCpan conventions applied to a '
            'different predictor whose reference distribution is its own.',
            'Rank normalisation assumes each allele presents a comparable share of the '
            'peptide universe. Real repertoire breadth differs between alleles, so part '
            'of the nM offset may be biological rather than a scale artefact.',
            'No single scale is uniquely allele-neutral. The geometric mean used here is '
            'a choice, and the resulting order is not invariant to it.',
            'HLA-C down-weighting uses published ratios from normal lymphocytes of a '
            'different haplotype. Allele-specific surface density was not measured for '
            'Yash, and the bundled vendor expression is gene-level RNA, neither '
            'allele-resolved nor a surface-protein measurement.',
            'HLA loss of heterozygosity was not assayed. Haplotype-loss rows enumerate '
            'the four possible A/B/C phasings of the recorded alleles as hypotheses.',
            'Crossing a gate removes a region from a count, not its peptides. Next-best '
            'affinities and percentiles are in single-allele-dependence.tsv.',
            'Target peptide sets, flanks, the 25-aa windows and the expression proxies '
            'are the recorded ones; only the per-allele nM scale and the HLA-C weight '
            'are varied.'],
    )

    outputs = {
        'allele-calibration.tsv': calib_rows,
        'locus-coverage.tsv': loci,
        'locus-background.tsv': background,
        'coverage-matrix.tsv': coverage,
        'escape-robustness.tsv': escape,
        'single-allele-dependence.tsv': single,
        'region-rescoring.tsv': rescoring,
        'reference-scale-sensitivity.tsv': references,
        'splice-window-rescoring.tsv': splice,
        'hla-c-sensitivity.tsv': sensitivity,
    }
    manifest = {relative: dict(sha256=sha256(REPO / relative),
                               bytes=(REPO / relative).stat().st_size)
                for relative in sorted(set(CALIBRATION_SOURCES)
                                       | {REGION_ANNOTATIONS, PEPTIDE_ANNOTATIONS,
                                          SELECTED, HLA_CLASS_I, SCORE_CONFIG,
                                          CONTRIBUTIONS})}
    return outputs, summary, manifest


def write(outputs, summary, manifest):
    for name, rows in outputs.items():
        write_tsv(RESULTS / name, rows)
    write_json(RESULTS / 'inputs.manifest.json', manifest)
    write_json(RESULTS / 'summary.json', summary)


def verify(outputs, summary, manifest):
    problems = []
    for name, rows in outputs.items():
        path = RESULTS / name
        if not path.is_file():
            problems.append('missing ' + name)
            continue
        buffer = io.StringIO()
        fields = list(dict.fromkeys(key for row in rows for key in row))
        writer = csv.DictWriter(buffer, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
        with path.open(newline='') as stream:
            recorded_text = stream.read()
        if buffer.getvalue() != recorded_text:
            problems.append('changed ' + name)
    for name, value in (('inputs.manifest.json', manifest), ('summary.json', summary)):
        path = RESULTS / name
        if not path.is_file():
            problems.append('missing ' + name)
        elif json.loads(path.read_text()) != json.loads(json.dumps(value, default=str)):
            problems.append('changed ' + name)
    if problems:
        raise ValueError('recorded results disagree with a fresh offline run: '
                         + ', '.join(problems))
    return dict(status='VERIFIED', files=len(outputs) + 2,
                regions=summary['regions_audited'],
                new_inference=False, shared_evidence_changed=False)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--write', action='store_true',
                        help='refresh audit/ instead of checking it')
    args = parser.parse_args(argv)
    outputs, summary, manifest = build()
    if args.write:
        write(outputs, summary, manifest)
        print(json.dumps(dict(status='WRITTEN', files=len(outputs) + 2,
                              regions=summary['regions_audited'])))
    else:
        print(json.dumps(verify(outputs, summary, manifest)))


if __name__ == '__main__':
    main()
