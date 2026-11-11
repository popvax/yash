"""Calculate segmented DNA eligibility from native calls and population measurements.

This local stage does not repeat Mutect2, population lookups, RNA reconstruction
or model inference.
"""
import collections
import json

from .io import sha256, write_json, write_tsv


def assess(repo, out):
    evidence = repo / 'data/mutation-evidence'
    policy_path = evidence / 'policy.json'
    policy = json.loads(policy_path.read_text())
    if ('minimum_alternate_median_mapping_quality' in policy['DNA_screen']
            or 'minimum_alternate_median_base_quality' in policy['DNA_screen']
            or 'minimum_group_AN' in policy['population_screen']):
        raise ValueError('Current policy must use native quality filtering and ancestry AF without a sampling minimum')
    source_path = evidence / 'annotate/variant-evidence.json'
    source = json.loads(source_path.read_text())
    rule = policy['DNA_screen']
    pop_rule = policy['population_screen']
    decisions = []
    eligible = {}
    native_counts = collections.Counter()
    eligible_counts = collections.Counter()
    for record in source:
        population_failures = set()
        for cohort in record['population']:
            for observation in cohort['records']:
                if (observation['FILTER'] != ['PASS'] or observation['AF'] is None
                        or not observation['AN']):
                    continue
                suffix = cohort['cohort'].upper()
                if observation['AF'] >= pop_rule['global_AF_exclusion_at_or_above']:
                    population_failures.add('COMMON_GLOBAL_' + suffix)
                if (observation['AF_grpmax'] is not None and observation['AF_grpmax']
                        >= pop_rule['ancestry_group_AF_exclusion_at_or_above']):
                    population_failures.add('COMMON_ANCESTRY_GROUP_' + suffix)
        context = sorted(set(record['manual_review_flags']))
        models = {}
        for name, model in sorted(record['model_evidence'].items()):
            ad, dp = model.get('AD'), model.get('DP')
            af = model.get('AF')
            if isinstance(af, list):
                af = af[0] if af else None
            germq = model.get('GERMQ')
            failures = []
            for label, fails in [
                ('NATIVE_FILTERED', model['effective_filter'] != 'PASS'),
                ('LOW_DEPTH', dp is None or dp < rule['minimum_sample_depth']),
                ('LOW_ALT_AD', not ad or len(ad) < 2 or ad[1] < rule['minimum_alternate_AD']),
                ('LOW_TUMOUR_AF', af is None or af < rule['minimum_model_tumour_AF']),
                ('LOW_GERMQ', germq is None or germq < rule['minimum_GERMQ']),
            ]:
                if fails:
                    failures.append(label)
            admitted = not (failures or population_failures or context)
            models[name] = admitted
            native_counts[name] += model['effective_filter'] == 'PASS'
            eligible_counts[name] += admitted
            decisions.append(dict(
                variant_id=record['variant_id'], gene=record['gene'],
                protein_change=record['protein_change'], model=name,
                native_FILTER=model['effective_filter'], DNA_DP=dp,
                DNA_alternate_AD=ad[1] if ad and len(ad) > 1 else None,
                DNA_model_AF=af, GERMQ=germq,
                DNA_screen_failures=';'.join(failures),
                population_exclusions=';'.join(sorted(population_failures)),
                manual_review_flags=';'.join(context), primary_DNA_eligible=admitted))
        eligible[record['variant_id']] = models
    decisions.sort(key=lambda row: (row['variant_id'], row['model']))
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / 'decisions.tsv', decisions)
    report = dict(
        status='VERIFIED', scope='Segmented DNA eligibility from native calls and population observations',
        source_protein_changing_variants=len(source),
        native_PASS_protein_changing_by_model=dict(sorted(native_counts.items())),
        primary_DNA_eligible_by_model=dict(sorted(eligible_counts.items())),
        extra_mapping_or_base_quality_cutoff=False,
        ancestry_frequency_minimum_sample_size=None,
        ancestry_frequency_cutoff=pop_rule['ancestry_group_AF_exclusion_at_or_above'],
        policy_sha256=sha256(policy_path),
        source_allele_evidence_sha256=sha256(source_path),
        source_sha256=sha256(__file__), decisions_sha256=sha256(out / 'decisions.tsv'),
        variant_calling_rerun=False, population_lookup_rerun=False,
        source_measurements_modified=False)
    write_json(out / 'verification.json', report)
    return eligible, report


def verify(repo, out):
    eligible, report = assess(repo, out)
    reference = repo / 'data/mutation-evidence/eligibility'
    for name in ('decisions.tsv', 'verification.json'):
        if sha256(out / name) != sha256(reference / name):
            raise ValueError('Current mutation-policy reassessment differs: ' + name)
    return eligible, report


if __name__ == '__main__':
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Choose a new output directory')
    _, result = assess(Path(__file__).resolve().parents[2], args.output.resolve())
    print(json.dumps(result, indent=2))
