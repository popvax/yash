"""Export the segmented mutation evidence used by the initial design.

Read completed stage outputs; retain native-PASS coding calls, the primary DNA
pool's RNA outcomes, and all scored windows/predictions for its usable contexts.
No variant calling, RNA reconstruction or model inference is performed here.
"""
import argparse
import json
import shutil
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'src'))
from yash_vaccine.io import rows, sha256, write_json, write_tsv


def table(path, records):
    write_tsv(path, [{k: json.dumps(v) if isinstance(v, (dict, list)) else v
                      for k, v in row.items()} for row in records])


def export(source, output):
    if output.exists():
        raise ValueError('Choose a new output directory')
    output.mkdir(parents=True)
    load = lambda name: json.loads((source / name).read_text())
    policy = {key: value for key, value in json.loads((REPO / 'data/mutation-evidence/policy.json').read_text()).items() if key in {
        'scope', 'native_filtering', 'DNA_screen', 'population_screen', 'manual_review_flags',
        'clinical_annotations', 'RNA_screen', 'downstream', 'context_review'}}
    write_json(output / 'policy.json', policy)
    evidence = []
    for row in load('annotate/variant-evidence.json'):
        model = row['model_evidence']['segmented']
        if model['effective_filter'] != 'PASS':
            continue
        row = {key: value for key, value in row.items() if key in {
            'variant_id', 'normalized_allele', 'gene', 'gene_id', 'protein_change', 'consequence',
            'transcript_id', 'other_effects', 'native_FILTER', 'effective_allele_filter',
            'DNA_DP', 'DNA_AD', 'DNA_model_AF', 'GERMQ', 'MBQ', 'MMQ', 'native_population_AF_prior',
            'Caris_exact_matches', 'Caris_flags', 'model_evidence', 'population', 'DNA_screen_failures',
            'population_exclusions', 'manual_review_flags', 'population_status', 'somatic_status',
            'primary_DNA_eligible', 'primary_DNA_eligible_by_model', 'QC_PASS_combined_global_AF',
            'RNA_review_eligible', 'evidence_tier'}}
        row['model_evidence'] = {'segmented': model}
        row['DNA_DP'], row['DNA_AD'] = model['DP'], model['AD']
        row['DNA_model_AF'], row['GERMQ'] = model['AF'][0], model['GERMQ']
        row['manual_review_flags'] = [v for v in row['manual_review_flags'] if v in policy['manual_review_flags']]
        good = [(c['cohort'].upper(), x) for c in row['population'] for x in c['records']
                if x['FILTER'] == ['PASS'] and x['AF'] is not None and x['AN']]
        population = set()
        for cohort, x in good:
            if x['AF'] >= policy['population_screen']['global_AF_exclusion_at_or_above']:
                population.add('COMMON_GLOBAL_' + cohort)
            if x['AF_grpmax'] is not None and x['AF_grpmax'] >= policy['population_screen']['ancestry_group_AF_exclusion_at_or_above']:
                population.add('COMMON_ANCESTRY_GROUP_' + cohort)
        row['population_exclusions'] = sorted(population)
        rule = policy['DNA_screen']
        failures = [name for name, bad in (
            ('LOW_DEPTH', row['DNA_DP'] is None or row['DNA_DP'] < rule['minimum_sample_depth']),
            ('LOW_ALT_AD', not row['DNA_AD'] or row['DNA_AD'][1] < rule['minimum_alternate_AD']),
            ('LOW_TUMOUR_AF', row['DNA_model_AF'] is None or row['DNA_model_AF'] < rule['minimum_model_tumour_AF']),
            ('LOW_GERMQ', row['GERMQ'] is None or row['GERMQ'] < rule['minimum_GERMQ'])) if bad]
        row['DNA_screen_failures'] = failures
        admitted = not (failures or population or row['manual_review_flags'])
        row['primary_DNA_eligible'] = admitted
        row['primary_DNA_eligible_by_model'] = {'segmented': admitted}
        row['RNA_review_eligible'] = admitted
        row['evidence_tier'] = ('EXCLUDED_POPULATION_COMMON' if population else 'REVIEW_DNA_QUALITY'
                                if failures else 'REVIEW_CONTEXT' if row['manual_review_flags'] else 'PUTATIVE_TUMOUR_ONLY')
        evidence.append(row)
    write_json(output / 'annotate/variant-evidence.json', evidence)
    table(output / 'annotate/variant-evidence.tsv', evidence)
    admitted = {r['variant_id'] for r in evidence if r['primary_DNA_eligible']}
    rna = [r for r in load('rna/rna-evidence.json') if r['variant_id'] in admitted]
    if {r['variant_id'] for r in rna} != admitted:
        raise ValueError('Some primary DNA variants lack RNA outcomes')
    proteins = [r for r in load('rna/primary-proteins.json') if r['variant_id'] in admitted]
    for r in proteins:
        r['primary_DNA_eligible'] = True
    protein_ids = {r['protein_id'] for r in proteins}
    for name, records in [('rna-evidence', rna), ('primary-proteins', proteins)]:
        write_json(output / ('rna/' + name + '.json'), records)
        table(output / ('rna/' + name + '.tsv'), records)
    (output / 'rna/primary-proteins.fasta').write_text(''.join(
        '>' + p['protein_id'] + ' ' + p['gene'] + '\n' + p['sequence'] + '\n' for p in proteins))
    for name in ('fragment-accounting.json', 'indel-representation-audit.json'):
        write_json(output / 'rna' / name, [r for r in load('rna/' + name) if r['variant_id'] in admitted])
    occurrences = [r for r in rows(source / 'rank/peptide-occurrences.tsv') if r['protein_id'] in protein_ids]
    peptides = {r['peptide'] for r in occurrences}
    table(output / 'rank/peptide-occurrences.tsv', occurrences)
    for name in ('peptide-HLA-affinity.tsv', 'affinity-contributions.tsv'):
        table(output / 'rank' / name, [r for r in rows(source / 'rank' / name) if r['peptide'] in peptides])
    table(output / 'rank/all-scored-windows.tsv', [r for r in rows(source / 'rank/all-scored-windows.tsv') if r['protein_id'] in protein_ids])
    write_json(output / 'rank/reference-matches.json', {p: r for p, r in load('rank/reference-matches.json').items() if p in peptides})
    rank_fields = {
        'variant_id', 'protein_id', 'gene', 'DNA_annotation_protein_change', 'DNA_annotation_transcript_id',
        'RNA_supported_transcript_ids', 'sequence', 'window_start_0based', 'window_end_exclusive',
        'target_start_in_window', 'target_end_in_window', 'affinity_contribution_sum', 'best_target_affinity_nm',
        'passes_500nm_gate', 'supplied_gene_value', 'rna_alternate_fraction', 'expression_proxy',
        'weighted_score', 'primary_DNA_eligible', 'gene_assignment_unambiguous', 'DNA_evidence_tier',
        'somatic_status', 'model_rank', 'filter_model', 'DNA_annotation_gene', 'DNA_DP', 'DNA_AD',
        'DNA_model_AF', 'model_GERMQ', 'updated_population_prior', 'population_status', 'primary_DNA_eligible_by_model'}
    rank_path = source / 'summarize/segmented-primary-ranked-regions.json'
    if not rank_path.is_file():
        rank_path = source / 'rank/segmented-primary-ranked-regions.json'
    ranked = [{k: v for k, v in r.items() if k in rank_fields} for r in json.loads(rank_path.read_text())]
    for r in ranked:
        r['primary_DNA_eligible_by_model'] = {'segmented': True}
    write_json(output / 'rank/segmented-primary-ranked-regions.json', ranked)
    table(output / 'rank/segmented-primary-ranked-regions.tsv', ranked)
    copied = ['case.json', 'analysis-genome-config.json', 'analysis-environment-verified.json',
              'exome-source-probe.json', 'filter/segmented/commands.json',
              'filter/segmented/filtering-stats.tsv', 'filter/segmented/complete.json',
              'annotate/population-source-metadata.json', 'rna/settings.json',
              'rank/score-config.json', 'rank/environment-packages.json',
              'rank/verified-affinity-model-files.json', 'rank/native-scoring-source-hashes.json']
    copied += [str(p.relative_to(source)) for p in sorted((source / 'population').glob('chr*/*.json'))]
    for name in copied:
        target = output / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source / name, target)
    plan = load('calling-plan.json')
    plan.pop('population_review', None)
    write_json(output / 'calling-plan.json', plan)
    write_json(output / 'source.json', dict(
        volume='research-lens-workspace',
        result_root='/data/results/TN25-112103/mutation-evidence-v3-20261102',
        native_filter_model='segmented', native_filter_receipt_sha256=sha256(source / 'filter/segmented/complete.json'),
        exported_native_PASS_coding_calls=len(evidence), primary_DNA_candidates=len(admitted),
        RNA_outcomes=len(rna), usable_RNA_protein_contexts=len(proteins), ranked_mutation_regions=len(ranked),
        source_annotation_sha256=sha256(source / 'annotate/variant-evidence.json'),
        source_RNA_sha256=sha256(source / 'rna/rna-evidence.json'),
        source_affinities_sha256=sha256(source / 'rank/peptide-HLA-affinity.tsv'),
        scope='Measurements and predictions exported for segmented eligibility and the primary RNA-supported mutation pool',
        variant_calling_rerun=False, model_inference_rerun=False))
    return dict(DNA=len(admitted), RNA_contexts=len(proteins), ranked_regions=len(ranked))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(export(args.source, args.output), indent=2))
