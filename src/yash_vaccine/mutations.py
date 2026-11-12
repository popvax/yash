"""Replay segmented tumour-only mutation windows from recorded native evidence."""
from yash_vaccine.paths import resolve_input
import collections
import json
import math
from .io import rows, write_tsv, write_json, sha256


def build(repo, out):
    evidence = repo / 'data/mutation-evidence'
    from .mutation_policy import verify as verify_policy
    dna_eligible, policy_report = verify_policy(repo, out / 'mutation-policy')
    dna = {r['variant_id']: r for r in json.loads((resolve_input(evidence / 'annotate/variant-evidence.json')).read_text())}
    rna = {r['variant_id']: r for r in json.loads((resolve_input(evidence / 'rna/rna-evidence.json')).read_text())}
    missing_rna = sorted(vid for vid, models in dna_eligible.items()
                         if models.get('segmented') and vid not in rna)
    if missing_rna:
        raise ValueError('Newly eligible DNA variants need RNA assessment: ' + ', '.join(missing_rna))
    proteins = {r['protein_id']: r for r in json.loads((resolve_input(evidence / 'rna/primary-proteins.json')).read_text())}
    matches = json.loads((resolve_input(evidence / 'rank/reference-matches.json')).read_text())
    predictions = {(r['peptide'], r['allele']): float(r['affinity_nm'])
                   for r in rows(evidence / 'rank/peptide-HLA-affinity.tsv')}
    alleles = json.loads((resolve_input(repo / 'config/case.json')).read_text())['hla_class_I']
    occurrences = collections.defaultdict(list)
    for r in rows(evidence / 'rank/peptide-occurrences.tsv'):
        occurrences[r['protein_id']].append(r)
    by_protein = collections.defaultdict(list)
    checked = 0
    def contribution(nm):
        value = 0 if nm >= 5000 else (1 + math.exp(-350/150)) / (1 + math.exp((nm-350)/150))
        return value if value >= 1e-5 else 0.
    for w in rows(evidence / 'rank/all-scored-windows.tsv'):
        p = proteins[w['protein_id']]
        a, b = int(w['window_start_0based']), int(w['window_end_exclusive'])
        assert p['sequence'][a:b] == w['sequence'] and b-a == 25
        assert a < p['mutation_end_0based_exclusive'] and b > p['mutation_start_0based']
        inside = [o for o in occurrences[p['protein_id']] if a <= int(o['start']) and int(o['end']) <= b
                  and not matches[o['peptide']]['occurs']]
        values = [predictions[o['peptide'], h] for o in inside for h in alleles]
        score = sum(contribution(v) for v in values)
        assert math.isclose(score, float(w['affinity_contribution_sum']), rel_tol=1e-8, abs_tol=1e-8)
        best = min(values, default=None)
        assert best is None and not w['best_target_affinity_nm'] or best == float(w['best_target_affinity_nm'])
        by_protein[p['protein_id']].append(w)
        checked += 1
    result = []
    for pid, windows in by_protein.items():
        # File order retains native Vaxrank centring order for exact score ties.
        w = max(windows, key=lambda r: float(r['affinity_contribution_sum']))
        p = proteins[pid]; d = dna[w['variant_id']]; e = rna[w['variant_id']]
        if (not dna_eligible[w['variant_id']].get('segmented') or not p['gene_assignment_unambiguous']
                or not w['best_target_affinity_nm'] or float(w['best_target_affinity_nm']) > 500
                or not w['supplied_gene_value'] or float(w['supplied_gene_value']) <= 0):
            continue
        assert e['outcome'] == 'RNA_SUPPORTED_PROTEIN' and p['usable_mutant_context']
        assert e['rna_alt_fraction'] is not None and e['rna_alt_fraction'] > 0
        a, b = int(w['window_start_0based']), int(w['window_end_exclusive'])
        inside = [o for o in occurrences[pid] if a <= int(o['start']) and int(o['end']) <= b
                  and not matches[o['peptide']]['occurs']]
        best_o, best_h = min(((o,h) for o in inside for h in alleles), key=lambda x: predictions[x[0]['peptide'],x[1]])
        value = float(w['supplied_gene_value']); score = float(w['affinity_contribution_sum'])
        proxy = value * e['rna_alt_fraction']; model = d['model_evidence']['segmented']
        result.append(dict(variant_id=w['variant_id'], protein_id=pid, gene=w['gene'],
            protein_change=w['DNA_annotation_protein_change'].removeprefix('p.'),
            transcript_id=w['DNA_annotation_transcript_id'], normalized_allele=json.dumps(d['normalized_allele']),
            window_start_0based=a, window_end_0based_exclusive=b, sequence=w['sequence'], window_length=25,
            target_start_in_window=w['target_start_in_window'], target_end_in_window=w['target_end_in_window'],
            post_self_affinity_sum=score, weighted_score=score * math.sqrt(proxy),
            supplied_transformed_tpm=value, rna_alternate_fraction=e['rna_alt_fraction'],
            alternate_supporting_rna_fragments=e['rna_alt_fragments'], total_rna_fragments_union=e['rna_total_fragments'],
            expression_proxy=proxy, rna_multiplier=math.sqrt(proxy),
            best_affinity_nm=predictions[best_o['peptide'],best_h], best_affinity_peptide=best_o['peptide'], best_affinity_allele=best_h,
            DNA_filter_model='segmented', DNA_FILTER=model['effective_filter'], DNA_DP=model['DP'],
            DNA_reference_AD=model['AD'][0], DNA_alternate_AD=model['AD'][1], DNA_model_AF=model['AF'][0], GERMQ=model['GERMQ'],
            Caris_flags=','.join(d['Caris_flags']), somatic_status='UNCONFIRMED_NO_MATCHED_NORMAL',
            reference_proteome_status='CHECKED_ENSEMBL103_EXACT_MATCH',
            wt_affinity_policy='Exact reference-peptide exclusion; no wild-type affinity weighting'))
    result.sort(key=lambda r: (-r['weighted_score'], -r['post_self_affinity_sum'], r['variant_id']))
    assert len({r['variant_id'] for r in result}) == len(result)
    for i, r in enumerate(result, 1): r['rank'] = i
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / 'mutation-regions.tsv', result)
    write_json(out / 'mutation-score-verification.json', dict(status='VERIFIED', windows_checked=checked,
        model='segmented', ranked_mutation_regions=len(result),
        DNA_eligibility_recomputed_from_measurements=True,
        DNA_policy_sha256=policy_report['policy_sha256'],
        primary_DNA_candidates=policy_report['primary_DNA_eligible_by_model']['segmented'],
        exact_reference_peptides_excluded=True, native_affinity_sums_independently_recounted=True,
        weighting='binding contribution sum * sqrt(supplied transformed gene value * alternate RNA fragment fraction)',
        RNA_supported_sequences_verified=True, model_inference_rerun=False))
    return result


def verify(repo, out):
    build(repo, out)
    if sha256(out / 'mutation-regions.tsv') != sha256(repo / 'data/checkpoints/mutation-regions.tsv'):
        raise ValueError('Mutation checkpoint differs from the segmented DNA/RNA/affinity evidence')
