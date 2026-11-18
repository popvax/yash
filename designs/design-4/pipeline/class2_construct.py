"""Resolve target Class-II predictions on the active construct's exact flanks."""
from pathlib import Path
from Bio import SeqIO
from curate import HERE
from dataio import rows, load
from protein import assembled_class2_context


def target_class2(target_root):
    root = Path(target_root)
    selected = {r['gene']: r for r in load(root / 'selection/selected.json')}
    components = {r['gene']: r for r in rows(root / 'construct/protein/components.tsv') if r['gene']}
    sequence = str(SeqIO.read(root / 'construct/protein/cassette.fasta', 'fasta').seq)
    context = {(r['Peptide'], r['Context']): r for r in rows(root / 'immunogenicity/mixmhc2-context.tsv')}
    no_context = {r['Peptide']: r for r in rows(root / 'immunogenicity/mixmhc2-no-context.tsv')}
    result = []
    for original in rows(HERE / 'immunogenicity/class2-target-scores.tsv'):
        if original['scenario'] != 'encoded' or original['gene'] not in selected:
            continue
        r = dict(original)
        st = int(components[r['gene']]['start_aa']) + int(r['start_in_region'])
        if sequence[st:st+15] != r['peptide']:
            raise ValueError('Class-II placement differs from selected sequence')
        ctx = assembled_class2_context(sequence, st, st+15)
        prediction = context[r['peptide'], ctx]  # Missing exact requests must fail.
        allele = r['allele']; p1 = int(prediction['CoreP1_'+allele])-1
        local = int(r['start_in_region'])+p1
        r.update(cassette_start_aa=st, context=ctx,
                 context_rank=float(prediction['%Rank_'+allele]),
                 no_context_rank=float(no_context[r['peptide']]['%Rank_'+allele]),
                 core_overlaps_recorded_target=(local < int(r['target_end_in_region']) and local+9 > int(r['target_start_in_region'])) if prediction['SubSpec_'+allele] != '-1' else None,
                 binding_subspecificity=prediction['SubSpec_'+allele],
                 core_coordinate_qualification='reverse_binding_interval_not_independently_interpreted' if prediction['SubSpec_'+allele]=='-1' else 'model_forward_CoreP1')
        result.append(r)
    reference = load(HERE / 'references/peptide-matches.json')
    summaries = []
    for gene, region in selected.items():
        eligible = [r for r in result if r['gene']==gene and (region['source']=='CTA/SELF' or not reference[r['peptide']]['occurs'])]
        if not eligible: continue
        best = min(eligible, key=lambda r: r['context_rank'])
        scored = context[best['peptide'], best['context']]; core = scored['Core_best']
        places = {i for candidate in (core,core[::-1]) for i in range(7) if best['peptide'].startswith(candidate,i)}
        overlaps = {int(best['start_in_region'])+i < int(region['target_end_in_window']) and int(best['start_in_region'])+i+9 > int(region['target_start_in_window']) for i in places}
        summaries.append(dict(gene=gene,class2_best_context_rank=float(scored['%Rank_best']),class2_best_core=core,
                              class2_best_core_overlaps_recorded_target=next(iter(overlaps)) if len(overlaps)==1 else None))
    return result, summaries


def verify_class2_annotations(root):
    """Check cached comparison annotations against exact construct predictions."""
    from dataio import cell
    from collections import defaultdict
    target, summaries = target_class2(root)
    eligible = defaultdict(list)
    for r in target:
        if r['context_rank'] <= 2 or r['no_context_rank'] <= 2:
            eligible[r['gene'], r['peptide']].append(r)
    for r in rows(Path(root)/'safety/wetlab-class2-wt-controls.tsv'):
        best = min(eligible[r['vaccine_gene'],r['mutant_15mer']],key=lambda x:(x['context_rank'],x['no_context_rank']))
        for field, source in [('best_class2_allele','allele'),('best_context_rank_pct','context_rank'),('best_no_context_rank_pct','no_context_rank'),('core_overlaps_recorded_target','core_overlaps_recorded_target')]:
            expected = cell(best[source])
            actual = r[field]
            if field.endswith('rank_pct'):
                matches = float(actual) == float(expected)
            else: matches = actual == expected
            if not matches: raise ValueError(f'Stale Class-II annotation: {r["vaccine_gene"]} {field}')
    meta = {r['gene']:r for r in summaries}
    for r in rows(Path(root)/'safety/class2-core-homologs.tsv'):
        if r['core_origin'] == 'TARGET_CLASS2_CORE':
            m=meta[r['vaccine_gene']]
            if r['query_core_9aa'] != m['class2_best_core'] or float(r['context_rank_pct']) != m['class2_best_context_rank']:
                raise ValueError('Stale construct-specific Class-II core annotation')


def eligible_core_rank(rank):
    """A valid zero percentile meets the reporting cutoff."""
    return rank is not None and 0 <= rank <= 2.0
