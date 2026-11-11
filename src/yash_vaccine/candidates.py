"""Rebuild compact regions and expression weights from frozen window evidence.

Native sequence recovery, reference checks and Vaxrank window scores were
computed upstream. Their executed sources and raw-replay entrypoints are under
workflow/stages. This module provides the small, offline checkpoint replay.
"""
import collections
import hashlib
import math
import shutil
from .io import rows, write_tsv, write_json
from .splice_evidence import verify as verify_splice
from .splice_checks import verify as verify_normal_screen


def patient_splice_windows(checkpoints, out):
    evidence = checkpoints.parent / 'splice-evidence'
    report = verify_splice(evidence, out / 'splice-rna-verification.json')
    normal = evidence / 'normal-screen'
    checked = verify_normal_screen(evidence, normal)
    assert checked['status'] == 'PASS'
    summary = __import__('json').loads((normal / 'summary.json').read_text())
    assert summary['normal_detection_required'] == 'NOT_DETECTED_IN_SEARCHED_REFERENCE'
    write_json(out / 'splice-frame-summary.json', dict(report, normal_screen=summary))
    shutil.copyfile(normal / 'gene-ranking.tsv', out / 'splice-gene-ranking.tsv')
    shutil.copyfile(normal / 'event-ranking.tsv', out / 'splice-event-ranking.tsv')
    retained = []
    for original in rows(normal / 'all-windows.tsv.gz'):
        if original['primary_window_eligible'] != 'True':
            continue
        row = dict(original)
        assert row['normal_detection'] == 'NOT_DETECTED_IN_SEARCHED_REFERENCE'
        value = row['supplied_transformed_gene_tpm']
        row.update(context_id=row['window_id'], window_start=0, window_end=25,
                   expression_proxy_left=value, expression_proxy_right=value,
                   weighted_score_right=row['weighted_score_left'],
                   expression_weighting_method='supplied_primary_gene_only; native_affinity_sum * sqrt(value)',
                   splice_fraction_used_for_weighting=False,
                   annotated_upstream_frame=True, patient_RNA_sequence_confirmed=True,
                   upstream_observed_stop_codons_excluded=True)
        retained.append(row)
    return retained


def inframe_fusion(checkpoints):
    """Explicit driver-junction inclusion, with no fabricated expression score."""
    import json
    from .io import sha256
    evidence = checkpoints.parent / 'fusion-evidence'
    admission = json.loads((evidence / 'inframe-selection.json').read_text())
    request = json.loads((evidence / 'inframe-request.json').read_text())
    review = json.loads((evidence / 'inframe-binding-review.json').read_text())
    assert admission['source_request_sha256'] == sha256(evidence / 'inframe-request.json')
    assert admission['source_prediction_sha256'] == sha256(evidence / 'inframe-binding-review.json')
    assert review['request_sha256'] == admission['source_request_sha256']
    assert review['model_files_verified'] and review['mhcflurry_version'] == '2.3.14'
    a, b = admission['source_protein_window_0based']
    sequence = request['source_sequence'][a:b]
    assert len(sequence) == 25
    best = review['best_nonself']
    assert best['peptide'] in sequence
    assert not any(r['exact_reference_match'] for r in review['predictions']
                   if r['peptide'] == best['peptide'])
    event = 'CDH1::NRG1:chr16:68819425:+--chr8:32727949:+:INFRAME'
    return dict(source='FUSION', event_id=event, gene='CDH1::NRG1',
        region_id='REG_' + hashlib.sha256((event+':'+str(a)+':'+sequence).encode()).hexdigest()[:20],
        sequence=sequence, context_id='CDH1_NRG1_INFRAME', window_start=a, window_end=b,
        native_affinity_sum='', weighted_score_left='', weighted_score_right='',
        expression_proxy_left='', expression_proxy_right='',
        manual_inclusion='True', admission_exception='DRIVER_JUNCTION_ABOVE_500NM',
        fusion_product='in-frame CDH1::NRG1 junction',
        best_target_affinity_nm=best['affinity_nm'], best_target_peptide=best['peptide'],
        best_target_allele=best['allele'], best_target_percentile_rank=best['percentile_rank'],
        rna_junction_reads=admission['original_junction_read_count'],
        rna_spanning_fragments=admission['original_spanning_fragment_count'],
        target_start_in_window=12, target_end_in_window=13, target_mode='junction_codon',
        expression_notes=admission['expression'])

def build_pool(checkpoints,out):
    best={};windows=0;window_fields=None
    for w in rows(checkpoints/'lens-scored-windows.tsv.gz'):
        if window_fields is None:window_fields=list(w)
        windows+=1
        if w['primary_window_eligible']!='True' or w['source']=='SPLICE':continue
        key=(w['source'],w['event_id'])
        if key not in best or (-float(w['native_affinity_sum']),w['window_id'])<(-float(best[key]['native_affinity_sum']),best[key]['window_id']):best[key]=w
    compatible_splice=patient_splice_windows(checkpoints,out)
    for w in compatible_splice:
        key=(w['source'],w['event_id'])
        if key not in best or (-float(w['native_affinity_sum']),w['window_id'])<(-float(best[key]['native_affinity_sum']),best[key]['window_id']):best[key]=w
    lenses=[]
    for key,w in sorted(best.items()):
        w=dict(w,current_construct_region='False')
        if w['source']=='SPLICE':
            w.update(weight_gene=w['gene'],splice_fraction_used_for_weighting='False')
        lenses.append(w)
    # Normalise mutation evidence directly instead of shipping a second table
    # containing the same sequences and scores. Keep the common export schema.
    mutation_fields=list(dict.fromkeys([
        'source','event_id','gene','region_id','sequence','native_affinity_sum',
        'weighted_score_left','weighted_score_right','current_construct_region',
        'prior_rank','research_only',*(window_fields or []),
        'coverage_order_within_event','new_target_peptides','new_target_count']))
    mutations=[]
    for m in rows(checkpoints/'mutation-regions.tsv'):
        score=float(m['post_self_affinity_sum'])
        expression=float(m['supplied_transformed_tpm'])*float(m['rna_alternate_fraction'])
        assert math.isfinite(score) and score>=0 and math.isfinite(expression) and expression>=0
        assert math.isclose(float(m['weighted_score']),score*math.sqrt(expression),rel_tol=1e-12,abs_tol=1e-12)
        row=dict.fromkeys(mutation_fields,'')
        row.update(source='MUTATION',event_id=m['variant_id'],gene=m['gene'],
            region_id='REG_'+hashlib.sha256(m['sequence'].encode()).hexdigest()[:20],
            sequence=m['sequence'],native_affinity_sum=m['post_self_affinity_sum'],
            weighted_score_left=m['weighted_score'],weighted_score_right=m['weighted_score'],
            current_construct_region='True',prior_rank=m['rank'],research_only='True',
            **{k:v for k,v in m.items() if k not in {'gene','sequence','variant_id','rank','weighted_score'}})
        mutations.append(row)
    pool=lenses+mutations+[inframe_fusion(checkpoints)]
    # The input files are frozen and hash-checked by the caller; selection quotas
    # and identity filtering do not change this upstream candidate evidence.
    assert all(len(r['sequence'])==25 for r in pool)
    write_tsv(out/'compact-pool.tsv',pool)
    write_json(out/'pool-summary.json',dict(rows=len(pool),windows_examined=windows,
        source_counts=dict(collections.Counter(r['source'] for r in pool)),
        splice_weighted=sum(r['source']=='SPLICE' and r['weighted_score_left']!='' for r in pool),
        splice_frame_policy='patient_RNA_annotated_frame',sequence_identity_filter=False,
        unique_sequences=len({r['sequence'] for r in pool})))
    return pool
