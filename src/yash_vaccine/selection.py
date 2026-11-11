"""Source quotas, with one highest-ranked retained splice region per gene."""
import collections
from .io import rows,ranked,write_tsv,write_json
from .hydrophobicity import screen_pool, POLICY_ID

def select(pool,settings,checkpoints,out):
    if settings['hydrophobicity_policy'] != POLICY_ID:raise ValueError('Hydrophobicity policy mismatch')
    hydropathy=screen_pool(pool,out)
    mutations={r['variant_id']:r for r in rows(checkpoints/'mutation-regions.tsv')}
    selected=[];audit=[];sequences=set()
    for source,quota in settings['quotas'].items():
        seen=set();n=0
        for rank,row in enumerate(ranked([r for r in pool if r['source']==source]),1):
            if hydropathy[row['region_id']]['status'] != 'PASS':
                audit.append(dict(source=source,event_id=row['event_id'],gene=row['gene'],action='SKIP_HYDROPHOBICITY'))
                continue
            manual = source == 'FUSION' and row.get('manual_inclusion') == 'True'
            if row['weighted_score_left']=='' and not manual:
                audit.append(dict(source=source,event_id=row['event_id'],gene=row['gene'],action='SKIP_MISSING_EXPRESSION'))
                continue
            if n>=quota:
                audit.append(dict(source=source,event_id=row['event_id'],gene=row['gene'],action='SKIP_QUOTA_FILLED'))
                continue
            if source=='SPLICE' and settings['splice_one_region_per_gene'] and (not row['gene'] or row['gene'] in seen):
                audit.append(dict(source=source,event_id=row['event_id'],gene=row['gene'],action='SKIP_GENE_ALREADY_SELECTED' if row['gene'] else 'SKIP_MISSING_GENE'))
                continue
            if row['sequence'] in sequences:
                audit.append(dict(source=source,event_id=row['event_id'],gene=row['gene'],action='SKIP_DUPLICATE_SEQUENCE'))
                continue
            sequences.add(row['sequence'])
            n+=1;seen.add(row['gene']);m=mutations.get(row['event_id'],{})
            selected.append(dict(row,antigen_id=f'AG{len(selected)+1:02d}',source_rank=rank,source_slot=n,
                protein_change=m.get('protein_change',''),Caris_flags=m.get('Caris_flags',''),
                somatic_status=m.get('somatic_status',''),
                best_affinity_nm=row.get('best_target_affinity_nm','') or m.get('best_affinity_nm',''),
                best_affinity_allele=row.get('best_target_allele','') or m.get('best_affinity_allele','')))
            audit.append(dict(source=source,event_id=row['event_id'],gene=row['gene'],action='SELECTED'))
        if n!=quota:raise ValueError(f'{source}: requested {quota} regions but only {n} satisfy hydrophobicity, expression, uniqueness and gene constraints; inspect the pool or lower the quota')
    assert len(selected)==sum(settings['quotas'].values()) and len(sequences)==len(selected)
    splice=[r for r in selected if r['source']=='SPLICE']
    if settings['splice_one_region_per_gene']:assert len({r['gene'] for r in splice})==len(splice)
    write_tsv(out/'selected.tsv',selected)
    write_tsv(out/'selection-audit.tsv',audit)
    with (out/'selected.fasta').open('w') as f:
        for r in selected:f.write('>'+r['antigen_id']+'|'+r['source']+'|'+r['gene']+'|'+r['event_id']+'\n'+r['sequence']+'\n')
    write_json(out/'selection-summary.json',dict(total_regions=len(selected),antigen_residues=sum(len(r['sequence']) for r in selected),
        sources=dict(collections.Counter(r['source'] for r in selected)),
        splice_genes=[r['gene'] for r in splice],one_splice_region_per_gene=settings['splice_one_region_per_gene'],
        all_sequences_25aa=all(len(r['sequence'])==25 for r in selected),distinct_sequences=len(sequences),ordering_by_source_and_rank=True,
        hydrophobicity_policy=POLICY_ID,all_selected_pass_hydrophobicity=True))
    return selected
