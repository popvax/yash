"""Build all junction costs, jointly optimize, independently verify exports."""
from yash_vaccine.paths import resolve_input
import collections,gzip,json
from pathlib import Path
import numpy as np
from ..io import rows,write_tsv,write_json,sha256
from .boundaries import junction_peptides,start_peptides,all_queries,assemble,scan_boundaries
from .optimize import solve,objective
from .inputs import inputs
from ..settings import frozen_recipe
from ..hydrophobicity import summarize, write_protein_screen

def run(repo,out,seconds=120,*,regions_path=None,cache_path=None,linkers_path=None):
    repo=Path(repo);out=Path(out);frozen=not any(x is not None for x in [regions_path,cache_path,linkers_path])
    regions,linkers,context,region_file,linker_file=inputs(repo,regions_path,linkers_path)
    cache=resolve_input(cache_path or repo/'data/checkpoints/construct-predictions.npz')
    if cache==resolve_input(repo/'data/checkpoints/construct-predictions.npz'):
        expected=json.loads((resolve_input(repo/'data/checkpoints/manifest.json')).read_text())[cache.name]['sha256']
        if sha256(cache)!=expected:raise ValueError('Frozen prediction cache checksum mismatch')
    if frozen and sha256(linker_file)!=frozen_recipe(repo)['linkers_sha256']:
        raise ValueError('Frozen replay expects the recorded linkers. Use explicit --regions, --linkers and a new --output for an experiment.')
    saved=np.load(resolve_input(cache));peptides=saved['peptides'].tolist();alleles=saved['alleles'].tolist()
    required=all_queries(regions,linkers);missing=sorted(set(required)-set(peptides))
    if missing:raise ValueError(f'Prediction cache is missing {len(missing)} required peptides. Run prepare-design and tools/launch_predictions.py for these regions/linkers, then supply --prediction-cache. Examples: {missing[:3]}')
    # A verified cache may be a superset after a construct-compatibility exclusion.
    # Coverage is checked above; required and cached counts remain distinct below.
    if alleles!=context['hla_class_I']:raise ValueError('Prediction-cache alleles do not match the recorded patient HLA-I types')
    if len(set(peptides))!=len(peptides):raise ValueError('Duplicate peptide keys in prediction cache')
    affinity=saved['affinity_nm'];rank=saved['affinity_percentile']
    assert affinity.shape==rank.shape==(len(peptides),len(alleles))
    assert np.isfinite(affinity).all() and (affinity>0).all() and np.isfinite(rank).all() and ((rank>=0)&(rank<=100)).all()
    prediction_counts=dict(required_peptides=len(required),required_peptide_hla_pairs=len(required)*len(alleles),
                           cache_peptides=len(peptides),cache_peptide_hla_pairs=len(peptides)*len(alleles))
    out.mkdir(parents=True,exist_ok=True)
    write_json(out/'design-inputs.json',dict(mode='frozen_replay' if frozen else 'new_design',
        regions=str(region_file.resolve()),regions_sha256=sha256(region_file),
        linkers=str(linker_file.resolve()),linkers_sha256=sha256(linker_file),
        prediction_cache=str(cache.resolve()),prediction_cache_sha256=sha256(cache),
        **prediction_counts,alleles=alleles))
    index={p:i for i,p in enumerate(peptides)};n=len(regions)
    def cost(queries,added):
        ix=[index[p] for p in queries]
        a=affinity[ix];r=rank[ix]
        return dict(strong_rank=int((r<=.5).sum()),mild_rank=int((r<=2).sum()),
                    affinity_500=int((a<=500).sum()),affinity_50=int((a<=50).sum()),
                    min_rank=float(r.min()) if ix else 100.,
                    min_affinity_nm=float(a.min()) if ix else '',
                    added_aa=added,peptide_occurrences=len(ix))
    starts=[cost(start_peptides(x['sequence']),int(not x['sequence'].startswith('M'))) for x in regions]
    for row, start in zip(regions, starts):
        start['hydrophobicity_pass']=summarize(('' if row['sequence'].startswith('M') else 'M')+row['sequence'])['status']=='PASS'
    edges=[]
    for i,x in enumerate(regions):
        for j,y in enumerate(regions):
            if i!=j:
                for linker in linkers:
                    edges.append(dict(left=i,right=j,left_id=x['antigen_id'],right_id=y['antigen_id'],
                                      linker=linker['name'],sequence=linker['sequence'],
                                      **cost(junction_peptides(x['sequence'],linker['sequence'],y['sequence']),len(linker['sequence']))))
    assert len(edges)==n*(n-1)*len(linkers)
    all_edges=len(edges)
    for edge in edges:
        screen=summarize(regions[edge['left']]['sequence']+edge['sequence']+regions[edge['right']]['sequence'])
        edge['hydrophobicity_pass']=screen['status']=='PASS'
    edges=[edge for edge in edges if edge['hydrophobicity_pass']]
    def export(label,order,chosen,meta):
        chosen_linkers=[dict(name=x['linker'],sequence=x['sequence']) for x in chosen]
        sequence,parts=assemble(regions,order,chosen_linkers)
        assert sorted(order)==list(range(n))
        assert [p['aa'] for p in parts if p['kind']=='antigen']==[regions[i]['sequence'] for i in order]
        assert all(sequence[p['start_aa']:p['end_aa']]==p['aa'] for p in parts)
        hydrophobicity=write_protein_screen(out,sequence,parts)
        pred=[]
        for row in scan_boundaries(sequence,parts):
            ix=index[row['peptide']]
            pred.extend(dict(**row,allele=h,affinity_nm=float(affinity[ix,j]),affinity_percentile=float(rank[ix,j])) for j,h in enumerate(alleles))
        totals={k:starts[order[0]][k]+sum(e[k] for e in chosen) for k in ['strong_rank','mild_rank','affinity_500','affinity_50','added_aa','peptide_occurrences']}
        assert sum(x['affinity_nm']<=500 for x in pred)==totals['affinity_500']
        assert sum(x['affinity_percentile']<=.5 for x in pred)==totals['strong_rank']
        assert sum(x['affinity_percentile']<=2 for x in pred)==totals['mild_rank']
        assert len(pred)==totals['peptide_occurrences']*len(alleles)
        assert len(sequence)==sum(len(x['sequence']) for x in regions)+totals['added_aa']
        assert min(x['affinity_percentile'] for x in pred)==-objective(order,chosen,starts)[2]
        (out/(label+'.fasta')).write_text('>'+label+'_RESEARCH_DRAFT\n'+'\n'.join(sequence[p:p+80] for p in range(0,len(sequence),80))+'\n')
        write_tsv(out/(label+'-components.tsv'),parts)
        write_tsv(out/(label+'-junctions.tsv'),[dict(junction_index=i+1,**x) for i,x in enumerate(chosen)])
        # The uncompressed table stays small and is directly inspectable.
        write_tsv(out/(label+'-boundary-predictions.tsv'),pred)
        result=dict(design=label,regions=n,antigen_residues=sum(len(x['sequence']) for x in regions),protein_length_aa=len(sequence),
                    order=[regions[i]['antigen_id'] for i in order],linkers=[x['linker'] for x in chosen],
                    linker_counts=dict(collections.Counter(x['linker'] for x in chosen)),
                    boundary_totals_including_start_M=totals,
                    objective_tuple=objective(order,chosen,starts),
                    hydrophobicity=hydrophobicity,
                    all_antigens_exactly_once=True,independent_full_chain_scan_verified=True,**meta)
        write_json(out/(label+'.json'),result)
        return result
    order,chosen,meta=solve(n,edges,starts,seconds)
    joint=export('joint',order,chosen,dict(optimization='joint antigen order and linker choice',solver=meta))
    write_json(out/'verification.json',dict(status='VERIFIED',antigens=n,directed_pairs=n*(n-1),linker_choices=len(linkers),
              all_pair_linker_combinations=all_edges,hydrophobicity_excluded_pair_linker_choices=all_edges-len(edges),
              hydrophobicity_excluded_starts=sum(not start['hydrophobicity_pass'] for start in starts),
              hydrophobicity_pass=True,**prediction_counts,
              cache_sha256=sha256(cache),joint_full_chain_verified=True,
              jointly_optimized_all_phases_proven_optimal=joint['solver']['all_phases_proven_optimal'],
              antigen_sequences_unchanged=True,scope='Protein antigen cassette; no CDS/UTR/LNP or cleavage/presentation optimization'))
    return joint
