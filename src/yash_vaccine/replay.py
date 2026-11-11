"""Frozen replay and configurable selection share the same verified input pool."""
import json,shutil
from pathlib import Path
from .candidates import build_pool
from .selection import select
from .io import sha256,write_json
from .verify_scores import verify
from .settings import checkpoint_context,frozen_recipe,load_selection,validate_selection


def execute(repo,out,settings,*,frozen=False):
    repo=Path(repo).resolve();out=Path(out).resolve()
    context=checkpoint_context(repo);validate_selection(settings)
    checkpoints=repo/'data/checkpoints';manifest=json.loads((checkpoints/'manifest.json').read_text())
    for name,meta in manifest.items():
        if sha256(checkpoints/name)!=meta['sha256']:raise ValueError('Checkpoint checksum mismatch: '+name)
    out.mkdir(parents=True,exist_ok=True)
    verify(checkpoints,out)
    from .mutations import verify as verify_mutations
    verify_mutations(repo,out)
    pool=build_pool(checkpoints,out)
    selected=select(pool,settings,checkpoints,out)
    write_json(out/'selection-settings.json',settings)
    write_json(out/'selection-provenance.json',dict(mode='frozen_replay' if frozen else 'new_selection',
        selection_settings=settings,checkpoint_context_sha256=sha256(repo/'config/case.json'),
        checkpoint_manifest_sha256=sha256(checkpoints/'manifest.json'),
        checkpoint_hashes_verified=True,region_length_aa=context['region_length_aa'],
        selected_regions_sha256=sha256(out/'selected.tsv'),model_inference_rerun=False))
    report=dict(status='VERIFIED',compact_pool=len(pool),sequence_identity_filter_applied=False,
        final_regions=len(selected),splice_genes=[r['gene'] for r in selected if r['source']=='SPLICE'],
        checkpoint_hashes_verified=True,scope='Selection from verified frozen sequence/score evidence; no raw-read analysis or new model inference')
    if frozen:
        recipe=frozen_recipe(repo);expected=recipe['expected']
        assert len(pool)==expected['compact_pool']
        assert len(selected)==expected['selected_regions'] and report['splice_genes']==expected['splice_genes']
        assert sum(r['source']=='SPLICE' for r in pool)==expected['retained_splice']
        aliases={'selected.tsv':f'selected-{len(selected)}.tsv','selected.fasta':f'selected-{len(selected)}.fasta',
                 'compact-pool.tsv':'candidate-pool.tsv'}
        for source,dest in aliases.items():shutil.copyfile(out/source,out/dest)
        for name,expected_hash in recipe['candidate_artifact_sha256'].items():
            if sha256(out/name)!=expected_hash:raise ValueError('Frozen artifact differs: '+name)
        report['frozen_artifacts_match']=True
    write_json(out/('replay-verification.json' if frozen else 'selection-verification.json'),report)
    return selected


def run(repo,out):
    """Always replay the immutable recipe, regardless of editable selection.json."""
    return execute(repo,out,frozen_recipe(repo)['selection'],frozen=True)


def experiment(repo,out,settings_path=None):
    out=Path(out)
    if out.exists():raise ValueError('Choose a new selection output directory: '+str(out))
    settings=load_selection(settings_path) if settings_path else json.loads((Path(repo)/'designs/initial/selection/design.json').read_text())['selection']
    return execute(repo,out,settings)
