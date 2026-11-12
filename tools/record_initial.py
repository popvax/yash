"""Record verified initial-design outputs and their current source hashes.

This is a maintainer recording step, not a constructor or a substitute for
reproduction. It requires completed selection, junction and nucleotide checks.
Use independent design directories for alternative approaches.
"""
from yash_vaccine.paths import resolve_input
import argparse
import collections
import json
import os
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--execute', action='store_true', help='Write the current recording manifests')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    import sys
    sys.path.insert(0, str(root / 'src'))
    from yash_vaccine.io import rows, sha256, write_json
    initial = root / 'designs/initial'
    selected = list(rows(initial / 'selection/selected.tsv'))
    selection = json.loads((resolve_input(initial / 'selection/verification.json')).read_text())
    protein = json.loads((resolve_input(initial / 'construct/protein/construct.json')).read_text())
    junction = json.loads((resolve_input(initial / 'construct/final-check.json')).read_text())
    nucleotide = json.loads((resolve_input(initial / 'construct/mrna/verification.json')).read_text())
    assert selection['status'] == 'VERIFIED'
    assert junction['status'] == 'FRESH_FULL_JUNCTION_CHECK_PASSED'
    assert protein['solver']['all_phases_proven_optimal']
    assert protein['hydrophobicity']['status']=='PASS' and junction['hydrophobicity']==protein['hydrophobicity']
    assert nucleotide['all_constraints_pass'] and nucleotide['polyA_status'] == 'unspecified'
    assert nucleotide['genbank_annotation']['round_trip_annotations_identical']
    assert len(selected) == protein['antigen_regions'] == junction['antigen_regions']
    assert junction['protein_sha256'] == sha256(initial / 'construct/protein/cassette.fasta')
    assert nucleotide['input_sha256']['../protein/cassette.fasta'] == junction['protein_sha256']
    assert all(r['normal_detection'] == 'NOT_DETECTED_IN_SEARCHED_REFERENCE'
               for r in selected if r['source'] == 'SPLICE')
    print(json.dumps(dict(selected_regions=len(selected), protein_aa=protein['protein_length_aa'],
        sources=dict(collections.Counter(r['source'] for r in selected)), execute=args.execute), indent=2))
    if not args.execute:
        return
    frozen_path = root / 'config/frozen-replay.json'
    frozen = json.loads(resolve_input(frozen_path).read_text())
    recipe = json.loads((resolve_input(initial / 'selection/design.json')).read_text())
    pool = list(rows(initial / 'selection/candidate-pool.tsv'))
    frozen['selection'] = recipe['selection']
    frozen['expected'] = dict(compact_pool=len(pool),
        retained_splice=sum(r['source'] == 'SPLICE' for r in pool), selected_regions=len(selected),
        splice_genes=[r['gene'] for r in selected if r['source'] == 'SPLICE'])
    frozen['candidate_artifact_sha256'] = {
        name: sha256(initial / 'selection' / name) for name in frozen['candidate_artifact_sha256']}
    write_json(frozen_path, frozen)
    nucleotide_outputs = ['construct/mrna/' + name for name in (
        'cds.fasta', 'components.tsv', 'construct.json', 'mrna-segments.fasta', 'optimizer.log',
        'rna-layout.txt', 'template-segments.fasta', 'template-segments.gb', 'verification.json')]
    files = [p for p in resolve_input(initial).rglob('*') if resolve_input(p).is_file()
             and p.name != 'candidate.json' and '__pycache__' not in p.parts]
    for directory in ('src', 'workflow'):
        files.extend(p for p in (resolve_input(root / directory)).rglob('*')
                     if resolve_input(p).is_file() and p.suffix in {'.py', '.sh', '.json'} and '__pycache__' not in p.parts)
    files.extend(p for p in (resolve_input(root / 'tools')).glob('*.py'))
    files.extend((resolve_input(root / 'config')).glob('*.json'))
    files.extend((resolve_input(root / 'data/checkpoints/provenance')).glob('*.json'))
    files.extend([root / 'data/checkpoints/manifest.json', root / 'pyproject.toml',
                  root / 'data/mutation-evidence/policy.json',
                  root / 'data/mutation-evidence/eligibility/decisions.tsv',
                  root / 'data/mutation-evidence/eligibility/verification.json'])
    manifest = dict(schema_version=2, name='initial',
        description=f'Initial vaccine design: {len(selected)} regions, hydropathy-screened mutation/splice/CTA/fusion regions, ND-only splice admission, patient-supported routing, joint junction optimization and PopRNA encoding.',
        selected_regions=len(selected), protein_length_aa=protein['protein_length_aa'],
        rna_length_nt=None, template_context_length_bp=None, entrypoint='run.py',
        reproduce_command='yash-vaccine reproduce', install_command="python -m pip install -e '.[mrna]'",
        stages=dict(selection='selection/', normal_reference_admission='../../data/splice-evidence/normal-screen/',
                    protein_assembly='construct/', junction_verification='construct/final-check.json', RNA_encoding='construct/mrna/'),
        final_outputs=dict(selected_antigens='selection/selected.tsv', protein='construct/protein/cassette.fasta',
                           coding_DNA='construct/mrna/cds.fasta', RNA_segments='construct/mrna/mrna-segments.fasta',
                           annotated_DNA_segments='construct/mrna/template-segments.gb'),
        nucleotide_outputs=nucleotide_outputs,
        artifact_sha256={os.path.relpath(p, initial): sha256(p) for p in sorted(set(files))},
        reproduction_scope='Local rebuild from recorded measurements and predictions; no raw-read analysis, new inference, downloads or cloud jobs.',
        polyA=dict(name='PolyA', status='unspecified', length_nt=None, sequence_dna=None, display='AAAAA…AAAAA'),
        complete_RNA_template_sequence=False)
    write_json(initial / 'candidate.json', manifest)


if __name__ == '__main__':
    main()
