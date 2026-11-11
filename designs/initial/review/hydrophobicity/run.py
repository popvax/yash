"""Recount the construct screen and inspect alternate recorded mutation windows.

Run from the repository root. This reads included evidence only, does not change
the candidate, and does not add a second DNA-eligibility rule.
"""
import argparse
from pathlib import Path
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.ProtParamData import kd
from yash_vaccine.hydrophobicity import KD_TENTHS, POLICY, scan, summarize
from yash_vaccine.io import rows, write_json, write_tsv, sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Choose a new output directory')
    repo = Path(__file__).resolve().parents[4]
    pool_path = repo/'designs/initial/selection/candidate-pool.tsv'
    windows_path = repo/'data/mutation-evidence/rank/all-scored-windows.tsv'
    protein_path = repo/'designs/initial/construct/protein/cassette.fasta'
    full_windows_path = repo/'designs/initial/construct/protein/hydrophobicity-windows.tsv'
    pool = list(rows(pool_path))
    # Independent scale table and GRAVY calculation for every complete window.
    assert set(kd) == set(KD_TENTHS)
    assert all(round(kd[a]*10) == KD_TENTHS[a] for a in kd)
    checked = 0
    excluded = []
    for row in pool:
        for window in scan(row['sequence']):
            independent = ProteinAnalysis(window['sequence']).gravy()
            assert abs(independent-window['mean_kd']) < 1e-12
            checked += 1
        if summarize(row['sequence'])['status'] != 'PASS':
            excluded.append(row)
    excluded_variants = {r['event_id'] for r in excluded if r['source']=='MUTATION'}
    alternatives = [dict(variant_id=r['variant_id'],gene=r['gene'],sequence=r['sequence'],
                         window_start_in_context=r['window_start_0based'],**summarize(r['sequence']))
                    for r in rows(windows_path) if r['variant_id'] in excluded_variants]
    protein = ''.join(line for line in protein_path.read_text().splitlines() if not line.startswith('>'))
    full_windows = list(rows(full_windows_path))
    assert len(full_windows) == len(protein)-18
    for i,window in enumerate(full_windows):
        independent = ProteinAnalysis(protein[i:i+19]).gravy()
        assert int(window['start_aa']) == i and int(window['end_aa']) == i+19
        assert window['sequence'] == protein[i:i+19]
        assert abs(independent-float(window['mean_kd'])) < 1e-12
    assert not any(w['above_threshold']=='True' for w in full_windows
                   if w['classification']=='ANTIGEN_LINKER_CASSETTE')
    args.output.mkdir(parents=True)
    write_tsv(args.output/'excluded-mutation-windows.tsv', alternatives)
    write_json(args.output/'verification.json', dict(status='VERIFIED',policy=POLICY,
        candidate_pool_regions=len(pool),independently_checked_windows=checked,
        independently_checked_full_protein_windows=len(full_windows),
        full_protein_sha256=sha256(protein_path),
        independent_implementation='Biopython ProteinAnalysis.gravy and ProtParamData.kd',
        excluded_regions=[dict(gene=r['gene'],region_id=r['region_id']) for r in excluded],
        alternate_recorded_windows=len(alternatives),
        alternate_windows_passing=sum(r['status']=='PASS' for r in alternatives),
        automatic_window_rescue=False,
        input_sha256={str(p.relative_to(repo)):sha256(p) for p in [pool_path,windows_path,protein_path,full_windows_path]},
        code_sha256=sha256(Path(__file__))))


if __name__ == '__main__':
    main()
