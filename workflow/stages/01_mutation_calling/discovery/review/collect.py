"""Export readable completed results; large VCFs, references and reads stay on Modal."""
from pathlib import Path
import argparse
import concurrent.futures
import hashlib
import json

HERE = Path(__file__).resolve().parent
REMOTE = 'results/TN25-112103/mutation-evidence-v3-20261102'
FILES = {
    'setup': ['commands.json'],
    'filter/segmented': ['commands.json', 'filtering-stats.tsv'],
    'gather': [],
    'annotate': ['variant-evidence.json', 'variant-evidence.tsv', 'rna-input.json', 'population-source-metadata.json', 'annotation-errors.json'],
    'rna': ['rna-evidence.json', 'rna-evidence.tsv', 'primary-proteins.json', 'primary-proteins.tsv', 'primary-proteins.fasta', 'indel-representation-audit.json', 'fragment-accounting.json', 'settings.json'],
    'rank': ['all-ranked-regions.tsv', 'primary-ranked-regions.tsv', 'primary-ranked-regions.json', 'primary-regions.fasta',
             'all-scored-windows.tsv', 'coverage-exclusions.tsv', 'best-peptide-per-HLA.tsv', 'peptide-HLA-affinity.tsv',
             'peptide-occurrences.tsv', 'reference-matches.json', 'affinity-contributions.tsv', 'score-config.json',
             'environment-packages.json', 'verified-affinity-model-files.json', 'native-scoring-source-hashes.json'],
    'summarize': ['segmented-primary-ranked-regions.json', 'segmented-primary-ranked-regions.tsv', 'segmented-regions.fasta', 'summary.json'],
    'audit': ['checks.json']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--execute', action='store_true')
    parser.add_argument('--output', type=Path, default=HERE / 'exports')
    args = parser.parse_args()
    print(json.dumps(dict(execute=args.execute, remote_root=REMOTE, output=str(args.output), stages=list(FILES)), indent=2), flush=True)
    if not args.execute:
        return
    import modal
    volume = modal.Volume.from_name('research-lens-workspace')
    jobs = []
    for stage, names in FILES.items():
        data = b''.join(volume.read_file(REMOTE + '/' + stage + '/complete.json'))
        completed = json.loads(data)
        assert completed['status'] == 'COMPLETED'
        jobs.append((stage + '/complete.json', None))
        for name in names:
            key = 'setup/' + name if stage == 'setup' else name
            jobs.append((stage + '/' + name, completed.get('outputs', {}).get(key, {}).get('sha256')))
    for c in list(range(1, 23)) + ['X', 'Y']:
        jobs.append(('population/chr' + str(c) + '/complete.json', None))
        jobs.append(('population/chr' + str(c) + '/exomes-source.json', None))
    for name in ['policy.json', 'case.json', 'calling-plan.json', 'analysis-genome-config.json', 'analysis-environment-verified.json', 'exome-source-probe.json']:
        jobs.append((name, None))

    def download(job):
        name, expected = job
        body = b''.join(volume.read_file(REMOTE + '/' + name))
        digest = hashlib.sha256(body).hexdigest()
        if expected:
            assert digest == expected, name
        path = args.output / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(body)
        return dict(file=name, bytes=len(body), sha256=digest, receipt_digest_checked=bool(expected), remote_path='/data/' + REMOTE + '/' + name)

    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        manifest = list(pool.map(download, jobs))
    (args.output / 'manifest.json').write_text(json.dumps(dict(status='EXPORTED', volume='research-lens-workspace',
        remote_root='/data/' + REMOTE, files=manifest, complete_raw_data_export=False,
        large_files='Full native VCFs, all-allele population observations, source references and patient reads remain on the persistent volume'), indent=2) + '\n')
    print('EXPORTED', len(manifest), 'files;', sum(x['bytes'] for x in manifest), 'bytes', flush=True)


if __name__ == '__main__':
    main()
