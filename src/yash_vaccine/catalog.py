"""The reader-facing file catalog and its locked public download identities."""
from pathlib import Path
from urllib.parse import quote
import csv
import hashlib
import json
import re

PRIMARY_DATASETS = ('tumour-dna', 'tumour-rna', 'normal-dna', 'expression', 'hla', 'clinical-reports')


def load(repo, *, check_local=False):
    repo = Path(repo).resolve()
    with (repo / 'data/files.tsv').open() as stream:
        records = list(csv.DictReader(stream, delimiter='\t'))
    ids, keys = set(), set()
    for row in records:
        if row['file_id'] in ids:
            raise ValueError('Duplicate catalog file ID: ' + row['file_id'])
        ids.add(row['file_id'])
        row['bytes'] = int(row['bytes'])
        if row['bytes'] < 0 or not re.fullmatch(r'[0-9a-f]{64}', row['sha256']):
            raise ValueError('Invalid size/checksum: ' + row['file_id'])
        if row['storage_key']:
            if row['storage_key'] in keys:
                raise ValueError('Duplicate download key: ' + row['storage_key'])
            keys.add(row['storage_key'])
        if row['repository_path']:
            path = (repo / row['repository_path']).resolve()
            if not path.is_relative_to(repo) or not path.is_file():
                raise ValueError('Missing/unsafe catalog path: ' + row['repository_path'])
            if check_local and (path.stat().st_size != row['bytes'] or hashlib.sha256(path.read_bytes()).hexdigest() != row['sha256']):
                raise ValueError('Catalog file changed: ' + row['repository_path'])
    return records


def public_manifest(repo):
    """Compile download entries from the catalog without changing the release."""
    repo = Path(repo).resolve()
    release = json.loads((repo / 'provenance/data-catalog/published-release.json').read_text())
    native_bytes = (repo / 'provenance/data-catalog/original-object-manifest.json').read_bytes()
    if hashlib.sha256(native_bytes).hexdigest() != release['manifest_metadata']['native_manifest_sha256']:
        raise ValueError('Original download identities changed')
    rows = {r['storage_key']: r for r in load(repo) if r['storage_key']}
    key_hash = hashlib.sha256(('\n'.join(sorted(rows)) + '\n').encode()).hexdigest()
    if len(rows) != release['public_objects'] or key_hash != release['public_keys_sha256']:
        raise ValueError('Public download set changed; review the release before publishing')
    access = json.loads((repo / 'data/catalog/access.json').read_text())
    objects = []
    for original in json.loads(native_bytes)['objects']:
        row = rows.get(original['key'])
        if row is None:
            continue
        if row['bytes'] != original['bytes'] or row['sha256'] != original['sha256']:
            raise ValueError('Published data identity changed: ' + original['key'])
        if original['category'] not in ('raw', 'analysis'):
            raise ValueError('Catalog must not mirror general references or software')
        expected_url = access['public_base_url'].rstrip('/') + '/' + quote(original['key'], safe='/')
        if row['download_url'] != expected_url:
            raise ValueError('Download URL differs from the published endpoint: ' + original['key'])
        objects.append(original)
    manifest = dict(release['manifest_metadata'], objects=objects)
    if hashlib.sha256(encoded(manifest)).hexdigest() != release['manifest_sha256']:
        raise ValueError('Generated download manifest differs from the published release')
    return manifest


def encoded(value):
    return (json.dumps(value, indent=2) + '\n').encode()


def select(records, *, datasets=None, files=None):
    if datasets:
        missing = set(datasets) - {r['dataset'] for r in records}
        if missing:
            raise ValueError('Unknown dataset: ' + ', '.join(sorted(missing)))
        records = [r for r in records if r['dataset'] in datasets]
    if files:
        missing = set(files) - {r['file_id'] for r in records}
        if missing:
            raise ValueError('Unknown file ID in the selected datasets: ' + ', '.join(sorted(missing)))
        records = [r for r in records if r['file_id'] in files]
    return records


def download_keys(repo, *, datasets=None, files=None):
    records = select(load(repo), datasets=datasets, files=files)
    return {r['storage_key'] for r in records if r['storage_key']}


def size(value):
    for unit, scale in [('GB', 10**9), ('MB', 10**6), ('kB', 10**3)]:
        if value >= scale:
            return f'{value / scale:.2f} {unit}'
    return f'{value} B'
