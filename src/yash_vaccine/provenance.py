"""Verify bundled historical project code and separately report upstream identities.

Third-party RNA-tool source is obtained from its providers, not redistributed.
Offline design replay does not execute those tools or pretend their absent bytes
were verified. Scientific measurements remain subject to their original hashes.
"""
from yash_vaccine.paths import resolve_input
from pathlib import Path
import hashlib
import json
import os
import tarfile

EXTERNAL_SOURCE_FILES = {
    'workflow/stages/02_lens/metadata.annotation-header-fix.py',
    'workflow/stages/02_lens/metadata.splice-review.py',
    'workflow/stages/02_lens/msnaf-core.8-11.py',
}

def digest(body):
    return hashlib.sha256(body).hexdigest()

class RecordedSource:
    def __init__(self, repo):
        self.repo = Path(repo).resolve()
        directory = self.repo / 'provenance/snapshots'
        self.record = json.loads((resolve_input(directory / 'recorded-source-2026-11-09.json')).read_text())
        archive = directory / self.record['archive']
        if digest(resolve_input(archive).read_bytes()) != self.record['archive_sha256']:
            raise ValueError('Recorded-source archive checksum mismatch')
        self.blobs = {}
        self.external = self.record.get('external_files', {})
        self.byte_verified = set()
        self.metadata_only = set()
        if set(self.external) - EXTERNAL_SOURCE_FILES:
            raise ValueError('Only named upstream RNA modules may be external source references')
        for name, meta in self.external.items():
            if any(meta[k] != self.record['files'][name][k] for k in ['bytes', 'sha256']):
                raise ValueError('External source identity differs from original record: '+name)
        with tarfile.open(resolve_input(archive), 'r:gz') as source:
            if {m.name for m in source.getmembers()} != set(self.record['files']) - set(self.external):
                raise ValueError('Recorded-source archive inventory mismatch')
            for member in source.getmembers():
                if not member.isfile() or Path(member.name).is_absolute() or '..' in Path(member.name).parts:
                    raise ValueError('Unsafe recorded-source entry')
                body = source.extractfile(member).read()
                meta = self.record['files'][member.name]
                if len(body) != meta['bytes'] or digest(body) != meta['sha256']:
                    raise ValueError('Recorded-source member checksum mismatch: ' + member.name)
                self.blobs[member.name] = body

    def external_path(self, name):
        directory = Path(os.environ.get('YASH_UPSTREAM_SOURCES_DIR', self.repo/'work/upstream-patched'))
        return directory / Path(name).name

    def read(self, name, expected, size=None, *, original_ledger=False):
        """Return exact bytes; absent external source is an explicit error."""
        name = str(name)
        snapshot_name = self.record.get('original_ledger_aliases', {}).get(name, name) if original_ledger else name
        if snapshot_name in self.blobs:
            body = self.blobs[snapshot_name]
        elif snapshot_name in self.external:
            path = self.external_path(snapshot_name)
            if not resolve_input(path).is_file():
                raise ValueError('Upstream source is not bundled: '+name+'. Obtain provider source and apply workflow/stages/02_lens/prepare_upstream.py; see docs/external-resources.md. Offline design replay does not execute this source.')
            body = resolve_input(path).read_bytes()
        else:
            path = (self.repo / name).resolve()
            if not path.is_relative_to(self.repo):
                raise ValueError('Historical input escapes repository: ' + name)
            body = resolve_input(path).read_bytes()
        if digest(body) != expected or (size is not None and len(body) != size):
            raise ValueError('Historical scientific input changed: ' + name)
        self.byte_verified.add(name)
        self.metadata_only.discard(name)
        return body

    def verify_files(self, mapping, *, original_ledger=False):
        """Verify bundled bytes; distinguish unbundled upstream hash records."""
        checked = 0
        for name, meta in mapping.items():
            snapshot_name = self.record.get('original_ledger_aliases', {}).get(name, name) if original_ledger else name
            if snapshot_name in self.external:
                recorded = self.external[snapshot_name]
                if meta['sha256'] != recorded['sha256'] or meta.get('bytes', recorded['bytes']) != recorded['bytes']:
                    raise ValueError('External source identity changed: '+name)
                if not resolve_input(self.external_path(snapshot_name)).is_file():
                    self.metadata_only.add(name)
                    continue
            self.read(name, meta['sha256'], meta.get('bytes'), original_ledger=original_ledger)
            checked += 1
        return checked

    def verification_scope(self):
        return dict(byte_verified_files=len(self.byte_verified),
                    external_source_hash_records_only=sorted(self.metadata_only),
                    all_recorded_source_bytes_verified=not self.metadata_only,
                    third_party_source_bundled=False)

    def materialize(self, destination, mapping, *, allow_external_references=False):
        """Materialize requested code; offline callers may explicitly omit RNA tools."""
        destination = Path(destination).resolve()
        destination.mkdir(parents=True, exist_ok=True)
        for name, meta in mapping.items():
            if name in self.external and allow_external_references:
                self.verify_files({name:meta})
                continue
            target = (destination / name).resolve()
            if not target.is_relative_to(destination):
                raise ValueError('Historical input escapes destination: ' + name)
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(self.read(name, meta['sha256'], meta.get('bytes')))
