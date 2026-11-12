"""Resolve recorded input identities to their current repository locations.

Original recipes and execution receipts retain their paths and hashes. Only
reads are mapped; outputs and remote /data locations are never redirected.
"""
from functools import lru_cache
from pathlib import Path
import json


@lru_cache(maxsize=512)
def _repository(parent):
    parent = Path(parent)
    for candidate in (parent, *parent.parents):
        if (candidate / 'provenance/path-map.json').is_file() and (candidate / 'src/yash_vaccine').is_dir():
            return candidate
    return None


@lru_cache(maxsize=8)
def _mapping(repo):
    record = json.loads((repo / 'provenance/path-map.json').read_text())
    if record['schema_version'] != 1:
        raise ValueError('Unsupported repository path-map version')
    return record


def resolve_input(path, repo=None):
    """Map an original local input path; leave unrelated paths unchanged."""
    original = Path(path).resolve()
    root = Path(repo).resolve() if repo is not None else _repository(original.parent)
    if root is None or not original.is_relative_to(root):
        return original
    name = original.relative_to(root).as_posix()
    record = _mapping(root)
    target = record['files'].get(name)
    if target is None:
        for before, after in sorted(record['prefixes'].items(), key=lambda item: len(item[0]), reverse=True):
            if name == before or name.startswith(before + '/'):
                target = after + name[len(before):]
                break
    if target is None:
        return original
    resolved = (root / target).resolve()
    if not resolved.is_relative_to(root):
        raise ValueError('Mapped input escapes repository: ' + name)
    return resolved
