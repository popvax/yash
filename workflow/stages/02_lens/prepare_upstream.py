"""Apply minimal patches to provider-obtained RNA source; no network or inference.

Default: print a plan. --execute requires originals extracted from the named
upstream images and writes a new private/local output directory. Both original
and patched bytes must match the recorded execution hashes.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re

HERE = Path(__file__).resolve().parent


def sha(body):
    return hashlib.sha256(body).hexdigest()


def apply_patch(base, patch):
    original = base.decode().splitlines(keepends=True)
    lines = patch.decode().splitlines(keepends=True)
    result = []; cursor = 0; i = 2
    while i < len(lines):
        match = re.match(r'@@ -(\d+)(?:,(\d+))? \+(\d+)(?:,(\d+))? @@', lines[i])
        if not match:
            raise ValueError('Unsupported patch header')
        start = int(match[1]); old_count = int(match[2] or 1); new_count = int(match[4] or 1)
        position = start if old_count == 0 else start - 1
        before = []; after = []; i += 1
        while i < len(lines) and not lines[i].startswith('@@ '):
            line = lines[i]
            if line[0] in '- ':before.append(line[1:])
            if line[0] in '+ ':after.append(line[1:])
            if line[0] not in '-+ ':raise ValueError('Unsupported patch line')
            i += 1
        if len(before) != old_count or len(after) != new_count:
            raise ValueError('Patch hunk size mismatch')
        if position < cursor or original[position:position+old_count] != before:
            raise ValueError('Patch does not match exact upstream source')
        result.extend(original[cursor:position]); result.extend(after)
        cursor = position + old_count
    result.extend(original[cursor:])
    return ''.join(result).encode()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sources-dir', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    sources = json.loads((HERE/'upstream-sources.json').read_text())['sources']
    if not args.execute:
        print(json.dumps(dict(execute=False, sources=sources,
                              instruction='Obtain original files from their providers under their terms, then supply --sources-dir, --output and --execute.'), indent=2))
        return
    if args.sources_dir is None or args.output is None:
        parser.error('--execute requires --sources-dir and --output')
    if args.output.exists():
        parser.error('Choose a new output directory; source files are never overwritten')
    prepared = {}
    for record in sources:
        base = (args.sources_dir/record['base_filename']).read_bytes()
        patch = (HERE/record['patch']).read_bytes()
        if len(base) != record['base_bytes'] or sha(base) != record['base_sha256']:
            raise ValueError('Upstream source does not match recorded version: '+record['base_filename'])
        if sha(patch) != record['patch_sha256']:
            raise ValueError('Patch checksum changed: '+record['patch'])
        body = apply_patch(base, patch)
        if len(body) != record['bytes'] or sha(body) != record['sha256']:
            raise ValueError('Patched output does not match executed bytes: '+record['path'])
        prepared[Path(record['path']).name] = body
    args.output.mkdir(parents=True)
    for name, body in prepared.items():(args.output/name).write_bytes(body)
    receipt = dict(status='VERIFIED', source='provider_obtained_originals',
                   target_sha256={name:sha(body) for name,body in prepared.items()},
                   model_inference_rerun=False, rights_granted=False)
    (args.output/'upstream-preparation.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
