"""Small shared utilities for the segmented discovery workflow; import has no side effects."""
from pathlib import Path
import datetime
import hashlib
import json
import subprocess

DISCOVERY = Path('/data/results/TN25-112103/mutation-evidence-v2-20261102')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(value, indent=2, default=str) + '\n')
    tmp.replace(path)


def run(directory, name, args):
    args = list(map(str, args))
    save(directory / 'status.json', dict(phase=name, updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    print(name, flush=True)
    with (directory / (name + '.log')).open('w') as log:
        subprocess.run(args, stdout=log, stderr=log, check=True)
    return args


def receipt(directory, files, **kwargs):
    result = dict(status='COMPLETED', **kwargs,
                  completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  outputs={str(p.relative_to(directory)): dict(bytes=p.stat().st_size, sha256=sha(p)) for p in files})
    save(directory / 'complete.json', result)
    print(json.dumps({k: v for k, v in result.items() if k != 'outputs'}), flush=True)
    return result


def scalar(v, index=0):
    if isinstance(v, (list, tuple)):
        return v[index] if len(v) > index else None
    return v
