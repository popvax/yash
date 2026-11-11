"""Plan or explicitly execute the supplementary selected-SNV check on Modal."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--execute', action='store_true')
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    remote = '/data/results/TN25-112103/selected-mutation-review-' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    state = dict(execute=args.execute, image='im-hX1h6HknDMvasm9Z6q530H',
                 volume='research-lens-workspace', remote=remote, cpu=2, memory_mib=4096,
                 terminated=False, changes_eligibility=False)
    print(json.dumps(state, indent=2), flush=True)
    if not args.execute:
        return
    if args.output.exists():
        parser.error('Choose a new output directory')
    args.output.mkdir(parents=True)
    def save():
        (args.output / 'modal-state.json').write_text(json.dumps(state, indent=2) + '\n')
    import modal
    volume = modal.Volume.from_name(state['volume'])
    sandbox = None
    try:
        sandbox = modal.Sandbox.create(app=modal.App.lookup('research-yash-vaccine', create_if_missing=True),
            image=modal.Image.from_id(state['image']), cpu=2, memory=4096, timeout=600,
            volumes={'/data': volume})
        state['sandbox_id'] = sandbox.object_id
        save()
        for name in ('check.py', 'request.json'):
            sandbox.filesystem.write_bytes((here / name).read_bytes(), remote + '/' + name)
        process = sandbox.exec('python', '-u', remote + '/check.py', '--request', remote + '/request.json',
                               '--output', remote + '/result.json', timeout=500)
        for line in process.stdout:
            print(line, end='', flush=True)
        (args.output / 'stderr.log').write_text(process.stderr.read())
        state['exit_code'] = process.wait()
        if state['exit_code']:
            raise RuntimeError('Read check failed; inspect stderr.log')
        (args.output / 'result.json').write_bytes(sandbox.filesystem.read_bytes(remote + '/result.json'))
        result = json.loads((args.output / 'result.json').read_text())
        assert result['request_sha256'] == hashlib.sha256((here / 'request.json').read_bytes()).hexdigest()
        assert result['worker_sha256'] == hashlib.sha256((here / 'check.py').read_bytes()).hexdigest()
        state['completed'] = True
    finally:
        if sandbox:
            sandbox.terminate(wait=True)
            state['terminated'] = True
        state['finished_utc'] = datetime.now(timezone.utc).isoformat()
        save()
    stored = b''.join(volume.read_file(remote.removeprefix('/data/') + '/result.json'))
    state['volume_persistence_verified'] = stored == (args.output / 'result.json').read_bytes()
    save()
    assert state['volume_persistence_verified']


if __name__ == '__main__':
    main()
