"""Execute segmented DNA filtering and RNA/affinity stages on Modal; explicit execution required."""
from pathlib import Path
import argparse
import concurrent.futures
import datetime
import hashlib
import json
import shutil

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p / 'pyproject.toml').is_file())
REMOTE = '/data/results/TN25-112103/mutation-evidence-v3-20261102'
CHROMS = ['chr' + str(c) for c in list(range(1, 23)) + ['X', 'Y']]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=['prepare', 'finish', 'all', 'setup', 'population', 'segmented', 'gather', 'annotate', 'rna', 'rank', 'summarize', 'audit'])
    parser.add_argument('--execute', action='store_true')
    parser.add_argument('--workers', type=int, default=8, choices=range(1, 13))
    parser.add_argument('--output', type=Path, default=REPO / 'work/mutation-evidence-v3-20261102')
    parser.add_argument('--chrom', choices=CHROMS)
    args = parser.parse_args()
    print(json.dumps(dict(phase=args.phase, execute=args.execute, remote=REMOTE, volume='research-lens-workspace',
                          population_workers=args.workers, DNA_discovery_repeated=False,
                          primary_filter='segmented'), indent=2), flush=True)
    if not args.execute:
        return
    import modal
    volume = modal.Volume.from_name('research-lens-workspace')
    app = modal.App.lookup('yash-mutation-evidence-v3', create_if_missing=True)
    args.output.mkdir(parents=True, exist_ok=True)
    common_sha = digest(HERE / 'common.py')
    probe = args.output / 'exome-source-probe.json'
    assert probe.is_file(), 'Prepare the documented exome-source HEAD metadata first'
    with volume.batch_upload(force=True) as upload:
        upload.put_file(str(HERE / 'common.py'), REMOTE.removeprefix('/data/') + '/code/common.py')
        upload.put_file(str(probe), REMOTE.removeprefix('/data/') + '/exome-source-probe.json')

    def run(stage, chrom=None):
        if stage == 'segmented':
            worker = HERE / 'filter.py'
            name, receipt_path = 'filter-' + stage, 'filter/' + stage + '/complete.json'
            parameters = []
        elif stage == 'population':
            worker = HERE / 'population.py'
            assert chrom in CHROMS
            name, receipt_path = 'population-' + chrom, 'population/' + chrom + '/complete.json'
            parameters = [chrom]
        else:
            worker = HERE / (stage + '.py')
            name, receipt_path, parameters = stage, stage + '/complete.json', []
        assert worker.is_file(), worker
        image = 'im-hX1h6HknDMvasm9Z6q530H' if stage in ['annotate', 'rna', 'rank', 'summarize', 'audit'] else 'im-sxSffXmGfIo7J1b56AqMFq'
        memory = 32768 if stage in ['annotate', 'rna', 'rank', 'summarize', 'audit'] else 16384
        state_path = args.output / (name + '-state.json')
        worker_sha = digest(worker)
        sandbox = None
        if state_path.exists():
            old = json.loads(state_path.read_text())
            existing = modal.Sandbox.from_id(old['sandbox_id'])
            rc = existing.poll()
            if rc is None or rc == 0:
                assert old['worker_sha256'] == worker_sha and old['common_sha256'] == common_sha, 'Code changed for a live/successful stage; preserve its worker'
                sandbox = existing
            else:
                archive = args.output / 'attempts' / name / old['sandbox_id']
                archive.mkdir(parents=True, exist_ok=True)
                for previous in [state_path, args.output / (name + '.stdout.log'), args.output / (name + '.stderr.log')]:
                    if previous.exists():
                        shutil.copyfile(previous, archive / previous.name)
        if sandbox is None:
            with volume.batch_upload(force=True) as upload:
                upload.put_file(str(worker), REMOTE.removeprefix('/data/') + '/code/' + worker.name)
            sandbox = modal.Sandbox.create('python', '-u', REMOTE + '/code/' + worker.name, REMOTE, *parameters,
                app=app, image=modal.Image.from_id(image), cpu=4, memory=memory, timeout=21600,
                volumes={'/data': volume})
            state_path.write_text(json.dumps(dict(sandbox_id=sandbox.object_id, stage=stage, chromosome=chrom,
                worker_sha256=worker_sha, common_sha256=common_sha, image=image, memory_mib=memory,
                started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()), indent=2) + '\n')
        print(name, sandbox.object_id, flush=True)
        stdout, stderr = sandbox.stdout.read(), sandbox.stderr.read()
        sandbox.wait()
        (args.output / (name + '.stdout.log')).write_text(stdout)
        (args.output / (name + '.stderr.log')).write_text(stderr)
        if sandbox.poll() != 0:
            raise RuntimeError((name, sandbox.poll(), stdout[-3000:], stderr[-3000:]))
        data = b''.join(volume.read_file(REMOTE.removeprefix('/data/') + '/' + receipt_path))
        report = json.loads(data)
        assert report['status'] == 'COMPLETED'
        if 'worker_sha256' in report:
            assert report['worker_sha256'] == worker_sha
        (args.output / (name + '.complete.json')).write_bytes(data)
        brief = {k: v for k, v in report.items() if k in ['all_emitted_alleles', 'normalized_alleles', 'counts', 'tiers', 'input_variants', 'usable_primary_proteins', 'primary_ranked_regions']}
        if 'counts' in brief:
            brief['counts'] = {k: v for k, v in brief['counts'].items() if not k.startswith('native_') or k == 'native_PASS'}
        print('COMPLETED', name, json.dumps(brief), flush=True)
        return report

    if args.phase in ['setup', 'prepare', 'all']:
        run('setup')
    if args.phase in ['population', 'prepare', 'all']:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
            list(pool.map(lambda c: run('population', c), [args.chrom] if args.chrom else CHROMS))
    if args.phase in ['finish', 'all']:
        # Native filtering covers the entire emitted set. Completed, verified
        # chromosome receipts are prerequisites; no prior shortlist is read.
        for c in CHROMS:
            report = json.loads(b''.join(volume.read_file(REMOTE.removeprefix('/data/') + '/population/' + c + '/complete.json')))
            assert report['status'] == 'COMPLETED'
        run('segmented')
        for stage in ['gather', 'annotate', 'rna', 'rank', 'summarize', 'audit']:
            run(stage)
    elif args.phase not in ['setup', 'population', 'prepare']:
        run(args.phase)


if __name__ == '__main__':
    main()
