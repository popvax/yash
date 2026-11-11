"""Plan or explicitly run the recorded initial RNA encoder on Modal."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
INPUTS = [
    'run.py', 'optimizer.py', 'genbank.py', 'backbone.json', 'requirements.txt',
    'annotation-sources.json', 'human-codon-usage.json',
    '../../selection/selected.tsv', '../architecture.json',
    '../protein/cassette.fasta', '../protein/components.tsv', '../final-check.json',
]
OUTPUTS = [
    'cds.fasta', 'mrna-segments.fasta', 'template-segments.fasta', 'template-segments.gb',
    'components.tsv', 'construct.json', 'verification.json', 'rna-layout.txt', 'optimizer.log',
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--execute', action='store_true')
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--volume', default='research-lens-workspace')
    parser.add_argument('--record-new-encoding', action='store_true',
        help='Explicitly generate and pin a new coding-sequence hash before the checked export')
    args = parser.parse_args()
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    remote = '/data/results/TN25-112103/poprna-final-' + stamp
    inputs = [(HERE / rel).resolve() for rel in INPUTS]
    state = dict(volume=args.volume, remote=remote, cpu=2, memory_mib=8192, gpu=None,
                 execute=args.execute, timeout_seconds=2400 if args.record_new_encoding else 1200,
                 record_new_encoding=args.record_new_encoding,
                 input_sha256={str(path.relative_to(REPO)): hashlib.sha256(path.read_bytes()).hexdigest()
                               for path in inputs})
    print(json.dumps(state, indent=2), flush=True)
    if not args.execute:
        return
    if args.output.exists():
        parser.error('Choose a new output directory')
    args.output.mkdir(parents=True)
    path = args.output / 'modal-state.json'
    def save():
        path.write_text(json.dumps(state, indent=2) + '\n')
    save()
    import modal
    image = modal.Image.debian_slim(python_version='3.11').pip_install(
        'ViennaRNA==2.7.2', 'python-codon-tables==0.1.18', 'biopython==1.86', 'numpy==1.26.4')
    volume = modal.Volume.from_name(args.volume)
    sandbox = None
    try:
        with modal.enable_output():
            sandbox = modal.Sandbox.create(
                app=modal.App.lookup('research-yash-vaccine', create_if_missing=True),
                image=image, cpu=2, memory=8192, timeout=state['timeout_seconds'], volumes={'/data': volume})
        state.update(sandbox_id=sandbox.object_id, image_id=image.object_id,
                     started_utc=datetime.now(timezone.utc).isoformat())
        save()
        for source in inputs:
            sandbox.filesystem.write_bytes(source.read_bytes(), remote + '/' + str(source.relative_to(REPO)))
        if args.record_new_encoding:
            cfg_path = HERE / 'backbone.json'
            cfg = json.loads(cfg_path.read_text())
            settings = cfg['optimization']
            protein = ''.join(line.strip() for line in (HERE / '../protein/cassette.fasta').read_text().splitlines()
                              if not line.startswith('>'))
            coding_path = remote + '/prepared-coding.txt'
            preparation = sandbox.exec('python', remote + '/designs/initial/construct/mrna/optimizer.py',
                '--protein', '-i', protein, '--output', coding_path,
                '--organism', settings['species_table'], '--cutoff', str(settings['codon_cutoff']),
                '--seed', str(settings['seed']), '--fold-iterations', str(settings['structure_iterations_per_pass']),
                '--sequences', *settings['avoid_sequences_rna'], timeout=500)
            (args.output / 'preparation-stdout.log').write_text(preparation.stdout.read())
            (args.output / 'preparation-stderr.log').write_text(preparation.stderr.read())
            if preparation.wait():
                raise RuntimeError('New PopRNA preparation failed; inspect preparation logs')
            coding = sandbox.filesystem.read_bytes(coding_path).decode().strip()
            from Bio.Seq import Seq
            if set(coding) - set('ACGT') or str(Seq(coding).translate()) != protein:
                raise RuntimeError('Prepared coding sequence does not encode the selected protein')
            cfg['optimization']['recorded_coding_sequence_sha256'] = hashlib.sha256(coding.encode()).hexdigest()
            cfg_path.write_text(json.dumps(cfg, indent=2) + '\n')
            state['input_sha256'][str(cfg_path.relative_to(REPO))] = hashlib.sha256(cfg_path.read_bytes()).hexdigest()
            state['prepared_coding_sequence_sha256'] = cfg['optimization']['recorded_coding_sequence_sha256']
            sandbox.filesystem.write_bytes(cfg_path.read_bytes(), remote + '/' + str(cfg_path.relative_to(REPO)))
            save()
        process = sandbox.exec('python', remote + '/designs/initial/construct/mrna/run.py',
                               '--output', remote + '/output', timeout=1100)
        stdout, stderr = process.stdout.read(), process.stderr.read()
        (args.output / 'stdout.log').write_text(stdout)
        (args.output / 'stderr.log').write_text(stderr)
        state['exit_code'] = process.wait()
        if state['exit_code']:
            raise RuntimeError('Modal encoder failed: ' + stderr)
        for name in OUTPUTS:
            (args.output / name).write_bytes(sandbox.filesystem.read_bytes(remote + '/output/' + name))
        state['completed'] = True
    finally:
        if sandbox is not None:
            sandbox.terminate(wait=True)
            state['terminated'] = True
        state['finished_utc'] = datetime.now(timezone.utc).isoformat()
        save()
    receipt = b''.join(volume.read_file(remote.removeprefix('/data/') + '/output/verification.json'))
    state['volume_persistence_verified'] = receipt == (args.output / 'verification.json').read_bytes()
    state['output_sha256'] = {name: hashlib.sha256((args.output / name).read_bytes()).hexdigest()
                              for name in OUTPUTS}
    save()
    if not state['volume_persistence_verified']:
        raise RuntimeError('Durable Modal verification receipt differs')
    print('Modal RNA encoding and persistent-volume verification passed.', flush=True)


if __name__ == '__main__':
    main()
