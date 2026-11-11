"""Compare selected mutation windows with compatible GENCODE reference contexts."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference-proteins', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[4]
    selected = [r for r in csv.DictReader((repo / 'designs/initial/selection/selected.tsv').open(), delimiter='\t') if r['source'] == 'MUTATION']
    proteins = {r['protein_id']: r for r in json.loads((repo / 'data/mutation-evidence/rna/primary-proteins.json').read_text())}
    wanted = {tx for r in selected for tx in proteins[r['protein_id']]['transcript_ids']}
    references = {}
    with gzip.open(args.reference_proteins, 'rt') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            tx = record.id.split('|')[1].split('.')[0]
            if tx in wanted:
                references[tx] = str(record.seq)
    regions = {}
    for row in selected:
        sequence = row['sequence']
        matches = []
        for tx in proteins[row['protein_id']]['transcript_ids']:
            if tx not in references:
                continue
            reference = references[tx]
            for start in range(len(reference) - 24):
                differences = [i for i in range(25) if reference[start+i] != sequence[i]]
                if len(differences) <= 2:
                    matches.append(dict(transcript=tx, reference_start_aa_0based=start,
                        reference=reference[start:start+25], candidate=sequence,
                        differences=[dict(index=i, reference=reference[start+i], candidate=sequence[i]) for i in differences]))
        if not matches:
            raise ValueError('No compatible annotated context within two substitutions: ' + row['gene'])
        best = min(len(m['differences']) for m in matches)
        regions[row['event_id']] = [m for m in matches if len(m['differences']) == best]
    output = dict(source='GENCODE 37 translated proteins; annotated supporting transcripts',
        reference_url='https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.pc_translations.fa.gz',
        reference_sha256=hashlib.sha256(args.reference_proteins.read_bytes()).hexdigest(),
        selection_sha256=hashlib.sha256((repo / 'designs/initial/selection/selected.tsv').read_bytes()).hexdigest(),
        scope='Ungapped 25-aa context comparison allowing at most two substitutions; reference annotations are not normal tissue measurements. This review does not change eligibility.',
        regions=regions)
    args.output.write_text(json.dumps(output, indent=2) + '\n')
    print('Compared', len(regions), 'mutation regions with compatible annotated contexts.')


if __name__ == '__main__':
    main()
