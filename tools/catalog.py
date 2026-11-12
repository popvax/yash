"""Check the file catalog and regenerate its download manifest/README tables.

Run with --write to update local generated files. This never uploads, publishes,
downloads data or changes access permissions.
"""
from pathlib import Path
import argparse
import os
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from yash_vaccine import catalog


def table(rows, headers):
    return '\n'.join(['| ' + ' | '.join(headers) + ' |',
                      '|' + '|'.join('---' for _ in headers) + '|',
                      *['| ' + ' | '.join(row) + ' |' for row in rows]])


def tables(records):
    files = {r['file_id']: r for r in records}
    def link(file_id, label, *, data_page=False):
        row = files[file_id]
        target = os.path.relpath(ROOT / row['repository_path'], ROOT / 'data' if data_page else ROOT) if row['repository_path'] else row['download_url']
        return '[' + label + '](' + target + ')'
    join = lambda values: ' · '.join(values)
    root_rows = [
        ['Tumour DNA', join([link('tumour-dna-bam','BAM'),link('tumour-dna-index','index'),link('caris-variants','Caris VCF')]), catalog.size(files['tumour-dna-bam']['bytes']) + ' BAM'],
        ['Tumour RNA', join([link('tumour-rna-r1','FASTQ R1'),link('tumour-rna-r2','FASTQ R2'),link('tumour-rna-bam','BAM'),link('tumour-rna-index','index')]), catalog.size(files['tumour-rna-r1']['bytes'] + files['tumour-rna-r2']['bytes']) + ' FASTQ pair; ' + catalog.size(files['tumour-rna-bam']['bytes']) + ' BAM'],
        ['Matched-normal DNA', join([link('normal-dna-bam','BAM'),link('normal-dna-index','index')]), catalog.size(files['normal-dna-bam']['bytes']) + ' BAM'],
        ['Gene expression',link('supplied-expression','Supplied expression values'),catalog.size(files['supplied-expression']['bytes'])],
        ['HLA typing',link('hla-class-I','Six class-I alleles used in our analyses'),'Small TSV'],
        ['Clinical reports','[Stanford, Caris and Guardant PDFs](data/README.md#clinical-reports)','Four PDFs'],
        ['Processed data','[Paired variant calls, RNA evidence and candidate tables](data/README.md#processed-data-and-candidate-evidence)','Listed by dataset']]
    data_rows = []
    for label, ids, description in [
        ('Tumour DNA',[('tumour-dna-bam','BAM'),('tumour-dna-index','index')],'Caris-aligned tumour DNA reads.'),
        ('Caris variant calls',[('caris-variants','VCF')],'Original vendor calls; separate from our paired tumour–normal calls below.'),
        ('Tumour RNA reads',[('tumour-rna-r1','FASTQ R1'),('tumour-rna-r2','FASTQ R2')],'Original paired RNA reads.'),
        ('Tumour RNA alignment',[('tumour-rna-bam','BAM'),('tumour-rna-index','index')],'Supplied tumour RNA alignment.'),
        ('Matched-normal DNA',[('normal-dna-bam','BAM'),('normal-dna-index','index')],'White-blood-cell genomic DNA; sorted, duplicates marked and base qualities recalibrated.')]:
        sizes = ' + '.join(catalog.size(files[ident]['bytes']) for ident,_ in ids)
        if label=='Matched-normal DNA':sizes='**'+sizes+'**'
        data_rows.append([label,join([link(ident,name,data_page=True) for ident,name in ids]),sizes,description])
    return {'README.md': table(root_rows,['Dataset','Files','Size']),
            'data/README.md': table(data_rows,['Dataset','Files','Size','What it contains'])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true', help='Regenerate local JSON and marked README tables')
    args = parser.parse_args()
    records = catalog.load(ROOT, check_local=True)
    outputs = {'data/catalog/public-object-manifest.json': catalog.encoded(catalog.public_manifest(ROOT))}
    for name, generated in tables(records).items():
        original = (ROOT / name).read_text()
        begin, end = '<!-- datasets:start -->', '<!-- datasets:end -->'
        if original.count(begin)!=1 or original.count(end)!=1:
            raise ValueError('Expected one generated dataset table in ' + name)
        prefix, rest = original.split(begin,1)
        _, suffix = rest.split(end,1)
        outputs[name] = (prefix + begin + '\n' + generated + '\n' + end + suffix).encode()
    stale = [name for name, body in outputs.items() if (ROOT / name).read_bytes()!=body]
    if args.write:
        for name in stale:(ROOT / name).write_bytes(outputs[name])
    elif stale:
        raise ValueError('Generated files are stale: ' + ', '.join(stale) + '. Run python tools/catalog.py --write')
    print(f'CATALOG VERIFIED: {len(records)} entries; 145 published objects unchanged; ' + ('local views regenerated.' if args.write else 'generated views current.'))


if __name__ == '__main__':
    try: main()
    except (ValueError, KeyError, FileNotFoundError) as error: raise SystemExit(str(error))
