"""Inspect official exome source metadata; --execute performs read-only HEAD requests."""
from pathlib import Path
import argparse
import concurrent.futures
import json
import urllib.request

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p/'pyproject.toml').is_file())
BASE = 'https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/gnomad/v4.1/exomes/'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--execute', action='store_true')
    parser.add_argument('--output', type=Path, default=REPO / 'work/mutation-evidence-v3-20261102/exome-source-probe.json')
    args = parser.parse_args()
    print(json.dumps(dict(execute=args.execute, source=BASE, chromosomes=24, output=str(args.output)), indent=2))
    if not args.execute:
        return
    def probe(chrom):
        url = BASE + 'gnomad.exomes.v4.1.sites.' + chrom + '_trimmed.vcf.bgz'
        with urllib.request.urlopen(urllib.request.Request(url, method='HEAD'), timeout=45) as response:
            return dict(chrom=chrom, bytes=int(response.headers['Content-Length']), etag=response.headers.get('ETag'),
                        last_modified=response.headers.get('Last-Modified'), url=url)
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        results = list(pool.map(probe, ['chr' + str(c) for c in list(range(1, 23)) + ['X', 'Y']]))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results, indent=2) + '\n')
    print('SOURCE_BYTES', sum(r['bytes'] for r in results))


if __name__ == '__main__':
    main()
