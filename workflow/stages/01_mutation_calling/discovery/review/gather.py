"""Collect unique native-PASS alleles from segmented FilterMutectCalls."""
from pathlib import Path
import gzip
import json
import sys
import pysam
from common import sha, run, receipt

ROOT = Path(sys.argv[1])
OUT = ROOT / 'gather'
OUT.mkdir(parents=True, exist_ok=True)
source = ROOT / 'filter/segmented/normalized.vcf.gz'
completed = json.loads((source.parent / 'complete.json').read_text())
assert completed['status'] == 'COMPLETED' and sha(source) == completed['outputs'][source.name]['sha256']
records = {}
with pysam.VariantFile(source) as reader:
    header = reader.header.copy()
    for r in reader:
        if r.info['V3_EFFECTIVE_FILTER'][0] != 'PASS':
            continue
        key = r.contig, r.pos, r.ref, r.alts[0]
        if key in records:
            signature = lambda record: json.dumps(dict(samples={k: dict(v.items()) for k, v in record.samples.items()}, GERMQ=record.info.get('GERMQ')), sort_keys=True)
            assert signature(records[key]) == signature(r), ('Conflicting duplicate', key)
        else:
            records[key] = r.copy()
order = {name: index for index, name in enumerate(header.contigs)}
output = OUT / 'pass.vcf.gz'
with pysam.VariantFile(output, 'wz', header=header) as writer, gzip.open(OUT / 'model-evidence.jsonl.gz', 'wt') as metadata:
    for key in sorted(records, key=lambda k: (order[k[0]], k[1], k[2], k[3])):
        r = records[key]
        raw_id = r.id
        r.id = 'v3_' + '_'.join(map(str, key))
        writer.write(r)
        sample = next(iter(r.samples.values()))
        model = dict(effective_filter='PASS', GERMQ=r.info.get('GERMQ'), AD=sample.get('AD'),
                     DP=sample.get('DP'), AF=sample.get('AF'), raw_id=raw_id,
                     raw_alt_index=r.info['V3_ALT_INDEX'][0])
        metadata.write(json.dumps(dict(variant_id=r.id, normalized_allele=key, models={'segmented': model})) + '\n')
run(OUT, 'index-PASS', ['bcftools', 'index', '-t', '-f', output])
receipt(OUT, [output, Path(str(output) + '.tbi'), OUT / 'model-evidence.jsonl.gz'],
        counts=dict(segmented_PASS=len(records)), primary_model='segmented',
        Caris_nominations=False, worker_sha256=sha(__file__))
