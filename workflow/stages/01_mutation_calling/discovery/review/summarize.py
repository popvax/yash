"""Record the segmented RNA-supported mutation ranking."""
from pathlib import Path
import csv
import json
import sys
from common import sha, save, receipt

ROOT = Path(sys.argv[1])
OUT = ROOT / 'summarize'
OUT.mkdir(parents=True, exist_ok=True)
ranked = json.loads((ROOT / 'rank/primary-ranked-regions.json').read_text())
for rank, row in enumerate(ranked, 1):
    row.update(model_rank=rank, filter_model='segmented')
save(OUT / 'segmented-primary-ranked-regions.json', ranked)
fields = list(dict.fromkeys(k for r in ranked for k in r)) or ['variant_id']
with (OUT / 'segmented-primary-ranked-regions.tsv').open('w') as handle:
    writer = csv.DictWriter(handle, delimiter='\t', fieldnames=fields, lineterminator='\n')
    writer.writeheader()
    writer.writerows({k: json.dumps(v) if isinstance(v, (dict, list)) else v for k, v in r.items()} for r in ranked)
(OUT / 'segmented-regions.fasta').write_text(''.join('>' + r['variant_id'] + ' ' + r['gene'] + '\n' + r['sequence'] + '\n' for r in ranked))
summary = dict(primary_model='segmented', ranked_regions_by_model={'segmented': len(ranked)},
               ranked_genes=len({r['gene'] for r in ranked}), somatic_status='UNCONFIRMED_NO_MATCHED_NORMAL')
save(OUT / 'summary.json', summary)
receipt(OUT, [OUT / name for name in ['segmented-primary-ranked-regions.json', 'segmented-primary-ranked-regions.tsv',
        'segmented-regions.fasta', 'summary.json']], ranked_regions_by_model=summary['ranked_regions_by_model'],
        worker_sha256=sha(__file__))
