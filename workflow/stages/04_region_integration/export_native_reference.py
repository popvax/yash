"""Freeze the exact native Vaxrank annotated reference inventory for fast scans."""
import dataclasses,json,sys
from pathlib import Path
from pyensembl import Genome
from vaxrank.reference_proteome import _protein_source_snapshot
R=Path(sys.argv[1]);S=Path('/data/results/TN25-112103')
cfg=json.loads((S/'vaxrank-wildtype-selfcheck-20261030/genome-config.json').read_text())
cfg['cache_directory_path']='/work/combined-regions/reference-cache'
g=Genome(**cfg);g.index();snapshot=_protein_source_snapshot(g)
with (R/'native-reference-snapshot.jsonl').open('w') as f:
    for sequence,sources in snapshot:f.write(json.dumps(dict(sequence=sequence,sources=[dataclasses.asdict(s) for s in sources]))+'\n')
print('Native reference inventory:',len(snapshot),'sequences;',sum(len(s) for _,s in snapshot),'sources',flush=True)
