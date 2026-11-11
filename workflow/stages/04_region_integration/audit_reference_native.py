"""Native Vaxrank reference parity anchor for all admitted event seed peptides."""
import dataclasses,json,shutil,sys
from pathlib import Path
from pyensembl import Genome
from vaxrank.reference_proteome import self_reference_matches
from vaxrank.vaccine_antigen import VaccineAntigen,TargetableMask,AminoAcidInterval,TumorSpecificityAttestation
R=Path(sys.argv[1]);S=Path('/data/results/TN25-112103')
contexts=[json.loads(l) for l in (R/'seed-sequence-contexts.jsonl').read_text().splitlines()]
queries=sorted({t['peptide'] for c in contexts if c['source']!='CTA/SELF' for t in c['targets']})
cfg=json.loads((S/'vaxrank-wildtype-selfcheck-20261030/genome-config.json').read_text())
scratch=Path('/work/combined-regions/reference-cache')
if not scratch.exists():shutil.copytree(cfg['cache_directory_path'],scratch)
cfg['cache_directory_path']=str(scratch);genome=Genome(**cfg);genome.index()
antigen=VaccineAntigen(kind='splice',amino_acids='A',targetable_mask=TargetableMask((AminoAcidInterval(0,1),)),
    tumor_specificity=TumorSpecificityAttestation(status='overridden',evidence_kind='RNA_only_research_region_catalogue',
    evidence_source='Frozen LENS RNA review',requires_review=True,override_reason='Reference comparison only'))
matches=self_reference_matches(queries,antigen,genome)
assert all(m.source_provenance_complete for m in matches.values())
(R/'event-peptide-self-reference.json').write_text(json.dumps({p:dataclasses.asdict(m) for p,m in matches.items()},indent=2)+'\n')
print('Native complete reference comparison:',len(matches),'peptides;',sum(m.occurs for m in matches.values()),'matches',flush=True)
