"""Exact peptide search in explicitly scoped normal genotype/HLA models."""
from pathlib import Path
import collections
import json


def patient_self_matches(root,peptides):
    root=Path(root);directory=root/'normal-self'
    done=json.loads((directory/'complete.json').read_text());assert done['status']=='COMPLETED'
    contexts=json.loads((directory/'normal-sequence-contexts.json').read_text())
    wanted=set(peptides);lengths=sorted({len(p) for p in wanted});hits=collections.defaultdict(list)
    for row in contexts:
        sequence=row['sequence'];found=set()
        for length in lengths:
            found.update(sequence[i:i+length] for i in range(len(sequence)-length+1) if sequence[i:i+length] in wanted)
        for peptide in sorted(found):
            hits[peptide].append(dict(context_id=row['context_id'],gene=row['gene'],transcript_id=row.get('transcript_id'),
                source_interpretation=row['interpretation'],sequence_offsets=[i for i in range(len(sequence)-len(peptide)+1) if sequence[i:i+len(peptide)]==peptide]))
    return {peptide:dict(matches=hits[peptide],occurs_in_searched_models=bool(hits[peptide]),
                interpretation='POTENTIAL_PATIENT_SELF_EXCLUDE_FROM_ALTERED_TARGET_SCORING' if hits[peptide] else 'NO_HIT_IN_SEARCHED_NORMAL_CODING_AND_TYPED_HLA_MODELS',
                proves_normal_tissue_absence=False) for peptide in sorted(wanted)}
