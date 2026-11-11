"""Supported local controls and the separately recorded checkpoint context."""
import json,math
from pathlib import Path
from .io import sha256

SOURCES=('MUTATION','SPLICE','CTA/SELF','FUSION')
SELECTION_KEYS={'quotas','splice_one_region_per_gene','splice_frame_policy','hydrophobicity_policy'}

def load_selection(path):
    path=Path(path);value=json.loads(path.read_text())
    if not isinstance(value,dict):raise ValueError(str(path)+': selection settings must be a JSON object')
    if set(value)!=SELECTION_KEYS:
        raise ValueError(f'{path}: expected exactly {sorted(SELECTION_KEYS)}; unsupported settings must not be silently ignored')
    validate_selection(value)
    return value

def validate_selection(value):
    from .hydrophobicity import POLICY_ID
    if value['hydrophobicity_policy'] != POLICY_ID:raise ValueError('Unsupported hydrophobicity policy; see src/yash_vaccine/hydrophobicity.py')
    quotas=value['quotas']
    if not isinstance(quotas,dict) or set(quotas)!=set(SOURCES):raise ValueError('quotas must specify MUTATION, SPLICE, CTA/SELF and FUSION; use 0 to omit a source')
    if any(type(n) is not int or n<0 for n in quotas.values()) or not sum(quotas.values()):raise ValueError('Quotas must be nonnegative integers with at least one selected region')
    if value['splice_frame_policy']!='patient_RNA_annotated_frame':raise ValueError('This helper requires patient_RNA_annotated_frame; use an independent design for another frame policy')
    for key in ['splice_one_region_per_gene']:
        if type(value[key]) is not bool:raise ValueError(key+' must be true or false')

def frozen_recipe(repo):return json.loads((Path(repo)/'config/frozen-replay.json').read_text())

def checkpoint_context(repo):
    repo=Path(repo);path=repo/'config/case.json';recipe=frozen_recipe(repo)
    if sha256(path)!=recipe['checkpoint_context_sha256']:
        raise ValueError('config/case.json describes the frozen upstream evidence. Restore it for this checkpoint run; create a named design for selection experiments. New HLA types, window lengths or expression models require rebuilt evidence and predictions.')
    return json.loads(path.read_text())
