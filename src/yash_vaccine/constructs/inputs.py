"""Shared selection/linker validation and explicit prediction requests."""
import json
from pathlib import Path
from ..io import rows,sha256,write_json
from ..settings import checkpoint_context
from .boundaries import all_queries,LENGTHS

AMINO_ACIDS=set('ACDEFGHIKLMNPQRSTVWY')

def inputs(repo,regions_path=None,linkers_path=None):
    repo=Path(repo);context=checkpoint_context(repo)
    region_file=Path(regions_path or repo/'designs/initial/selection/selected.tsv')
    linker_file=Path(linkers_path or repo/'config/linkers.json')
    regions=list(rows(region_file));linker_config=json.loads(linker_file.read_text())
    if not isinstance(linker_config,dict) or not isinstance(linker_config.get('choices'),list):raise ValueError('Linker file must be a JSON object with a choices list')
    linkers=linker_config['choices']
    if len(regions)<2:raise ValueError('Joint construct design requires at least two regions')
    if any(not x.get('antigen_id') or len(x.get('sequence',''))!=context['region_length_aa'] or set(x['sequence'])-AMINO_ACIDS for x in regions):
        raise ValueError('Regions must have antigen_id and a canonical 25-aa sequence; new window lengths require rebuilt upstream evidence')
    if len({x['antigen_id'] for x in regions})!=len(regions) or len({x['sequence'] for x in regions})!=len(regions):raise ValueError('Antigen IDs and region sequences must each be unique')
    if not linkers or any(not isinstance(x,dict) or not x.get('name') or not isinstance(x.get('sequence'),str) or set(x['sequence'])-AMINO_ACIDS for x in linkers):raise ValueError('Linker choices need a name and a canonical amino-acid sequence; empty sequence means a direct join')
    if len({x['name'] for x in linkers})!=len(linkers):raise ValueError('Linker names must be unique')
    return regions,linkers,context,region_file,linker_file

def request(repo,regions_path=None,linkers_path=None):
    regions,linkers,context,rf,lf=inputs(repo,regions_path,linkers_path)
    peptides=all_queries(regions,linkers)
    metadata=dict(regions=str(rf.resolve()),regions_sha256=sha256(rf),region_count=len(regions),
        linkers=str(lf.resolve()),linkers_sha256=sha256(lf),linker_choices=len(linkers),
        alleles=context['hla_class_I'],peptide_lengths=list(LENGTHS),unique_peptides=len(peptides),
        peptide_hla_pairs=len(peptides)*len(context['hla_class_I']),model_inference_rerun=False)
    return peptides,metadata

def prepare(repo,out,regions_path,linkers_path=None):
    out=Path(out)
    if out.exists():raise ValueError('Choose a new prediction-request directory: '+str(out))
    peptides,metadata=request(repo,regions_path,linkers_path);out.mkdir(parents=True)
    (out/'peptides.txt').write_text(''.join(p+'\n' for p in peptides))
    write_json(out/'alleles.json',metadata['alleles']);write_json(out/'prediction-request.json',metadata)
    return metadata
