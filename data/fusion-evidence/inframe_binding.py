"""Complete the in-frame fusion junction comparison without changing selection."""
import os
os.environ['OMP_NUM_THREADS']='2'
os.environ['MKL_NUM_THREADS']='2'
from pathlib import Path
import hashlib
import importlib.metadata
import json
import sys
import torch
from mhctools import MHCflurry

root=Path(sys.argv[1]);request=json.loads((root/'inframe-request.json').read_text())
models=Path('/data/results/TN25-112103/openvax-mutations-20261018/references/mhcflurry/2.3.0/models_class1_presentation/models')
manifest=json.loads((root/'model-files.json').read_text())
for r in manifest:
    p=models/r['path']
    assert p.stat().st_size==r['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'],str(p)
torch.set_num_threads(2)
predictor=MHCflurry(alleles=request['alleles'],default_peptide_lengths=[8,9,10,11],models_path=str(models))
frame=predictor.predict_peptides(request['peptides']).to_dataframe()
assert len(frame)==len(request['peptides'])*len(request['alleles'])
rows=[]
for x in frame.to_dict('records'):
    rows.append(dict(peptide=x['peptide'],allele=x['allele'],affinity_nm=float(x['affinity']),
                     percentile_rank=float(x['percentile_rank']),exact_reference_match=x['peptide'] in request['reference_matching_peptides']))
nonself=[r for r in rows if not r['exact_reference_match']]
out=dict(status='COMPLETED',worker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         request_sha256=hashlib.sha256((root/'inframe-request.json').read_bytes()).hexdigest(),
         model_files_verified=True,mhcflurry_version=importlib.metadata.version('mhcflurry'),
         mhctools_version=importlib.metadata.version('mhctools'),predictions=rows,
         best_nonself=min(nonself,key=lambda r:r['affinity_nm']),nonself_rows_le500=sum(r['affinity_nm']<=500 for r in nonself),
         changes_selection=False)
(root/'inframe-binding-review.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='predictions'}),flush=True)
