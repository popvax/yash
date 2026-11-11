"""One native MHCflurry affinity prediction per unique peptide and patient allele.
No RNA abundance enters the model. Raw predictions are retained before any
binding threshold or event-expression weighting. Run in the pinned CUDA image.
"""
import os
os.environ['OMP_NUM_THREADS']='2';os.environ['MKL_NUM_THREADS']='2'
import csv,json,time,hashlib,importlib.metadata
from pathlib import Path
import numpy as np
import torch
from mhcflurry import Class1PresentationPredictor
ROOT=Path(os.environ.get('LENS_RUN_ROOT','/data/results/TN25-112103/lens-expression-aware-20261030'))
MODELS=Path('/data/results/TN25-112103/openvax-mutations-20261018/references/mhcflurry/2.3.0/models_class1_presentation/models')
ALLELES=['HLA-A*01:01','HLA-A*33:03','HLA-B*08:01','HLA-B*15:08','HLA-C*01:02','HLA-C*07:02']
torch.set_num_threads(2);torch.set_float32_matmul_precision('highest')
torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False
assert torch.cuda.is_available(), 'GPU expected'
assert importlib.metadata.version('mhcflurry')=='2.3.14'
peptides=(ROOT/'prediction/unique-peptides.txt').read_text().splitlines()
assert len(peptides)==len(set(peptides)) and all(8<=len(p)<=11 for p in peptides)
m=Class1PresentationPredictor.load(str(MODELS)).affinity_predictor
start=time.time();count=0
with (ROOT/'prediction/affinity-all-unique.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['peptide','allele','mhcflurry_2.3.14.aff','mhcflurry_2.3.14.aff_perc'])
 for i in range(0,len(peptides),131072):
  batch=peptides[i:i+131072]
  for allele in ALLELES:
   affinities=m.predict(peptides=batch,allele=allele,model_kwargs={'batch_size':16384})
   percentiles=m.percentile_ranks(affinities,allele=allele)
   assert np.isfinite(affinities).all() and (affinities>0).all()
   w.writerows(zip(batch,[allele]*len(batch),affinities,percentiles));count+=len(batch)
  f.flush();print(json.dumps({'completed_peptides':i+len(batch),'prediction_rows':count}),flush=True)
summary={'mhcflurry':importlib.metadata.version('mhcflurry'),'torch':torch.__version__,'model_bundle':'2.3.0','models_path':str(MODELS),'alleles':ALLELES,'lengths':[8,9,10,11],'unique_peptides':len(peptides),'prediction_rows':count,'seconds':time.time()-start,'affinity_output_sha256':hashlib.sha256((ROOT/'prediction/affinity-all-unique.tsv').read_bytes()).hexdigest(),'input_sha256':hashlib.sha256((ROOT/'prediction/unique-peptides.txt').read_bytes()).hexdigest(),'model_files_sha256':{str(p.relative_to(MODELS)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(MODELS.rglob('*')) if p.is_file()},'note':'Affinity-only; expression weighting is separate. No presentation or immunogenicity claim.'}
(ROOT/'prediction/model-run.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k!='model_files_sha256'}),flush=True)
