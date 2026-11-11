"""Check model identity and CPU/cache agreement for the three added queries."""
import csv,hashlib,importlib.metadata,json,math,sys
from pathlib import Path
import torch
from mhcflurry import Class1PresentationPredictor
R=Path(sys.argv[1]);L=Path('/data/results/TN25-112103/lens-expression-aware-20261030')
parent=json.loads((L/'prediction/model-run.json').read_text());models=Path(parent['models_path'])
assert importlib.metadata.version('mhcflurry')=='2.3.14'
for relative,expected in parent['model_files_sha256'].items():assert hashlib.sha256((models/relative).read_bytes()).hexdigest()==expected,relative
torch.set_num_threads(4);torch.set_float32_matmul_precision('highest')
model=Class1PresentationPredictor.load(str(models)).affinity_predictor
samples=[]
with (R/'scoring-affinities.tsv').open() as f:
    for row in csv.DictReader(f,delimiter='\t'):
        if row['peptide'] not in samples:samples.append(row['peptide'])
        if len(samples)==20:break
expected={}
with (R/'scoring-affinities.tsv').open() as f:
    for row in csv.DictReader(f,delimiter='\t'):
        if row['peptide'] in samples:expected[(row['peptide'],row['allele'])]=float(row['mhcflurry_2.3.14.aff'])
max_relative=0
for allele in parent['alleles']:
    observed=model.predict(peptides=samples,allele=allele,model_kwargs={'batch_size':4096})
    for p,value in zip(samples,observed):
        ref=expected[(p,allele)];delta=abs(float(value)-ref)/ref;max_relative=max(max_relative,delta)
        assert math.isclose(float(value),ref,rel_tol=1e-4,abs_tol=1e-3),(p,allele,ref,value)
result=dict(status='VERIFIED',model_files_checked=len(parent['model_files_sha256']),
            cache_CPU_pairs_checked=len(expected),max_relative_affinity_difference=max_relative,
            relative_tolerance=1e-4,absolute_tolerance_nm=1e-3,
            mhcflurry=importlib.metadata.version('mhcflurry'),torch=torch.__version__,models_path=str(models))
(R/'prediction-runtime-verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
