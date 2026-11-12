"""Run inside the recorded CUDA environment; all queries are newly predicted."""
import os
os.environ['OMP_NUM_THREADS']='4'
os.environ['MKL_NUM_THREADS']='4'
import argparse, hashlib, importlib.metadata, inspect, json, sys, time
from pathlib import Path
import numpy as np

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',type=Path,required=True)
    p.add_argument('--models',type=Path,required=True)
    p.add_argument('--regions',type=Path,required=True)
    p.add_argument('--linkers',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();root=args.repo
    if args.output.exists():raise ValueError('Choose a new prediction output directory')
    sys.path.insert(0,str(root/'src'))
    from yash_vaccine.paths import resolve_input
    from yash_vaccine.io import rows,write_json,sha256
    from yash_vaccine.constructs.boundaries import all_queries,junction_peptides,LENGTHS
    import torch
    torch.set_num_threads(4);torch.set_float32_matmul_precision('highest')
    torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False
    assert torch.cuda.is_available()
    from mhcflurry import Class1AffinityPredictor
    from mhctools.mhcflurry import MHCflurry_Affinity
    from mhctools.pred import Kind
    import vaxrank.junction_swap as native_js
    from vaxrank.vaccine_library import JUNCTION_SWAP_CANDIDATES,get_linker
    from yash_vaccine.constructs.inputs import inputs
    regions,linkers,context,region_file,linker_file=inputs(root,args.regions,args.linkers)
    for linker in linkers:
        if linker.get('origin')=='native_default':assert get_linker(linker['name']).amino_acids==linker['sequence']
        if linker['sequence']:
            assert junction_peptides(regions[0]['sequence'],linker['sequence'],regions[1]['sequence'])==native_js.junction_kmers(regions[0]['sequence'],linker['sequence'],regions[1]['sequence'],LENGTHS)
    alleles=context['hla_class_I']
    queries=all_queries(regions,linkers);out=args.output;out.mkdir(parents=True)
    model_manifest=root/'data/checkpoints/provenance/construct-model-files.json'
    weights=[dict(path=str(path.relative_to(args.models)),bytes=resolve_input(path).stat().st_size,sha256=sha256(path))
             for path in sorted(resolve_input(args.models).rglob('*')) if resolve_input(path).is_file()]
    assert weights==json.loads(resolve_input(model_manifest).read_text()),'Model files differ from the pinned construct model bundle'
    start=time.perf_counter();model=Class1AffinityPredictor.load(str(args.models/'affinity_predictor'))
    affinity=np.empty((len(queries),len(alleles)));ranks=np.empty_like(affinity)
    for h,allele in enumerate(alleles):
        df=model.predict_to_dataframe(queries,allele=allele,model_kwargs={'batch_size':16384},
                                     include_percentile_ranks=True,include_confidence_intervals=False)
        assert df.peptide.tolist()==queries
        affinity[:,h]=df.prediction.to_numpy();ranks[:,h]=df.prediction_percentile.to_numpy()
        print('PREDICTED',allele,len(queries),flush=True)
    assert np.isfinite(affinity).all() and (affinity>0).all()
    assert np.isfinite(ranks).all() and ((ranks>=0)&(ranks<=100)).all()
    # Wrapper parity checks cover every allele, all peptide lengths and source blocks.
    check=sorted(set(queries[::max(1,len(queries)//96)]))
    wrapper=MHCflurry_Affinity(alleles=alleles,default_peptide_lengths=LENGTHS,models_path=str(args.models/'affinity_predictor'))
    ix={s:i for i,s in enumerate(queries)};hix={a:i for i,a in enumerate(alleles)}
    diffs=[];rdiffs=[];count=0
    for pep,result in zip(check,wrapper.predict(check)):
        for pred in result.preds:
            if pred.kind==Kind.pMHC_affinity:
                i,h=ix[pep],hix[pred.allele]
                diffs.append(abs(float(pred.value)/affinity[i,h]-1))
                rdiffs.append(abs(float(pred.percentile_rank)-ranks[i,h]));count+=1
    assert count==len(check)*len(alleles)
    assert max(diffs)<1e-4 and max(rdiffs)<1e-4,(max(diffs),max(rdiffs))
    np.savez_compressed(out/'construct-predictions.npz',peptides=np.array(queries),alleles=np.array(alleles),
                        affinity_nm=affinity,affinity_percentile=ranks)
    write_json(out/'construct-model-files.json',weights)
    write_json(out/'construct-inference.json',dict(
        regions_sha256=sha256(region_file),linkers_sha256=sha256(linker_file),regions=len(regions),alleles=alleles,
        gpu=torch.cuda.get_device_name(),cuda_allocated_peak_bytes=torch.cuda.max_memory_allocated(),
        elapsed_seconds=time.perf_counter()-start,queries=len(queries),peptide_hla_pairs=len(queries)*len(alleles),
        all_queries_freshly_predicted=True,processing_or_presentation_model_used=False,
        all_model_files_match_reference=True,prediction_worker_sha256=sha256(Path(__file__)),
        model_manifest_input_sha256=sha256(model_manifest),mhctools_adapter='MHCflurry_Affinity',
        versions={x:importlib.metadata.version(x) for x in ['mhcflurry','mhctools','vaxrank','torch','numpy']},
        native_linker_defaults_verified=True,nonempty_enumerator_matches_native=True,
        mhctools_parity_pairs=count,max_relative_affinity_difference=max(diffs),max_absolute_rank_difference=max(rdiffs),
        model_files_manifest='construct-model-files.json',cache_sha256=sha256(out/'construct-predictions.npz')))
    print('VERIFIED',len(queries),count,flush=True)
if __name__=='__main__':main()
