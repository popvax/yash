"""Execute a matched-normal stage explicitly on Modal, retaining immutable code."""
from pathlib import Path
import argparse
import hashlib
import json

CODE = Path(__file__).resolve().parent
HERE = CODE.parent
REPO=HERE.parents[1]
REMOTE="/data/results/TN25-112103/matched-normal-20261107"
VOLUME="research-lens-workspace"
BASE_IMAGE="im-sxSffXmGfIo7J1b56AqMFq"
ANALYSIS_IMAGE="im-hX1h6HknDMvasm9Z6q530H"
ANALYSIS_STAGES={'annotate','rna','normal_self','rank','integrate','report'}
DEPENDENCIES={'identity_precal':['fingerprint_panel.py'],'identity':['fingerprint_panel.py'],
              'annotate':['normal_evidence.py'],'rank':['self_search.py'],'integrate':['self_search.py'],'report':['normal_evidence.py']}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage",choices=["resources","normal_metadata","identity_precal","identity","audit","annotate","rna","normal_self","rank","integrate","report"])
    parser.add_argument("--execute",action="store_true")
    parser.add_argument("--retry-failed",action="store_true")
    args=parser.parse_args()
    worker=CODE/(args.stage+".py")
    assert worker.is_file(),"Stage implementation is not ready: "+args.stage
    common=CODE / "common.py"
    files=[worker,common]+[CODE/name for name in DEPENDENCIES.get(args.stage,[])]
    digest=hashlib.sha256(b''.join(path.read_bytes() for path in files)).hexdigest()
    image=ANALYSIS_IMAGE if args.stage in ANALYSIS_STAGES else BASE_IMAGE
    plan=dict(stage=args.stage,remote=REMOTE,volume=VOLUME,worker_bundle_sha256=digest,
              cpu=4,memory_mib=16384,image_id=image,execute=args.execute)
    print(json.dumps(plan,indent=2),flush=True)
    if not args.execute:return
    import modal
    output=REPO/"work/matched-normal-20261107";output.mkdir(parents=True,exist_ok=True)
    sf=output/(args.stage+"-state.json")
    if sf.exists():
        old=json.loads(sf.read_text());sandbox=modal.Sandbox.from_id(old["sandbox_id"])
        rc=sandbox.poll();print(json.dumps(dict(existing_sandbox=sandbox.object_id,returncode=rc)),flush=True)
        if rc is None or rc==0:
            assert old["worker_bundle_sha256"]==digest,"Preserve executed code for a live/successful worker"
            return
        if not args.retry_failed:
            raise RuntimeError("Prior stage is terminal; inspect logs before retry")
        archive=output/"attempts"/old["sandbox_id"];archive.mkdir(parents=True,exist_ok=True)
        for path in [sf,output/(args.stage+".stdout.log"),output/(args.stage+".stderr.log")]:
            if path.exists():path.rename(archive/path.name)
    directory=REMOTE+"/code/"+args.stage+"-"+digest[:12]
    volume=modal.Volume.from_name(VOLUME)
    input_files=[(HERE/'selection'/name if name=='paired-policy.json' else CODE/'inputs'/name) for name in {'normal_metadata':['platform-default-verification.json'],'annotate':['paired-policy.json'],'normal_self':['patient-hla-proteins.json']}.get(args.stage,[])]
    if args.stage=='integrate':
        from frozen_inputs import package
        input_files.extend(package(REPO,output))
    for path in input_files:
        key=(REMOTE+'/inputs/'+path.name).removeprefix('/data/')
        try:
            with volume.batch_upload() as upload:upload.put_file(str(path),key)
        except FileExistsError:assert b''.join(volume.read_file(key))==path.read_bytes()
    for path in files:
        key=directory.removeprefix("/data/")+"/"+path.name
        try:
            with volume.batch_upload() as upload:upload.put_file(str(path),key)
        except FileExistsError:assert b''.join(volume.read_file(key))==path.read_bytes(),"Immutable worker source differs"
    sandbox=modal.Sandbox.create("python","-u",directory+"/"+worker.name,REMOTE,
        app=modal.App.lookup("research-yash-matched-normal",create_if_missing=True),
        image=modal.Image.from_id(image),cpu=4,memory=16384,timeout=21600,
        volumes={"/data":volume})
    sf.write_text(json.dumps(dict(plan,sandbox_id=sandbox.object_id,remote_code=directory),indent=2)+"\n")
    print(json.dumps(dict(sandbox_id=sandbox.object_id)),flush=True)


if __name__=="__main__":main()
