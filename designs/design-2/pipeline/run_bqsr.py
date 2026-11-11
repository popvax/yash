"""Wait for the verified alignment, then run native scatter/gather BQSR."""
from pathlib import Path
import argparse
import concurrent.futures
import datetime as dt
import hashlib
import json
import subprocess
import sys
import threading
import time

CODE = Path(__file__).resolve().parent
HERE = CODE.parent;REPO=HERE.parents[1]
REMOTE="/data/results/TN25-112103/matched-normal-20261107"


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--execute",action="store_true")
    p.add_argument("--wait-alignment",action="store_true")
    p.add_argument("--workers",type=int,default=24,choices=range(1,33))
    p.add_argument("--application",choices=["scatter","whole"],default="scatter")
    args=p.parse_args();print(json.dumps(dict(stage="BQSR",execute=args.execute,workers=args.workers)),flush=True)
    if not args.execute:return
    import modal
    local=REPO/"work/matched-normal-20261107";local.mkdir(parents=True,exist_ok=True)
    volume=modal.Volume.from_name("research-lens-workspace")
    source=json.loads((local/"alignment-state.json").read_text())
    while True:
        sb=modal.Sandbox.from_id(source["sandbox_id"]);code=sb.poll()
        if code==0:break
        if code is not None:raise RuntimeError("Normal alignment is terminal with exit "+str(code))
        if not args.wait_alignment:raise RuntimeError("Alignment is still running; use --wait-alignment")
        print(json.dumps(dict(waiting_for_alignment=sb.object_id,verified_live=True,utc=dt.datetime.now(dt.timezone.utc).isoformat())),flush=True)
        time.sleep(30)
    alignment=json.loads(b"".join(volume.read_file(REMOTE.removeprefix("/data/")+"/alignment/complete.json")))
    assert alignment["status"]=="COMPLETED"
    (local/"alignment.complete.json").write_text(json.dumps(alignment,indent=2)+"\n")
    # Audit metadata; use a verified native default without rewriting the BAM.
    subprocess.run([sys.executable,str(CODE / "run_stage.py"),"normal_metadata","--execute"],cwd=REPO,check=True)
    metadata_state=json.loads((local/"normal_metadata-state.json").read_text())
    metadata_sandbox=modal.Sandbox.from_id(metadata_state["sandbox_id"])
    stdout=metadata_sandbox.stdout.read();stderr=metadata_sandbox.stderr.read();metadata_sandbox.wait()
    (local/"normal_metadata.stdout.log").write_text(stdout);(local/"normal_metadata.stderr.log").write_text(stderr)
    assert metadata_sandbox.poll()==0,"Platform metadata audit failed; preserve its source and logs"
    metadata=json.loads(b"".join(volume.read_file(REMOTE.removeprefix("/data/")+"/normal_metadata/complete.json")))
    assert metadata["status"]=="COMPLETED" and metadata["platform_metadata"]=="UNKNOWN"
    (local/"normal_metadata.complete.json").write_text(json.dumps(metadata,indent=2)+"\n")
    dictionary=json.loads((CODE/"inputs/reference-contigs.json").read_text())
    assert dictionary["reference_sha256"]==alignment["reference_sha256"]
    lengths={r["name"]:r["length"] for r in dictionary["contigs"]}
    shards=[]
    for chrom in ["chr"+str(c) for c in list(range(1,23))+["X","Y"]]:
        for i,start in enumerate(range(1,lengths[chrom]+1,20000000),1):
            end=min(lengths[chrom],start+20000000-1)
            shards.append(dict(id=f"{chrom}-{i:03d}",chrom=chrom,interval=f"{chrom}:{start}-{end}"))
    plan=dict(normal_bam=metadata["normal_bam"],normal_bam_sha256=metadata["normal_bam_sha256"],
              normal_bam_bytes=metadata["normal_bam_bytes"],source_alignment_sha256=metadata["source_alignment_sha256"],
              platform_metadata="UNKNOWN",actual_instrument_unconfirmed=True,
              recalibration_default_platform="UNKNOWN",original_alignment_header_unchanged=True,
              reference=alignment["reference"],reference_sha256=alignment["reference_sha256"],
              known_sites=[REMOTE+"/resources/"+n for n in ["Homo_sapiens_assembly38.dbsnp138.vcf.gz","Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"]],
              shards=shards,scope="Train across primary nuclear chromosomes; apply once to the complete BAM without interval scattering",
              boundary_note="Native interval-scattered BQSR training can count reads overlapping shard boundaries in both reports; application is unscattered and record counts must remain identical")
    plan_path=local/"bqsr-plan.json";plan_path.write_text(json.dumps(plan,indent=2)+"\n")
    worker=CODE / "bqsr_worker.py";common=CODE / "common.py"
    digest=hashlib.sha256(worker.read_bytes()+common.read_bytes()+plan_path.read_bytes()).hexdigest()
    code_root=REMOTE+"/code/bqsr-"+digest[:12]
    # Upload once. Retries reuse immutable code after digest verification.
    for path,dest in [(worker,code_root+"/bqsr_worker.py"),(common,code_root+"/common.py"),(plan_path,REMOTE+"/bqsr-plan.json")]:
        try:
            with volume.batch_upload() as up:up.put_file(str(path),dest.removeprefix("/data/"))
        except FileExistsError:
            assert b"".join(volume.read_file(dest.removeprefix("/data/")))==path.read_bytes()
    state_path=local/"bqsr-state.json"
    state=json.loads(state_path.read_text()) if state_path.exists() else dict(bundle_sha256=digest,workers={})
    assert state["bundle_sha256"]==digest
    lock=threading.Lock()
    def save_state():
        temporary=state_path.with_suffix(".tmp");temporary.write_text(json.dumps(state,indent=2)+"\n");temporary.replace(state_path)
    def run(action,shard=None):
        name=action+("-"+shard["id"] if shard else "")
        with lock:old=state["workers"].get(name)
        if old:
            sb=modal.Sandbox.from_id(old["sandbox_id"])
            if sb.poll() not in [None,0]:raise RuntimeError(name+" is terminal; inspect before retry")
        else:
            command=["python","-u",code_root+"/bqsr_worker.py",REMOTE,action]+([shard["id"]] if shard else [])
            sb=modal.Sandbox.create(*command,app=modal.App.lookup("research-yash-matched-normal",create_if_missing=True),
                image=modal.Image.from_id("im-sxSffXmGfIo7J1b56AqMFq"),cpu=2 if action=="train" else 4,
                memory=8192 if action=="train" else 16384,timeout=43200,volumes={"/data":volume})
            with lock:
                state["workers"][name]=dict(sandbox_id=sb.object_id,status="running");save_state()
            print("START",name,sb.object_id,flush=True)
        stdout=sb.stdout.read();stderr=sb.stderr.read();sb.wait();rc=sb.poll()
        (local/(name+".stdout.log")).write_text(stdout);(local/(name+".stderr.log")).write_text(stderr)
        with lock:
            state["workers"][name].update(status="completed" if rc==0 else "failed",exit_code=rc);save_state()
        if rc:raise RuntimeError(name+" failed: "+stderr[-2000:])
        receipt_path=("bqsr/training/"+shard["id"] if action=="train" else "bqsr/"+action)+"/complete.json"
        result=json.loads(b"".join(volume.read_file(REMOTE.removeprefix("/data/")+"/"+receipt_path)))
        assert result["status"]=="COMPLETED"
        print("DONE",name,flush=True)
        return result
    # One verified pilot catches environment/input failures before a broad scatter.
    pilot=next(s for s in shards if s["id"]=="chr22-002")
    run("train",pilot)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        list(pool.map(lambda s:run("train",s),[s for s in shards if s!=pilot]))
    run("gather")
    if args.application=="whole":run("apply")
    else:
        subprocess.run([sys.executable,str(CODE / "run_bqsr_scatter.py"),"--execute"],cwd=REPO,check=True)
        state=json.loads(state_path.read_text())
    state["status"]="COMPLETED";save_state()
    print("BQSR_COMPLETE",flush=True)


if __name__=="__main__":main()
