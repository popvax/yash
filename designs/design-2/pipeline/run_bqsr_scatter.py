"""Apply a gathered BQSR table on disjoint chromosomes after a native parity probe."""
from pathlib import Path
import argparse
import concurrent.futures
import hashlib
import json
import threading

CODE = Path(__file__).resolve().parent
HERE = CODE.parent;REPO=HERE.parents[1]
ROOT="/data/results/TN25-112103/matched-normal-20261107"
LOCAL=REPO/"work/matched-normal-20261107"


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execute",action="store_true")
    parser.add_argument("--workers",type=int,default=27,choices=range(1,33))
    args=parser.parse_args()
    print(json.dumps(dict(stage="BQSR chromosome-disjoint application",execute=args.execute,workers=args.workers)),flush=True)
    if not args.execute:return
    import modal
    volume=modal.Volume.from_name("research-lens-workspace")
    files=[CODE / "bqsr_scatter.py",CODE / "common.py"]
    digest=hashlib.sha256(b"".join(p.read_bytes() for p in files)).hexdigest()
    code_root=ROOT+"/code/bqsr-scatter-"+digest[:12]
    for path in files:
        dest=(code_root+"/"+path.name).removeprefix("/data/")
        try:
            with volume.batch_upload() as upload:upload.put_file(str(path),dest)
        except FileExistsError:assert b"".join(volume.read_file(dest))==path.read_bytes()
    state_path=LOCAL/"bqsr-scatter-state.json"
    state=json.loads(state_path.read_text()) if state_path.exists() else dict(bundle_sha256=digest,workers={})
    assert state["bundle_sha256"]==digest
    lock=threading.Lock()
    def save_state():
        temp=state_path.with_suffix(".tmp");temp.write_text(json.dumps(state,indent=2)+"\n");temp.replace(state_path)
    def launch(action,shard=None):
        name=action+("-"+shard["id"] if shard else "")
        with lock:old=state["workers"].get(name)
        # Reattach the already launched probe when compatible.
        if name=="probe" and old is None:
            early_path=LOCAL/"bqsr-scatter-probe-state.json"
            if early_path.exists():
                early=json.loads(early_path.read_text())
                assert early["bundle_sha256"]==digest
                old=dict(sandbox_id=early["sandbox_id"],status="running")
        if old:
            sb=modal.Sandbox.from_id(old["sandbox_id"])
            assert sb.poll() in [None,0],name+" failed; inspect before retrying"
        else:
            sb=modal.Sandbox.create("python","-u",code_root+"/bqsr_scatter.py",ROOT,action,
                  *([shard["id"]] if shard else []),
                  app=modal.App.lookup("research-yash-matched-normal",create_if_missing=True),
                  image=modal.Image.from_id("im-sxSffXmGfIo7J1b56AqMFq"),cpu=8 if action=="gather" else 4,
                  memory=16384 if action=="gather" else 8192,timeout=43200,volumes={"/data":volume})
        with lock:
            state["workers"][name]=dict(sandbox_id=sb.object_id,status="running");save_state()
        if action=="gather":
            # The existing continuation controller observes this final native
            # gather/validation worker as the calibrated-BAM completion handle.
            base_path=LOCAL/"bqsr-state.json";base=json.loads(base_path.read_text())
            base["application_strategy"]="chromosome-disjoint native ApplyBQSR"
            base["application_bundle_sha256"]=digest
            base["workers"]["apply"]=dict(sandbox_id=sb.object_id,status="running",strategy="scatter gather and validation")
            temp=base_path.with_suffix(".tmp");temp.write_text(json.dumps(base,indent=2)+"\n");temp.replace(base_path)
        print("START",name,sb.object_id,flush=True)
        stdout=sb.stdout.read();stderr=sb.stderr.read();sb.wait();rc=sb.poll()
        (LOCAL/("scatter-"+name+".stdout.log")).write_text(stdout)
        (LOCAL/("scatter-"+name+".stderr.log")).write_text(stderr)
        with lock:
            state["workers"][name].update(status="completed" if rc==0 else "failed",exit_code=rc);save_state()
        if rc:raise RuntimeError(name+" failed: "+stderr[-4000:])
        remote=("bqsr/scatter-probe" if action=="probe" else "bqsr/scatter-apply/"+shard["id"] if action=="apply" else "bqsr/apply")
        rec=json.loads(b"".join(volume.read_file(ROOT.removeprefix("/data/")+"/"+remote+"/complete.json")))
        assert rec["status"]=="COMPLETED"
        (LOCAL/("scatter-"+name+".complete.json")).write_text(json.dumps(rec,indent=2)+"\n")
        print("DONE",name,flush=True)
        return rec
    probe=launch("probe")
    proof=b"".join(volume.read_file(ROOT.removeprefix("/data/")+"/bqsr/scatter-probe/equivalence.json"))
    (LOCAL/"BQSR-scatter-equivalence.json").write_bytes(proof)
    assert json.loads(proof)["identical_recalibrated_records"]
    # Read the exact original BAM dictionary order, not an assumed reference order.
    dictionary=[]
    for line in (LOCAL/"normal-idxstats.txt").read_text().splitlines():
        chrom,length,mapped,unmapped=line.split("\t")
        if chrom!="*":dictionary.append(chrom)
    primary={"chrM","chrX","chrY",*["chr"+str(i) for i in range(1,23)]}
    shards=[dict(id=chrom,intervals=[chrom]) for chrom in dictionary if chrom in primary]
    others=[c for c in dictionary if c not in primary]
    if others:shards.append(dict(id="other-contigs",intervals=others))
    shards.append(dict(id="no-coordinate",intervals=["unmapped"]))
    assert [c for s in shards[:-1] for c in s["intervals"]]==dictionary
    plan=dict(shards=shards,all_contigs_once=True,no_coordinate_tail_separate=True,
              probe_sha256=hashlib.sha256(proof).hexdigest(),bundle_sha256=digest,
              final_count_requirement="Identical mapped/unmapped counts on every contig and identical no-coordinate record count")
    plan_path=LOCAL/"bqsr-scatter-plan.json";plan_path.write_text(json.dumps(plan,indent=2)+"\n")
    dest=ROOT.removeprefix("/data/")+"/bqsr-scatter-plan.json"
    try:
        with volume.batch_upload() as up:up.put_file(str(plan_path),dest)
    except FileExistsError:assert b"".join(volume.read_file(dest))==plan_path.read_bytes()
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        list(pool.map(lambda s:launch("apply",s),shards))
    final=launch("gather")
    assert final["all_original_records_preserved"]
    state["status"]="COMPLETED";save_state()
    base_path=LOCAL/"bqsr-state.json";base=json.loads(base_path.read_text())
    base["workers"]["apply"].update(status="completed",exit_code=0)
    base["status"]="COMPLETED";base_path.write_text(json.dumps(base,indent=2)+"\n")
    print("BQSR_SCATTER_APPLICATION_COMPLETE",flush=True)


if __name__=="__main__":main()
