"""Prepare a normal alignment on Modal; --execute starts compute."""
from pathlib import Path
import argparse
import datetime as dt
import hashlib
import json

CODE = Path(__file__).resolve().parent
HERE = CODE.parent
REPO = HERE.parents[1]
REMOTE = "/data/results/TN25-112103/matched-normal-20261107"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--retry-failed", action="store_true")
    args = parser.parse_args()
    plan = dict(stage="normal_alignment", remote=REMOTE,
                volume="research-lens-workspace", cpu=32, memory_mib=131072,
                bwa_mem2="2.3", fastp_mode="QC-only, no read trimming/filtering",
                execute=args.execute)
    print(json.dumps(plan, indent=2), flush=True)
    if not args.execute:
        return
    import modal
    output = REPO / "work/matched-normal-20261107"
    output.mkdir(parents=True, exist_ok=True)
    sf = output / "alignment-state.json"
    worker = CODE / "normal_alignment.py"
    digest = hashlib.sha256(worker.read_bytes()).hexdigest()
    if sf.exists():
        old = json.loads(sf.read_text())
        sb = modal.Sandbox.from_id(old["sandbox_id"])
        rc = sb.poll()
        print(json.dumps({"existing_sandbox": sb.object_id, "returncode": rc}), flush=True)
        if rc is None or rc == 0:
            if old["worker_sha256"] != digest:
                raise RuntimeError("A live/successful worker differs; preserve its executed source")
            return
        if not args.retry_failed:
            raise RuntimeError("Prior alignment is terminal; inspect its logs before a new attempt")
        archive = output / "attempts" / old["sandbox_id"]
        archive.mkdir(parents=True, exist_ok=True)
        for path in [sf, output / "normal_alignment.executed.py", output / "alignment-attempt-1-stop.json"]:
            if path.exists():
                path.rename(archive / path.name)
    image = (modal.Image.from_id("im-sxSffXmGfIo7J1b56AqMFq")
             .apt_install("fastp", "pigz")
             .run_commands(
                 "curl -fL https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.3/bwa-mem2-2.3_x64-linux.tar.bz2 -o /tmp/bwa-mem2.tar.bz2",
                 "echo '112f3a3ebf3f8c2377f52d61446ead392c4a98ecf89e7629617d3ed16f4e73cb  /tmp/bwa-mem2.tar.bz2' | sha256sum -c -",
                 "tar -xjf /tmp/bwa-mem2.tar.bz2 -C /opt",
             ))
    volume = modal.Volume.from_name("research-lens-workspace")
    remote_worker = REMOTE + "/code/normal_alignment-" + digest[:12] + ".py"
    with volume.batch_upload() as upload:
        upload.put_file(str(worker), remote_worker.removeprefix("/data/"))
    with modal.enable_output():
        sb = modal.Sandbox.create(
            "python", "-u", remote_worker, REMOTE,
            app=modal.App.lookup("research-yash-matched-normal", create_if_missing=True),
            image=image, cpu=32, memory=131072, timeout=86400,
            volumes={"/data": volume},
        )
    state = dict(plan, sandbox_id=sb.object_id, worker_sha256=digest,
                 image_id=image.object_id, remote_worker=remote_worker,
                 started_utc=dt.datetime.now(dt.timezone.utc).isoformat())
    sf.write_text(json.dumps(state, indent=2) + "\n")
    (output / "normal_alignment.executed.py").write_bytes(worker.read_bytes())
    print(json.dumps(state, indent=2), flush=True)


if __name__ == "__main__":
    main()
