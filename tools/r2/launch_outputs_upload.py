"""Explicitly launch/resume a verified Modal-to-R2 output transfer."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import re

from common import load_credentials

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--state", type=Path, default=REPO / "work/design-2-r2/transfer-state.json")
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--run-name", default="r2-export-20261107", help="Distinct persistent run directory for an immutable transfer plan")
    parser.add_argument("--retry-failed", action="store_true", help="Start another worker only after the prior handle has failed and its logs were reviewed")
    args = parser.parse_args()
    body = args.plan.read_bytes(); plan = json.loads(body)
    worker = HERE / "upload_modal_outputs.py"
    helper = HERE / "ranged_readback.py"
    digest = hashlib.sha256(worker.read_bytes()).hexdigest()
    plan_digest = hashlib.sha256(body).hexdigest()
    print(json.dumps(dict(execute=args.execute, objects=len(plan["objects"]),
                         bytes=sum(r["bytes"] for r in plan["objects"]),
                         original_normal_FASTQs_included=plan["original_normal_FASTQs_included"],
                         plan_sha256=plan_digest, worker_sha256=digest), indent=2), flush=True)
    if not args.execute:
        return
    import modal
    if args.state.exists():
        state = json.loads(args.state.read_text())
        assert state["plan_sha256"] == plan_digest
        sandbox = modal.Sandbox.from_id(state["sandbox_id"])
        rc = sandbox.poll()
        print(json.dumps(dict(sandbox_id=sandbox.object_id, returncode=rc, resumed=True)), flush=True)
        if rc in [None, 0]:
            assert state["worker_sha256"] == digest, "Executed worker source changed; use its preserved state and code for observation"
            return
        if not args.retry_failed:
            raise RuntimeError("Prior transfer failed; inspect its logs before using --retry-failed")
        archive = args.state.parent / "attempts"; archive.mkdir(exist_ok=True)
        (archive / (state["sandbox_id"] + ".json")).write_text(json.dumps(state, indent=2) + "\n")
    load_credentials()
    assert re.fullmatch(r"r2-[a-z0-9-]+", args.run_name)
    remote = "/data/results/TN25-112103/matched-normal-20261107/" + args.run_name
    code = remote + "/code-" + digest[:12]
    volume = modal.Volume.from_name(plan["volume"])
    for path, key in [(worker, code + "/worker.py"), (helper, code + "/ranged_readback.py"), (args.plan, remote + "/plan.json")]:
        key = key.removeprefix("/data/")
        try:
            with volume.batch_upload() as upload:
                upload.put_file(str(path), key)
        except FileExistsError:
            assert b"".join(volume.read_file(key)) == path.read_bytes()
    secret = modal.Secret.from_dict({name: os.environ[name] for name in
        ["R2_ACCESS_KEY_ID", "R2_SECRET_ACCESS_KEY", "R2_ENDPOINT", "R2_BUCKET"]})
    image = modal.Image.debian_slim(python_version="3.11").pip_install("boto3==1.43.123")
    with modal.enable_output():
        sandbox = modal.Sandbox.create("python", "-u", code + "/worker.py", remote + "/plan.json", remote,
            app=modal.App.lookup("research-yash-data-export", create_if_missing=True), image=image,
            cpu=2, memory=8192, timeout=21600, volumes={"/data": volume}, secrets=[secret])
    args.state.parent.mkdir(parents=True, exist_ok=True)
    state = dict(sandbox_id=sandbox.object_id, worker_sha256=digest, plan_sha256=plan_digest,
                 remote=remote, volume=plan["volume"], image_id=image.object_id, credentials_on_volume=False,
                 ranged_readback_helper_sha256=hashlib.sha256(helper.read_bytes()).hexdigest())
    args.state.write_text(json.dumps(state, indent=2) + "\n")
    print(json.dumps(state, indent=2), flush=True)


if __name__ == "__main__":
    main()
