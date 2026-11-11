"""Wait for actual successful transfers, collect proofs and publish their catalog."""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent


def digest(body):
    return hashlib.sha256(body).hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--transfer-directory", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    if not args.execute:
        print(json.dumps(dict(execute=False, transfers=list(map(str, args.transfer_directory)), output=str(args.output))))
        return
    import modal
    args.output.mkdir(parents=True, exist_ok=True)
    objects = []; proofs = {}; sources = []
    for number, directory in enumerate(args.transfer_directory):
        state = json.loads((directory / "transfer-state.json").read_text())
        body = (directory / "plan.json").read_bytes(); plan = json.loads(body)
        assert digest(body) == state["plan_sha256"]
        sandbox = modal.Sandbox.from_id(state["sandbox_id"])
        while True:
            try:
                rc = sandbox.poll()
            except Exception as error:
                print("OBSERVATION_ERROR_RETRY_SAME_HANDLE", type(error).__name__, sandbox.object_id, flush=True)
                time.sleep(30); continue
            if rc is not None:
                break
            print("WAITING_FOR_TRANSFER", sandbox.object_id, flush=True); time.sleep(30)
        assert rc == 0, "Transfer failed; inspect the same native handle before retry: " + sandbox.object_id
        volume = modal.Volume.from_name(state["volume"])
        receipt_body = b"".join(volume.read_file(state["remote"].removeprefix("/data/") + "/verification.json"))
        receipt = json.loads(receipt_body)
        assert receipt["status"] == "VERIFIED_BY_FULL_GET_SHA256" and receipt["plan_sha256"] == digest(body)
        assert set(receipt["objects"]) == {r["key"] for r in plan["objects"]}
        (directory / "verification.json").write_bytes(receipt_body)
        (args.output / f"transfer-{number}-verification.json").write_bytes(receipt_body)
        assert not set(proofs).intersection(receipt["objects"])
        objects.extend(plan["objects"]); proofs.update(receipt["objects"])
        sources.append(dict(sandbox_id=sandbox.object_id, native_exit_code=rc, plan_sha256=digest(body),
                            verification_sha256=digest(receipt_body), worker_sha256=state["worker_sha256"],
                            ranged_readback_helper_sha256=state.get("ranged_readback_helper_sha256")))
    for name in ["before-index.json", "before-remote-manifest.json", "before-bundled-manifest.json"]:
        shutil.copyfile(args.transfer_directory[0] / name, args.output / name)
    plan = dict(schema_version=1, bucket="yash-rna-foundry", volume="research-lens-workspace", design="design-2",
                original_normal_FASTQs_included=False, objects=objects)
    save(args.output / "plan.json", plan)
    receipt = dict(status="VERIFIED_BY_FULL_GET_SHA256", objects=proofs, source_transfer_receipts=sources,
                   original_normal_FASTQs_included=False, plan_sha256=digest((args.output / "plan.json").read_bytes()))
    save(args.output / "verification.json", receipt)
    subprocess.run([sys.executable, str(HERE / "publish_outputs_catalog.py"), "--work", str(args.output), "--execute"], check=True)
    print("TRANSFER_VERIFICATION_AND_CATALOG_PUBLICATION_COMPLETE", len(objects), flush=True)


if __name__ == "__main__":
    main()
