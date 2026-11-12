"""Add coordinate-aware normal and limited strand checks to existing splice evidence."""
import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["measure", "status", "collect", "apply", "verify"])
    parser.add_argument("--source", type=Path, required=True, help="Existing complete local splice-evidence directory")
    parser.add_argument("--output", type=Path, required=True, help="New analysis directory; source is never modified")
    parser.add_argument("--remote-source", default="/data/results/TN25-112103/patient-rna-splices-20261103-final")
    parser.add_argument("--remote-root", default="/data/results/TN25-112103/splice-normal-strand-20261104")
    parser.add_argument("--checks-policy", type=Path, default=HERE / "checks-policy.json")
    parser.add_argument("--retry-failed", action="store_true")
    parser.add_argument("--execute", action="store_true", help="Required for cloud work/downloads; apply and verify are local")
    args = parser.parse_args()
    if args.source.resolve() == args.output.resolve():
        parser.error("Use a new output directory; preserve source evidence")
    for path in (args.remote_root, args.remote_source):
        if not path.startswith("/data/results/TN25-112103/") or ".." in path.split("/"):
            parser.error("Use case-specific paths under /data/results/TN25-112103/")
    if args.remote_root == args.remote_source:
        parser.error("The remote source must remain unchanged")
    print(json.dumps(dict(stage=args.stage, source=str(args.source.resolve()), output=str(args.output.resolve()),
        remote_source=args.remote_source, remote_root=args.remote_root,
        new_model_inference=False, distinct_pattern_filter_added=False,
        substitution_filter_added=False, execute=args.execute), indent=2), flush=True)
    if args.stage in ("apply", "verify"):
        from yash_vaccine.splice_checks import apply, verify
        result = {"apply": apply, "verify": verify}[args.stage](args.source, args.output)
        print(json.dumps(result, indent=2), flush=True)
        return
    if not args.execute:
        return
    import modal
    volume = modal.Volume.from_name("research-lens-workspace")
    args.output.mkdir(parents=True, exist_ok=True)
    if args.stage in ("status", "collect"):
        for mode in ("normal", "strand"):
            state_path = args.output / (mode + "-state.json")
            if not state_path.is_file():
                raise RuntimeError("No recorded launch for " + mode)
            state = json.loads(state_path.read_text())
            sandbox = modal.Sandbox.from_id(state["sandbox_id"])
            code = sandbox.poll()
            print(json.dumps(dict(mode=mode, sandbox_id=state["sandbox_id"], exit_code=code)), flush=True)
            if code is None:
                continue
            (args.output / (mode + "-stdout.log")).write_text(sandbox.stdout.read())
            (args.output / (mode + "-stderr.log")).write_text(sandbox.stderr.read())
            if code != 0:
                raise RuntimeError(mode + " job failed; inspect its saved stderr")
            if args.stage == "collect":
                prefix = state["remote_root"].removeprefix("/data/") + "/"
                receipt_name = mode + "-receipt.json"
                record = b"".join(volume.read_file(prefix + receipt_name))
                receipt = json.loads(record)
                assert receipt["worker_sha256"] == state["uploaded_sha256"]["check_measurements.py"]
                assert receipt["checks_policy_sha256"] == state["uploaded_sha256"]["checks-policy.json"]
                assert receipt["request_sha256"] == state["uploaded_sha256"]["request.json"]
                for name, sha in receipt["output_sha256"].items():
                    payload = b"".join(volume.read_file(prefix + name))
                    if hashlib.sha256(payload).hexdigest() != sha:
                        raise ValueError("Downloaded checksum mismatch: " + name)
                    (args.output / name).write_bytes(payload)
                (args.output / receipt_name).write_bytes(record)
        return
    request = dict(source_reconstruction_sha256=digest(args.source / "reconstruction/patient-windows.jsonl.gz"),
        source_events_sha256=digest(args.source / "reconstruction/input-events.tsv"),
        source_policy_sha256=digest(args.source / "policy.json"))
    policy_bytes = args.checks_policy.read_bytes()
    for name, data in [("request.json", (json.dumps(request, indent=2) + "\n").encode()),
                       ("checks-policy.json", policy_bytes)]:
        target = args.output / name
        if target.exists() and target.read_bytes() != data:
            raise RuntimeError("An existing launch has different inputs: use a new output/remote root")
        target.write_bytes(data)
    uploads = {"check_measurements.py": HERE / "check_measurements.py",
               "request.json": args.output / "request.json", "checks-policy.json": args.output / "checks-policy.json"}
    # Check every existing handle before uploading or starting any new work.
    launch = []
    for mode in ("normal", "strand"):
        state_path = args.output / (mode + "-state.json")
        if state_path.exists():
            state = json.loads(state_path.read_text())
            code = modal.Sandbox.from_id(state["sandbox_id"]).poll()
            if not args.retry_failed or code in (None, 0):
                raise RuntimeError(mode + " already launched; status/collect it. Retry is only for confirmed failures.")
            state_path.rename(args.output / (mode + "-failed-" + state["sandbox_id"] + ".json"))
        launch.append(mode)
    with volume.batch_upload(force=True) as upload:
        for name, path in uploads.items():
            upload.put_file(str(path), args.remote_root.removeprefix("/data/") + "/" + name)
    app = modal.App.lookup("yash-splice-rna-checks", create_if_missing=True)
    for mode in launch:
        image = modal.Image.from_id("im-hX1h6HknDMvasm9Z6q530H")
        if mode == "normal":
            image = image.pip_install("h5py==3.12.1")
        sandbox = modal.Sandbox.create("python", "-u", args.remote_root + "/check_measurements.py", mode,
            args.remote_root, args.remote_source, app=app, image=image,
            cpu=4, memory=24576 if mode == "normal" else 16384, timeout=7200,
            volumes={"/data": volume})
        state = dict(mode=mode, sandbox_id=sandbox.object_id, remote_root=args.remote_root,
            remote_source=args.remote_source, controller_sha256=digest(__file__),
            uploaded_sha256={name: digest(path) for name, path in uploads.items()}, modifies_original_inputs=False)
        save(args.output / (mode + "-state.json"), state)
        print("LAUNCHED", mode, sandbox.object_id, flush=True)


if __name__ == "__main__":
    main()
