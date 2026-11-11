"""Plan/launch patient-RNA splice reconstruction and fresh HLA-I scoring."""
import argparse
import hashlib
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["reconstruct", "merge", "review", "reference", "predict-score", "verify", "status", "collect"])
    parser.add_argument("--remote-root", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--proteins", type=Path)
    parser.add_argument("--shards", type=int, default=1)
    parser.add_argument("--scoring-name", default="scoring", help="New scoring attempt directory; reconstruction is reused without modification")
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    if not args.scoring_name.startswith("scoring") or not all(c.isalnum() or c in "-_" for c in args.scoring_name):
        parser.error("--scoring-name must start with scoring and contain only letters, digits, - or _")
    if not 1 <= args.shards <= 16:
        parser.error("--shards must be between 1 and 16")
    if not args.remote_root.startswith("/data/results/TN25-112103/"):
        parser.error("Use a new case result directory under /data/results/TN25-112103/")
    repo = Path(__file__).resolve().parents[4]
    here = Path(__file__).resolve().parent
    plan = dict(stage=args.stage, remote_root=args.remote_root,
                output=str(args.output.resolve()), execute=args.execute,
                volume="research-lens-workspace", image="im-hX1h6HknDMvasm9Z6q530H",
                scoring_directory=args.scoring_name)
    if args.stage == "reconstruct":
        plan["runtime_additions"] = ["cyvcf2==0.31.4"]
    if args.stage == "reference":
        plan.update(image="im-DBMghmAQfB3Oz7McdNZeic", runtime_additions=["pyahocorasick==2.1.0"])
    if args.stage == "predict-score":
        plan.update(image="im-DBMghmAQfB3Oz7McdNZeic", gpu="A10")
    print(json.dumps(plan, indent=2), flush=True)
    if not args.execute:
        return
    import modal
    volume = modal.Volume.from_name(plan["volume"])
    args.output.mkdir(parents=True, exist_ok=True)
    if args.stage == "status":
        for path in sorted(args.output.glob("*-state.json")):
            state = json.loads(path.read_text())
            for run in state.get("runs", [state]):
                sandbox = modal.Sandbox.from_id(run["sandbox_id"])
                code = sandbox.poll()
                print(json.dumps(dict(stage=run["stage"], shard_index=run.get("shard_index", 0),
                      sandbox_id=run["sandbox_id"], exit_code=code,
                      status="RUNNING" if code is None else "COMPLETED" if code == 0 else "FAILED")), flush=True)
                if code is not None:
                    (args.output / (run["sandbox_id"] + ".stdout.log")).write_text(sandbox.stdout.read())
                    (args.output / (run["sandbox_id"] + ".stderr.log")).write_text(sandbox.stderr.read())
        return
    if args.stage == "collect":
        from concurrent.futures import ThreadPoolExecutor
        remote = args.remote_root.removeprefix("/data/")
        def fetch(folder, name, expected=None):
            target = args.output / folder / name
            target.parent.mkdir(parents=True, exist_ok=True)
            if expected and target.is_file() and hashlib.sha256(target.read_bytes()).hexdigest() == expected:
                return
            remote_folder = args.scoring_name if folder == "scoring" else folder
            remote_path = "/".join(part for part in (remote, remote_folder, name) if part)
            data = b"".join(volume.read_file(remote_path))
            if expected and hashlib.sha256(data).hexdigest() != expected:
                raise ValueError("Downloaded checksum mismatch: " + folder + "/" + name)
            temporary = target.with_name(target.name + ".download")
            temporary.write_bytes(data)
            temporary.replace(target)
        jobs = [("", name, None) for name in ["policy.json", "model-files.json"]]
        folders = {
            "reconstruction": ["input-events.tsv", "reference-transcripts.tsv", "patient-windows.jsonl.gz",
                               "native-event-evidence.json.gz", "event-audit.tsv"],
            "review": ["peptide-review.json.gz"],
            "reference-check": ["matches.json.gz"],
            "scoring": ["affinities.tsv.gz", "reference-matches.json.gz", "all-windows.tsv.gz",
                        "event-ranking.tsv", "gene-ranking.tsv", "selected.tsv", "selected.fasta",
                        "native-score-config.json", "verification.json"]}
        for folder, names in folders.items():
            fetch(folder, "receipt.json")
            receipt = json.loads((args.output / folder / "receipt.json").read_text())
            hashes = receipt.get("output_sha256", {})
            if isinstance(hashes, str):
                hashes = {names[0]: hashes}
            jobs.extend((folder, name, hashes.get(name)) for name in names)
        jobs.append(("read-verification", "verification.json", None))
        with ThreadPoolExecutor(max_workers=6) as pool:
            list(pool.map(lambda job: fetch(*job), jobs))
        print("COLLECTED", len(jobs), "verified files", flush=True)
        return
    state_name = args.stage + ("-" + args.scoring_name if args.stage == "predict-score" and args.scoring_name != "scoring" else "") + "-state.json"
    state_path = args.output / state_name
    if state_path.exists():
        previous = json.loads(state_path.read_text())
        codes = [modal.Sandbox.from_id(run["sandbox_id"]).poll()
                 for run in previous.get("runs", [previous])]
        if any(code is None for code in codes):
            raise RuntimeError("This stage has a live sandbox. Inspect its handle; do not launch a duplicate.")
        if codes and all(code == 0 for code in codes):
            raise RuntimeError("This stage already completed. Collect its results or choose a new result directory.")
    local = {"policy.json": here / "policy.json", "splice_rna.py": repo / "src/yash_vaccine/splice_rna.py",
             "msnaf-core.8-11.py": here.parent.parent / "02_lens/msnaf-core.8-11.py"}
    if args.stage == "reconstruct":
        if not args.proteins or not args.proteins.is_file():
            parser.error("--proteins must name the complete pinned GENCODE37 protein FASTA.gz")
        local.update({"reconstruct.py": here / "reconstruct.py", "gencode37.proteins.fasta.gz": args.proteins})
    elif args.stage == "review":
        local.update({"review.py": here / "review.py"})
    elif args.stage == "verify":
        local.update({"verify.py": here / "verify.py"})
    elif args.stage == "merge":
        local.update({"merge.py": here / "merge.py"})
    elif args.stage == "reference":
        local.update({"reference_exact.py": here / "reference_exact.py"})
    else:
        local.update({"predict_score.py": here / "predict_score.py", "model-files.json": repo / "designs/initial/review/model-files.json"})
    with volume.batch_upload(force=True) as upload:
        for name, path in local.items():
            upload.put_file(str(path), args.remote_root.removeprefix("/data/") + "/" + name)
    program = {"reconstruct": "reconstruct.py", "merge": "merge.py", "review": "review.py", "reference": "reference_exact.py", "predict-score": "predict_score.py", "verify": "verify.py"}[args.stage]
    image = modal.Image.from_id(plan["image"])
    if plan.get("runtime_additions"):
        image = image.pip_install(*plan["runtime_additions"])
    (args.output / "attempts").mkdir(exist_ok=True)
    states = []
    for index in range(args.shards if args.stage == "reconstruct" else 1):
        command = ["python", "-u", args.remote_root + "/" + program, args.remote_root]
        if args.stage == "reconstruct":
            command.extend([str(index), str(args.shards)])
        elif args.stage == "merge":
            command.append(str(args.shards))
        elif args.stage == "predict-score":
            command.append(args.scoring_name)
        sandbox = modal.Sandbox.create(*command,
            app=modal.App.lookup("yash-splice-rna", create_if_missing=True),
            image=image, cpu=4, memory=16384,
            timeout=7200, volumes={"/data": volume},
            **({"gpu": "A10"} if args.stage == "predict-score" else {}))
        state = dict(plan, sandbox_id=sandbox.object_id, shard_index=index, shards=args.shards,
                     controller_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                     uploaded_sha256={name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in local.items()})
        (args.output / "attempts" / (sandbox.object_id + ".json")).write_text(json.dumps(state, indent=2) + "\n")
        states.append(state)
        state_path.write_text(json.dumps(dict(runs=states), indent=2) + "\n")
        print("LAUNCHED", index, sandbox.object_id, flush=True)


if __name__ == "__main__":
    main()
