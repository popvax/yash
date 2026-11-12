"""Plan or run explicit affinity requests on Modal; never starts compute implicitly."""
from yash_vaccine.paths import resolve_input
import argparse
import datetime
from pathlib import Path
from common import HERE, REPO, read_json, write_json, terminal_request, load_predictions


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--request", type=Path, help="Omit to predict all terminal choices")
    p.add_argument("--output", type=Path, default=HERE / "terminal-inference")
    p.add_argument("--execute", action="store_true")
    p.add_argument("--image", default="im-DBMghmAQfB3Oz7McdNZeic")
    p.add_argument("--volume", default="research-lens-workspace")
    p.add_argument("--models", default="/data/results/TN25-112103/openvax-mutations-20261018/references/mhcflurry/2.3.0/models_class1_presentation/models")
    args = p.parse_args()
    req = read_json(args.request) if args.request else terminal_request()
    state = dict(image=args.image, volume=args.volume, models=args.models, gpu="A10", cpu=4,
        memory_mib=16384, unique_peptides=len(req["peptides"]),
        peptide_hla_pairs=len(req["peptides"])*len(req["alleles"]), junction_options=len(req["joins"]),
        output=str(args.output.resolve()), execute=args.execute, terminated=False)
    print(__import__("json").dumps(state, indent=2), flush=True)
    if not args.execute:
        return
    if resolve_input(args.output).exists():
        p.error("Choose a new output directory")
    import modal
    args.output.mkdir(parents=True)
    write_json(args.output / "inputs.json", req)
    sandbox = None
    remote = "/tmp/routed-joint"
    try:
        sandbox = modal.Sandbox.create(app=modal.App.lookup("research-yash-vaccine", create_if_missing=True),
            image=modal.Image.from_id(args.image), gpu="A10", cpu=4, memory=16384, timeout=1200,
            volumes={"/data": modal.Volume.from_name(args.volume)})
        state.update(sandbox_id=sandbox.object_id, started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
        write_json(args.output / "modal-state.json", state)
        for local, name in [(HERE / "predict.py", "predict.py"), (args.output / "inputs.json", "inputs.json"),
                            (REPO / "data/checkpoints/provenance/construct-model-files.json", "models.json")]:
            sandbox.filesystem.write_bytes(resolve_input(local).read_bytes(), remote + "/" + name)
        proc = sandbox.exec("python", remote + "/predict.py", "--inputs", remote + "/inputs.json",
            "--model-manifest", remote + "/models.json", "--models", args.models,
            "--output", remote + "/output", timeout=900)
        for line in proc.stdout:
            print(line, end="", flush=True)
        (args.output / "stderr.log").write_text(proc.stderr.read())
        state["exit_code"] = proc.wait()
        if state["exit_code"] != 0:
            raise RuntimeError("Inference failed; inspect stderr.log")
        for name in ["predictions.tsv", "inference.json"]:
            (args.output / name).write_bytes(sandbox.filesystem.read_bytes(remote + "/output/" + name))
        load_predictions(args.output, req)
        state["completed"] = True
    finally:
        if sandbox:
            sandbox.terminate(wait=True)
            state["terminated"] = True
        state["finished_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        write_json(args.output / "modal-state.json", state)


if __name__ == "__main__":
    main()
