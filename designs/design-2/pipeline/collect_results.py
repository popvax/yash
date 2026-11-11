"""Collect completed Modal results without launching scientific computation.

Preserve native receipts and bytes. Large files remain on the persistent volume
and are listed with checksums instead of being silently omitted from the record.
"""
from pathlib import Path
import argparse
import hashlib
import json

CODE = Path(__file__).resolve().parent
HERE = CODE.parent
PATHS = json.loads((HERE / 'execution/local-paths.json').read_text())['paths']
REMOTE = "results/TN25-112103/matched-normal-20261107"
VOLUME = "research-lens-workspace"
STAGES = {
    "alignment": "execution/alignment",
    "normal_metadata": "execution/normal-metadata",
    "bqsr": "execution/bqsr",
    "identity": "execution/identity",
    "normal-coding/gather": "execution/normal-coding",
    "paired/gather": "execution/paired-gather",
    "paired/filter": "execution/paired-filter",
    "normal-self": "execution/normal-self",
    "annotate": "execution/annotation",
    "rna": "execution/rna",
    "rank": "execution/rank",
    "integrate": "execution/integration",
    "construct/prepare": "construct/prepare",
    "construct/all-inference": "construct/all-inference",
    "construct/protein": "construct/protein",
    "construct/selected-inference": "construct/selected-inference",
    "construct/verify": "construct/verify",
    "construct/mrna": "construct/mrna",
    "report": "review/final",
}


def digest(body):
    return hashlib.sha256(body).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    print(json.dumps(dict(execute=args.execute, volume=VOLUME, remote=REMOTE,
                          compute_started=False, stages=STAGES), indent=2))
    if not args.execute:
        return
    import modal
    volume = modal.Volume.from_name(VOLUME)
    records = []

    def fetch(remote, local=None, expected=None):
        body = b"".join(volume.read_file(REMOTE + "/" + remote))
        if expected:
            assert len(body) == expected["bytes"]
            assert digest(body) == expected["sha256"], remote
        if local:
            local = PATHS.get(local, local)
            path = HERE / local
            path.parent.mkdir(parents=True, exist_ok=True)
            if path.exists():
                assert path.read_bytes() == body, "Refusing to overwrite different output: " + local
            else:
                path.write_bytes(body)
        records.append(dict(volume=VOLUME, remote="/data/" + REMOTE + "/" + remote,
                            local=local, bytes=len(body), sha256=digest(body)))
        return body

    receipts = {}
    # Native BAM/VCF processing outputs are intentionally not downloaded here.
    full_small_stages = {s for s in STAGES if s.startswith("construct/")} | {"report"}
    for stage, local in STAGES.items():
        body = fetch(stage + "/complete.json", local + "/complete.json")
        receipt = json.loads(body)
        assert receipt["status"] == "COMPLETED", stage
        receipts[stage] = receipt
        for name, metadata in receipt["outputs"].items():
            # The bqsr/complete.json receipt aliases the gathered apply stage;
            # its output files live beside the normal_bam, under bqsr/apply/.
            remote_directory = str(Path(receipt["normal_bam"]).parent.relative_to(Path("/data") / REMOTE)) if stage == "bqsr" else stage
            destination = local + "/" + name if stage in full_small_stages and metadata["bytes"] <= 10_000_000 else None
            if destination:
                fetch(remote_directory + "/" + name, destination, metadata)
            else:
                records.append(dict(volume=VOLUME, remote="/data/" + REMOTE + "/" + remote_directory + "/" + name,
                                    local=None, **metadata))

    for name in ["selected.tsv", "selected.json", "selected.fasta"]:
        fetch("integrate/" + name, "selection/" + name, receipts["integrate"]["outputs"][name])
    for name in ["run.py", "optimizer.py", "genbank.py", "backbone.json", "annotation-sources.json",
                 "human-codon-usage.json", "requirements.txt"]:
        fetch("rebuild/construct/mrna/" + name, "construct/mrna/" + name)
    fetch("construct/architecture.json", "construct/architecture.json")
    fetch("construct/final-check.json", "construct/final-check.json")
    execution = HERE / "execution"
    (execution / "artifacts.json").write_text(json.dumps(records, indent=2) + "\n")
    (execution / "stage-receipts.json").write_text(json.dumps(receipts, indent=2) + "\n")
    print(json.dumps(dict(status="COLLECTED", records=len(records),
                          local_files=sum(bool(r["local"]) for r in records), compute_started=False)))


if __name__ == "__main__":
    main()
