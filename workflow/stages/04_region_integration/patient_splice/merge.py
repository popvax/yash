"""Merge disjoint splice shards and prove complete event accounting."""
from collections import Counter
import csv
import datetime
import gzip
import hashlib
import json
from pathlib import Path
import re
import sys


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows(path):
    return list(csv.DictReader(Path(path).open(), delimiter="\t"))


def table(path, records):
    fields = list(dict.fromkeys(k for record in records for k in record))
    with Path(path).open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)


def main(root, shards):
    receipts, events, audit, windows, evidence, transcripts = [], [], [], [], {}, {}
    for index in range(shards):
        folder = root / "reconstruction-shards" / str(index)
        receipt = json.loads((folder / "receipt.json").read_text())
        assert receipt["status"] == "COMPLETED" and receipt["shard_index"] == index and receipt["shards"] == shards
        if receipts:
            assert receipt["input_sha256"] == receipts[0]["input_sha256"] and receipt["policy"] == receipts[0]["policy"]
        receipts.append(receipt)
        events.extend(rows(folder / "input-events.tsv"))
        audit.extend(rows(folder / "event-audit.tsv"))
        for row in rows(folder / "reference-transcripts.tsv"):
            key = row["transcript_id"]
            assert key not in transcripts or row == transcripts[key]
            transcripts[key] = row
        with gzip.open(folder / "patient-windows.jsonl.gz", "rt") as handle:
            windows.extend(json.loads(line) for line in handle)
        with gzip.open(folder / "native-event-evidence.json.gz", "rt") as handle:
            child = json.load(handle)
            assert not (set(evidence) & set(child))
            evidence.update(child)
    total = receipts[0]["total_pool_events"]
    assert len(events) == len({r["event_id"] for r in events}) == total
    assert len(audit) == len({r["event_id"] for r in audit}) == total
    assert {r["event_id"] for r in events} == {r["event_id"] for r in audit}
    assert len(windows) == len({w["window_id"] for w in windows})
    # Surface known annotation and non-linear coordinate geometry explicitly.
    # These are annotations, not extra filters or claims of tumour specificity.
    wanted = {event["gene_id"] for event in events}
    exons, metadata = {}, {}
    with Path(receipts[0]["policy"]["gtf"]).open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            attrs = dict(re.findall(r'(\w+) "([^\"]*)"', fields[8]))
            gene = attrs.get("gene_id", "").split(".")[0]
            if gene not in wanted:
                continue
            tx = attrs["transcript_id"]
            exons.setdefault(tx, []).append((int(fields[3]), int(fields[4])))
            metadata[tx] = (gene, fields[0], fields[6])
    known = set()
    for tx, parts in exons.items():
        gene, chrom, strand = metadata[tx]
        parts.sort()
        for left, right in zip(parts, parts[1:]):
            known.add((gene, chrom, left[1], right[0], strand))
    event_metadata = {}
    for event in events:
        match = re.fullmatch(r"([^:]+):(\d+)-(\d+)\(([+-])\)", event["coordinates"])
        annotated, topology = False, "UNRESOLVED_COORDINATES"
        if match:
            chrom, low, high, strand = match[1], int(match[2]), int(match[3]), match[4]
            annotated = (event["gene_id"], chrom, low, high, strand) in known
            topology = "NON_LINEAR_BOUNDARY_ORDER" if high <= low else "CONTIGUOUS_INTRON_BOUNDARY" if high == low + 1 else "LINEAR_SPLICE_GAP"
        event_metadata[event["event_id"]] = dict(junction_is_annotated=annotated, coordinate_topology=topology)
    for record in [*events, *audit, *windows]:
        record.update(event_metadata[record["event_id"]])
    out = root / "reconstruction"
    out.mkdir(exist_ok=True)
    table(out / "input-events.tsv", sorted(events, key=lambda r: r["event_id"]))
    table(out / "event-audit.tsv", sorted(audit, key=lambda r: r["event_id"]))
    table(out / "reference-transcripts.tsv", [transcripts[k] for k in sorted(transcripts)])
    with gzip.open(out / "patient-windows.jsonl.gz", "wt") as handle:
        for window in sorted(windows, key=lambda w: w["window_id"]):
            handle.write(json.dumps(window) + "\n")
    with gzip.open(out / "native-event-evidence.json.gz", "wt") as handle:
        json.dump(evidence, handle)
    receipt = dict(status="COMPLETED", completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        policy=receipts[0]["policy"], input_sha256=receipts[0]["input_sha256"],
        merged_shards=shards, complete_disjoint_event_partition=True, input_events=total,
        event_counts=dict(Counter(r["status"] for r in audit)), patient_windows=len(windows),
        patient_window_events=len({w["event_id"] for w in windows}),
        native_PASS_window_events=len({w["event_id"] for w in windows if w["native_event_review_tier"] == "PASS"}),
        input_coordinate_topology_counts=dict(Counter(event["coordinate_topology"] for event in events)),
        annotated_junction_patient_window_events=len({w["event_id"] for w in windows if w["junction_is_annotated"]}),
        old_MSNAF_coordinates_checked=receipts[0]["old_MSNAF_coordinates_checked"],
        verified_transcripts=sum(r["reference_translation_exact"] == "True" for r in transcripts.values()),
        shard_receipts=receipts, merge_worker_sha256=digest(Path(__file__)),
        output_sha256={p.name: digest(p) for p in sorted(out.iterdir()) if p.name != "receipt.json"})
    (out / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({k: v for k, v in receipt.items() if k not in {"input_sha256", "policy", "shard_receipts"}}, indent=2), flush=True)


if __name__ == "__main__":
    main(Path(sys.argv[1]), int(sys.argv[2]))
