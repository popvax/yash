"""Apply native LENS peptide flags to actual patient coding sequences.

The reference-only peptide-consistency test is superseded by complete patient
window reconstruction. All other native flag rules, including its coding-kmer
support and cysteine penalties, are applied without changing their thresholds.
"""
import ast
from collections import defaultdict
import csv
import datetime
import gzip
import hashlib
import json
from pathlib import Path
import re
import sys

from Bio.Seq import Seq
import pysam
import yaml

from splice_rna import target_peptides


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def key(event, coding):
    return hashlib.sha256((event + "|" + coding).encode()).hexdigest()


def main(root):
    policy = json.loads((root / "policy.json").read_text())
    lens = Path(policy["lens_root"])
    source = lens / "metadata.splice-review.py"
    names = {"_extract_junctions_from_read", "_junction_spanning_read_seq", "_kmer_set",
             "_safe_fraction", "_splice_read_coding_support", "_flags_from_evidence"}
    syntax = ast.parse(source.read_text())
    functions = [node for node in syntax.body if isinstance(node, ast.FunctionDef) and node.name in names]
    assert {node.name for node in functions} == names
    native = {"Seq": Seq}
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(source), "exec"), native)
    cfg = yaml.safe_load((lens / "review_scoring_defaults.expression-aware.yaml").read_text())
    with gzip.open(root / "reconstruction/patient-windows.jsonl.gz", "rt") as handle:
        windows = [json.loads(line) for line in handle]
    with gzip.open(root / "reconstruction/native-event-evidence.json.gz", "rt") as handle:
        events = json.load(handle)
    targets, coordinates = defaultdict(dict), {}
    for window in windows:
        if window["native_event_review_tier"] != "PASS":
            continue
        event = window["event_id"]
        coordinates[event] = window["coordinates"]
        for target in target_peptides(window):
            coding = window["coding_sequence"][3 * target["start"]:3 * target["end"]]
            targets[event][coding] = target["peptide"]
    rna = pysam.AlignmentFile(policy["rna_bam"], "rb")
    results = {}
    for i, (event, coding_targets) in enumerate(sorted(targets.items())):
        match = re.fullmatch(r"([^:]+):(\d+)-(\d+)\(([+-])\)", coordinates[event])
        chrom, low, high = match[1], int(match[2]), int(match[3])
        bridges = defaultdict(list)
        for read in rna.fetch(chrom, max(0, low - 201), low + 200):
            # These are deliberately the native review criteria. The separate
            # patient-window evidence uses stricter primary/NH/BQ criteria.
            if read.is_unmapped or read.is_secondary or read.mapping_quality < cfg["extract"]["min_mapq"]:
                continue
            for junction in native["_extract_junctions_from_read"](read):
                if (junction["start"], junction["end"]) == (low + 1, high - 1):
                    seq = native["_junction_spanning_read_seq"](read, junction)
                    if seq:
                        bridges[read.query_name].append(seq)
        bridge_sequences = [sorted(values, key=len, reverse=True)[0] for values in bridges.values()]
        ev = events[event]
        assert len(bridge_sequences) == ev["tumor_rna"]["support_reads"], (event, len(bridge_sequences), ev["tumor_rna"]["support_reads"])
        for coding, peptide in sorted(coding_targets.items()):
            support = native["_splice_read_coding_support"](bridge_sequences, coding, k=9)
            evidence = dict(ev)
            evidence["event_specific"] = dict(ev["event_specific"], peptide_read_support=support)
            flags, hard = native["_flags_from_evidence"](dict(source_type="SPLICE", peptide=peptide), evidence, cfg)
            score = max(0, min(100, 100 - sum(cfg["scoring"]["penalties"].get(f, 0) for f in flags)))
            thresholds = cfg["scoring"]["tier_thresholds"]
            tier = "FAIL" if hard or score < thresholds["review_min"] else "REVIEW" if score < thresholds["pass_min"] else "PASS"
            results[key(event, coding)] = dict(event_id=event, peptide=peptide, coding_sequence=coding,
                  native_review_tier=tier, native_review_score=score, native_review_flags=flags,
                  coding_read_support=support)
        if i % 50 == 0:
            print("REVIEWED", i + 1, "/", len(targets), "events", len(results), "coding targets", flush=True)
    out = root / "review"
    out.mkdir(exist_ok=True)
    with gzip.open(out / "peptide-review.json.gz", "wt") as handle:
        json.dump(results, handle)
    receipt = dict(status="COMPLETED", completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        native_review_source_sha256=digest(source), worker_sha256=digest(Path(__file__)),
        reconstruction_sha256=digest(root / "reconstruction/patient-windows.jsonl.gz"),
        native_event_evidence_sha256=digest(root / "reconstruction/native-event-evidence.json.gz"),
        policy_sha256=digest(root / "policy.json"), review_events=len(targets),
        coding_targets=len(results), PASS_coding_targets=sum(r["native_review_tier"] == "PASS" for r in results.values()),
        output_sha256=digest(out / "peptide-review.json.gz"),
        replaced_check="Reference-only peptide consistency replaced by complete read-supported patient sequence and annotated-CDS anchor; other native flags unchanged")
    (out / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2), flush=True)


if __name__ == "__main__":
    main(Path(sys.argv[1]))
