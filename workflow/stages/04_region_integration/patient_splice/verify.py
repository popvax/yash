"""Independently recheck every recorded window against BAM query bases.

Uses pysam's query/reference-position arrays rather than the reconstruction
worker's CIGAR path builder. Checks named supporting fragments, including
overlapping mates, without using the genomic reference to supply any base.
"""
from collections import defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import sys

from Bio.Seq import Seq
import pysam


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def eligible(read):
    return (not any((read.is_unmapped, read.is_secondary, read.is_supplementary,
                     read.is_duplicate, read.is_qcfail)) and read.mapping_quality >= 20
            and (not read.has_tag("NH") or read.get_tag("NH") == 1))


def exact_event(read, low, high):
    if high == low + 1:
        positions = read.get_reference_positions(full_length=True)
        return any(a == low - 1 and b == high - 1 for a, b in zip(positions, positions[1:]))
    reference = read.reference_start
    for operation, length in read.cigartuples or ():
        if operation == 3 and reference == low and reference + length == high - 1:
            return True
        if operation in (0, 2, 3, 7, 8):
            reference += length
    return False


def fragment_matches(reads, window, low, high, strand):
    expected = window["genomic_positions_0based"]
    indexes = {p: i for i, p in enumerate(expected) if p is not None}
    observed = defaultdict(set)
    if not any(exact_event(r, low, high) for r in reads):
        return False
    for read in reads:
        first, last = read.query_alignment_start, read.query_alignment_end
        positions = read.get_reference_positions(full_length=True)[first:last]
        sequence = read.query_sequence[first:last]
        qualities = list(read.query_qualities[first:last]) if read.query_qualities is not None else []
        if strand == "-":
            positions = list(reversed(positions))
            sequence = str(Seq(sequence).reverse_complement())
            qualities.reverse()
        offsets = {indexes[p] - q for q, p in enumerate(positions) if p in indexes}
        valid = []
        for offset in offsets:
            pairs = [(q, q + offset) for q in range(len(positions)) if 0 <= q + offset < len(expected)]
            if pairs and all(positions[q] == expected[i] for q, i in pairs):
                valid.append(pairs)
        if len(valid) != 1:
            continue
        for q, i in valid[0]:
            if qualities and qualities[q] >= 20:
                observed[i].add(sequence[q])
    return all(observed[i] == {base} for i, base in enumerate(window["coding_sequence"]))


def main(root):
    policy = json.loads((root / "policy.json").read_text())
    with gzip.open(root / "reconstruction/patient-windows.jsonl.gz", "rt") as handle:
        windows = [json.loads(line) for line in handle]
    grouped = defaultdict(list)
    for window in windows:
        grouped[window["event_id"]].append(window)
    rna = pysam.AlignmentFile(policy["rna_bam"], "rb")
    failures, verified, fragments = [], 0, 0
    for event, group in sorted(grouped.items()):
        anchor_match_cache = {}
        match = re.fullmatch(r"([^:]+):(\d+)-(\d+)\(([+-])\)", group[0]["coordinates"])
        chrom, low, high, strand = match[1], int(match[2]), int(match[3]), match[4]
        reads, seen = defaultdict(list), set()
        for a, b in [(max(0, low - 201), low + 200), (max(0, high - 201), high + 200)]:
            for read in rna.fetch(chrom, a, b):
                key = (read.query_name, read.flag, read.reference_start, read.cigarstring)
                if eligible(read) and key not in seen:
                    reads[read.query_name].append(read)
                    seen.add(key)
        for window in group:
            assert str(Seq(window["coding_sequence"]).translate()) == window["sequence"]
            assert len(window["coding_sequence"]) == 75 and len(window["sequence"]) == 25
            assert window["fragment_count"] == len(set(window["fragment_names"])) >= policy["min_full_window_fragments"]
            assert all((a["upstream_CDS_index_0based"] - a["upstream_query_offset"] + window["start_nt_relative_to_join"]) % 3 == 0
                       for a in window["upstream_frame_anchors"])
            confirmed, failed_witnesses, recorded_names = set(), [], set()
            for witness in window["genomic_path_support"]:
                local = dict(window, genomic_positions_0based=witness["genomic_positions_0based"])
                recorded_names.update(witness["fragment_names"])
                for name in witness["fragment_names"]:
                    anchors = [a for a in witness["upstream_frame_anchors"] if name in a["fragment_names"]]
                    frame_confirmed = False
                    for anchor in anchors:
                        prefix = anchor["upstream_coding_sequence"]
                        assert len(prefix) % 3 == 0 and "*" not in str(Seq(prefix).translate())
                        anchor_window = dict(coding_sequence=prefix,
                            genomic_positions_0based=anchor["upstream_coding_positions_0based"])
                        key = (name, prefix, tuple(anchor["upstream_coding_positions_0based"]))
                        if key not in anchor_match_cache:
                            anchor_match_cache[key] = fragment_matches(reads[name], anchor_window, low, high, strand)
                        if anchor_match_cache[key]:
                            frame_confirmed = True
                            break
                    if frame_confirmed and fragment_matches(reads[name], local, low, high, strand):
                        confirmed.add(name)
                    else:
                        failed_witnesses.append(name)
            assert recorded_names == set(window["fragment_names"])
            if len(confirmed) != window["fragment_count"] or failed_witnesses:
                failures.append(dict(window_id=window["window_id"], event_id=event,
                                     recorded_fragments=window["fragment_count"], verified_fragments=len(confirmed),
                                     failed_path_witnesses=len(failed_witnesses)))
            else:
                verified += 1
                fragments += len(confirmed)
        if verified % 500 < len(group):
            print("VERIFIED", verified, "/", len(windows), "windows", flush=True)
    result = dict(status="PASS" if not failures else "FAIL", windows_examined=len(windows),
                  windows_independently_confirmed=verified, named_fragment_window_observations_confirmed=fragments,
                  events_examined=len(grouped), failures=failures,
                  reconstruction_sha256=digest(root / "reconstruction/patient-windows.jsonl.gz"),
                  worker_sha256=digest(Path(__file__)), policy_sha256=digest(root / "policy.json"))
    out = root / "read-verification"
    out.mkdir(exist_ok=True)
    (out / "verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "failures"}, indent=2), flush=True)
    assert not failures, failures[:5]


if __name__ == "__main__":
    main(Path(sys.argv[1]))
