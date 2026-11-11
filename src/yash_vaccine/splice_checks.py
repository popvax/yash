"""Coordinate-normal and limited strand gates over existing splice evidence.

This module does not change sequences, fragment-name requirements, native QC,
reference matching, affinity predictions or expression scoring. Frozen inputs
remain immutable. The recorded measurements support a fully offline replay.
"""
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import re

from .splice_evidence import verify as verify_source


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load(path):
    opener = gzip.open if Path(path).suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return json.load(handle)


def rows(path):
    opener = gzip.open if Path(path).suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def encoded_table(records, fields=None):
    fields = fields or list(dict.fromkeys(k for record in records for k in record))
    assert fields, "Even an empty result must retain its column schema"
    handle = io.StringIO(newline="")
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(records)
    return handle.getvalue().encode()


def encoded_json(value):
    return (json.dumps(value, indent=2) + "\n").encode()


def coordinate(value):
    match = re.fullmatch(r"([^:]+):(\d+)-(\d+)\(([+-])\)", value)
    if not match:
        return None
    chrom, a, b, strand = match.groups()
    return (chrom.removeprefix("chr"), int(a), int(b)), strand


def measurement_inputs(source, checks):
    source, checks = Path(source), Path(checks)
    policy, request = load(checks / "checks-policy.json"), load(checks / "request.json")
    assert policy["schema_version"] == 1
    assert policy["selection"]["added_distinct_pattern_filter"] is False
    assert policy["selection"]["added_substitution_filter"] is False
    assert policy["strand"]["motif_disagreement"].startswith("REPORT_ONLY")
    assert request["source_reconstruction_sha256"] == digest(source / "reconstruction/patient-windows.jsonl.gz")
    assert request["source_events_sha256"] == digest(source / "reconstruction/input-events.tsv")
    assert request["source_policy_sha256"] == digest(source / "policy.json")
    reconstruction_inputs = load(source / "reconstruction/receipt.json")["input_sha256"]
    for mode in ("normal", "strand"):
        receipt = load(checks / (mode + "-receipt.json"))
        assert receipt["status"] == "COMPLETED"
        assert receipt["request_sha256"] == digest(checks / "request.json")
        assert receipt["checks_policy_sha256"] == digest(checks / "checks-policy.json")
        assert receipt["source_reconstruction_sha256"] == request["source_reconstruction_sha256"]
        assert receipt["source_events_sha256"] == request["source_events_sha256"]
        for path, sha in receipt["input_sha256"].items():
            if path.endswith("/policy.json"):
                assert sha == request["source_policy_sha256"]
            elif path in reconstruction_inputs:
                assert sha == reconstruction_inputs[path], "Measurement used a different reconstruction input: " + path
        for name, sha in receipt["output_sha256"].items():
            assert digest(checks / name) == sha, name
        state = checks / (mode + "-state.json")
        if state.is_file():
            uploaded = load(state)["uploaded_sha256"]
            assert receipt["worker_sha256"] == uploaded["check_measurements.py"]
    normal, strand = load(checks / "normal-counts.json.gz"), load(checks / "strand-counts.json.gz")
    events = rows(source / "reconstruction/input-events.tsv")
    assert normal["source_events"] == len(events)
    assert set(normal["cohorts"]) == {c["name"] for c in policy["normal"]["cohorts"]}
    return policy, events, normal, strand


def normal_decisions(events, normal, policy):
    decisions, tissue_rows = {}, []
    for event in events:
        event_id, parsed = event["event_id"], coordinate(event["coordinates"])
        reasons, match_count, positive_opposite, positive_same, worst_tissue_mean = [], 0, 0, 0, None
        if not parsed:
            reasons.append("NORMAL_COORDINATES_UNRESOLVED")
        for name, cohort in normal["cohorts"].items():
            assert len(cohort["tissues"]) == len(cohort["sample_names"]) == cohort["samples"]
            aliases = cohort["event_matches"].get(event_id, [])
            assert len(aliases) == len(set(aliases))
            groups = {"same_strand": [], "opposite_strand": []}
            for uid in aliases:
                match = cohort["rows"][uid]
                resolved = coordinate(match["coordinates"])
                assert parsed and resolved and resolved[0] == parsed[0]
                assert resolved[1] == match["strand"]
                group = "same_strand" if resolved[1] == parsed[1] else "opposite_strand"
                groups[group].append(uid)
            match_count += len(aliases)
            for group, uids in groups.items():
                if not uids:
                    continue
                # Alternate UID annotations must not be added together.
                counts = [0.0] * cohort["samples"]
                for uid in uids:
                    for column, count in cohort["rows"][uid]["nonzero_counts"]:
                        assert 0 <= column < len(counts) and count > 0
                        counts[column] = max(counts[column], count)
                by_tissue = defaultdict(list)
                for tissue, count in zip(cohort["tissues"], counts, strict=True):
                    by_tissue[tissue].append(count)
                for tissue, values in sorted(by_tissue.items()):
                    mean = sum(values) / len(values)
                    positive = sum(v > 0 for v in values)
                    row = dict(event_id=event_id, gene_id=event["gene_id"], coordinates=event["coordinates"],
                        cohort=name, strand_relation=group, tissue=tissue, sample_count=len(values),
                        positive_samples=positive, mean_count=mean, max_count=max(values),
                        matched_UID_count=len(uids), control_UIDs=";".join(sorted(uids)))
                    tissue_rows.append(row)
                    if group == "opposite_strand":
                        positive_opposite += positive
                    else:
                        positive_same += positive
                        worst_tissue_mean = max(worst_tissue_mean or 0.0, mean)
                        if mean >= policy["normal"]["same_strand_max_tissue_mean_exclusive"]:
                            reasons.append("NORMAL_TISSUE_MEAN_AT_OR_ABOVE_THRESHOLD:" + name + ":" + tissue)
                        if positive and tissue.strip().lower() in {"unknown", "none", "nan", "", "na"}:
                            reasons.append("NORMAL_TISSUE_METADATA_UNKNOWN:" + name)
        if positive_opposite:
            reasons.append("POSITIVE_OPPOSITE_STRAND_NORMAL_MATCH_REQUIRES_REVIEW")
        status = "REVIEW" if reasons else "PASS"
        detection = ("COORDINATES_UNRESOLVED" if not parsed else "COORDINATE_MATCH_IN_CATALOG"
                     if match_count else "NOT_DETECTED_IN_SEARCHED_REFERENCE")
        decisions[event_id] = dict(normal_status=status, normal_detection=detection,
            normal_reasons=";".join(sorted(set(reasons))), normal_matched_UIDs=match_count,
            normal_positive_same_sample_observations=positive_same,
            normal_positive_opposite_sample_observations=positive_opposite,
            highest_same_strand_tissue_mean=worst_tissue_mean)
    return decisions, tissue_rows


def library_decision(strand, policy):
    lib, rule = strand["library"], policy["strand"]
    counts = lib["counts"]
    total = sum(counts.get(k, 0) for k in ("sense", "antisense"))
    assert total == lib["fragments"]
    dominant = max(("sense", "antisense"), key=lambda k: counts.get(k, 0))
    fraction = counts.get(dominant, 0) / total if total else 0.0
    established = (total >= rule["library_min_fragments"] and lib["genes"] >= rule["library_min_genes"]
                   and fraction >= rule["library_min_dominant_fraction"])
    return dict(status="PASS" if established else "REVIEW", read1_orientation=dominant if established else "unresolved",
        dominant_fraction=fraction, fragments=total, genes=lib["genes"])


def strand_decision(window_measurements, library):
    observations = window_measurements["fragment_orientations"]
    counts = Counter(observations.values())
    assert dict(counts) == window_measurements["counts"]
    reasons = []
    if library["status"] != "PASS":
        reasons.append("LIBRARY_ORIENTATION_UNRESOLVED")
    expected = library["read1_orientation"]
    opposite = "antisense" if expected == "sense" else "sense"
    agreeing, opposing = counts.get(expected, 0), counts.get(opposite, 0)
    unresolved = sum(n for label, n in counts.items() if label not in ("sense", "antisense"))
    if unresolved:
        reasons.append("FULL_WINDOW_WITNESS_ORIENTATION_UNRESOLVED")
    if agreeing <= opposing:
        reasons.append("NO_MAJORITY_GENE_DIRECTION_SUPPORT")
    return dict(strand_status="REVIEW" if reasons else "PASS", strand_reasons=";".join(reasons),
        strand_agreeing_fragments=agreeing, strand_opposing_fragments=opposing,
        strand_unresolved_fragments=unresolved)


def decisions(source, checks):
    policy, events, normal, strand = measurement_inputs(source, checks)
    normal_checks, tissue_rows = normal_decisions(events, normal, policy)
    # Measurement receipts retain the policy used when the counts were collected.
    # A separately hashed admission policy can tighten selection without claiming
    # that the original measurements were rerun under a different policy.
    selection_path = Path(checks) / "selection-policy.json"
    if selection_path.is_file():
        selection = load(selection_path)
        assert selection["schema_version"] == 1
        required = selection["normal_detection_required"]
        assert required == "NOT_DETECTED_IN_SEARCHED_REFERENCE"
        policy["selection"]["normal_detection_required"] = required
        for check in normal_checks.values():
            if check["normal_detection"] != required:
                reasons = set(filter(None, check["normal_reasons"].split(";")))
                reasons.add("NORMAL_REFERENCE_NOT_ND")
                check["normal_status"] = "REVIEW"
                check["normal_reasons"] = ";".join(sorted(reasons))
    library = library_decision(strand, policy)
    window_checks, event_checks, per_event = {}, [], defaultdict(Counter)
    with gzip.open(Path(source) / "reconstruction/patient-windows.jsonl.gz", "rt") as handle:
        for line in handle:
            window = json.loads(line)
            wid, event_id = window["window_id"], window["event_id"]
            observed = strand["windows"][wid]
            assert observed["event_id"] == event_id
            assert set(observed["fragment_orientations"]) == set(window["fragment_names"])
            orientation = strand_decision(observed, library)
            event_orientation = strand["events"][event_id]
            assert event_orientation["motif_and_XS_are_report_only"] is True
            status = "PASS" if normal_checks[event_id]["normal_status"] == orientation["strand_status"] == "PASS" else "REVIEW"
            check = dict(window_id=wid, event_id=event_id, gene=window["gene"], additional_status=status,
                **normal_checks[event_id], **orientation,
                genomic_plus_motif=event_orientation["genomic_plus_motif"],
                motif_agrees_with_gene=event_orientation["motif_agrees_with_gene"])
            window_checks[wid] = check
            per_event[event_id][status] += 1
    assert set(window_checks) == set(strand["windows"])
    assert set(strand["events"]) == set(per_event)
    for event in events:
        event_id = event["event_id"]
        measured = strand["events"].get(event_id, {})
        event_checks.append(dict(event_id=event_id, gene_id=event["gene_id"], coordinates=event["coordinates"],
            **normal_checks[event_id], windows_PASS=per_event[event_id]["PASS"],
            windows_REVIEW=per_event[event_id]["REVIEW"],
            genomic_plus_motif=measured.get("genomic_plus_motif"),
            motif_agrees_with_gene=measured.get("motif_agrees_with_gene"),
            junction_fragment_orientations=json.dumps(measured.get("junction_fragment_orientations", {}), sort_keys=True)))
    return policy, window_checks, event_checks, tissue_rows, library


def selection_products(source, checks):
    policy, checks_by_window, event_checks, tissue_rows, library = decisions(source, checks)
    original = rows(Path(source) / "scoring/all-windows.tsv.gz")
    assert {r["window_id"] for r in original} == set(checks_by_window)
    updated = []
    old_eligible, new_eligible = [], []
    for row in original:
        check = checks_by_window[row["window_id"]]
        assert row["event_id"] == check["event_id"]
        before = row["primary_window_eligible"] == "True"
        if before:
            old_eligible.append(row)
        after = before and check["additional_status"] == "PASS"
        result = dict(row, base_primary_window_eligible=row["primary_window_eligible"], **{
            k: v for k, v in check.items() if k not in {"window_id", "event_id", "gene"}})
        result["primary_window_eligible"] = str(after)
        if before and not after:
            result["rejection_reason"] = ";".join(x for x in (check["normal_reasons"], check["strand_reasons"]) if x)
        updated.append(result)
        if after:
            new_eligible.append(result)
    def rank_key(row):
        return (row["weighted_score_left"] == "", -float(row["weighted_score_left"] or 0),
                -float(row["native_affinity_sum"]), row["event_id"], row["window_id"])
    event_best = {}
    for row in sorted(new_eligible, key=rank_key):
        event_best.setdefault(row["event_id"], row)
    event_ranking = [dict(event_rank=i, **r) for i, r in enumerate(sorted(event_best.values(), key=rank_key), 1)]
    gene_best = {}
    for row in sorted(event_best.values(), key=rank_key):
        gene_best.setdefault(row["gene_id"], row)
    gene_ranking = [dict(gene_rank=i, **r) for i, r in enumerate(sorted(gene_best.values(), key=rank_key), 1)]
    selected = [r for r in gene_ranking if r["weighted_score_left"] != ""][:policy["selection"]["splice_gene_quota"]]
    previous_selected = rows(Path(source) / "scoring/selected.tsv")
    comparison = []
    for row in previous_selected:
        check = checks_by_window[row["window_id"]]
        new_gene = gene_best.get(row["gene_id"])
        comparison.append(dict(gene=row["gene"], original_window_id=row["window_id"],
            original_sequence=row["sequence"], **{k: v for k, v in check.items() if k not in {"gene", "window_id"}},
            best_surviving_gene_window=new_gene["window_id"] if new_gene else "",
            best_surviving_gene_sequence=new_gene["sequence"] if new_gene else "",
            original_window_selected_after_checks=any(r["window_id"] == row["window_id"] for r in selected),
            gene_selected_after_checks=any(r["gene_id"] == row["gene_id"] for r in selected)))
    summary = dict(source_input_events=len(event_checks), source_windows=len(original),
        source_primary_eligible_windows=len(old_eligible), primary_eligible_windows=len(new_eligible),
        primary_eligible_events=len(event_ranking), primary_eligible_genes=len(gene_ranking),
        positive_expression_genes=sum(r["weighted_score_left"] != "" for r in gene_ranking),
        excluded_previously_eligible_windows=dict(Counter(
            ("NORMAL_AND_STRAND" if r["normal_status"] != "PASS" and r["strand_status"] != "PASS" else
             "NORMAL" if r["normal_status"] != "PASS" else "STRAND")
            for r in updated if r["base_primary_window_eligible"] == "True" and r["primary_window_eligible"] == "False")),
        normal_event_status_counts=dict(Counter(r["normal_status"] for r in event_checks)),
        normal_detection_counts=dict(Counter(r["normal_detection"] for r in event_checks)),
        library_orientation=library, selected_genes=[r["gene"] for r in selected],
        selected_regions=len(selected), new_inference_predictions=0,
        distinct_pattern_filter_added=False, substitution_filter_added=False,
        motifs_are_report_only=True, source_sequences_and_scores_unchanged=True,
        scope="Existing normal-filtered 4,929-event pool only; no recovery of events excluded upstream, no construct update")
    if "normal_detection_required" in policy["selection"]:
        summary["normal_detection_required"] = policy["selection"]["normal_detection_required"]
    products = {
        "event-checks.tsv": encoded_table(event_checks),
        "normal-tissues.tsv.gz": gzip.compress(encoded_table(tissue_rows, [
            "event_id", "gene_id", "coordinates", "cohort", "strand_relation", "tissue", "sample_count",
            "positive_samples", "mean_count", "max_count", "matched_UID_count", "control_UIDs"]), mtime=0),
        "window-checks.tsv.gz": gzip.compress(encoded_table(list(checks_by_window.values())), mtime=0),
        "all-windows.tsv.gz": gzip.compress(encoded_table(updated), mtime=0),
        "event-ranking.tsv": encoded_table(event_ranking, ["event_rank"] + list(updated[0])),
        "gene-ranking.tsv": encoded_table(gene_ranking, ["gene_rank"] + list(updated[0])),
        "selected.tsv": encoded_table(selected, ["gene_rank"] + list(updated[0])),
        "previous-selection-checks.tsv": encoded_table(comparison),
        "selected.fasta": "".join(f">{r['window_id']} gene={r['gene']} event={r['event_id']}\n{r['sequence']}\n" for r in selected).encode(),
        "summary.json": encoded_json(summary)}
    return products, summary


def input_hashes(source, checks):
    result = {"source/" + name: digest(Path(source) / name) for name in (
        "policy.json", "reconstruction/patient-windows.jsonl.gz", "reconstruction/input-events.tsv",
        "scoring/all-windows.tsv.gz", "scoring/selected.tsv", "scoring/affinities.tsv.gz")} | {
        "checks/" + name: digest(Path(checks) / name) for name in (
            "checks-policy.json", "request.json", "normal-counts.json.gz", "strand-counts.json.gz",
            "normal-receipt.json", "strand-receipt.json")}
    selection_path = Path(checks) / "selection-policy.json"
    if selection_path.is_file():
        result["checks/selection-policy.json"] = digest(selection_path)
    return result


def apply(source, output):
    source, output = Path(source), Path(output)
    assert source.resolve() != output.resolve()
    # The original numerical/HLA/reference/read-verification contracts must pass
    # before saved scores are reused. This writes only to the new directory.
    baseline = verify_source(source, output / "source-verification.json")
    assert baseline["status"] == "VERIFIED"
    before = input_hashes(source, output)
    products, summary = selection_products(source, output)
    for name, content in products.items():
        (output / name).write_bytes(content)
    assert before == input_hashes(source, output)
    record = dict(status="PASS", mode="Offline selection using recorded normal/strand measurements and cached predictions",
        module_sha256=digest(__file__), input_sha256=before,
        output_sha256={name: digest(output / name) for name in products},
        source_verification_sha256=digest(output / "source-verification.json"), summary=summary)
    (output / "verification.json").write_bytes(encoded_json(record))
    return record


def verify(source, output):
    source, output = Path(source), Path(output)
    recorded = load(output / "verification.json")
    assert recorded["status"] == "PASS" and recorded["input_sha256"] == input_hashes(source, output)
    assert recorded["source_verification_sha256"] == digest(output / "source-verification.json")
    assert load(output / "source-verification.json")["status"] == "VERIFIED"
    products, summary = selection_products(source, output)
    for name, content in products.items():
        actual = (output / name).read_bytes()
        assert actual == content, "Offline replay mismatch: " + name
        assert digest(output / name) == recorded["output_sha256"][name]
    assert summary == recorded["summary"]
    return dict(status="PASS", offline_replayed_files=len(products), every_window_checked=summary["source_windows"],
        every_input_event_checked=summary["source_input_events"], unchanged_source_inputs=True,
        new_model_inference=False, **{k: summary[k] for k in ("selected_genes", "selected_regions")})
