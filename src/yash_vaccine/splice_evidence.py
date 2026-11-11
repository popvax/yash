"""Offline audit of all patient splice windows and the resulting gene ranking.

Independently recounts native affinity contributions from the saved fresh model
outputs. Does not start cloud jobs, run models or claim to reread patient BAMs.
"""
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

from Bio.Seq import Seq


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    opener = gzip.open if Path(path).suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return json.load(handle)


def rows(path):
    opener = gzip.open if Path(path).suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def close(left, right):
    return math.isclose(float(left), float(right), rel_tol=1e-8, abs_tol=1e-8)


def verify(results, output):
    policy = load(results / "policy.json")
    receipt = load(results / "reconstruction/receipt.json")
    assert receipt["status"] == "COMPLETED" and receipt["complete_disjoint_event_partition"]
    assert receipt["policy"] == policy
    for name, sha in receipt["output_sha256"].items():
        assert digest(results / "reconstruction" / name) == sha, name
    def patient_windows():
        with gzip.open(results / "reconstruction/patient-windows.jsonl.gz", "rt") as handle:
            for line in handle:
                yield json.loads(line)
    window_ids = set()
    checked_prefixes = set()
    events = rows(results / "reconstruction/input-events.tsv")
    audit = rows(results / "reconstruction/event-audit.tsv")
    assert len(events) == len({r["event_id"] for r in events}) == receipt["input_events"]
    assert {r["event_id"] for r in events} == {r["event_id"] for r in audit}
    assert all(r["cond"] == r["cond_add_tcga_control"] == "True" for r in events)
    assert dict(Counter(r["status"] for r in audit)) == receipt["event_counts"]
    for window in patient_windows():
        assert window["window_id"] not in window_ids
        window_ids.add(window["window_id"])
        assert len(window["coding_sequence"]) == 75 and len(window["sequence"]) == 25
        assert str(Seq(window["coding_sequence"]).translate()) == window["sequence"]
        names = set().union(*(set(p["fragment_names"]) for p in window["genomic_path_support"]))
        assert names == set(window["fragment_names"])
        assert len(names) == window["fragment_count"] >= policy["min_full_window_fragments"]
        assert window["genomic_path_count"] == len(window["genomic_path_support"])
        for path in window["genomic_path_support"]:
            assert len(path["genomic_positions_0based"]) == 75 and path["upstream_frame_anchors"]
            anchored_names = set()
            for anchor in path["upstream_frame_anchors"]:
                prefix = anchor["upstream_coding_sequence"]
                if prefix not in checked_prefixes:
                    assert len(prefix) % 3 == 0 and "*" not in str(Seq(prefix).translate())
                    checked_prefixes.add(prefix)
                assert len(anchor["upstream_coding_positions_0based"]) == len(prefix)
                assert (anchor["upstream_CDS_index_0based"] - anchor["upstream_query_offset"]
                        + anchor["first_complete_codon_query_offset"]) % 3 == 0
                anchored_names.update(anchor["fragment_names"])
            assert anchored_names == set(path["fragment_names"])
            assert all((anchor["upstream_CDS_index_0based"] - anchor["upstream_query_offset"]
                        + window["start_nt_relative_to_join"]) % 3 == 0 for anchor in path["upstream_frame_anchors"])
    assert len(window_ids) == receipt["patient_windows"]
    read_check = load(results / "read-verification/verification.json")
    assert read_check["status"] == "PASS" and not read_check["failures"]
    assert read_check["windows_examined"] == read_check["windows_independently_confirmed"] == len(window_ids)
    assert read_check["reconstruction_sha256"] == digest(results / "reconstruction/patient-windows.jsonl.gz")
    reviews = load(results / "review/peptide-review.json.gz")
    review_receipt = load(results / "review/receipt.json")
    assert review_receipt["status"] == "COMPLETED"
    assert review_receipt["output_sha256"] == digest(results / "review/peptide-review.json.gz")
    assert review_receipt["reconstruction_sha256"] == read_check["reconstruction_sha256"]
    inference = load(results / "scoring/receipt.json")
    assert inference["status"] == "COMPLETED" and inference["model_files_verified"]
    assert inference["reconstruction_sha256"] == read_check["reconstruction_sha256"]
    assert inference["patient_peptide_review_sha256"] == review_receipt["output_sha256"]
    assert inference["policy_sha256"] == digest(results / "policy.json")
    for name, sha in inference["output_sha256"].items():
        assert digest(results / "scoring" / name) == sha, name
    affinities, peptide_alleles = {}, defaultdict(set)
    for row in rows(results / "scoring/affinities.tsv.gz"):
        key = row["peptide"], row["allele"]
        assert key not in affinities
        affinities[key] = float(row["affinity_nm"])
        peptide_alleles[key[0]].add(key[1])
    assert len(affinities) == inference["prediction_rows"] == len(peptide_alleles) * len(policy["hla_alleles"])
    assert all(alleles == set(policy["hla_alleles"]) for alleles in peptide_alleles.values())
    matches = load(results / "scoring/reference-matches.json.gz")
    assert set(matches) == set(peptide_alleles)
    assert all(m["source_provenance_complete"] for m in matches.values())
    reference_check = load(results / "reference-check/receipt.json")
    assert reference_check["status"] == "PASS" and reference_check["all_reference_proteins_checked_by_two_exact_matchers"]
    assert reference_check["current_annotation_associations_verified"] and reference_check["current_protein_sequences_verified"]
    assert reference_check["reconstruction_sha256"] == read_check["reconstruction_sha256"]
    assert reference_check["output_sha256"] == digest(results / "reference-check/matches.json.gz")
    assert matches == load(results / "reference-check/matches.json.gz")
    score_config = load(results / "scoring/native-score-config.json")["epitope_config"]
    assert score_config["scoring_mode"] == "affinity" and score_config["filter_expr"] is None and score_config["score_expr"] is None
    midpoint = score_config["logistic_epitope_score_midpoint"]
    width = score_config["logistic_epitope_score_width"]
    cutoff = score_config["binding_affinity_cutoff"]
    minimum = score_config["min_epitope_score"]
    def contribution(nm):
        value = 0.0 if nm >= cutoff else (1 + math.exp(-midpoint / width)) / (1 + math.exp((nm - midpoint) / width))
        return value if value >= minimum else 0.0
    scored = {r["window_id"]: r for r in rows(results / "scoring/all-windows.tsv.gz")}
    assert set(scored) == window_ids
    eligible, checked, queried = [], 0, set()
    for window in patient_windows():
        wid = window["window_id"]
        row = scored[wid]
        assert row["sequence"] == window["sequence"] and int(row["fragment_count"]) == window["fragment_count"]
        if window["native_event_review_tier"] != "PASS":
            assert row["primary_window_eligible"] == "False"
            continue
        score, best, admitted = 0.0, None, False
        for length in (8, 9, 10, 11):
            for start in range(25 - length + 1):
                # Independent nucleotide-level crossing rule, rather than the
                # production amino-acid target-interval helper.
                first = window["start_nt_relative_to_join"] + 3 * start
                last = first + 3 * length
                if not first < 0 < last:
                    continue
                peptide = window["sequence"][start:start + length]
                coding = window["coding_sequence"][3 * start:3 * (start + length)]
                queried.add(peptide)
                if matches[peptide]["occurs"]:
                    continue
                values = [affinities[peptide, allele] for allele in policy["hla_alleles"]]
                score += sum(contribution(value) for value in values)
                best = min(values + ([best] if best is not None else []))
                review = reviews[hashlib.sha256((window["event_id"] + "|" + coding).encode()).hexdigest()]
                assert review["peptide"] == peptide
                admitted |= min(values) <= policy["admission_max_affinity_nm"] and review["native_review_tier"] == "PASS"
        assert close(score, row["native_affinity_sum"]), (wid, score, row["native_affinity_sum"])
        assert best is None and not row["best_target_affinity_nm"] or best == float(row["best_target_affinity_nm"])
        expression = window["supplied_transformed_gene_tpm"]
        if expression is not None and expression > 0:
            assert close(score * math.sqrt(expression), row["weighted_score_left"])
        else:
            assert row["weighted_score_left"] == ""
        assert (row["primary_window_eligible"] == "True") == (admitted and score > 0)
        if admitted and score > 0:
            eligible.append(row)
        checked += 1
    assert queried == set(peptide_alleles)
    def rank_key(row):
        return (row["weighted_score_left"] == "", -float(row["weighted_score_left"] or 0),
                -float(row["native_affinity_sum"]), row["event_id"], row["window_id"])
    event_best, gene_best = {}, {}
    for row in sorted(eligible, key=rank_key):
        event_best.setdefault(row["event_id"], row)
    for row in sorted(event_best.values(), key=rank_key):
        gene_best.setdefault(row["gene_id"], row)
    genes = sorted(gene_best.values(), key=rank_key)
    selection = [r for r in genes if r["weighted_score_left"] != ""][:5]
    assert [r["window_id"] for r in rows(results / "scoring/gene-ranking.tsv")] == [r["window_id"] for r in genes]
    assert [r["window_id"] for r in rows(results / "scoring/selected.tsv")] == [r["window_id"] for r in selection]
    report = dict(status="VERIFIED", input_events=len(events), patient_windows=len(window_ids),
         fresh_affinity_rows=len(affinities), scored_windows_independently_recounted=checked,
         eligible_windows=len(eligible), eligible_events=len(event_best), eligible_genes=len(genes),
         selected_genes=[r["gene"] for r in selection],
         upstream_observed_stop_codons_excluded=True,
         verification_source_sha256=digest(Path(__file__)),
         raw_read_verification_receipt_checked=True, all_artifact_hashes_checked=True,
         scope="Offline replay from recorded read evidence and fresh inference outputs; no BAM rereading or new inference")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n")
    return report
