"""Verify the collected matched-normal result bundle without network or new jobs."""
from yash_vaccine.paths import resolve_input
from pathlib import Path
from collections import Counter
import argparse
import csv
import datetime
import hashlib
import importlib.util
import json

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
LOCAL_RECORDS = json.loads((resolve_input(HERE / 'execution/local-paths.json')).read_text())
PATHS = LOCAL_RECORDS['paths']


def resolve(name):
    return resolve_input(HERE / PATHS.get(name, name))


def read(name):
    return json.loads(resolve_input(resolve(name)).read_text())


def digest(path):
    return hashlib.sha256(resolve_input(path).read_bytes()).hexdigest()


def rows(name):
    with resolve_input(resolve(name)).open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Optionally save this check; verification is read-only by default")
    args = parser.parse_args()
    artifacts = read("execution/artifacts.json")
    for record in artifacts:
        if record["local"]:
            path = resolve(record["local"])
            assert resolve_input(path).stat().st_size == record["bytes"] and digest(path) == record["sha256"], record["local"]
    for name, origin in LOCAL_RECORDS['export_sources'].items():
        receipt = read(origin['receipt'])
        if origin['stage'] is not None:receipt = receipt[origin['stage']]
        expected = receipt['outputs'][origin['output']]
        path = resolve(name)
        assert resolve_input(path).stat().st_size == expected['bytes'] and digest(path) == expected['sha256'], name
    receipts = read("execution/stage-receipts.json")
    assert all(r["status"] == "COMPLETED" for r in receipts.values())
    execution = read("execution/native-execution.json")
    assert len(execution["workers"]) == 552
    assert all(r.get("observed_exit_code") == 0 for r in execution["workers"])
    frozen = read("execution/frozen-inputs.json")
    from yash_vaccine.provenance import RecordedSource
    source = RecordedSource(REPO)
    source.verify_files(frozen['files'])

    audit = read("review/final/current-17-final-audit.json")
    outcomes = Counter(r["final_reassessment"] for r in audit)
    assert len(audit) == 17 and outcomes == {"SUPPORTED_PAIRED_SOMATIC": 10, "PRESENT_IN_NORMAL": 7}
    assert all(r["new_normal_evidence"]["adequate_qualifying_normal_fragments"] for r in audit)
    assert receipts["identity"]["identity"] == "MATCH"
    assert receipts["annotate"]["primary_DNA_eligible"] == 99
    assert receipts["rna"]["input_variants"] == 101 and receipts["rna"]["all_inputs_accounted_for"]
    assert receipts["rank"]["primary_ranked_regions"] == 12
    assert receipts["rank"]["new_peptide_HLA_predictions"] == 3624

    selected = read("selection/selected.json")
    sources = Counter(r["source"] for r in selected)
    assert sources == {"MUTATION": 12, "SPLICE": 5, "CTA/SELF": 5, "FUSION": 2}
    assert len({r["region_id"] for r in selected}) == len(selected) == 24
    assert all(len(r["sequence"]) == 25 and r["hydrophobicity"]["status"] == "PASS" for r in selected)
    mutations = [r for r in selected if r["source"] == "MUTATION"]
    assert all(r["primary_DNA_eligible"] and r["somatic_status"] == "SUPPORTED_PAIRED_SOMATIC" for r in mutations)
    old = rows("../initial/selection/selected.tsv")
    old_sequences = {r["sequence"] for r in old}; new_sequences = {r["sequence"] for r in selected}
    assert (len(old_sequences & new_sequences), len(old_sequences - new_sequences), len(new_sequences - old_sequences)) == (22, 7, 2)
    assert {(r["gene"], r["sequence"]) for r in old if r["source"] != "MUTATION"} == {
        (r["gene"], r["sequence"]) for r in selected if r["source"] != "MUTATION"}

    protein = "".join(line for line in (resolve_input(HERE / "construct/protein/cassette.fasta")).read_text().splitlines() if not line.startswith(">"))
    parts = rows("construct/protein/components.tsv")
    assert len(protein) == 724 and "".join(p["aa"] for p in parts) == protein
    assert {(p["region_id"], p["aa"]) for p in parts if p["kind"] == "antigen"} == {(r["region_id"], r["sequence"]) for r in selected}
    check = read("construct/final-check.json")
    assert check["status"] == "FRESH_FULL_JUNCTION_CHECK_PASSED" and check["all_antigens_exactly_once"]
    assert check["protein_sha256"] == digest(HERE / "construct/protein/cassette.fasta")
    assert check["components_sha256"] == digest(HERE / "construct/protein/components.tsv")
    assert check["hydrophobicity"]["unexpected_high_windows"] == 0
    assert read("construct/protein/construct.json")["solver"]["all_phases_proven_optimal"]

    nucleotide = read("construct/mrna/verification.json")
    for name, expected in nucleotide["input_sha256"].items():
        assert digest(HERE / "construct/mrna" / name) == expected, name
    for name, expected in nucleotide["output_sha256"].items():
        assert digest(HERE / "construct/mrna" / name) == expected, name
    assert nucleotide["all_constraints_pass"] and nucleotide["exact_parent_protein_translation"]
    assert not nucleotide["complete_construct_verified"] and nucleotide["polyA_status"] == "unspecified"
    assert not nucleotide["folding_heuristic_target_retained"] and nucleotide["final_paired_bracket_run_nt"] == 21

    # Apply the maintained strict GenBank verifier again to the exported record.
    from Bio import SeqIO
    spec = importlib.util.spec_from_file_location("matched_normal_genbank", HERE / "construct/mrna/genbank.py")
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    mrna = HERE / "construct/mrna"
    cfg = read("construct/mrna/backbone.json")
    sources_record, selected_by_id, architecture, _ = module.load_inputs(mrna)
    components = rows("construct/mrna/components.tsv")
    for index, part in enumerate(components):
        for key in ["rna_start_0based", "rna_end_exclusive", "template_start_0based", "template_end_exclusive", "length_nt"]:
            if part[key] != "":
                part[key] = int(part[key])
        if part["kind"] == "polyA_unspecified":
            # Reconstruct the encoder's internal caret at the known boundary,
            # solely for annotation verification. The public table correctly
            # leaves the unknown tail's coordinates/length blank.
            previous = components[index - 1]
            part.update(rna_start_0based=previous["rna_end_exclusive"],
                        rna_end_exclusive=previous["rna_end_exclusive"],
                        template_start_0based=previous["template_end_exclusive"],
                        template_end_exclusive=previous["template_end_exclusive"], length_nt=0)
    dna_segment = str(SeqIO.read(mrna / "template-segments.fasta", "fasta").seq)
    rna_segment = str(SeqIO.read(mrna / "mrna-segments.fasta", "fasta").seq)
    records = module.build_record(dna_segment, protein, components, cfg, sources_record, selected_by_id, architecture)
    genbank = module.verify_export(mrna / "template-segments.gb", records, components, protein, [rna_segment], selected_by_id)
    assert genbank == nucleotide["genbank_annotation"]
    assert genbank["feature_count"] == 64 and genbank["antigen_regions"] == 24
    assert read("execution/local-replay.json")["all_recorded_output_hashes_identical"]
    review = read("review/antigen-review.json")
    assert {r["region_id"] for r in review} == {r["region_id"] for r in selected}

    # Verify the original files directly so this also works in a one-commit
    # clone or exported tree. Old revision IDs remain provenance, not runtime
    # dependencies on discarded Git history.
    original = read("execution/original-scientific-files.json")
    assert original["verification_requires_git_history"] is False
    original_bytes_verified = source.verify_files(original['files'], original_ledger=True)
    requirements = [
        dict(requirement="Matched-normal alignment and calibrated record preservation", status="PASS", evidence=["execution/alignment/complete.json", "execution/bqsr/complete.json", "qc/normal-read-accounting.json"]),
        dict(requirement="Sample identity, contamination and scoped normal callability", status="PASS", evidence=["qc/fingerprint-verdict.json", "qc/normal-callability-summary.json"]),
        dict(requirement="All 17 original mutations reassessed independently of rediscovery", status="PASS", evidence=["review/final/current-17-final-audit.json"]),
        dict(requirement="Independent primary-nuclear paired calling and native/population filtering", status="PASS", evidence=["review/paired-native-command-audit.json", "review/paired-discovery-interval-audit.json", "execution/paired-filter/complete.json"]),
        dict(requirement="Fresh RNA, class-I affinity ranking and personal-normal self search", status="PASS", evidence=["exports/RNA-reconstruction-completion.json", "exports/class-I-ranking-completion.json", "review/normal-self-scope.json"]),
        dict(requirement="Integrated selection, joint assembly, junction/hydropathy and nucleotide checks", status="PASS_WITH_RECORDED_LIMITATIONS", evidence=["selection/selected.tsv", "construct/final-check.json", "construct/mrna/verification.json", "execution/local-replay.json"]),
        dict(requirement="Every-antigen evidence review, methods and reproducible comparison", status="PASS", evidence=["review/antigen-review.md", "methods.md", "exports/integrated-comparison.tsv"]),
    ]
    result = dict(status="COMPLETED", design_id="design-2", matched_normal_DNA_used=True,
                  observed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), requirements=requirements,
                  verified_local_files=sum(bool(r["local"]) for r in artifacts), native_workers_with_exit_zero=len(execution["workers"]),
                  selected_regions=24, source_counts=dict(Counter(r["source"] for r in selected)), protein_aa=724,
                  original_mutation_outcomes=dict(outcomes),
                  vaccine_sequences_and_frozen_computational_inputs_unchanged=True,
                  original_files_verified=original_bytes_verified, verification_requires_git_history=False,
                  source_snapshot_verification=source.verification_scope(),
                  scientific_limitations=["RNA counts may be sparse; gene-expression calibration unresolved", "Normal self search is bounded, not proof of healthy-tissue absence", "Recorded RNA-source translation, extra-base and strand uncertainties remain", "In-frame driver region explicitly included above 500 nM", "Residual junction binders remain; no cleavage/presentation/immunogenicity validation", "Folding heuristic 21 versus target 14; hard sequence checks pass", "PolyA/downstream extension unspecified; no complete molecular sequence", "Instrument unknown; compatible analysis reference lacks original vendor masks"],
                  website_modified=False, pushed=False)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
