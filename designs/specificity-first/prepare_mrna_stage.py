"""Stage the RNA-encoding step for this design.

Copies the recorded PopRNA stage from designs/design-2 (genbank.py, optimizer.py,
human-codon-usage.json, requirements.txt, run.py), copies the recorded routing
architecture.json unchanged, extends construct/components.tsv with the
region_id/event_id/gene annotation columns the GenBank exporter requires, and
writes this design's backbone.json and annotation-sources.json with fresh
checksums. The recorded coding-sequence pin is left null until the first
deterministic run, then pinned by pin_coding_sequence.py.
"""
import csv
import hashlib
import json
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
LEAD = HERE.parent / "design-2" / "construct"
LEAD_MRNA = LEAD / "mrna"
MRNA = HERE / "construct" / "mrna"


def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    MRNA.mkdir(exist_ok=True)
    for name in ["genbank.py", "optimizer.py", "human-codon-usage.json",
                 "requirements.txt"]:
        shutil.copy2(LEAD_MRNA / name, MRNA / name)
    shutil.copy2(LEAD / "architecture.json", HERE / "construct" / "architecture.json")

    # run.py copied from the lead, with one design-specific change: the
    # recorded coding-sequence pin is optional (null until this design's own
    # first deterministic run pins it), because the pin is protein-specific.
    run_text = (LEAD_MRNA / "run.py").read_text()
    patched = run_text.replace(
        '    require(sequence_hash(coding) == settings["recorded_coding_sequence_sha256"],\n'
        '            "PopRNA output differs from the recorded design.")',
        '    # Design-specific pin: null until the first deterministic run of THIS\n'
        '    # protein, then frozen by pin_coding_sequence.py and re-checked on reruns.\n'
        '    if settings.get("recorded_coding_sequence_sha256"):\n'
        '        require(sequence_hash(coding) == settings["recorded_coding_sequence_sha256"],\n'
        '                "PopRNA output differs from the recorded design.")')
    require_patched = patched != run_text
    assert require_patched, "Lead run.py no longer matches the expected patch site."
    (MRNA / "run.py").write_text(patched)

    # Extend components.tsv with annotation columns keyed by antigen_id.
    selected = {row["antigen_id"]: row for row in
                csv.DictReader((HERE / "selected.tsv").open(newline=""), delimiter="\t")}
    parts_path = HERE / "construct" / "components.tsv"
    parts = list(csv.DictReader(parts_path.open(newline=""), delimiter="\t"))
    fieldnames = ["kind", "name", "aa", "start_aa", "end_aa", "region_id", "event_id", "gene"]
    with parts_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        for part in parts:
            row = dict(part)
            if part["kind"] == "antigen":
                record = selected[part["name"]]
                row.update(region_id=record["region_id"], event_id=record["event_id"],
                           gene=record["gene"])
                assert record["sequence"] == part["aa"]
            else:
                row.update(region_id="", event_id="", gene="")
            writer.writerow(row)

    backbone = json.loads((LEAD_MRNA / "backbone.json").read_text())
    backbone["name"] = "specificity-first-poprna-mrna"
    backbone["protein_design"] = ".."
    backbone["protein_file_sha256"] = file_hash(HERE / "construct" / "cassette.fasta")
    backbone["protein_components_sha256"] = file_hash(parts_path)
    backbone["protein_length_aa"] = 748
    backbone["junction_verification"] = "../final-check.json"
    backbone["optimization"]["recorded_coding_sequence_sha256"] = None
    # Preserve an existing pin from a previous deterministic run of this design.
    existing = MRNA / "backbone.json"
    if existing.exists():
        previous = json.loads(existing.read_text())
        if previous.get("protein_file_sha256") == backbone["protein_file_sha256"]:
            backbone["optimization"]["recorded_coding_sequence_sha256"] = \
                previous["optimization"].get("recorded_coding_sequence_sha256")
    (MRNA / "backbone.json").write_text(json.dumps(backbone, indent=2) + "\n")

    sources = json.loads((LEAD_MRNA / "annotation-sources.json").read_text())
    sources["selected_regions"] = {
        "path": "../../selected.tsv",
        "sha256": file_hash(HERE / "selected.tsv")}
    sources["routing_architecture"] = {
        "path": "../architecture.json",
        "sha256": file_hash(HERE / "construct" / "architecture.json")}
    sources["record_date"] = "12-NOV-2026"
    (MRNA / "annotation-sources.json").write_text(json.dumps(sources, indent=2) + "\n")
    print(json.dumps({"staged": sorted(p.name for p in MRNA.iterdir()),
                      "protein_sha256": backbone["protein_file_sha256"],
                      "components_sha256": backbone["protein_components_sha256"],
                      "architecture_sha256": sources["routing_architecture"]["sha256"],
                      "selected_sha256": sources["selected_regions"]["sha256"]}, indent=2))


if __name__ == "__main__":
    main()
