"""Fixed inputs, position-aware enumeration and affinity cache checks."""
import csv
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
LENGTHS = [8, 9, 10, 11]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + "\n")


def rows(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_rows(path, values):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(values[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(values)


def fasta(path, name, sequence):
    Path(path).write_text(">" + name + "\n" + "".join(sequence[i:i+80] + "\n" for i in range(0, len(sequence), 80)))


def read_fasta(path):
    lines = Path(path).read_text().splitlines()
    assert sum(x.startswith(">") for x in lines) == 1
    return "".join(x for x in lines if not x.startswith(">"))


def junction_peptides(left, linker, right):
    sequence = left[-10:] + linker + right[:10]
    begin, end = len(left[-10:]), len(left[-10:]) + len(linker)
    return [sequence[start:start+k] for k in LENGTHS for start in range(len(sequence)-k+1)
            if (start < end and start+k > begin if linker else start < begin < start+k)]


def inputs():
    initial = REPO / "designs/initial/selection"
    expected = read_json(initial / "verification.json")["artifact_sha256"]
    for name in ["selected.tsv"]:
        if sha(initial / name) != expected[name]:
            raise ValueError(f"Selection input differs: {name}")
    regions = rows(initial / "selected.tsv")
    linkers = read_json(REPO / "config/linkers.json")["choices"]
    recipe = read_json(initial / "design.json")
    if linkers != recipe["linkers"]["choices"]:
        raise ValueError("Linker choices differ from the selection recipe")
    architecture = read_json(HERE / "architecture.json")
    sec, mitd = architecture["elements"]["SEC"]["aa"], architecture["elements"]["MITD"]["aa"]
    source_architecture = read_json(HERE / "terminal-inference/source-architecture.json")
    if architecture["elements"] != source_architecture["elements"]:
        raise ValueError("Routing elements differ from the recorded inference inputs")
    case = read_json(REPO / "config/case.json")
    evidence = read_json(HERE / architecture["sequence_evidence"])
    typing = read_json(HERE / architecture["typing_evidence"])
    for element in architecture["elements"].values():
        allele = element["patient_allele"]
        if allele not in case["hla_class_I"]:
            raise ValueError("Routing element does not match a recorded patient HLA allele")
        locus = allele.removeprefix("HLA-").split("*")[0]
        for dataset in ("dna", "rna"):
            calls = typing["loci"][locus][dataset]["alleles"]
            if not any("HLA-" + c["allele"] == allele and c["passes_tool_quality"]
                       for c in calls if c["allele"]):
                raise ValueError("Patient HLA allele lacks DNA/RNA typing support")
            target = evidence["datasets"][dataset]["targets"][element["evidence_target"]]
            if (target["protein"] != element["aa"]
                    or target["protein_interval_1based"] != element["protein_interval_1based_inclusive"]
                    or not target["every_nucleotide_majority_matches_reference"]
                    or target["min_expected_nucleotide_fragments"] <= 0
                    or (dataset == "rna" and target["min_expected_codon_fragments"] <= 0)):
                raise ValueError("Routing-element sequence lacks the recorded read support")
    if sha(REPO / "config/case.json") != recipe["evidence"]["context_sha256"]:
        raise ValueError("Case context changed")
    assert len(regions) == len({r["region_id"] for r in regions}) == sum(recipe["selection"]["quotas"].values())
    assert len(linkers) == 11 and len(sec) == 26 and len(mitd) == 55
    assert all(len(r["sequence"]) == 25 for r in regions)
    sources = ["designs/initial/selection/selected.tsv", "config/linkers.json", "config/case.json",
               "designs/initial/construct/architecture.json",
               "designs/initial/construct/allele-evidence/read-evidence.json",
               "designs/initial/construct/allele-evidence/typing-comparison.json",
               "data/checkpoints/provenance/construct-inference.json",
               "data/checkpoints/provenance/construct-model-files.json"]
    # Canonical repository paths and checksums identify the inputs to this inference.
    locations = read_json(HERE / "input-locations.json")
    if set(locations) != set(sources):
        raise ValueError("Input-location manifest differs from the recorded request identifiers")
    return regions, linkers, sec, mitd, case["hla_class_I"], {s: sha(REPO / locations[s]) for s in sources}


def request(joins, scope, input_hashes):
    _, _, _, _, alleles, _ = inputs()
    peptides = sorted({p for j in joins for p in junction_peptides(j["left"], j["linker"], j["right"])})
    return dict(scope=scope, expected_versions=read_json(
        REPO / "data/checkpoints/provenance/construct-inference.json")["versions"],
        lengths=LENGTHS, alleles=alleles, input_sha256=input_hashes, joins=joins, peptides=peptides)


def terminal_request():
    regions, linkers, sec, mitd, _, hashes = inputs()
    joins = []
    for r in regions:
        for linker in linkers:
            joins.extend([
                dict(name="SEC_to_" + r["region_id"] + ":" + linker["name"],
                     left=sec, linker=linker["sequence"], right=r["sequence"]),
                dict(name=r["region_id"] + "_to_MITD:" + linker["name"],
                     left=r["sequence"], linker=linker["sequence"], right=mitd)])
    return request(joins, "All terminal antigen/linker options, with SEC and MITD fixed", hashes)


def load_predictions(directory, expected_request):
    directory = Path(directory)
    saved = read_json(directory / "inputs.json")
    record = read_json(directory / "inference.json")
    if saved != expected_request:
        raise ValueError("Prediction request differs from current inputs")
    if (record["input_sha256"] != sha(directory / "inputs.json") or
        record["script_sha256"] != sha(HERE / "predict.py") or
        record["predictions_sha256"] != sha(directory / "predictions.tsv") or
        record["versions"] != saved["expected_versions"] or
        record["model_manifest_sha256"] != sha(REPO / "data/checkpoints/provenance/construct-model-files.json") or
        not record["all_model_files_match_reference"] or
        not record["native_junction_enumeration_matches"] or
        record["processing_or_presentation_model_used"]):
        raise ValueError("Prediction provenance mismatch")
    values = {}
    for row in rows(directory / "predictions.tsv"):
        key = (row["peptide"], row["allele"])
        a, rank = float(row["affinity_nm"]), float(row["affinity_percentile"])
        if key in values or not (math.isfinite(a) and a > 0 and math.isfinite(rank) and 0 <= rank <= 100):
            raise ValueError("Invalid prediction")
        values[key] = (a, rank)
    expected = {(p, a) for p in saved["peptides"] for a in saved["alleles"]}
    if set(values) != expected or record["peptide_hla_pairs"] != len(expected):
        raise ValueError("Prediction coverage mismatch")
    return values


def full_scan(sequence, joins):
    """Independently scan full protein against placed linker intervals/direct joins."""
    found = []
    for length in LENGTHS:
        for start in range(len(sequence)-length+1):
            end = start+length
            touched = [j for j in joins if
                (start < j["linker_end_aa"] and end > j["linker_start_aa"]
                 if j["linker_end_aa"] > j["linker_start_aa"]
                 else start < j["linker_start_aa"] < end)]
            if len(touched) > 1:
                raise ValueError("One peptide spans multiple joins; pairwise objective is invalid")
            if touched:
                j = touched[0]
                found.append(dict(junction_index=j["junction_index"], category=j["category"],
                    start_aa=start, end_aa=end, peptide=sequence[start:end]))
    return found
