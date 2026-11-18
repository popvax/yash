"""Offline replay of curation, optimal objectives, canonical assembly and synonymous RNA."""
import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "pipeline"))

from curate import write
from dataio import iter_rows, load, read_text_auto, save, sha
from protein import HERE, emit_assembly, solve
from replay_helpers import (
    canonical_path,
    compare_json_with_tolerance,
    compare_optimizer_log,
    objective_values,
    read_edges,
    read_triplets,
)

EXACT_MRNA_SEQUENCE_AND_ANNOTATION_OUTPUTS = (
    "cds.fasta",
    "mrna-segments.fasta",
    "template-segments.fasta",
    "components.tsv",
    "template-segments.gb",
    "rna-layout.txt",
)
DIAGNOSTIC_MRNA_JSON_OUTPUTS = (
    "construct.json",
    "optimizer-audit.json",
)
FLOAT_TOL = 1e-12


def reproduce(output):
    output = Path(output).resolve()
    if output.exists():
        raise ValueError("Choose a new, nonexistent replay directory")
    subprocess.run([sys.executable, str(HERE / "verify_results.py")], check=True, capture_output=True, text=True)
    output.mkdir(parents=True)
    for name in ("exports", "references", "alphagenome2", "safety"):
        shutil.copytree(HERE / name, output / name, ignore=shutil.ignore_patterns("__pycache__"))
    write(output, panel="core")
    regions = load(output / "selection/selected.json")
    edges = read_edges(HERE / "construct/protein/edge-costs.tsv")
    triplets = read_triplets(HERE / "construct/protein/processing-triplets.tsv")
    _, solved = solve(len(regions), edges, triplets, seconds=600)
    recorded = load(HERE / "construct/protein/construct.json")
    established = {phase["objective"]: phase["value"] for phase in solved["phases"]}
    expected = {phase["objective"]: phase["value"] for phase in recorded["solver"]["phases"]}
    if established != expected:
        raise ValueError("Independently replayed objective values differ")
    path = canonical_path(
        len(regions), edges, load(HERE / "construct/consensus-receipt.json")["chosen_edge_ids"]
    )
    if objective_values(edges, path, triplets) != expected:
        raise ValueError("The canonical recorded path is not a replayed optimum")
    alleles = load(HERE / "exports/source-case.json")["hla_class_I"]
    values = {
        (row["peptide"], row["allele"]): (float(row["affinity_nm"]), float(row["affinity_percentile"]))
        for row in iter_rows(HERE / "construct/all-junction-affinity.tsv.gz")
    }
    emit_assembly(output, regions, path, recorded["solver"], recorded["SEC"], recorded["MITD"], alleles, values)
    for name in ("cassette.fasta", "components.tsv", "junctions.tsv", "antigen-order.tsv",
                 "boundary-predictions.tsv", "hydrophobicity.json", "hydrophobicity-windows.tsv"):
        if sha(output / "construct/protein" / name) != sha(HERE / "construct/protein" / name):
            raise ValueError("Canonical protein replay changed: " + name)
    for name in ("architecture.json", "final-check.json"):
        shutil.copyfile(HERE / "construct" / name, output / "construct" / name)
    mrna = output / "construct/mrna"
    mrna.mkdir()
    input_names = [
        "run.py", "optimizer.py", "genbank.py", "accessibility.py", "backbone.json",
        "annotation-sources.json", "human-codon-usage.json", "requirements.txt",
    ]
    if (HERE / "construct/mrna/structure_rescue.py").is_file():
        input_names.append("structure_rescue.py")
    for name in input_names:
        shutil.copyfile(HERE / "construct/mrna" / name, mrna / name)
    subprocess.run([sys.executable, str(mrna / "run.py"), "--output", str(mrna)], check=True,
                   capture_output=True, text=True)
    nucleotide = load(HERE / "construct/mrna/verification.json")
    for name in EXACT_MRNA_SEQUENCE_AND_ANNOTATION_OUTPUTS:
        if sha(mrna / name) != nucleotide["output_sha256"][name]:
            raise ValueError("Nucleotide sequence/annotation replay changed: " + name)
    max_float_diff = 0.0
    diagnostic_exact_sha = {}
    for name in DIAGNOSTIC_MRNA_JSON_OUTPUTS:
        diagnostic_exact_sha[name] = sha(mrna / name) == nucleotide["output_sha256"][name]
        diff = compare_json_with_tolerance(
            load(mrna / name), load(HERE / "construct/mrna" / name),
            rel_tol=FLOAT_TOL, abs_tol=FLOAT_TOL, path=name,
        )
        max_float_diff = max(max_float_diff, diff)
    diagnostic_exact_sha["optimizer.log"] = (
        sha(mrna / "optimizer.log") == nucleotide["output_sha256"]["optimizer.log"]
    )
    log_diff = compare_optimizer_log(
        read_text_auto(mrna / "optimizer.log"),
        read_text_auto(HERE / "construct/mrna/optimizer.log"),
        rel_tol=FLOAT_TOL, abs_tol=FLOAT_TOL,
    )
    max_float_diff = max(max_float_diff, log_diff)

    # Also replay the 15-region KCTD19-excluded alternative protein + mRNA
    strict_out = output / "strand-strict-15"
    write(strict_out, panel="strand-strict-core")
    if (HERE / "strand-strict-15/safety").is_dir():
        shutil.copytree(HERE / "strand-strict-15/safety", strict_out / "safety", ignore=shutil.ignore_patterns("__pycache__"))
    regions15 = load(strict_out / "selection/selected.json")
    kctd_idx = next(i for i, r in enumerate(regions) if r["gene"] == "KCTD19")
    def remap15(idx):
        if idx == kctd_idx:
            return None
        if idx == len(regions):
            return len(regions15)
        return idx if idx < kctd_idx else idx - 1
    edges15 = [
        dict(e, left_index=remap15(e["left_index"]), right_index=remap15(e["right_index"]))
        for e in edges if remap15(e["left_index"]) is not None and remap15(e["right_index"]) is not None
    ]
    kept15 = {e["edge_id"] for e in edges15}
    triplets15 = [
        dict(t, region_index=remap15(t["region_index"]))
        for t in triplets if remap15(t["region_index"]) is not None
        and t["incoming_edge_id"] in kept15 and t["outgoing_edge_id"] in kept15
    ]
    _, solved15 = solve(len(regions15), edges15, triplets15, seconds=180)
    rec15_receipt = load(HERE / "strand-strict-15/contingency-receipt.json")
    est15 = {ph["objective"]: ph["value"] for ph in solved15["phases"]}
    if est15 != rec15_receipt["solver_objectives"]:
        raise ValueError("Independently replayed strand-strict-15 objective values differ")
    path15 = canonical_path(len(regions15), edges15, rec15_receipt["chosen_edge_ids"])
    rec15_prot = load(HERE / "strand-strict-15/construct/protein/construct.json")
    emit_assembly(strict_out, regions15, path15, rec15_prot["solver"], rec15_prot["SEC"], rec15_prot["MITD"], alleles, values)
    for name in ("cassette.fasta", "components.tsv", "junctions.tsv", "antigen-order.tsv",
                 "boundary-predictions.tsv", "hydrophobicity.json", "hydrophobicity-windows.tsv"):
        if sha(strict_out / "construct/protein" / name) != sha(HERE / "strand-strict-15/construct/protein" / name):
            raise ValueError("Canonical strand-strict-15 protein replay changed: " + name)
    for name in ("architecture.json", "final-check.json"):
        shutil.copyfile(HERE / "strand-strict-15/construct" / name, strict_out / "construct" / name)
    mrna15 = strict_out / "construct/mrna"
    mrna15.mkdir()
    for name in input_names:
        shutil.copyfile(HERE / "strand-strict-15/construct/mrna" / name, mrna15 / name)
    subprocess.run([sys.executable, str(mrna15 / "run.py"), "--output", str(mrna15)], check=True,
                   capture_output=True, text=True)
    nuc15 = load(HERE / "strand-strict-15/construct/mrna/verification.json")
    for name in EXACT_MRNA_SEQUENCE_AND_ANNOTATION_OUTPUTS:
        if sha(mrna15 / name) != nuc15["output_sha256"][name]:
            raise ValueError("Fallback-15 nucleotide sequence/annotation replay changed: " + name)
    diagnostic_exact_sha_15 = {}
    for name in DIAGNOSTIC_MRNA_JSON_OUTPUTS:
        diagnostic_exact_sha_15[name] = sha(mrna15 / name) == nuc15["output_sha256"][name]
        diff15 = compare_json_with_tolerance(
            load(mrna15 / name), load(HERE / "strand-strict-15/construct/mrna" / name),
            rel_tol=FLOAT_TOL, abs_tol=FLOAT_TOL, path="strand-strict-15/" + name,
        )
        max_float_diff = max(max_float_diff, diff15)
    diagnostic_exact_sha_15["optimizer.log"] = (
        sha(mrna15 / "optimizer.log") == nuc15["output_sha256"]["optimizer.log"]
    )
    max_float_diff = max(max_float_diff, compare_optimizer_log(
        read_text_auto(mrna15 / "optimizer.log"),
        read_text_auto(HERE / "strand-strict-15/construct/mrna/optimizer.log"),
        rel_tol=FLOAT_TOL, abs_tol=FLOAT_TOL,
    ))

    result = {
        "status": "OFFLINE_REPLAY_MATCHED",
        "selected_regions": len(regions), "canonical_recorded_optimum_replayed": True,
        "solver_objectives_independently_reestablished": established,
        "strand_strict_15_replayed": True,
        "strand_strict_15_solver_objectives_independently_reestablished": est15,
        "nucleotide_sequence_and_annotation_outputs_byte_identical": True,
        "nucleotide_outputs_byte_identical": all(diagnostic_exact_sha.values()) and all(diagnostic_exact_sha_15.values()),
        "diagnostic_outputs_exact_sha256": diagnostic_exact_sha,
        "strand_strict_15_diagnostic_outputs_exact_sha256": diagnostic_exact_sha_15,
        "diagnostic_outputs_within_numerical_tolerance": True,
        "diagnostic_numerical_tolerance": {"rel_tol": FLOAT_TOL, "abs_tol": FLOAT_TOL},
        "diagnostic_max_abs_float_diff": max_float_diff,
        "fresh_model_inference": False, "raw_read_alignment_or_variant_calling": False,
        "new_cloud_jobs_or_uploads": False, "clinical_validation": False,
    }
    save(output / "replay-result.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    print(json.dumps(reproduce(parser.parse_args().output), indent=2))
