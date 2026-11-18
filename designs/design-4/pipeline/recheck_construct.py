"""Independent fresh re-inference of the adopted constructs; no biological validation."""
import argparse
import os
from importlib.metadata import version
from pathlib import Path

from dataio import load, rows, save, sha, table
from hydrophobicity import protein_screen
from protein import HERE, full_scan, infer_affinity, summarize_predictions


def recheck_target(target_root, seq_name, affinity_predictor, presentation, models_dir):
    directory = target_root / "construct/protein"
    sequence = "".join(line for line in (directory / "cassette.fasta").read_text().splitlines()
                       if not line.startswith(">"))
    parts = rows(directory / "components.tsv")
    joins = rows(directory / "junctions.tsv")
    for row in parts:
        for key in ("start_aa", "end_aa"):
            row[key] = int(row[key])
    for row in joins:
        for key in ("junction_index", "linker_start_aa", "linker_end_aa"):
            row[key] = int(row[key])
    regions = load(target_root / "selection/selected.json")
    expected = {(r["region_id"], r["sequence"]) for r in regions}
    actual = [(r["region_id"], r["aa"]) for r in parts if r["kind"] == "antigen"]
    if len(actual) != len(expected) or set(actual) != expected:
        raise ValueError("An antigen was lost or changed")
    _, hydro = protein_screen(sequence, parts)
    if hydro["status"] != "PASS":
        raise ValueError("Assembled hydropathy heuristic failed")
    occurrences = full_scan(sequence, joins)
    peptides = sorted({row["peptide"] for row in occurrences})
    alleles = load(HERE / "exports/source-case.json")["hla_class_I"]
    values = infer_affinity(affinity_predictor, peptides, alleles)
    fresh = [{
        **row, "allele": allele, "affinity_nm": values[row["peptide"], allele][0],
        "affinity_percentile": values[row["peptide"], allele][1],
    } for row in occurrences for allele in alleles]
    prior = rows(directory / "boundary-predictions.tsv")
    if len(prior) != len(fresh):
        raise ValueError("Fresh junction enumeration differs")
    differences = []
    for before, after in zip(prior, fresh):
        for key in ("peptide", "allele", "junction_index", "start_aa", "end_aa"):
            if str(before[key]) != str(after[key]):
                raise ValueError("Fresh placed peptide identity differs")
        nm, rank = float(before["affinity_nm"]), float(before["affinity_percentile"])
        if [nm <= 50, nm <= 500, rank <= 0.5, rank <= 2] != [
            after["affinity_nm"] <= 50, after["affinity_nm"] <= 500,
            after["affinity_percentile"] <= 0.5, after["affinity_percentile"] <= 2]:
            raise ValueError("Fresh threshold classification differs")
        differences.append(abs(after["affinity_nm"] / nm - 1))
    if max(differences) >= 1e-4:
        raise ValueError("Fresh model numerical difference exceeds tolerance")
    table(target_root / "construct/selected-inference/predictions.tsv", fresh)
    full = presentation.predict_sequences(
        sequences={seq_name: sequence}, alleles={"Yash": alleles}, result="all"
    )
    table(target_root / "immunogenicity/cassette-class1-presentation.tsv", full.to_dict("records"))
    by_position = {(int(row["pos"]), row["peptide"]): row for row in full.to_dict("records")}
    positions = {row["region_id"]: row for row in parts if row["kind"] == "antigen"}
    lead_path = directory / "lead-targets.json"
    if not lead_path.is_file():
        kept_rids = {r["region_id"]: idx for idx, r in enumerate(regions)}
        base_leads = [
            dict(item, region_index=kept_rids[item["region_id"]])
            for item in load(HERE / "construct/protein/lead-targets.json")
            if item["region_id"] in kept_rids
        ]
        save(lead_path, base_leads)
    placed, score = [], 0
    for lead in load(lead_path):
        start = positions[lead["region_id"]]["start_aa"] + lead["start_in_region"]
        prediction = by_position[start, lead["peptide"]]
        placed.append({
            **lead, "cassette_start_aa": start,
            "processing_score": float(prediction["processing_score"]),
            "presentation_score": float(prediction["presentation_score"]),
            "presentation_percentile": float(prediction["presentation_percentile"]),
            "presentation_best_allele": prediction["best_allele"],
        })
        if lead["needs_N_neighbor"] or lead["needs_C_neighbor"]:
            score += round(float(prediction["processing_score"]) * 10000)
    solver = load(directory / "construct.json")["solver"]
    if abs(score - solver["phases"][-1]["value"]) > 1:
        raise ValueError("Two-neighbor processing objective disagrees with the placed cassette")
    table(target_root / "immunogenicity/placed-class1-leads.tsv", placed)
    save(target_root / "construct/final-check.json", {
        "status": "FRESH_FULL_JUNCTION_CHECK_PASSED",
        "fresh_repeat_model_inference_performed": True,
        "protein_sha256": sha(directory / "cassette.fasta"),
        "components_sha256": sha(directory / "components.tsv"),
        "protein_length_aa": len(sequence), "antigen_regions": len(regions),
        "all_antigens_exactly_once": True, "totals": summarize_predictions(fresh, len(alleles)),
        "hydrophobicity": hydro, "max_relative_affinity_difference": max(differences),
        "threshold_classifications_identical": True,
        "processing_objective_matches_full_cassette": True,
        "scope": "Fresh local model inference and sequence checks, not measured processing, presentation, safety, efficacy or a complete therapeutic product.",
        "clinical_release_ready": False,
    })
    receipt_path = target_root / "construct/inference-receipt.json"
    receipt = load(receipt_path) if receipt_path.is_file() else {
        "affinity_manifest_sha256": sha(Path(models_dir) / "affinity_predictor/manifest.csv"),
        "presentation_export_sha256": sha(Path(models_dir) / "inference_export.json"),
        "scope": "Fresh repeat MHCflurry Class-I junction and full-cassette presentation inference.",
    }
    receipt.update(
        fresh_selected_junction_queries=len(peptides),
        placed_junction_hla_rows=len(fresh),
        final_protein_sha256=sha(directory / "cassette.fasta"),
        library_version=version("mhcflurry"),
    )
    save(receipt_path, receipt)
    print(f"Fresh recheck ({seq_name}): {len(sequence)} aa; {len(occurrences)} placed junction peptides; max_rel_diff={max(differences):.3e}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--models", type=Path, required=True)
    parser.add_argument("--target", choices=["main", "strand-strict-15", "both"], default="main")
    args = parser.parse_args()
    os.environ["OMP_NUM_THREADS"] = "2"
    import torch
    from mhcflurry import Class1AffinityPredictor, Class1PresentationPredictor
    torch.set_num_threads(2)
    affinity_predictor = Class1AffinityPredictor.load(str(args.models / "affinity_predictor"))
    presentation = Class1PresentationPredictor.load(str(args.models))
    if args.target in ("main", "both"):
        recheck_target(HERE, "Design4_core_consensus", affinity_predictor, presentation, args.models)
    if args.target in ("strand-strict-15", "both"):
        recheck_target(HERE / "strand-strict-15", "Design4_strand_strict_15", affinity_predictor, presentation, args.models)


if __name__ == "__main__":
    main()
