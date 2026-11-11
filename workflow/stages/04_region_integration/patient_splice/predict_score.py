"""Fresh exact-self comparison, MHCflurry affinity prediction and native ranking."""
import os
os.environ["OMP_NUM_THREADS"] = "4"
os.environ["MKL_NUM_THREADS"] = "4"

from collections import Counter
import csv
import dataclasses
import datetime
import gzip
import hashlib
import importlib.metadata
import inspect
import json
import math
from pathlib import Path
import sys

import msgspec
import numpy as np
import torch
from mhcflurry import Class1AffinityPredictor
from mhctools.mhcflurry import MHCflurry_Affinity
from mhctools.pred import Kind, Prediction
from vaxrank.candidate_epitope import CandidateEpitope
from vaxrank.epitope_config import EpitopeConfig
from vaxrank.epitope_dsl import attach_per_allele_scores, default_score_expr, epitopes_for_ranking
from vaxrank.vaccine_antigen import (AminoAcidInterval, TargetableMask,
                                    TumorSpecificityAttestation, VaccineAntigen,
                                    SelfReferenceMatch, SelfReferenceSource)
from vaxrank.vaccine_peptide import VaccinePeptide

from splice_rna import target_peptides


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, data):
    Path(path).write_text(json.dumps(data, indent=2) + "\n")


def table(path, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    opener = gzip.open if Path(path).suffix == ".gz" else open
    with opener(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main(root, scoring_name="scoring"):
    out = root / scoring_name
    out.mkdir(exist_ok=False)
    policy = json.loads((root / "policy.json").read_text())
    expected_versions = {"vaxrank": "3.46.0", "mhcflurry": "2.3.14", "mhctools": "4.5.0"}
    assert {name: importlib.metadata.version(name) for name in expected_versions} == expected_versions
    reconstruction = root / "reconstruction/patient-windows.jsonl.gz"
    receipt = json.loads((root / "reconstruction/receipt.json").read_text())
    assert receipt["status"] == "COMPLETED" and receipt["policy"] == policy
    assert receipt["input_sha256"][str(root / "splice_rna.py")] == digest(root / "splice_rna.py")
    with gzip.open(reconstruction, "rt") as handle:
        windows = [json.loads(line) for line in handle]
    review_path = root / "review/peptide-review.json.gz"
    with gzip.open(review_path, "rt") as handle:
        reviews = json.load(handle)
    review_receipt = json.loads((root / "review/receipt.json").read_text())
    assert review_receipt["status"] == "COMPLETED" and review_receipt["reconstruction_sha256"] == digest(reconstruction)
    assert review_receipt["output_sha256"] == digest(review_path)
    candidates = [w for w in windows if w["native_event_review_tier"] == "PASS"]
    peptides = sorted({t["peptide"] for w in candidates for t in target_peptides(w)})
    assert peptides, "No RNA-derived windows pass native event review"
    alleles = policy["hla_alleles"]
    models = Path("/data/results/TN25-112103/openvax-mutations-20261018/references/mhcflurry/2.3.0/models_class1_presentation/models")
    model_manifest = json.loads((root / "model-files.json").read_text())
    for record in model_manifest:
        path = models / record["path"]
        assert path.stat().st_size == record["bytes"] and digest(path) == record["sha256"], str(path)
    print("FRESH_QUERIES", len(peptides), "peptides", len(peptides) * len(alleles), "peptide/HLA predictions", flush=True)
    save(root / "status.json", dict(stage="fresh_reference_comparison", peptides=len(peptides),
                                    native_PASS_patient_windows=len(candidates)))
    reference_path = root / "reference-check/matches.json.gz"
    reference_receipt = json.loads((root / "reference-check/receipt.json").read_text())
    assert reference_receipt["status"] == "PASS"
    assert reference_receipt["reconstruction_sha256"] == digest(reconstruction)
    assert reference_receipt["output_sha256"] == digest(reference_path)
    assert reference_receipt["all_reference_proteins_checked_by_two_exact_matchers"]
    config = dict(reference_name="GRCh38", annotation_name="Ensembl", annotation_version=103,
                  reference_inventory_sha256=reference_receipt["inventory_sha256"],
                  protein_fasta_sha256=reference_receipt["reference_protein_fasta_sha256"],
                  GTF_sha256=reference_receipt["reference_gtf_sha256"],
                  reference_matching="Complete native inventory; multi-pattern scan independently checked by exhaustive k-mer enumeration")
    attestation = TumorSpecificityAttestation(
        status="overridden", evidence_kind="patient_RNA_splice_reconstruction",
        evidence_source="LENS normal screens and native read review; actual patient query bases",
        requires_review=True, override_reason="RNA-only research hypotheses; normal patient RNA and protein presentation unmeasured")
    with gzip.open(reference_path, "rt") as handle:
        raw_matches = json.load(handle)
    matches = {peptide: SelfReferenceMatch(**dict(record,
                   excluded_gene_ids=tuple(record["excluded_gene_ids"]),
                   sources=tuple(SelfReferenceSource(**source) for source in record["sources"])))
               for peptide, record in raw_matches.items()}
    assert set(matches) == set(peptides) and all(m.source_provenance_complete for m in matches.values())
    with gzip.open(out / "reference-matches.json.gz", "wt") as handle:
        json.dump({p: dataclasses.asdict(m) for p, m in matches.items()}, handle)
    print("EXACT_SELF", sum(m.occurs for m in matches.values()), flush=True)
    torch.set_num_threads(4)
    assert torch.cuda.is_available(), "The recorded fresh-inference worker requires CUDA"
    torch.set_float32_matmul_precision("highest")
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    np.random.seed(0)
    torch.manual_seed(0)
    torch.cuda.manual_seed_all(0)
    predictor = Class1AffinityPredictor.load(str(models / "affinity_predictor"))
    predictions, by_peptide = [], {}
    for allele in alleles:
        frame = predictor.predict_to_dataframe(peptides, allele=allele, model_kwargs={"batch_size": 16384},
                    include_percentile_ranks=True, include_confidence_intervals=False)
        assert frame.peptide.tolist() == peptides
        for row in frame.itertuples():
            record = dict(peptide=row.peptide, allele=allele, affinity_nm=float(row.prediction),
                          percentile_rank=float(row.prediction_percentile))
            assert math.isfinite(record["affinity_nm"]) and record["affinity_nm"] > 0 and 0 <= record["percentile_rank"] <= 100
            predictions.append(record)
            by_peptide.setdefault(record["peptide"], []).append(record)
        print("PREDICTED", allele, len(peptides), flush=True)
    # Verify all direct model results against the affinity-only OpenVax adapter.
    wrapper = MHCflurry_Affinity(alleles=alleles, default_peptide_lengths=policy["peptide_lengths"],
                                models_path=str(models / "affinity_predictor"))
    seen, largest_relative_difference = set(), 0.0
    values = {(r["peptide"], r["allele"]): r["affinity_nm"] for r in predictions}
    for peptide, result in zip(peptides, wrapper.predict(peptides), strict=True):
        for prediction in result.preds:
            key = (peptide, prediction.allele)
            assert prediction.kind == Kind.pMHC_affinity and prediction.peptide == peptide and key not in seen
            seen.add(key)
            largest_relative_difference = max(largest_relative_difference, abs(float(prediction.value) / values[key] - 1))
    assert seen == set(values) and largest_relative_difference < 1e-4
    assert set(by_peptide) == set(peptides)
    assert all({r["allele"] for r in rows} == set(alleles) and len(rows) == len(alleles)
               for rows in by_peptide.values())
    table(out / "affinities.tsv.gz", predictions)
    epitope_config = EpitopeConfig()
    bases, unscored = {}, []
    for peptide, rows in sorted(by_peptide.items()):
        preds = tuple(Prediction(kind=Kind.pMHC_affinity, peptide=peptide, allele=row["allele"],
                                value=row["affinity_nm"], score=row["affinity_nm"],
                                percentile_rank=row["percentile_rank"], predictor_name="mhcflurry",
                                predictor_version="2.3.14+models2.3.0") for row in rows)
        unscored.append(CandidateEpitope(sequence=peptide, source_sequence=peptide, offset=0,
                       predictions=preds, overlaps_targetable=True, patient_alleles=tuple(alleles),
                       self_reference_match=matches[peptide], occurs_in_reference=matches[peptide].occurs))
    for offset in range(0, len(unscored), 5000):
        scored = epitopes_for_ranking(attach_per_allele_scores(unscored[offset:offset + 5000], epitope_config), epitope_config)
        bases.update({e.sequence: e for e in scored})
    save(out / "native-score-config.json", dict(vaxrank_version=importlib.metadata.version("vaxrank"),
         epitope_config=msgspec.to_builtins(epitope_config), resolved_expression=default_score_expr(epitope_config),
         affinity_cutoff_nm=policy["admission_max_affinity_nm"], expression_weight=policy["expression_weight"],
         reference_genome_config=config, query_source="Actual patient RNA; no original LENS peptide seed required"))
    results, checks = [], Counter()
    for window in windows:
        row = {k: v for k, v in window.items()
               if k not in {"fragment_names", "genomic_positions_0based", "genomic_path_support", "upstream_frame_anchors", "native_event_review_flags"}}
        row["native_event_review_flags"] = ";".join(window["native_event_review_flags"])
        row["region_id"] = "REG_" + hashlib.sha256(window["sequence"].encode()).hexdigest()[:20]
        row.update(primary_window_eligible=False, rejection_reason="NATIVE_EVENT_REVIEW_NOT_PASS")
        if window["native_event_review_tier"] != "PASS":
            results.append(row)
            continue
        targets = list(target_peptides(window))
        # All native flags are applied to the actual patient coding sequence.
        # Additional weak peptides contribute; admission needs a PASS seed.
        qualified = [t for t in targets if not matches[t["peptide"]].occurs
                     and min(p["affinity_nm"] for p in by_peptide[t["peptide"]]) <= policy["admission_max_affinity_nm"]
                     and reviews[hashlib.sha256((window["event_id"] + "|" + window["coding_sequence"][3*t["start"]:3*t["end"]]).encode()).hexdigest()]["native_review_tier"] == "PASS"]
        epitopes = [dataclasses.replace(bases[t["peptide"]], source_sequence=window["sequence"],
                                       source_name=window["window_id"], offset=t["start"])
                    for t in targets if t["peptide"] in bases]
        antigen = VaccineAntigen(kind="splice", amino_acids=window["sequence"],
                    targetable_mask=TargetableMask((AminoAcidInterval(window["target_start_in_window"], window["target_end_in_window"]),)),
                    tumor_specificity=attestation, gene_name=window["gene"], source_identifier=window["window_id"])
        vaccine = VaccinePeptide(antigen=antigen, epitopes=epitopes, num_target_epitopes_to_keep=0,
                                 combined_score_expr="target_epitope_score", ranking_rules=("target_epitope_score",))
        expected = sum(e.epitope_score for e in epitopes if not e.occurs_in_self_reference)
        assert math.isclose(expected, vaccine.target_epitope_score, rel_tol=1e-9, abs_tol=1e-9)
        checks["native_window_score_parity"] += 1
        nonself = [t for t in targets if not matches[t["peptide"]].occurs]
        best = min((r["affinity_nm"] for t in nonself for r in by_peptide[t["peptide"]]), default=None)
        best_seed = min((r["affinity_nm"] for t in qualified for r in by_peptide[t["peptide"]]), default=None)
        best_record = min((r for t in nonself for r in by_peptide[t["peptide"]]), key=lambda r: r["affinity_nm"], default=None)
        expression = window["supplied_transformed_gene_tpm"]
        score = vaccine.target_epitope_score
        row.update(native_affinity_sum=score, best_target_affinity_nm=best,
                   best_target_peptide=best_record["peptide"] if best_record else "",
                   best_target_allele=best_record["allele"] if best_record else "",
                   best_admitted_affinity_nm=best_seed,
                   exact_self_target_occurrences=sum(matches[t["peptide"]].occurs for t in targets),
                   nonself_target_occurrences=len(nonself), seed_peptides=";".join(sorted({t["peptide"] for t in qualified})),
                   target_peptides=";".join(sorted({e.sequence for e in vaccine.target_epitopes})),
                   weighted_score_left=score * math.sqrt(expression) if expression is not None and expression > 0 else "",
                   primary_window_eligible=bool(qualified) and score > 0,
                   rejection_reason="" if qualified and score > 0 else "NO_NONSELF_PASS_PEPTIDE_LE500NM")
        results.append(row)
    table(out / "all-windows.tsv.gz", results)
    eligible_rows = [r for r in results if r["primary_window_eligible"]]
    def rank_key(row):
        return (row["weighted_score_left"] == "", -float(row["weighted_score_left"] or 0),
                -float(row["native_affinity_sum"]), row["event_id"], row["window_id"])
    event_best = {}
    for row in sorted(eligible_rows, key=rank_key):
        event_best.setdefault(row["event_id"], row)
    event_ranking = [dict(event_rank=i, **r) for i, r in enumerate(sorted(event_best.values(), key=rank_key), 1)]
    gene_best = {}
    for row in sorted(event_best.values(), key=rank_key):
        gene_best.setdefault(row["gene_id"], row)
    gene_ranking = [dict(gene_rank=i, **r) for i, r in enumerate(sorted(gene_best.values(), key=rank_key), 1)]
    selected = [r for r in gene_ranking if r["weighted_score_left"] != ""][:5]
    table(out / "event-ranking.tsv", event_ranking)
    table(out / "gene-ranking.tsv", gene_ranking)
    table(out / "selected.tsv", selected)
    (out / "selected.fasta").write_text("".join(
        f">{r['window_id']} gene={r['gene']} event={r['event_id']} fragments={r['fragment_count']}\n{r['sequence']}\n"
        for r in selected))
    verify = dict(status="PASS", prediction_rows_complete=len(predictions) == len(peptides) * len(alleles),
                  native_window_scores_verified=checks["native_window_score_parity"],
                  selected_distinct_genes=len({r["gene_id"] for r in selected}) == len(selected),
                  all_selected_have_patient_window_support=all(r["fragment_count"] >= policy["min_full_window_fragments"] for r in selected),
                  all_selected_have_qualifying_binding=all(r["best_admitted_affinity_nm"] <= policy["admission_max_affinity_nm"] for r in selected),
                  all_selected_have_positive_expression=all(r["supplied_transformed_gene_tpm"] > 0 for r in selected),
                  selected_regions=len(selected),
                  identity_filter_applied=False, splice_fraction_used=False)
    save(out / "verification.json", verify)
    native_modules = ["vaxrank.reference_proteome", "vaxrank.epitope_dsl", "vaxrank.vaccine_peptide", "mhctools", "mhcflurry"]
    receipt = dict(status="COMPLETED", completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        reconstruction_sha256=digest(reconstruction), reconstruction_receipt_sha256=digest(root / "reconstruction/receipt.json"),
        patient_peptide_review_sha256=digest(review_path),
        exact_reference_check_receipt_sha256=digest(root / "reference-check/receipt.json"),
        worker_sha256=digest(Path(__file__)), policy_sha256=digest(root / "policy.json"),
        model_manifest_sha256=digest(root / "model-files.json"), model_files_verified=True,
        package_versions={name: importlib.metadata.version(name) for name in ["vaxrank", "mhcflurry", "mhctools", "pyensembl", "torch"]},
        gpu=torch.cuda.get_device_name(), mhctools_affinity_adapter_parity_pairs=len(seen),
        maximum_adapter_relative_affinity_difference=largest_relative_difference,
        native_source_sha256={name: digest(inspect.getfile(__import__(name, fromlist=[""]))) for name in native_modules},
        native_PASS_patient_windows=len(candidates), unique_predicted_peptides=len(peptides), prediction_rows=len(predictions),
        exact_self_peptides=sum(m.occurs for m in matches.values()),
        reference_inventory_sha256=sorted({m.reference_sha256 for m in matches.values()}),
        eligible_windows=len(eligible_rows), eligible_events=len(event_ranking), eligible_genes=len(gene_ranking),
        weighted_genes=sum(r["weighted_score_left"] != "" for r in gene_ranking),
        selected_genes=[r["gene"] for r in selected],
        output_sha256={p.name: digest(p) for p in sorted(out.iterdir()) if p.name != "receipt.json"})
    save(out / "receipt.json", receipt)
    save(root / "status.json", dict(stage="scoring_complete", **receipt))
    print(json.dumps({k: v for k, v in receipt.items() if k not in {"native_source_sha256", "output_sha256"}}, indent=2), flush=True)


if __name__ == "__main__":
    main(Path(sys.argv[1]), sys.argv[2] if len(sys.argv) > 2 else "scoring")
