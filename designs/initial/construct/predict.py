"""Affinity-only worker for explicit routed-junction requests; recorded runtime."""
import os
os.environ["OMP_NUM_THREADS"] = "4"
os.environ["MKL_NUM_THREADS"] = "4"

import argparse
import csv
import hashlib
import importlib.metadata
import inspect
import json
import time
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--models", type=Path, required=True)
    parser.add_argument("--model-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    inputs = json.loads(args.inputs.read_text())
    if args.output.exists():
        raise ValueError("Choose a new output directory")

    import numpy as np
    import torch
    from mhcflurry import Class1AffinityPredictor
    from mhctools.mhcflurry import MHCflurry_Affinity
    from mhctools.pred import Kind
    import vaxrank.junction_swap as native

    versions = {name: importlib.metadata.version(name) for name in inputs["expected_versions"]}
    if versions != inputs["expected_versions"]:
        raise ValueError(f"Runtime differs from the recorded predictions: {versions}")
    if not torch.cuda.is_available():
        raise ValueError("Recorded runtime requires CUDA")
    torch.set_num_threads(4)
    torch.set_float32_matmul_precision("highest")
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    np.random.seed(0)
    torch.manual_seed(0)
    torch.cuda.manual_seed_all(0)

    # Require the same complete model directory as the original inference.
    expected_files = json.loads(args.model_manifest.read_text())
    observed_files = [dict(path=str(p.relative_to(args.models)), bytes=p.stat().st_size,
                           sha256=digest(p)) for p in sorted(args.models.rglob("*")) if p.is_file()]
    if observed_files != expected_files:
        raise ValueError("Model files differ from the recorded construct model bundle")

    # Check every supplied junction/linker against native Vaxrank enumeration.
    enumerated = set()
    for join in inputs["joins"]:
        enumerated.update(native.junction_kmers(
            join["left"], join["linker"], join["right"], inputs["lengths"]))
    if enumerated != set(inputs["peptides"]):
        raise ValueError("Request differs from native Vaxrank junction enumeration")

    peptides, alleles = inputs["peptides"], inputs["alleles"]
    start = time.perf_counter()
    predictor = Class1AffinityPredictor.load(str(args.models / "affinity_predictor"))
    values = {}
    for allele in alleles:
        frame = predictor.predict_to_dataframe(
            peptides, allele=allele, model_kwargs={"batch_size": 16384},
            include_percentile_ranks=True, include_confidence_intervals=False)
        if frame.peptide.tolist() != peptides:
            raise ValueError("Prediction peptide order mismatch")
        for row in frame.itertuples():
            affinity, rank = float(row.prediction), float(row.prediction_percentile)
            if not (np.isfinite(affinity) and affinity > 0 and np.isfinite(rank) and 0 <= rank <= 100):
                raise ValueError("Invalid model output")
            values[row.peptide, allele] = (affinity, rank)
        print(f"PREDICTED {allele}: {len(peptides)} peptides", flush=True)

    # Check the complete request against the mhctools adapter used by OpenVax.
    wrapper = MHCflurry_Affinity(alleles=alleles, default_peptide_lengths=inputs["lengths"],
                                models_path=str(args.models / "affinity_predictor"))
    results = list(wrapper.predict(peptides))
    if len(results) != len(peptides):
        raise ValueError("Wrapper result count mismatch")
    seen, affinity_diffs, rank_diffs = set(), [], []
    for peptide, result in zip(peptides, results):
        for pred in result.preds:
            if pred.kind != Kind.pMHC_affinity:
                raise ValueError("Affinity-only adapter returned another prediction kind")
            key = (peptide, pred.allele)
            if pred.peptide != peptide or key in seen or key not in values:
                raise ValueError("Unexpected or duplicated wrapper prediction")
            seen.add(key)
            affinity, rank = values[key]
            affinity_diffs.append(abs(float(pred.value) / affinity - 1))
            rank_diffs.append(abs(float(pred.percentile_rank) - rank))
    if seen != set(values) or max(affinity_diffs) >= 1e-4 or max(rank_diffs) >= 1e-4:
        raise ValueError("OpenVax adapter parity check failed")

    args.output.mkdir(parents=True)
    with (args.output / "predictions.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["peptide", "allele", "affinity_nm", "affinity_percentile"])
        for peptide in peptides:
            for allele in alleles:
                writer.writerow([peptide, allele, *values[peptide, allele]])
    record = dict(
        endpoint="MHC-I binding affinity", all_queries_freshly_predicted=True,
        processing_or_presentation_model_used=False, unique_peptides=len(peptides),
        peptide_hla_pairs=len(values), alleles=alleles, peptide_lengths=inputs["lengths"],
        model_bundle="MHCflurry 2.3.0 models_class1_presentation: affinity submodel only",
        versions=versions, gpu=torch.cuda.get_device_name(), elapsed_seconds=time.perf_counter() - start,
        cuda_allocated_peak_bytes=torch.cuda.max_memory_allocated(),
        input_sha256=digest(args.inputs), script_sha256=digest(Path(__file__)),
        model_manifest_sha256=digest(args.model_manifest), all_model_files_match_reference=True,
        native_junction_enumeration_matches=True, junction_options_checked=len(inputs["joins"]),
        native_junction_source_sha256=digest(Path(inspect.getfile(native))),
        mhctools_adapter="MHCflurry_Affinity", mhctools_parity_pairs=len(seen), max_relative_affinity_difference=max(affinity_diffs),
        max_absolute_rank_difference=max(rank_diffs), predictions_sha256=digest(args.output / "predictions.tsv"))
    (args.output / "inference.json").write_text(json.dumps(record, indent=2) + "\n")
    print(f"VERIFIED {len(values)} unique peptide-HLA predictions", flush=True)


if __name__ == "__main__":
    main()
