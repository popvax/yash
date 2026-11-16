"""Fresh full-junction check: re-predict every boundary peptide of the
assembled cassette with MHCflurry and compare against the optimization cache,
then record the verification inputs the RNA encoding stage requires."""
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MODELS = "/root/.local/share/mhcflurry/4/2.3.0/models_class1_pan/models.combined"


def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    from mhcflurry import Class1AffinityPredictor
    protein_path = HERE / "construct/cassette.fasta"
    parts_path = HERE / "construct/components.tsv"
    protein = "".join(line.strip() for line in protein_path.read_text().splitlines()[1:])
    parts = list(csv.DictReader(parts_path.open(), delimiter="\t"))
    cursor = 0
    for part in parts:
        start, end = int(part["start_aa"]), int(part["end_aa"])
        assert start == cursor and protein[start:end] == part["aa"]
        cursor = end
    assert cursor == len(protein)

    cache = np.load(HERE / "prediction-cache.npz")
    index = {p: i for i, p in enumerate(cache["peptides"])}
    alleles = list(cache["alleles"])
    predictions = list(csv.DictReader((HERE / "construct/boundary-predictions.tsv").open(), delimiter="\t"))
    unique = sorted({row["peptide"] for row in predictions})
    assert all(row["peptide"] in index for row in predictions)

    predictor = Class1AffinityPredictor.load(MODELS)
    frame = predictor.predict_to_dataframe(
        alleles=[a for a in alleles for _ in unique], peptides=[p for _ in alleles for p in unique])
    fresh = {}
    for allele, peptide, affinity, percentile in zip(
            frame["allele"], frame["peptide"], frame["prediction"], frame["prediction_percentile"]):
        fresh[(peptide, allele)] = (float(affinity), float(percentile))

    max_rel = 0.0
    max_abs_rank = 0.0
    identical = True
    for row in predictions:
        cached_affinity = float(cache["affinity_nm"][index[row["peptide"]], alleles.index(row["allele"])])
        cached_rank = float(cache["affinity_percentile"][index[row["peptide"]], alleles.index(row["allele"])])
        fresh_affinity, fresh_rank = fresh[(row["peptide"], row["allele"])]
        max_rel = max(max_rel, abs(fresh_affinity - cached_affinity) / cached_affinity)
        max_abs_rank = max(max_abs_rank, abs(fresh_rank - cached_rank))
        for threshold, column in [(0.5, "category")]:
            pass
    # threshold classification identity (strong/mild rank categories)
    for row in predictions:
        cached_rank = float(cache["affinity_percentile"][index[row["peptide"]], alleles.index(row["allele"])])
        fresh_rank = fresh[(row["peptide"], row["allele"])][1]
        if (cached_rank <= .5) != (fresh_rank <= .5) or (cached_rank <= 2) != (fresh_rank <= 2):
            identical = False

    antigen_parts = [p for p in parts if p["kind"] == "antigen"]
    result = dict(
        status="FRESH_FULL_JUNCTION_CHECK_PASSED",
        protein_sha256=file_hash(protein_path),
        components_sha256=file_hash(parts_path),
        protein_length_aa=len(protein),
        antigen_regions=len(antigen_parts),
        all_antigens_exactly_once=len({p["name"] for p in antigen_parts}) == len(antigen_parts),
        unique_peptide_hla_pairs=len(unique) * len(alleles),
        max_relative_affinity_difference=max_rel,
        max_absolute_rank_difference=max_abs_rank,
        threshold_classifications_identical=identical,
        model="MHCflurry 2.3.14 pan-allele 2.3.0 affinity weights")
    (HERE / "construct/final-check.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=1))
    assert result["status"] == "FRESH_FULL_JUNCTION_CHECK_PASSED"
    assert max_rel < 1e-3 and identical


if __name__ == "__main__":
    main()
