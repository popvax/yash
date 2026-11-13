"""Generate the MHCflurry prediction cache for the specificity-first construct.

Queries follow src/yash_vaccine/constructs/boundaries.all_queries exactly:
every 8-11mer inside each region, every junction-overlapping peptide for each
ordered region pair and linker choice, and start-M peptides. Affinities and
percentile ranks come from MHCflurry 2.3.14 with the recorded 2.3.0 pan-allele
affinity weights, the same tool and weights as the lead design's fresh checks.
"""
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
from yash_vaccine.constructs.boundaries import all_queries, LENGTHS  # noqa: E402

MODELS = "/root/.local/share/mhcflurry/4/2.3.0/models_class1_pan/models.combined"


def main():
    from mhcflurry import Class1AffinityPredictor
    regions = list(csv.DictReader((HERE / "regions.tsv").open(), delimiter="\t"))
    linkers = json.loads((ROOT / "config/linkers.json").read_text())["choices"]
    alleles = json.loads((ROOT / "config/case.json").read_text())["hla_class_I"]
    queries = all_queries(regions, linkers)
    print(f"regions={len(regions)} linkers={len(linkers)} unique peptides={len(queries)} "
          f"pairs={len(queries) * len(alleles)}", flush=True)

    predictor = Class1AffinityPredictor.load(MODELS)
    peptide_list = list(queries)
    # predict in allele blocks; the pan model accepts arbitrary alleles
    affinity = np.empty((len(peptide_list), len(alleles)), dtype=np.float64)
    percentile = np.empty((len(peptide_list), len(alleles)), dtype=np.float64)
    t0 = time.time()
    for j, allele in enumerate(alleles):
        frame = predictor.predict_to_dataframe(
            alleles=[allele] * len(peptide_list), peptides=peptide_list)
        affinity[:, j] = frame["prediction"].to_numpy()
        percentile[:, j] = frame["prediction_percentile"].to_numpy()
        print(f"{allele}: done ({time.time() - t0:.0f}s elapsed)", flush=True)

    out = HERE / "prediction-cache.npz"
    np.savez(out, peptides=np.array(peptide_list), alleles=np.array(alleles),
             affinity_nm=affinity, affinity_percentile=percentile)
    print("wrote", out, f"in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
