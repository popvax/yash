"""Secondary junction-model agreement preferences; no clinical validation is implied."""
import argparse
import os
import shutil
from pathlib import Path

from dataio import iter_rows, load, rows, save, sha, table
from protein import HERE, create_edges, emit_assembly, solve


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--models", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=HERE / "consensus-search")
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("Use a new consensus output directory")
    os.environ["OMP_NUM_THREADS"] = "2"
    import torch
    from mhcflurry import Class1PresentationPredictor
    torch.set_num_threads(2)
    regions = load(HERE / "selection/selected.json")
    linkers = load(HERE / "construct/linkers.json")["choices"]
    architecture = load(HERE / "construct/architecture.json")
    sec, mitd = (architecture["elements"][name]["aa"] for name in ("SEC", "MITD"))
    alleles = load(HERE / "exports/source-case.json")["hla_class_I"]
    values = {
        (row["peptide"], row["allele"]): (float(row["affinity_nm"]), float(row["affinity_percentile"]))
        for row in iter_rows(HERE / "construct/all-junction-affinity.tsv.gz")
    }
    mix = {row["Peptide"]: row for row in rows(HERE / "immunogenicity/mixmhcpred3.tsv")}
    prime = {row["Peptide"]: row for row in rows(HERE / "immunogenicity/prime21.tsv")}
    c2_map = {
        (row["Peptide"], row["Context"]): row
        for row in iter_rows(HERE / "construct/all-junction-class2.tsv.gz")
    }
    predictor = Class1PresentationPredictor.load(str(args.models))
    edges, triplets, leads = create_edges(
        regions, linkers, sec, mitd, alleles, values, predictor.processing_predictor_with_flanks,
        load(HERE / "references/peptide-matches.json"), mix, prime, c2_map,
    )
    path, solver = solve(len(regions), edges, triplets, seconds=600)
    (args.output / "selection").mkdir(parents=True)
    shutil.copyfile(HERE / "selection/selected.tsv", args.output / "selection/selected.tsv")
    save(args.output / "selection/selected.json", regions)
    table(args.output / "construct/protein/edge-costs.tsv", edges)
    table(args.output / "construct/protein/processing-triplets.tsv", triplets,
          fields=["region_index", "incoming_edge_id", "outgoing_edge_id", "peptide",
                  "n_flank", "c_flank", "processing_score_x10000"])
    save(args.output / "construct/protein/lead-targets.json", leads)
    sequence, _, _, _, totals = emit_assembly(
        args.output, regions, path, solver, sec, mitd, alleles, values
    )
    selected_ids = {edge["edge_id"] for edge in path}
    save(args.output / "receipt.json", {
        "protein_sha256": sha(args.output / "construct/protein/cassette.fasta"),
        "source_affinity_cache_sha256": sha(HERE / "construct/all-junction-affinity.tsv.gz"),
        "source_class2_cache_sha256": sha(HERE / "construct/all-junction-class2.tsv.gz"),
        "source_mix_sha256": sha(HERE / "immunogenicity/mixmhcpred3.tsv"),
        "source_prime_sha256": sha(HERE / "immunogenicity/prime21.tsv"),
        "totals": totals,
        "model_scope": "Class I 8-11mers across six patient alleles plus Class II 15mer synthetic 9-aa cores across nine conditional chain combinations",
        "mix_strong_rows": sum(edge["mix_strong_rank"] for edge in path),
        "mixmhc2_strong_core_rows": sum(edge["mixmhc2_strong_core"] for edge in path),
        "mixmhc2_mild_core_rows": sum(edge["mixmhc2_mild_core"] for edge in path),
        "prime_low_rank_rows": sum(edge["prime_low_rank"] for edge in path),
        "source_processing_terms": len(triplets),
        "chosen_edge_ids": sorted(selected_ids),
        "clinical_validation": False,
    })
    print("Consensus preview:", len(sequence), "aa", solver["consensus_objective_values"], totals)


if __name__ == "__main__":
    main()
