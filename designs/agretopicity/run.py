"""Agretopicity for the substitution regions of Design #2.

Agretopicity follows LENS (lenstools/src/metadata.py, add_agreto_metadata, line
3542): the ratio is mutant affinity divided by wild-type affinity, both in nM,
so a LOWER value means the mutant is more distinct from the self peptide it came
from. That is the reciprocal of the DAI convention used in much of the
literature, where higher is better. Getting the direction wrong inverts every
conclusion, so the column is named agretopicity_mt_over_wt throughout.

Both sides are predicted fresh with the same model. No cached prediction is
relabelled, and no threshold is applied: this reports a quantity and does not
rank, filter or change the selection.

Run from the repository root:

    python designs/agretopicity/run.py --models "$(mhcflurry downloads path \
        models_class1_presentation --release 2.3.0)/models/affinity_predictor"
"""
import argparse, csv, importlib.metadata, json, re, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
LENGTHS = (8, 9, 10, 11)


def cassette_labels(repo):
    """Map stable region_id to the AG label Design #2 gave it in its cassette.

    The AG numbering belongs to that one selection, so it is carried only as a
    reading aid. region_id is the identifier this analysis keys on.
    """
    selected = json.loads((repo / "designs/design-2/selection/selected.json").read_text())
    headers = [line[1:].strip().split("|") for line in
               (repo / "designs/design-2/selection/selected.fasta").read_text().splitlines()
               if line.startswith(">")]
    if len(headers) != len(selected):
        raise ValueError("selected.fasta and selected.json disagree on region count")
    labels = {}
    for record, header in zip(selected, headers):
        label, source, gene = header
        if (source, gene) != (record["source"], record["gene"]):
            raise ValueError("selected.fasta and selected.json disagree at " + label)
        labels[record["region_id"]] = label
    return labels


def substitution_regions(repo):
    """Selected MUTATION regions that are single-residue substitutions."""
    pool = repo / "designs/design-2/inputs/integrated/candidate-pool.tsv"
    with pool.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    labels = cassette_labels(repo)
    regions, skipped = [], []
    for row in rows:
        if row["source"] != "MUTATION" or not row["region_id"]:
            continue
        match = re.fullmatch(r"p\.([A-Z])(\d+)([A-Z])", row["protein_change"])
        if not match:
            skipped.append(dict(region_id=row["region_id"], gene=row["gene"],
                                protein_change=row["protein_change"],
                                reason="not a single-residue substitution, so it has no "
                                       "wild-type counterpart to compare against"))
            continue
        reference, position, alternate = match.group(1), int(match.group(2)), match.group(3)
        window, index = row["sequence"], int(row["target_start_in_window"])
        if window[index] != alternate:
            raise ValueError("Alternate residue is not at target_start_in_window: " + row["region_id"])
        regions.append(dict(
            region_id=row["region_id"], ag_label_in_design_2=labels[row["region_id"]],
            gene=row["gene"], variant_id=row["variant_id"], protein_change=row["protein_change"],
            protein_position=position, reference_residue=reference, alternate_residue=alternate,
            mutant_window=window, index_in_window=index,
            wildtype_window=window[:index] + reference + window[index + 1:],
            recorded_best_affinity_nm=row["best_target_affinity_nm"]))
    return regions, skipped


def covering_peptides(window, index):
    """Every 8 to 11mer of the window that contains the substituted position."""
    for length in LENGTHS:
        for start in range(0, len(window) - length + 1):
            if start <= index < start + length:
                yield length, start


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--repo", type=Path, default=REPO)
    parser.add_argument("--models", type=Path, required=True,
                        help="affinity_predictor directory inside models_class1_presentation 2.3.0")
    parser.add_argument("--output", type=Path, default=HERE / "results")
    parser.add_argument("--model-manifest", type=Path,
                        default=REPO / "designs/initial/inputs/checkpoints/provenance/construct-model-files.json",
                        help="pinned bundle manifest the models are checked against")
    args = parser.parse_args()

    alleles = json.loads((args.repo / "config/case.json").read_text())["hla_class_I"]
    regions, skipped = substitution_regions(args.repo)

    pairs = []
    for region in regions:
        for length, start in covering_peptides(region["mutant_window"], region["index_in_window"]):
            mutant = region["mutant_window"][start:start + length]
            wildtype = region["wildtype_window"][start:start + length]
            if len(mutant) != length or mutant == wildtype:
                raise ValueError("Unexpected peptide pair in " + region["region_id"])
            pairs.append(dict(region, length=length, offset_in_window=start,
                              peptide_mutant=mutant, peptide_wildtype=wildtype))
    queries = sorted({p["peptide_mutant"] for p in pairs} | {p["peptide_wildtype"] for p in pairs})

    import hashlib
    bundle = args.models.parent
    files = [dict(path=str(path.relative_to(bundle)), bytes=path.stat().st_size,
                  sha256=hashlib.sha256(path.read_bytes()).hexdigest())
             for path in sorted(bundle.rglob("*")) if path.is_file()]
    pinned = json.loads(args.model_manifest.read_text())
    models_match_pinned_bundle = files == pinned

    from mhcflurry import Class1AffinityPredictor
    started = time.perf_counter()
    model = Class1AffinityPredictor.load(str(args.models))
    affinity = {}
    for allele in alleles:
        frame = model.predict_to_dataframe(queries, allele=allele,
                                           include_percentile_ranks=True,
                                           include_confidence_intervals=False)
        if frame.peptide.tolist() != queries:
            raise ValueError("Predictor returned peptides out of order")
        for peptide, value, rank in zip(frame.peptide, frame.prediction, frame.prediction_percentile):
            affinity[peptide, allele] = (float(value), float(rank))
    if not all(value > 0 for value, _ in affinity.values()):
        raise ValueError("Non-positive affinity returned")

    args.output.mkdir(parents=True, exist_ok=True)
    records = []
    for pair in pairs:
        for allele in alleles:
            mutant_nm, mutant_rank = affinity[pair["peptide_mutant"], allele]
            wildtype_nm, wildtype_rank = affinity[pair["peptide_wildtype"], allele]
            records.append(dict(
                region_id=pair["region_id"], ag_label_in_design_2=pair["ag_label_in_design_2"],
                gene=pair["gene"], variant_id=pair["variant_id"],
                protein_change=pair["protein_change"], allele=allele, length=pair["length"],
                offset_in_window=pair["offset_in_window"],
                peptide_mutant=pair["peptide_mutant"], peptide_wildtype=pair["peptide_wildtype"],
                mutant_affinity_nm=mutant_nm, wildtype_affinity_nm=wildtype_nm,
                mutant_affinity_percentile=mutant_rank, wildtype_affinity_percentile=wildtype_rank,
                agretopicity_mt_over_wt=mutant_nm / wildtype_nm,
                wildtype_source="reference_residue_restored"))

    with (args.output / "agretopicity.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)

    summary = []
    for region in regions:
        mine = [r for r in records if r["region_id"] == region["region_id"]]
        best = min(mine, key=lambda r: r["mutant_affinity_nm"])
        ratios = sorted(r["agretopicity_mt_over_wt"] for r in mine)
        summary.append(dict(
            region_id=region["region_id"], ag_label_in_design_2=region["ag_label_in_design_2"],
            gene=region["gene"], protein_change=region["protein_change"],
            peptide_allele_pairs=len(mine),
            best_mutant_peptide=best["peptide_mutant"], best_allele=best["allele"],
            best_mutant_affinity_nm=best["mutant_affinity_nm"],
            wildtype_counterpart_nm=best["wildtype_affinity_nm"],
            agretopicity_of_best_binder=best["agretopicity_mt_over_wt"],
            lowest_agretopicity=ratios[0], median_agretopicity=ratios[len(ratios) // 2],
            recorded_best_affinity_nm=region["recorded_best_affinity_nm"]))
    summary.sort(key=lambda r: r["agretopicity_of_best_binder"])
    with (args.output / "region-summary.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)

    recorded_path = args.repo / "designs/design-2/inputs/predictions/mutations/best-peptide-per-HLA.tsv"
    with recorded_path.open() as handle:
        recorded = list(csv.DictReader(handle, delimiter="\t"))
    variants = {r["variant_id"] for r in records}
    compared, differences, outside = 0, [], 0
    for row in recorded:
        if row["variant_id"] not in variants:
            continue
        key = (row["peptide"], row["allele"])
        if key not in affinity:
            outside += 1
            continue
        differences.append(abs(affinity[key][0] / float(row["affinity_nm"]) - 1))
        compared += 1

    (args.output / "agretopicity-run.json").write_text(json.dumps(dict(
        scope="Fresh class I affinity for the mutant and wild-type forms of every 8 to 11mer "
              "covering the substituted residue in the selected substitution regions",
        direction="LENS convention, mutant nM divided by wild-type nM, so lower is more distinct from self",
        regions_in_scope=len(regions), regions_skipped=skipped,
        wildtype_source="reference residue restored from the DNA annotation protein change; "
                        "no BLAST search against an external database",
        alleles=alleles, peptide_lengths=list(LENGTHS),
        unique_peptides_predicted=len(queries), peptide_allele_predictions=len(queries) * len(alleles),
        mutant_wildtype_pairs=len(pairs), rows=len(records),
        threshold_applied=None, used_to_rank_or_filter=False,
        new_inference=True, cached_predictions_relabelled=False,
        models_match_pinned_bundle=models_match_pinned_bundle,
        model_files=len(files), model_manifest=str(args.model_manifest.relative_to(args.repo)),
        recorded_affinity_parity=dict(
            table=str(recorded_path.relative_to(args.repo)), rows_compared=compared,
            rows_whose_peptide_is_outside_this_scan=outside,
            max_relative_affinity_difference=max(differences) if differences else None),
        elapsed_seconds=time.perf_counter() - started,
        versions={name: importlib.metadata.version(name) for name in ["mhcflurry", "numpy", "torch"]},
    ), indent=2) + "\n")
    print("WROTE", args.output, len(records), "rows from", len(queries), "peptides;",
          "parity rows", compared, "max relative difference",
          max(differences) if differences else None, flush=True)


if __name__ == "__main__":
    main()
