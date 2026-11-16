"""Re-rank the Design #2 candidate peptides by presentation instead of affinity.

Design #2 ranks and admits class I candidates on predicted binding affinity.
MHCflurry's presentation score combines that affinity with an antigen processing
score through a two-feature logistic regression whose coefficients ship inside
the same model bundle this repository already pins. This script computes the
processing and presentation columns for the recorded candidate peptides and
measures how much the ranking moves.

It changes no selection. It adds columns and reports rank movement.

Run from the repository root:

    python designs/presentation-reranking/run.py --models "$(mhcflurry downloads path \
        models_class1_presentation --release 2.3.0)/models"
"""
import argparse, csv, hashlib, importlib.metadata, json, re, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
LENGTHS = (8, 9, 10, 11)
FLANK = 5
GATE_NM = 500.0


def spearman(first, second):
    """Spearman rank correlation with average ranks for ties, numpy only."""
    import numpy as np

    def ranks(values):
        values = np.asarray(values, dtype=float)
        order = values.argsort(kind="mergesort")
        out = np.empty(len(values), dtype=float)
        position = 0
        while position < len(order):
            end = position
            while end + 1 < len(order) and values[order[end + 1]] == values[order[position]]:
                end += 1
            out[order[position:end + 1]] = 0.5 * (position + end) + 1.0
            position = end + 1
        return out

    a, b = ranks(first), ranks(second)
    a, b = a - a.mean(), b - b.mean()
    denominator = float(np.sqrt((a * a).sum() * (b * b).sum()))
    return float((a * b).sum() / denominator) if denominator else float("nan")


def protein_contexts(repo):
    """Isovar reconstructed protein contexts, which are longer than a 25mer window."""
    path = repo / "designs/design-2/inputs/rna/primary-proteins.tsv"
    with path.open() as handle:
        return [row for row in csv.DictReader(handle, delimiter="\t")]


def candidates(repo):
    """Every recorded candidate peptide, with the window and context it came from.

    CTA/SELF, SPLICE and FUSION rows carry their scored peptides in the
    target_peptides column, so those lists are used as recorded. MUTATION rows do
    not, so their peptides are enumerated as every 8 to 11mer covering the
    substituted residue. The parity check at the end confirms that enumeration
    reproduces the recorded per-allele affinities.
    """
    pool = repo / "designs/design-2/inputs/integrated/candidate-pool.tsv"
    with pool.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    contexts = {row["variant_id"]: row for row in protein_contexts(repo)}

    entries, windows, skipped_windows = [], [], []
    for row in rows:
        window = row["sequence"]
        context, context_offset, context_kind = window, 0, "window_only"
        if row["source"] == "MUTATION":
            record = contexts.get(row["variant_id"])
            if record and window in record["sequence"]:
                context = record["sequence"]
                context_offset = context.index(window)
                context_kind = "isovar_protein_context"
            peptides = []
            match = re.fullmatch(r"p\.([A-Z])(\d+)([A-Z])", row["protein_change"])
            index = int(row["target_start_in_window"]) if row["target_start_in_window"] else None
            if index is None:
                raise ValueError("Mutation row without a target position: " + row["region_id"])
            for length in LENGTHS:
                for start in range(0, len(window) - length + 1):
                    if start <= index < start + length:
                        peptides.append((window[start:start + length], start))
            enumerated = True
        else:
            recorded = [p for p in row["target_peptides"].split(";") if p]
            peptides, enumerated = [], False
            for peptide in recorded:
                start = window.find(peptide)
                if start < 0:
                    raise ValueError("Recorded peptide is not inside its window: " + peptide)
                peptides.append((peptide, start))
            if not recorded:
                # The in-frame driver junction row (AG24, CDH1::NRG1) carries no
                # target span and no peptide list, unlike every other row in the
                # pool. Scoring its whole window would include peptides lying
                # entirely on one side of the junction, which are self sequence
                # and not neoantigenic, so the region is left out rather than
                # scored on a span this file does not state. See the README.
                skipped_windows.append(dict(region_id=row["region_id"], source=row["source"],
                                            gene=row["gene"], window=window,
                                            recorded_best_affinity_nm=row["best_target_affinity_nm"],
                                            reason="pool row records no target span and no target "
                                                   "peptides, so junction-covering peptides cannot "
                                                   "be identified from it"))
                continue
        windows.append(dict(region_id=row["region_id"], source=row["source"], gene=row["gene"],
                            window=window, context_kind=context_kind,
                            peptides=len(peptides), peptides_enumerated=enumerated,
                            recorded_best_affinity_nm=row["best_target_affinity_nm"]))
        for peptide, start in peptides:
            position = context_offset + start
            n_flank = context[max(0, position - FLANK):position]
            c_flank = context[position + len(peptide):position + len(peptide) + FLANK]
            entries.append(dict(
                region_id=row["region_id"], source=row["source"], gene=row["gene"],
                variant_id=row["variant_id"], peptide=peptide, length=len(peptide),
                offset_in_window=start, flank_context=context_kind,
                n_flank=n_flank, c_flank=c_flank,
                full_5_5_flanks=len(n_flank) == FLANK and len(c_flank) == FLANK))
    return entries, windows, skipped_windows


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--repo", type=Path, default=REPO)
    parser.add_argument("--models", type=Path, required=True,
                        help="models directory of models_class1_presentation 2.3.0")
    parser.add_argument("--output", type=Path, default=HERE)
    parser.add_argument("--model-manifest", type=Path,
                        default=REPO / "designs/initial/inputs/checkpoints/provenance/construct-model-files.json")
    args = parser.parse_args()

    alleles = json.loads((args.repo / "config/case.json").read_text())["hla_class_I"]
    entries, windows, skipped_windows = candidates(args.repo)
    by_peptide = {}
    for entry in entries:
        by_peptide.setdefault(entry["peptide"], entry)
    peptides = sorted(by_peptide)

    files = [dict(path=str(path.relative_to(args.models)), bytes=path.stat().st_size,
                  sha256=hashlib.sha256(path.read_bytes()).hexdigest())
             for path in sorted(args.models.rglob("*")) if path.is_file()]
    models_match_pinned_bundle = files == json.loads(args.model_manifest.read_text())

    import numpy as np
    from mhcflurry import Class1PresentationPredictor
    started = time.perf_counter()
    predictor = Class1PresentationPredictor.load(str(args.models))
    if not predictor.supports_presentation_prediction:
        raise ValueError("Bundle cannot predict presentation")
    lengths = {k: int(v) for k, v in predictor.processing_predictor_with_flanks.sequence_lengths.items()}
    if lengths != dict(peptide=15, n_flank=5, c_flank=5):
        raise ValueError("Unexpected flank lengths in this bundle: " + repr(lengths))
    if predictor.processing_predictor_with_flanks.percent_rank_transform is not None:
        raise ValueError("This bundle has a processing percentile calibration; revisit the report")
    weights = {row: {column: float(value) for column, value in record.items()}
               for row, record in predictor.weights_dataframe.to_dict("index").items()}

    flanked = [p for p in peptides if by_peptide[p]["n_flank"] and by_peptide[p]["c_flank"]]
    bare = [p for p in peptides if p not in set(flanked)]

    scored = {}
    for allele in alleles:
        for group, use_flanks in ((flanked, True), (bare, False)):
            if not group:
                continue
            arguments = dict(peptides=group, alleles=[allele],
                             include_affinity_percentile=True, verbose=0)
            if use_flanks:
                arguments.update(n_flanks=[by_peptide[p]["n_flank"] for p in group],
                                 c_flanks=[by_peptide[p]["c_flank"] for p in group])
            frame = predictor.predict(**arguments)
            for record in frame.to_dict("records"):
                scored[record["peptide"], allele] = dict(
                    affinity_nm=float(record["affinity"]),
                    affinity_percentile=float(record["affinity_percentile"]),
                    processing_score=float(record["processing_score"]),
                    presentation_score=float(record["presentation_score"]),
                    presentation_percentile=float(record["presentation_percentile"]),
                    coefficient_row="with_flanks" if use_flanks else "without_flanks")
    # Sensitivity: the same peptides scored with no flanks at all, to measure how
    # much real context moves the processing and presentation columns.
    no_flank = {}
    for allele in alleles:
        frame = predictor.predict(peptides=peptides, alleles=[allele],
                                  include_affinity_percentile=True, verbose=0)
        for record in frame.to_dict("records"):
            no_flank[record["peptide"], allele] = (float(record["processing_score"]),
                                                   float(record["presentation_score"]))

    args.output.mkdir(parents=True, exist_ok=True)
    rows = []
    for entry in entries:
        for allele in alleles:
            value = scored[entry["peptide"], allele]
            without = no_flank[entry["peptide"], allele]
            rows.append(dict(
                region_id=entry["region_id"], source=entry["source"], gene=entry["gene"],
                peptide=entry["peptide"], length=entry["length"], allele=allele,
                offset_in_window=entry["offset_in_window"],
                flank_context=entry["flank_context"],
                n_flank=entry["n_flank"], c_flank=entry["c_flank"],
                full_5_5_flanks=entry["full_5_5_flanks"],
                coefficient_row=value["coefficient_row"],
                affinity_nm=value["affinity_nm"], affinity_percentile=value["affinity_percentile"],
                processing_score=value["processing_score"],
                presentation_score=value["presentation_score"],
                presentation_percentile=value["presentation_percentile"],
                processing_score_no_flanks=without[0],
                presentation_score_no_flanks=without[1]))
    with (args.output / "presentation.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    # Peptide level, best allele by each quantity.
    best = {}
    for row in rows:
        current = best.setdefault(row["peptide"], dict(peptide=row["peptide"], source=row["source"],
                                                       full_5_5_flanks=row["full_5_5_flanks"],
                                                       flank_context=row["flank_context"],
                                                       affinity_nm=row["affinity_nm"],
                                                       affinity_allele=row["allele"],
                                                       processing_score=row["processing_score"],
                                                       presentation_score=row["presentation_score"],
                                                       presentation_allele=row["allele"]))
        if row["affinity_nm"] < current["affinity_nm"]:
            current.update(affinity_nm=row["affinity_nm"], affinity_allele=row["allele"])
        if row["presentation_score"] > current["presentation_score"]:
            current.update(presentation_score=row["presentation_score"],
                           presentation_allele=row["allele"],
                           processing_score=row["processing_score"])
    best = list(best.values())
    admitted = [b for b in best if b["affinity_nm"] <= GATE_NM]
    rejected = [b for b in best if b["affinity_nm"] > GATE_NM]

    rho = spearman([b["affinity_nm"] for b in admitted],
                   [-b["presentation_score"] for b in admitted])
    top_affinity = {b["peptide"] for b in sorted(admitted, key=lambda b: b["affinity_nm"])[:50]}
    top_presentation = {b["peptide"] for b in
                        sorted(admitted, key=lambda b: -b["presentation_score"])[:50]}
    median_admitted = float(np.median([b["presentation_score"] for b in admitted]))
    rejected_above_median = sum(1 for b in rejected if b["presentation_score"] > median_admitted)

    # Region level: rank regions by their best peptide under each quantity. Built
    # from the per-row table rather than the per-peptide table, because the same
    # peptide can occur in more than one region and deduplicating by peptide
    # would attribute it to whichever region happened to come first.
    regions = {}
    for row in rows:
        current = regions.setdefault(row["region_id"], dict(
            region_id=row["region_id"], source=row["source"], gene=row["gene"],
            best_affinity_nm=row["affinity_nm"], best_presentation_score=row["presentation_score"],
            best_affinity_peptide=row["peptide"], best_presentation_peptide=row["peptide"],
            best_affinity_allele=row["allele"], best_presentation_allele=row["allele"]))
        if row["affinity_nm"] < current["best_affinity_nm"]:
            current.update(best_affinity_nm=row["affinity_nm"], best_affinity_peptide=row["peptide"],
                           best_affinity_allele=row["allele"])
        if row["presentation_score"] > current["best_presentation_score"]:
            current.update(best_presentation_score=row["presentation_score"],
                           best_presentation_peptide=row["peptide"],
                           best_presentation_allele=row["allele"])
    regions = list(regions.values())
    by_affinity = sorted(regions, key=lambda r: r["best_affinity_nm"])
    by_presentation = sorted(regions, key=lambda r: -r["best_presentation_score"])
    affinity_rank = {r["region_id"]: i + 1 for i, r in enumerate(by_affinity)}
    presentation_rank = {r["region_id"]: i + 1 for i, r in enumerate(by_presentation)}
    for record in regions:
        record["affinity_rank"] = affinity_rank[record["region_id"]]
        record["presentation_rank"] = presentation_rank[record["region_id"]]
        record["rank_change"] = record["affinity_rank"] - record["presentation_rank"]
    regions.sort(key=lambda r: r["affinity_rank"])
    with (args.output / "region-rank-comparison.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(regions[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(regions)

    region_rho = spearman([r["affinity_rank"] for r in regions],
                          [r["presentation_rank"] for r in regions])

    # Parity against recorded affinities, mutation and fusion tables.
    parity = {}
    for name, path, peptide_key, allele_key, affinity_key in [
        ("mutations", "designs/design-2/inputs/predictions/mutations/best-peptide-per-HLA.tsv",
         "peptide", "allele", "affinity_nm"),
        ("fusions", "designs/design-2/inputs/predictions/fusions/affinities.tsv",
         "peptide", "allele", "affinity")]:
        with (args.repo / path).open() as handle:
            recorded = list(csv.DictReader(handle, delimiter="\t"))
        differences, outside = [], 0
        for record in recorded:
            key = (record[peptide_key], record[allele_key])
            if key not in scored:
                outside += 1
                continue
            differences.append(abs(scored[key]["affinity_nm"] / float(record[affinity_key]) - 1))
        parity[name] = dict(table=path, rows_compared=len(differences),
                            rows_outside_this_scan=outside,
                            max_relative_affinity_difference=max(differences) if differences else None)

    flank_moves = [abs(scored[p, a]["presentation_score"] - no_flank[p, a][1])
                   for p in flanked for a in alleles]
    processing_moves = [abs(scored[p, a]["processing_score"] - no_flank[p, a][0])
                        for p in flanked for a in alleles]

    summary = dict(
        scope="Processing and presentation columns for the recorded Design #2 candidate peptides, "
              "with the ranking compared against affinity. No selection is changed.",
        candidate_windows=len(windows), windows_skipped=skipped_windows, windows_by_source={
            source: sum(1 for w in windows if w["source"] == source)
            for source in sorted({w["source"] for w in windows})},
        peptides=len(peptides), peptide_allele_rows=len(rows), alleles=alleles,
        admitted_at_500nm=len(admitted), rejected_at_500nm=len(rejected),
        flanks=dict(
            peptides_with_both_flanks=len(flanked), peptides_scored_without_flanks=len(bare),
            peptides_with_full_5_5=sum(1 for p in peptides if by_peptide[p]["full_5_5_flanks"]),
            full_5_5_fraction=round(sum(1 for p in peptides if by_peptide[p]["full_5_5_flanks"]) / len(peptides), 4),
            flank_context_counts={kind: sum(1 for p in peptides if by_peptide[p]["flank_context"] == kind)
                                  for kind in sorted({e["flank_context"] for e in entries})},
            coefficient_rows_used=sorted({v["coefficient_row"] for v in scored.values()}),
            shipped_coefficients=weights,
            measured_flank_lengths=lengths),
        ranking=dict(
            spearman_affinity_vs_presentation_admitted_peptides=round(rho, 4),
            top50_overlap_of_admitted_peptides=len(top_affinity & top_presentation),
            admitted_peptides_with_processing_below_0_1=sum(
                1 for b in admitted if b["processing_score"] < 0.1),
            median_admitted_presentation_score=median_admitted,
            rejected_peptides_above_median_admitted_presentation=rejected_above_median,
            rejected_peptides_total=len(rejected),
            region_spearman=round(region_rho, 4),
            regions=len(regions),
            regions_moving_more_than_20_places=sum(1 for r in regions if abs(r["rank_change"]) > 20),
            regions_moving_more_than_5_places=sum(1 for r in regions if abs(r["rank_change"]) > 5)),
        flank_sensitivity=dict(
            peptides_compared=len(flanked),
            median_absolute_presentation_change=float(np.median(flank_moves)) if flank_moves else None,
            max_absolute_presentation_change=float(np.max(flank_moves)) if flank_moves else None,
            median_absolute_processing_change=float(np.median(processing_moves)) if processing_moves else None,
            max_absolute_processing_change=float(np.max(processing_moves)) if processing_moves else None),
        threshold_applied=None, selection_changed=False,
        new_inference=True, cached_predictions_relabelled=False,
        models_match_pinned_bundle=models_match_pinned_bundle, model_files=len(files),
        model_manifest=str(args.model_manifest.relative_to(args.repo)),
        recorded_affinity_parity=parity,
        elapsed_seconds=time.perf_counter() - started,
        versions={name: importlib.metadata.version(name) for name in ["mhcflurry", "numpy", "torch"]})
    (args.output / "presentation-run.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("WROTE", args.output, len(rows), "rows;", len(peptides), "peptides;",
          "admitted", len(admitted), "rho", round(rho, 4), "top50 overlap",
          len(top_affinity & top_presentation), flush=True)


if __name__ == "__main__":
    main()
