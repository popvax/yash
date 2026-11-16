#!/usr/bin/env python3
"""Re-score Design #2's mutation candidates with MHCflurry's presentation model.

Design #2 ranks and admits class I regions on predicted *affinity* alone; its
methods state plainly that "class II, presentation and immunogenicity
predictions do not contribute". MHCflurry 2.x also ships a presentation
predictor that combines binding affinity with an antigen-processing model
trained on mass-spectrometry eluted ligands, which is a better predictor of
what actually reaches the cell surface than affinity by itself.

This script reuses Design #2's own evidence, so nothing is re-derived from raw
reads:

  inputs/rna/primary-proteins.tsv                 mutant protein contexts
  inputs/predictions/mutations/all-ranked-regions.tsv   their affinity ranking

For every 8-11mer overlapping the mutated residue(s) it predicts affinity,
processing and presentation with the patient's six class I alleles, using the
surrounding context as N/C flanks (the processing model needs them). It then
reports, per region, the best presentation score alongside Design #2's affinity,
and whether the affinity-chosen 25-aa window is also the presentation-chosen
one.

Usage:
    pip install 'mhcflurry==2.3.14'
    mhcflurry-downloads fetch models_class1_presentation
    python score_presentation.py --repo-root ../../.. --out results
"""
from __future__ import annotations

import argparse, csv, json, os

ALLELES = ["HLA-A*01:01", "HLA-A*33:03", "HLA-B*08:01",
           "HLA-B*15:08", "HLA-C*01:02", "HLA-C*07:02"]
PEPTIDE_LENGTHS = (8, 9, 10, 11)
WINDOW = 25
FLANK = 15          # MHCflurry uses up to 15 residues of flanking context


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def overlapping_peptides(seq, mut_start, mut_end):
    """All 8-11mers inside seq that overlap the mutated span, with flanks."""
    out = []
    for L in PEPTIDE_LENGTHS:
        for start in range(0, len(seq) - L + 1):
            end = start + L
            if end <= mut_start or start >= mut_end:
                continue                      # does not touch the mutation
            out.append({
                "peptide": seq[start:end],
                "start": start,
                "end": end,
                "n_flank": seq[max(0, start - FLANK):start],
                "c_flank": seq[end:end + FLANK],
            })
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", default="../../..")
    ap.add_argument("--out", default=".")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    d2 = os.path.join(args.repo_root, "designs", "design-2")

    contexts = read_tsv(os.path.join(d2, "inputs/rna/primary-proteins.tsv"))
    ranked = {r["gene"]: r for r in
              read_tsv(os.path.join(d2, "inputs/predictions/mutations/all-ranked-regions.tsv"))}
    selected = {r["gene"] for r in
                read_tsv(os.path.join(d2, "selection/selected.tsv")) if r.get("variant_id")}

    from mhcflurry import Class1PresentationPredictor
    predictor = Class1PresentationPredictor.load()

    rows, per_peptide = [], []
    for ctx in contexts:
        seq = ctx["sequence"]
        ms, me = int(ctx["mutation_start_0based"]), int(ctx["mutation_end_0based_exclusive"])
        peps = overlapping_peptides(seq, ms, me)
        if not peps:
            continue
        pred = predictor.predict(
            peptides=[p["peptide"] for p in peps],
            n_flanks=[p["n_flank"] for p in peps],
            c_flanks=[p["c_flank"] for p in peps],
            alleles={"patient": ALLELES},
            include_affinity_percentile=True,
            verbose=0,
        )
        pred = pred.reset_index(drop=True)
        for i, p in enumerate(peps):
            r = pred.iloc[i]
            per_peptide.append({
                "gene": ctx["gene"],
                "protein_change": ctx["DNA_annotation_protein_change"],
                "peptide": p["peptide"], "start": p["start"], "end": p["end"],
                "best_allele": r["best_allele"],
                "affinity_nm": round(float(r["affinity"]), 3),
                "processing_score": round(float(r["processing_score"]), 5),
                "presentation_score": round(float(r["presentation_score"]), 5),
            })

        best = max(per_peptide[-len(peps):], key=lambda d: d["presentation_score"])
        best_aff = min(per_peptide[-len(peps):], key=lambda d: d["affinity_nm"])

        # Which 25-aa window maximises summed presentation? Compare with the
        # window Design #2 chose on summed affinity contribution.
        win_scores = {}
        if len(seq) >= WINDOW:
            for w in range(0, len(seq) - WINDOW + 1):
                inside = [d for d in per_peptide[-len(peps):]
                          if d["start"] >= w and d["end"] <= w + WINDOW]
                if not inside:
                    continue
                win_scores[w] = sum(d["presentation_score"] for d in inside)
        best_window = max(win_scores, key=win_scores.get) if win_scores else None
        d2row = ranked.get(ctx["gene"])
        d2_window_seq = d2row["sequence"] if d2row else ""
        pres_window_seq = seq[best_window:best_window + WINDOW] if best_window is not None else ""

        rows.append({
            "gene": ctx["gene"],
            "protein_change": ctx["DNA_annotation_protein_change"],
            "in_design2": ctx["gene"] in selected,
            "d2_best_affinity_nm": round(float(d2row["best_target_affinity_nm"]), 2) if d2row else None,
            "d2_passes_500nm": d2row["passes_500nm_gate"] if d2row else "",
            "d2_primary_rank": d2row["primary_rank"] if d2row else "",
            "our_best_affinity_nm": best_aff["affinity_nm"],
            "best_presentation_score": best["presentation_score"],
            "best_presentation_peptide": best["peptide"],
            "best_presentation_allele": best["best_allele"],
            "best_processing_score": round(max(d["processing_score"] for d in per_peptide[-len(peps):]), 5),
            "presentation_window_matches_d2": (pres_window_seq == d2_window_seq) if d2_window_seq and pres_window_seq else None,
            "presentation_window": pres_window_seq,
            "d2_window": d2_window_seq,
        })

    rows.sort(key=lambda r: -r["best_presentation_score"])
    for i, r in enumerate(rows, 1):
        r["presentation_rank"] = i

    with open(os.path.join(args.out, "region-presentation.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    with open(os.path.join(args.out, "peptide-presentation.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(per_peptide[0]), delimiter="\t")
        w.writeheader(); w.writerows(per_peptide)

    print(f"{'rank':>4} {'gene':10} {'change':13} {'inD2':>5} "
          f"{'D2_nM':>9} {'500nM':>6} {'D2rank':>7} {'pres':>7} {'proc':>7} {'window=':>8}")
    for r in rows:
        print(f"{r['presentation_rank']:>4} {r['gene']:10} {r['protein_change']:13} "
              f"{str(r['in_design2']):>5} {r['d2_best_affinity_nm']!s:>9} "
              f"{r['d2_passes_500nm']:>6} {r['d2_primary_rank'] or '-':>7} "
              f"{r['best_presentation_score']:>7.4f} {r['best_processing_score']:>7.4f} "
              f"{str(r['presentation_window_matches_d2']):>8}")
    print(f"\nwrote {args.out}/region-presentation.tsv and peptide-presentation.tsv "
          f"({len(per_peptide)} peptide-level rows)")


if __name__ == "__main__":
    main()
