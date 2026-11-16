"""HLA-robustness analysis and conditional-candidate discrimination.

1. Per-region presenting-allele dependence: for each region's target peptides
   (overlapping the mutation codon, splice junction, or the full CTA/fusion
   window), the eligible alleles at 50/150/500/1000 nM and single-allele flags.
2. Lost-allele sensitivity: for all 64 hypothetical allele-loss subsets, how
   many regions retain an eligible prediction at each threshold.
3. H2BC12 E114K and FES R804W mutant-versus-focal-reference discrimination.
"""
import csv
import itertools
import sys
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
from yash_vaccine.io import write_tsv  # noqa: E402
THRESHOLDS = [50, 150, 500, 1000]


def overlapping(sequence, lo, hi, lengths=(8, 9, 10, 11)):
    return sorted({sequence[s:s + k] for k in lengths for s in range(len(sequence) - k + 1)
                   if s < hi and s + k > lo})


def main():
    from mhcflurry import Class1AffinityPredictor
    predictor = Class1AffinityPredictor.load(
        "/root/.local/share/mhcflurry/4/2.3.0/models_class1_pan/models.combined")
    alleles = json.loads((ROOT / "config/case.json").read_text())["hla_class_I"]
    regions = list(csv.DictReader((HERE / "regions.tsv").open(), delimiter="\t"))

    # ---- per-region allele dependence (fresh predictions on target peptides)
    rows = []
    per_region_peptides = {}
    for row in regions:
        lo, hi = int(row["target_start_in_window"]), int(row["target_end_in_window"])
        peptides = overlapping(row["sequence"], lo, hi)
        per_region_peptides[row["antigen_id"]] = peptides
        frame = predictor.predict_to_dataframe(
            alleles=[a for a in alleles for _ in peptides],
            peptides=[p for _ in alleles for p in peptides])
        frame["allele"] = [a for a in alleles for _ in peptides]
        best = frame.groupby("allele")["prediction"].min()
        eligible = {t: [a for a in alleles if best[a] <= t] for t in THRESHOLDS}
        rows.append(dict(
            antigen_id=row["antigen_id"], gene=row["gene"], source=row["source"],
            evidence_tier=row["evidence_tier"],
            target_peptides=len(peptides),
            **{f"eligible_alleles_{t}nm": ";".join(eligible[t]) or "NONE" for t in THRESHOLDS},
            **{f"n_eligible_{t}nm": len(eligible[t]) for t in THRESHOLDS},
            single_allele_dependence_500nm=("YES" if len(eligible[500]) == 1 else
                                            "none" if len(eligible[500]) > 1 else "no_eligible"),
            best_affinity_nm=float(frame["prediction"].min())))
    write_tsv(HERE / "hla-allele-dependence.tsv", rows)

    # ---- lost-allele sensitivity over all 64 subsets
    subsets = []
    lost_options = [[], alleles]  # per-locus: keep or lose
    loci = {"A": ["HLA-A*01:01", "HLA-A*33:03"], "B": ["HLA-B*08:01", "HLA-B*15:08"],
            "C": ["HLA-C*01:02", "HLA-C*07:02"]}
    scenarios = []
    for a_choice in [[], loci["A"]]:
        for b_choice in [[], loci["B"]]:
            for c_choice in [[], loci["C"]]:
                scenarios.append(a_choice + b_choice + c_choice)
    for scenario in scenarios:
        lost = set(scenario)
        retained = {t: 0 for t in THRESHOLDS}
        for row in rows:
            for t in THRESHOLDS:
                eligible = set(filter(None, row[f"eligible_alleles_{t}nm"].split(";")))
                if eligible and eligible - lost:
                    retained[t] += 1
        subsets.append(dict(
            lost_alleles=";".join(sorted(lost)) or "none",
            n_lost=len(lost),
            **{f"regions_retained_{t}nm": retained[t] for t in THRESHOLDS}))
    write_tsv(HERE / "lost-allele-sensitivity.tsv", subsets)

    # ---- conditional candidates: mutant vs focal-reference discrimination
    conditional = []
    candidates = [
        ("H2BC12_E114K", "ELAKHAVSK", "ELAKHAVSE", "HLA-A*33:03"),
        ("H2BC12_E114K", "AVSKGTKAV", "AVSEGTKAV", "HLA-C*01:02"),
        ("FES_R804W_short_context", "AYEPGQWPSF", None, "HLA-C*07:02"),
    ]
    for name, mutant, reference, allele in candidates:
        peptides = [p for p in [mutant, reference] if p]
        frame = predictor.predict_to_dataframe(
            alleles=[allele] * len(peptides), peptides=peptides)
        aff = dict(zip(frame["peptide"], frame["prediction"]))
        conditional.append(dict(
            candidate=name, allele=allele, mutant_peptide=mutant,
            mutant_affinity_nm=aff[mutant],
            reference_peptide=reference or "",
            reference_affinity_nm=aff.get(reference, ""),
            discrimination=(f"{aff[reference] / aff[mutant]:.2f}x mutant/reference gain" if reference else
                            "reference_peptide_outside_25aa_window"),
            note="FES 15-aa context EQCWAYEPGQWPSFS; the focal reference counterpart "
                 "spans the mutated base and is not a valid 10-mer comparison" if name.startswith("FES") else ""))
    write_tsv(HERE / "conditional-candidates.tsv", conditional)
    print(json.dumps(rows, indent=1)[:1500])
    print(json.dumps(conditional, indent=1))


if __name__ == "__main__":
    main()
