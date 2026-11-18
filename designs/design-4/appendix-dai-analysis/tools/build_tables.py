#!/usr/bin/env python3
"""Regenerate the three DAI tables and their receipt from Design 4's stored files.

Scope: rows cover MUTATION-SPANNING single-amino-acid wild-type comparisons only
(Design 4's `single_aa_WT_reversion` scenario). They are NOT every peptide in an
encoded window. Flanking peptides that do not span the mutation are excluded because
their wild-type counterpart is IDENTICAL to them, so the ratio is trivially 1 and
carries no information — not because a counterpart is missing.

Reads the supplied Design 4 inputs and writes derived tables to --out.
Runs no model inference.
Usage: build_tables.py --design4 <path to designs/design-4> --out <data dir>
"""
import argparse, csv, json, math, os
from collections import defaultdict
from validation import open_text_auto, residual_low_dai

GUARDANT_COVERED = {"RXRA", "TGFBR2", "CDH1::NRG1"}
"""Curated historical annotation, not derived from the five prediction inputs.
These three encoded alterations were detected in the 14 October 2025 cfDNA draw.
See PROVENANCE.md and the pinned Guardant360 report for its source and limits.
"""


def tsv(p):
    with open_text_auto(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--design4", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    sel = tsv(os.path.join(a.design4, "selection/selected.tsv"))
    cf = tsv(os.path.join(a.design4, "immunogenicity/class1-counterfactual-contrasts.tsv"))
    if not sel:
        raise ValueError("selection/selected.tsv contains no selected regions")
    if not cf:
        raise ValueError("class1-counterfactual-contrasts.tsv contains no comparison rows")

    window = {r["gene"]: r["sequence"] for r in sel}
    lead = {r["gene"]: (r["best_target_peptide"], r["best_target_allele"]) for r in sel}

    # per peptide-allele DAI, restricted to peptides inside the encoded window
    rows, bygene, low = [], defaultdict(list), []
    for r in cf:
        if r["scenario"] != "single_aa_WT_reversion":
            continue
        g = r["gene"]
        if g not in window or r["mutant_peptide"] not in window[g]:
            continue
        m = float(r["mutant_affinity_nm"])
        w = float(r["counterfactual_affinity_nm"])
        if not (math.isfinite(m) and math.isfinite(w) and m > 0 and w > 0):
            raise ValueError(f"invalid affinity in source row for {g}")
        d = w / m
        if not math.isfinite(d):
            raise ValueError(f"invalid DAI in source row for {g}")
        is_lead = (r["mutant_peptide"], r["allele"]) == lead[g]
        rows.append(dict(gene=g, mutant_peptide=r["mutant_peptide"],
                         wildtype_peptide=r["counterfactual_peptide"], allele=r["allele"],
                         mutant_affinity_nm=round(m, 4), wildtype_affinity_nm=round(w, 4),
                         dai_wt_over_mut=round(d, 4),
                         is_design4_lead=is_lead,
                         binds_le_500nM=m <= 500.0,
                         encoded_in_design4_window=True))
        bygene[g].append((d, m, r["mutant_peptide"], r["allele"]))
        if residual_low_dai(d, m):
            low.append((g, r["mutant_peptide"], r["allele"], round(m, 4), round(d, 4)))

    if not rows:
        raise ValueError("no encoded single-AA WT comparison rows were found")

    with open(os.path.join(a.out, "encoded-peptide-dai.tsv"), "w", newline="") as fh:
        w_ = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w_.writeheader(); w_.writerows(sorted(rows, key=lambda x: (x["gene"], x["dai_wt_over_mut"])))

    # per-region summary
    with open(os.path.join(a.out, "region-dai-summary.tsv"), "w", newline="") as fh:
        w_ = csv.writer(fh, delimiter="\t")
        w_.writerow(["gene", "window_aa", "design4_lead", "design4_lead_allele", "design4_lead_nm",
                     "design4_lead_dai", "best_dai_in_window_le_500nM",
                     "peptide_achieving_best_dai", "allele_achieving_best_dai",
                     "higher_dai_choice_exists", "choice_kind",
                     "n_mutation_spanning_pep_allele_dai_le_1_and_le500nM",
                     "guardant360_covered", "dai_defined"])
        for r in sel:
            g = r["gene"]; c = bygene.get(g, [])
            lp, la = lead[g]
            ld = next((d for d, m, p, al in c if p == lp and al == la), "")
            ok = [x for x in c if x[1] <= 500]
            best = max(ok, key=lambda x: x[0]) if ok else None
            nlow = sum(residual_low_dai(d, m) for d, m, p, al in c)
            w_.writerow([g, len(window[g]), lp, la, r["best_target_affinity_nm"],
                         round(ld, 4) if ld != "" else "",
                         round(best[0], 4) if best else "", best[2] if best else "",
                         best[3] if best else "",
                         bool(best and ld != "" and best[0] > ld + 1e-9),
                         ("" if not (best and ld != "" and best[0] > ld + 1e-9)
                          else ("allele_only" if best[2] == lp else "peptide_change")),
                         nlow, g in GUARDANT_COVERED, bool(c)])

    # the residual-peptide finding
    low = sorted(set(low), key=lambda x: x[4])
    with open(os.path.join(a.out, "residual-low-dai-peptides.tsv"), "w", newline="") as fh:
        w_ = csv.writer(fh, delimiter="\t")
        w_.writerow(["gene", "peptide", "allele", "mutant_affinity_nm", "dai",
                     "encoded_window_aa", "note"])
        for g, p, al, nm, d in low:
            w_.writerow([g, p, al, nm, d, len(window[g]),
                         f"inside the encoded {len(window[g])}-aa window; relabelling a "
                         "region's lead does NOT remove this peptide from the construct"])

    receipt = {"regions_examined": len(sel),
               "peptide_allele_rows": len(rows),
               "residual_low_dai_pep_allele_pairs": len(low),
               "source": "designs/design-4 (unmodified)",
               "model_inference_run": False}
    with open(os.path.join(a.out, "receipt.json"), "w") as fh:
        json.dump(receipt, fh, indent=1, sort_keys=True)
    print(f"{len(rows)} peptide-allele rows, {len(low)} residual low-DAI pairs, "
          f"{len(sel)} regions")


if __name__ == "__main__":
    main()
