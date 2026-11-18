#!/usr/bin/env python3
"""Corroborate every fresh-inference affinity cell against Design 4's OWN stored
values, using a named, reproducible lookup. No model inference.

Two Design 4 sources are consulted, in this order, keyed on exact (peptide, allele):

  1. `construct/all-junction-affinity.tsv.gz` — columns peptide, allele,
     affinity_nm, affinity_percentile. Direct key lookup.
  2. `safety/offtarget-homologs.tsv.gz` — the homolog's value on the shared allele
     is `shared_allele_homolog_nm`. The query's best-across-alleles value is not
     used as evidence for a particular shared allele.

A cell counts as corroborated when a source value exists for its exact key and
agrees within TOL_NM. Cells with no source key are reported as uncovered, not
silently dropped.

Usage: corroborate_cells.py --design4 <dir> --table <tier1-direction.tsv> [--out <json>]
"""
import argparse, csv, gzip, json, os
from validation import number, open_text_auto

TOL_NM = 0.05


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--design4", required=True)
    ap.add_argument("--table", required=True)
    ap.add_argument("--out")
    a = ap.parse_args()

    jp = os.path.join(a.design4, "construct/all-junction-affinity.tsv.gz")
    cache = {}
    with gzip.open(jp, "rt") as fh:
        for line, r in enumerate(csv.DictReader(fh, delimiter="\t"), 2):
            cache[(r["peptide"], r["allele"])] = number(
                r.get("affinity_nm"), f"{jp}:{line} affinity_nm", positive=True)

    hp = os.path.join(a.design4, "safety/offtarget-homologs.tsv.gz")
    hom = {}
    with gzip.open(hp, "rt") as fh:
        for line, r in enumerate(csv.DictReader(fh, delimiter="\t"), 2):
            if r["model_priority_category"] != "TIER1_1MM_PARALOG_SAME_HLA" or not r["shared_binding_allele"]:
                continue
            al = r["shared_binding_allele"]
            hom[(r["homolog_peptide"], al)] = number(
                r.get("shared_allele_homolog_nm"),
                f"{hp}:{line} shared_allele_homolog_nm", positive=True)

    with open_text_auto(a.table) as fh:
        tr = list(csv.DictReader(fh, delimiter="\t"))
    checked, ok, uncovered, maxdiff, worst = 0, 0, [], 0.0, None
    per_source = {"all-junction-affinity.tsv.gz": 0, "offtarget-homologs.tsv.gz": 0}
    for r in tr:
        al = r["shared_allele"]
        for side in ("query", "homolog"):
            pep = r[f"{side}_peptide"]
            got = number(r.get(f"mhcflurry_{side}_nm"),
                         f"{a.table} {pep} {al} {side} affinity", positive=True)
            checked += 1
            src, which = None, None
            if (pep, al) in cache:
                src, which = cache[(pep, al)], "all-junction-affinity.tsv.gz"
            elif (pep, al) in hom:
                src, which = hom[(pep, al)], "offtarget-homologs.tsv.gz"
            if src is None:
                uncovered.append((pep, al, side))
                continue
            d = abs(src - got)
            if d > maxdiff:
                maxdiff, worst = d, (pep, al, src, got)
            if d <= TOL_NM:
                ok += 1
                per_source[which] += 1
            else:
                uncovered.append((pep, al, side, f"differs by {d}"))

    res = {"cells_checked": checked, "cells_corroborated": ok,
           "cells_uncovered_or_disagreeing": len(uncovered),
           "uncovered_detail": uncovered[:20],
           "max_abs_difference_nm": round(maxdiff, 8),
           "largest_difference_cell": worst,
           "tolerance_nm": TOL_NM,
           "lookup_method": (
               "exact (peptide, allele) key; source 1 = "
               "construct/all-junction-affinity.tsv.gz, source 2 = "
               "shared_allele_homolog_nm from safety/offtarget-homologs.tsv.gz"),
           "per_source_corroborated": per_source,
           "model_inference_run": False}
    if a.out:
        with open(a.out, "w") as fh:
            json.dump(res, fh, indent=1, sort_keys=True)
    print(json.dumps({k: v for k, v in res.items() if k != "uncovered_detail"}, indent=1))


if __name__ == "__main__":
    main()
