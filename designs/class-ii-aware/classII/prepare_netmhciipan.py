#!/usr/bin/env python3
"""Turn the newly typed class II alleles into a ready-to-run NetMHCIIpan job.

Why this is a preparation script and not a prediction script: every class II
binding predictor I could find is licence-encumbered. NetMHCIIpan and the IEDB
class II tools require accepting an academic licence, and MixMHC2pred's licence
takes effect on download. Accepting one of those on the project's behalf, or
adding an academic-only dependency to a public repository that tracks
third-party notices, is not mine to decide. The maintainers almost certainly
hold a NetMHCIIpan licence already; this script does everything up to the point
where that binary is invoked.

It writes:
  alleles.txt     the patient's class II molecules in NetMHCIIpan allele syntax
  regions.fasta   Design #2's 24 antigen regions, one record each
  cassette.fasta  the full 724-aa construct, which also exposes the epitopes
                  that only exist across antigen junctions
  run.sh          the exact commands

Usage:
    python prepare_netmhciipan.py --repo-root ../../.. --out results
    cd results && sh run.sh          # needs netMHCIIpan on PATH
"""
from __future__ import annotations

import argparse, csv, json, os

# Class II molecules implied by the typing in ../hla/results/hla_calls.json.
# DR beta chains pair with the invariant DRA, so they are named alone. DQ and DP
# are alpha/beta heterodimers: the two cis (same-haplotype) pairs are the ones
# certainly expressed, and the two trans pairs do form for DQ and are included
# so they are not silently ignored.
DR = ["DRB1_0301", "DRB1_0802", "DRB3_0202"]
DQ_CIS = ["HLA-DQA10501-DQB10201", "HLA-DQA10401-DQB10402"]
DQ_TRANS = ["HLA-DQA10501-DQB10402", "HLA-DQA10401-DQB10201"]
DP_CIS = ["HLA-DPA10103-DPB10201", "HLA-DPA10201-DPB11301"]
DP_TRANS = ["HLA-DPA10103-DPB11301", "HLA-DPA10201-DPB10201"]

HAPLOTYPE_NOTE = """\
Haplotype 1 : DRB1*03:01 - DRB3*02:02 - DQA1*05:01 - DQB1*02:01   (DR3-DQ2, DR52)
Haplotype 2 : DRB1*08:02 -      none  - DQA1*04:01 - DQB1*04:02   (DR8-DQ4)
DRB4 and DRB5 absent, which is what this pair of haplotypes predicts.
"""


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", default="../../..")
    ap.add_argument("--out", default="job")
    ap.add_argument("--include-trans", action="store_true", default=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    d2 = os.path.join(args.repo_root, "designs", "design-2")

    alleles = DR + DQ_CIS + DP_CIS + (DQ_TRANS + DP_TRANS if args.include_trans else [])
    with open(os.path.join(args.out, "alleles.txt"), "w") as fh:
        fh.write("# " + HAPLOTYPE_NOTE.replace("\n", "\n# ").rstrip("# ") + "\n")
        fh.write("\n".join(alleles) + "\n")

    order = read_tsv(os.path.join(d2, "construct/protein/antigen-order.tsv"))
    with open(os.path.join(args.out, "regions.fasta"), "w") as fh:
        for r in order:
            fh.write(f">{r['antigen_id']}|{r['gene']}|pos{r['position']}\n{r['sequence']}\n")

    cassette = []
    with open(os.path.join(d2, "construct/protein/cassette.fasta")) as fh:
        for line in fh:
            if not line.startswith(">"):
                cassette.append(line.strip())
    cassette = "".join(cassette)
    with open(os.path.join(args.out, "cassette.fasta"), "w") as fh:
        fh.write(f">design2_cassette_{len(cassette)}aa\n")
        for i in range(0, len(cassette), 60):
            fh.write(cassette[i:i + 60] + "\n")

    with open(os.path.join(args.out, "run.sh"), "w") as fh:
        fh.write(
            "#!/bin/sh\n"
            "# Class II binding for Yash's typed alleles. Needs netMHCIIpan on PATH\n"
            "# (academic licence, obtained by the user, not bundled here).\n"
            "set -eu\n"
            f"ALLELES=$(grep -v '^#' alleles.txt | paste -sd, -)\n\n"
            "# Per-region: does each 25-aa antigen block contain a CD4 epitope at all?\n"
            "netMHCIIpan -f regions.fasta -inptype 0 -a \"$ALLELES\" -length 15 \\\n"
            "    -BA -xls -xlsfile regions.classII.xls > regions.classII.txt\n\n"
            "# Whole cassette: also covers epitopes formed across antigen junctions,\n"
            "# which per-region scoring cannot see.\n"
            "netMHCIIpan -f cassette.fasta -inptype 0 -a \"$ALLELES\" -length 15 \\\n"
            "    -BA -xls -xlsfile cassette.classII.xls > cassette.classII.txt\n"
        )

    meta = {"alleles": alleles, "n_regions": len(order), "cassette_aa": len(cassette),
            "haplotypes": HAPLOTYPE_NOTE.strip().splitlines()}
    with open(os.path.join(args.out, "job.json"), "w") as fh:
        json.dump(meta, fh, indent=2)

    print(HAPLOTYPE_NOTE)
    print(f"{len(alleles)} class II molecules:")
    for a in alleles:
        print("   ", a)
    print(f"\n{len(order)} regions and a {len(cassette)}-aa cassette written to {args.out}/")
    print(f"Run: cd {args.out} && sh run.sh   (requires netMHCIIpan)")


if __name__ == "__main__":
    main()
