#!/usr/bin/env python3
"""Held-out validation of the genotyper on two public donors.

Recovering the patient's own class I type is a weak check: the pair-selection
tie-breaks were settled while looking at that answer, so it is in-sample. These
two donors are not. Neither was used to choose any parameter except where noted.

NA12878 / HG001  (GIAB, 30x downsample, GRCh37, contig "6")
    Truth: Cellosaurus CVCL_7526, from PubMed 27792722 (Dilthey et al., HLA*PRG).
    Covers class I *and* class II, mostly at 1-field resolution.
    This sample was used to diagnose and fix a defect (see below), so after that
    fix it is a development set, not a clean test.

HG002 / NA24385  (GIAB, 60x, GRCh38, contig "chr6")
    Truth: clinical class I typing reported for NA24385 in the HLAminer paper
    (arXiv 2209.09155): A*01:01/26:01, B*35:08/38:01, C*04:01/12:03.
    Class I only, but at full 2-field resolution, and untouched by any tuning.

Defect this found
-----------------
With the original admission margin of 0.05, a divergent second allele that
loses reads to mismapping sits below `max_cov - margin` and never enters the
pair search, so the caller returns a false homozygote. On NA12878 this lost
DRB1*03:01 (cov 0.854 against a perfectly covered DRB1*01:01) and DQB1*02:01
(cov 0.886 against DQB1*05:01 at 1.000): 14/16. Both failures showed a low
`pair explains` (0.47 and 0.55) while every correct locus scored 0.875-1.000,
so the metric detects its own failure mode. Widening the margin to 0.15 or
more gives 16/16, and the patient's eight loci are unchanged at every margin
from 0.05 to 0.30, so the defect never touched the patient's calls.

Usage:
    pip install pysam
    python validate_held_out.py --out results
"""
from __future__ import annotations

import argparse, json, os, sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "hla"))
import type_hla as T  # noqa: E402

GIAB = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data"

SAMPLES = {
    "NA12878": {
        "bam": f"{GIAB}/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam",
        "contig": "6",
        # GRCh37 spans; DRB1 widened over the whole DRB region as in type_hla.py.
        "loci": {
            "A": (29909037, 29913661), "B": (31321649, 31324989),
            "C": (31236526, 31239913), "DRB1": (32520000, 32557613),
            "DQA1": (32595956, 32614839), "DQB1": (32627241, 32636160),
            "DPA1": (33032346, 33048555), "DPB1": (33043703, 33057473),
        },
        "truth": {
            "A": ["A*01", "A*11"], "B": ["B*08", "B*56"], "C": ["C*01", "C*07"],
            "DRB1": ["DRB1*01", "DRB1*03:01"],
            "DQA1": ["DQA1*01", "DQA1*05"], "DQB1": ["DQB1*02", "DQB1*05"],
            "DPA1": ["DPA1*01:03", "DPA1*02:01"],
            "DPB1": ["DPB1*04:01", "DPB1*14:01"],
        },
        "source": "Cellosaurus CVCL_7526 (PubMed 27792722)",
        "note": "used to diagnose the margin defect; development set after that",
    },
    "HG002": {
        "bam": f"{GIAB}/AshkenazimTrio/HG002_NA24385_son/"
               f"NIST_HiSeq_HG002_Homogeneity-10953946/"
               f"NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.60x.1.bam",
        "contig": "chr6",
        "loci": {k: (v[1], v[2]) for k, v in T.LOCI.items() if k in ("A", "B", "C")},
        "truth": {"A": ["A*01:01", "A*26:01"], "B": ["B*35:08", "B*38:01"],
                  "C": ["C*04:01", "C*12:03"]},
        "source": "clinical typing for NA24385 reported in arXiv 2209.09155",
        "note": "clean test: no parameter was chosen using this sample",
    },
}


def matches(call, truth):
    """Truth is given at 1 or 2 fields; compare at the resolution provided."""
    return call.startswith(truth + ":") or call == truth if ":" not in truth \
        else call == truth


def run_sample(name, spec, outdir, cache):
    import pysam
    bam = pysam.AlignmentFile(spec["bam"], "rb")
    os.makedirs(os.path.join(outdir, name), exist_ok=True)
    hits = total = 0
    rows = []
    for locus, (s, e) in spec["loci"].items():
        path = os.path.join(outdir, name, f"{locus}.reads.fasta")
        if not os.path.exists(path):
            seqs = {}
            for r in bam.fetch(spec["contig"], s - T.PAD, e + T.PAD):
                if r.is_duplicate or r.is_secondary or r.is_supplementary:
                    continue
                if r.query_sequence:
                    seqs[(r.query_name, r.is_read2)] = r.query_sequence
            with open(path, "w") as fh:
                for (n, r2), q in seqs.items():
                    fh.write(f">{n}/{2 if r2 else 1}\n{q}\n")
        depth = T.fasta_depth(path, e - s + 2 * T.PAD)

        reads = T.load_reads(path)
        ksets, cov, observed, support = T.score_alleles(reads, T.fetch_imgt(locus, cache))
        _, (g1, a1, g2, a2, expl, solo) = T.call_pair(ksets, cov, observed, support)

        called, truth = sorted([g1, g2]), sorted(spec["truth"][locus])
        pool, ok = list(called), []
        for t in truth:
            m = next((c for c in pool if matches(c, t)), None)
            ok.append(m is not None)
            if m:
                pool.remove(m)
        hits += sum(ok); total += len(truth)
        print(f"  {'OK  ' if all(ok) else 'MISS'} HLA-{locus:5} "
              f"called {'/'.join(called):26} truth {'/'.join(truth):24} "
              f"depth={depth:5.1f}x explains={expl:.3f}")
        rows.append({"locus": locus, "called": called, "truth": truth,
                     "matched": ok, "depth_x": depth, "pair_explains": round(expl, 4)})
    print(f"  -> {name}: {hits}/{total} alleles over {len(spec['loci'])} loci "
          f"({spec['note']})")
    return {"source": spec["source"], "note": spec["note"],
            "concordance": [hits, total], "loci": rows}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=".")
    ap.add_argument("--cache", default="../hla/imgt")
    ap.add_argument("--samples", default=",".join(SAMPLES))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    os.makedirs(args.cache, exist_ok=True)

    report = {}
    for name in args.samples.split(","):
        print(f"\n{name}  ({SAMPLES[name]['contig']}, {SAMPLES[name]['source']})")
        report[name] = run_sample(name, SAMPLES[name], args.out, args.cache)
    with open(os.path.join(args.out, "held-out.json"), "w") as fh:
        json.dump(report, fh, indent=2)
    print(f"\nwrote {args.out}/held-out.json")


if __name__ == "__main__":
    main()
