#!/usr/bin/env python3
"""Type Yash's HLA class II (and class I, as a control) from the matched-normal WGS.

The repository records six class I alleles and no class II typing, so no
class II-aware design is currently possible. The matched-normal white-blood-cell
BAM already published with this project is germline WGS at ~40-95x over the MHC,
which is enough to type class II directly. The BAM is served over HTTPS with
byte-range support, so only the MHC region is transferred (~17 s, a few hundred
MB), not the 143 GB file.

Method
------
For each locus, reads overlapping the gene are pulled from the BAM and reduced to
canonical 31-mers. Each IMGT/HLA genomic allele sequence is scored two ways:

  cov(a)   fraction of allele a's k-mers observed in the reads. An allele the
           donor carries has essentially all of them at this depth; one they do
           not carry has distinctive k-mers that are missing.
  expl(p)  fraction of the locus's observed k-mers that an allele *pair*
           explains, where the observed set is restricted to k-mers belonging to
           some allele of this locus (this drops intronic, flanking and
           paralogue k-mers). A real heterozygote explains strictly more than
           either allele alone, so zygosity is measured, not assumed.

Pairs are ranked by expl (coarsely, so a true heterozygote always beats a
homozygous call), then by total coverage, then by how well documented the
alleles are. That last tie-break matters because IMGT contains many rare alleles
differing from a common one by a few bases; without it a singleton deposit can
win by a hair.

Validation
----------
Run with --loci A,B,C the method must reproduce the six class I alleles that
Caris (config/case.json) and Guardant360 independently report. It does, 6/6.
That control is the reason to trust the class II calls, which have no
independent report to check against. DRB3/4/5 presence is also tested, since
those genes are absent from the GRCh38 primary assembly this BAM was aligned to
and their reads land on DRB1.

Usage
-----
    pip install pysam
    python type_hla.py --out results            # all loci
    python type_hla.py --loci A,B,C --out ctrl  # validation only
"""
from __future__ import annotations

import argparse, collections, itertools, json, os, sys, urllib.request

BAM_URL = ("https://yash-data.samarth-3f0.workers.dev/"
           "raw/TN25-112103/normal/normal.recal.bam")
IMGT_BASE = "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/fasta"

K = 31
MIN_ABUNDANCE = 3          # drops sequencing-error k-mers
MIN_LEN_FRAC = 0.6         # keep full-length genomic allele entries only
PAD = 2000

# GRCh38 gene spans. DRB1 is widened to the whole DRB region so that reads from
# a divergent DRB1 allele, or from DRB3/4/5, are collected rather than lost.
LOCI = {
    "A":    ("chr6", 29942532, 29945870),
    "B":    ("chr6", 31353875, 31357188),
    "C":    ("chr6", 31268749, 31272130),
    "DRB1": ("chr6", 32546000, 32589848),
    "DQA1": ("chr6", 32628179, 32647062),
    "DQB1": ("chr6", 32659467, 32668383),
    "DPA1": ("chr6", 33064569, 33080775),
    "DPB1": ("chr6", 33075926, 33089696),
}
# Scored against the DRB1 read set; these genes have no GRCh38 primary locus.
DRB_PARALOGUES = ["DRB3", "DRB4", "DRB5"]

COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


def kmers(seq: str, k: int = K):
    for i in range(len(seq) - k + 1):
        km = seq[i:i + k]
        if "N" in km:
            continue
        r = rc(km)
        yield km if km < r else r


def read_fasta(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf)
                name, buf = line[1:], []
            else:
                buf.append(line.upper())
    if name:
        yield name, "".join(buf)


def two_field(allele: str) -> str:
    gene, _, rest = allele.partition("*")
    f = rest.split(":")
    return f"{gene}*{f[0]}:{f[1]}" if len(f) >= 2 else allele


def fasta_depth(path, span):
    total = sum(len(s) for _, s in read_fasta(path))
    return round(total / span, 1)


def fetch_reads(locus, outdir):
    """Pull reads over one locus from the remote BAM into a FASTA."""
    import pysam
    path = os.path.join(outdir, f"{locus}.reads.fasta")
    chrom, start, end = LOCI[locus]
    span = end - start + 2 * PAD
    if os.path.exists(path):
        return path, fasta_depth(path, span)
    bam = pysam.AlignmentFile(BAM_URL, "rb")
    seqs = {}
    for r in bam.fetch(chrom, start - PAD, end + PAD):
        if r.is_duplicate or r.is_secondary or r.is_supplementary:
            continue
        if r.query_sequence:
            seqs[(r.query_name, r.is_read2)] = r.query_sequence
    with open(path, "w") as fh:
        for (name, r2), s in seqs.items():
            fh.write(f">{name}/{2 if r2 else 1}\n{s}\n")
    depth = round(sum(len(s) for s in seqs.values()) / span, 1)
    print(f"  {locus:5} reads={len(seqs):6} approx_depth={depth:5.1f}x", flush=True)
    return path, depth


def fetch_imgt(gene, cachedir):
    path = os.path.join(cachedir, f"{gene}_gen.fasta")
    if not os.path.exists(path):
        urllib.request.urlretrieve(f"{IMGT_BASE}/{gene}_gen.fasta", path)
    return path


def load_reads(path):
    cnt = collections.Counter()
    for _, s in read_fasta(path):
        cnt.update(kmers(s))
    return {k for k, v in cnt.items() if v >= MIN_ABUNDANCE}


def score_alleles(reads, imgt_path):
    alleles = [(h.split()[1], s) for h, s in read_fasta(imgt_path)]
    maxlen = max(len(s) for _, s in alleles)
    ksets, cov = {}, {}
    locus_kmers = set()
    for a, s in alleles:
        if len(s) < MIN_LEN_FRAC * maxlen:
            continue
        ks = set(kmers(s))
        if not ks:
            continue
        ksets[a], cov[a] = ks, len(ks & reads) / len(ks)
        locus_kmers |= ks
    support = collections.Counter(two_field(a) for a, _ in alleles)
    return ksets, cov, locus_kmers & reads, support


def call_pair(ksets, cov, observed, support, cov_margin=0.20):
    ranked_alleles = sorted(ksets, key=lambda a: -cov[a])
    best_rep = {}
    for a in ranked_alleles:
        best_rep.setdefault(two_field(a), a)
    ranked = sorted(best_rep.items(), key=lambda kv: -cov[kv[1]])

    # Candidate pool. Top-by-coverage alone lets one allele group (e.g. A*01)
    # fill every slot and hide the second allele, so also keep the best few
    # 2-field groups from each first field.
    per_first = collections.defaultdict(list)
    for g, a in ranked:
        per_first[g.split(":")[0]].append((g, a))
    cands, seen = [], set()
    for g, a in ranked[:20]:
        if g not in seen:
            seen.add(g); cands.append((g, a))
    for first in sorted(per_first, key=lambda f: -cov[per_first[f][0][1]]):
        for g, a in per_first[first][:3]:
            if g not in seen:
                seen.add(g); cands.append((g, a))
    # Coverage ceilings are locus-specific: class I alleles reach 1.00 whereas
    # DRB1 tops out near 0.85, because long introns leave gaps that every
    # allele at the locus shares. So admit alleles within a margin of the best
    # at *this* locus rather than applying one absolute threshold everywhere.
    #
    # The margin has to be generous. A divergent second allele loses reads to
    # mismapping and can sit 0.10-0.15 below a perfectly covered first allele;
    # a tight margin then excludes the true partner and the search returns a
    # false homozygote. That is exactly how a 0.05 margin lost DRB1*03:01 and
    # DQB1*02:01 on NA12878 (see validation/). Both failures announced
    # themselves through a low `pair explains` value, so treat anything below
    # ~0.8 there as "second allele possibly missed" rather than a firm call.
    cov_floor = max(cov[a] for _, a in ranked) - cov_margin
    pool = [(g, a) for g, a in cands if cov[a] >= cov_floor][:90] or ranked[:30]

    best = None
    for (g1, a1), (g2, a2) in itertools.combinations_with_replacement(pool, 2):
        expl = len((ksets[a1] | ksets[a2]) & observed) / len(observed)
        key = (round(expl, 2), round(cov[a1] + cov[a2], 4),
               support[g1] + support[g2], expl)
        if best is None or key > best[0]:
            best = (key, expl, g1, a1, g2, a2)
    _, expl, g1, a1, g2, a2 = best
    solo = max(len(ksets[a] & observed) / len(observed) for _, a in pool)
    return ranked, (g1, a1, g2, a2, expl, solo)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=".")
    ap.add_argument("--cache", default="imgt")
    ap.add_argument("--loci", default=",".join(LOCI))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    os.makedirs(args.cache, exist_ok=True)
    loci = args.loci.split(",")

    print("Fetching reads over the MHC from the matched-normal BAM:")
    fetched = {l: fetch_reads(l, args.out) for l in loci}
    read_paths = {l: v[0] for l, v in fetched.items()}
    depths = {l: v[1] for l, v in fetched.items()}

    out = {}
    for locus in loci:
        reads = load_reads(read_paths[locus])
        ksets, cov, observed, support = score_alleles(
            reads, fetch_imgt(locus, args.cache))
        ranked, (g1, a1, g2, a2, expl, solo) = call_pair(ksets, cov, observed, support)
        zyg = "homozygous" if g1 == g2 else "heterozygous"
        print(f"\nHLA-{locus}")
        for g, a in ranked[:6]:
            print(f"   cov={cov[a]:7.4f}  {g:16} ({a})")
        print(f"   -> {g1} / {g2}  [{zyg}]  pair explains {expl:.4f} "
              f"of locus k-mers (best single allele {solo:.4f})")
        out[locus] = {"call": [g1, g2], "representative": [a1, a2],
                      "coverage": [cov[a1], cov[a2]], "zygosity": zyg,
                      "approx_depth_x": depths[locus],
                      "pair_explains": expl, "best_single_allele": solo,
                      "top": [[cov[a], g, a] for g, a in ranked[:10]]}

        if locus == "DRB1":
            # DRB3/4/5 are absent from the GRCh38 primary assembly, so their
            # reads land here. Presence/absence is itself a check: a DRB1*03
            # haplotype carries DRB3, and DRB4/DRB5 should be clearly absent.
            out["DRB345"] = {}
            for para in DRB_PARALOGUES:
                pk, pc, pobs, psup = score_alleles(reads, fetch_imgt(para, args.cache))
                top = sorted(((pc[a], two_field(a)) for a in pk), reverse=True)[:1]
                c, g = top[0]
                present = c >= 0.70
                print(f"   {para}: best {g} cov={c:.4f} -> "
                      f"{'PRESENT' if present else 'absent'}")
                out["DRB345"][para] = {"best": g, "coverage": c, "present": present}

    with open(os.path.join(args.out, "hla_calls.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nwrote {args.out}/hla_calls.json")


if __name__ == "__main__":
    sys.exit(main())
