#!/usr/bin/env python3
"""Reproducibly audit the 13 nucleotide-resolved Guardant360 cfDNA loci across BAMs, VCF and MHCflurry.

With --execute, queries indexed BAMs (local paths or public read-only HTTPS with local .bai indices),
checks paired Mutect2 VCF status when supplied, and runs MHCflurry 2.3.14 Class-I affinity predictions
for every coding missense/frameshift variant against Yash's six Class-I alleles.
Without --execute, verifies the recorded TSVs and receipt offline with zero network or model calls.
"""
import argparse
import collections
import gzip
import json
import os
import platform
import sys
from importlib.metadata import version
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "pipeline"))

from dataio import load, resolve_gencode_fasta, rows, save, sha, table

DEFAULT_URLS = {
    "normal_bam": "https://yash-data.samarth-3f0.workers.dev/raw/TN25-112103/normal/normal.recal.bam",
    "tumour_dna_bam": "https://yash-data.samarth-3f0.workers.dev/raw/TN25-112103/DNA_TN25-112103.bam",
    "tumour_rna_bam": "https://yash-data.samarth-3f0.workers.dev/raw/TN25-112103/RNA_TN25-112103.bam",
}

def _derive_observation_status(locus, npile, dpile, rpile, vrec):
    """Derive bounded observation status strictly from computed pileup/VCF counts and report annotations."""
    if locus["gene"] == "DNMT3A":
        return (
            f"cfDNA_{locus['mar_2025_cfdna_pct']}_to_{locus['oct_2025_cfdna_pct']}_"
            f"normalWBC_alt_{npile['alt']}_liverWES_alt_{dpile['alt']}_liverWTS_alt_{rpile['alt']}; "
            "Guardant360_annotates_potential_CHIP_hematopoietic_vs_tumour_origin_unconfirmed"
        )
    parts = []
    if locus["mar_2025_cfdna_pct"] != "ND":
        parts.append("both_cfDNA_draws")
    else:
        parts.append("oct2025_cfDNA_only")
    if vrec.get("filter") == "PASS":
        parts.append("paired_Mutect2_PASS_liver_WES")
    elif dpile["alt"] > 0:
        parts.append(f"uncalled_trace_liver_WES_alt_{dpile['alt']}")
    else:
        parts.append("zero_alt_in_liver_WES")
    if rpile["coding_alt"] >= 2:
        parts.append(f"liver_WTS_coding_alt_{rpile['coding_alt']}")
    elif rpile["alt"] > 0:
        parts.append(f"liver_WTS_shallow_alt_{rpile['alt']}")
    else:
        parts.append("zero_alt_in_liver_WTS")
    if locus["encoded_in_design4_main"]:
        parts.append("encoded_in_Design4")
    qual = "clonality_and_lesion_origin_unestablished"
    if locus["variant_kind"] == "nonsense_snv" and rpile["alt"] == 0:
        qual = "zero_mutant_RNA_compatible_with_NMD_or_absent_mutant_transcription; clonality_unestablished"
    elif locus["variant_kind"] == "nonsense_snv":
        qual = "stop_gain_creates_no_novel_amino_acid_sequence; clonality_unestablished"
    elif locus["variant_kind"] == "promoter_snv":
        qual = "noncoding_promoter_variant; clonality_unestablished"
    return "_".join(parts) + "; " + qual
# Guardant360 plasma draws were collected on 2025-03-29 (epigenomic TF 1.4%) and
# 2025-10-14 (epigenomic TF 31.9%). BCOR L1632P (0.2% Mar -> ND Oct)
# is disclosed separately in clinical JSON and was not included in this executed 13-locus audit.
# Reference BCOR transcripts and an inferable SNV candidate exist; the report row supplies no transcript/nucleotide allele.
# CDH1::NRG1 fusion is a structural RNA/cfDNA rearrangement audited in exports/inframe-request.json.
CTDNA_LOCI = (
    {
        "variant_label": "ROS1_F1699F",
        "gene": "ROS1",
        "ensembl_id": "ENSG00000047936",
        "gene_strand": "-",
        "report_label": "p.F1699F (NM_002944.3:c.5097T>C)",
        "chrom": "chr6",
        "pos_1based": 117337323,
        "ref": "A",
        "alt": "G",
        "variant_kind": "synonymous_snv",
        "mar_2025_cfdna_pct": "1.6%",
        "oct_2025_cfdna_pct": "21.4%",
        "encoded_in_design4_main": False,
        "protein_spec": None,
    },
    {
        "variant_label": "TERT_c.-124C>T",
        "gene": "TERT",
        "ensembl_id": "ENSG00000164362",
        "gene_strand": "-",
        "report_label": "Promoter SNV (NM_198253.3:c.-124C>T)",
        "chrom": "chr5",
        "pos_1based": 1295113,
        "ref": "G",
        "alt": "A",
        "variant_kind": "promoter_snv",
        "mar_2025_cfdna_pct": "0.6%",
        "oct_2025_cfdna_pct": "20.5%",
        "encoded_in_design4_main": False,
        "protein_spec": None,
    },
    {
        "variant_label": "MERTK_E192*",
        "gene": "MERTK",
        "ensembl_id": "ENSG00000153208",
        "gene_strand": "+",
        "report_label": "p.E192* (NM_006343.3:c.574G>T)",
        "chrom": "chr2",
        "pos_1based": 111945051,
        "ref": "G",
        "alt": "T",
        "variant_kind": "nonsense_snv",
        "mar_2025_cfdna_pct": "1.0%",
        "oct_2025_cfdna_pct": "19.4%",
        "encoded_in_design4_main": False,
        "protein_spec": None,
    },
    {
        "variant_label": "JUN_R276*",
        "gene": "JUN",
        "ensembl_id": "ENSG00000177606",
        "gene_strand": "-",
        "report_label": "p.R276* (NM_002228.4:c.826A>T)",
        "chrom": "chr1",
        "pos_1based": 58782245,
        "ref": "T",
        "alt": "A",
        "variant_kind": "nonsense_snv",
        "mar_2025_cfdna_pct": "1.2%",
        "oct_2025_cfdna_pct": "14.8%",
        "encoded_in_design4_main": False,
        "protein_spec": None,
    },
    {
        "variant_label": "RXRA_A214S",
        "gene": "RXRA",
        "ensembl_id": "ENSG00000186350",
        "gene_strand": "+",
        "report_label": "p.A214S (Guardant report) / p.A241S (NM_002957.6:c.721G>T; ENST00000481739.2:c.721G>T)",
        "chrom": "chr9",
        "pos_1based": 134417268,
        "ref": "G",
        "alt": "T",
        "variant_kind": "missense_snv",
        "mar_2025_cfdna_pct": "0.5%",
        "oct_2025_cfdna_pct": "12.1%",
        "encoded_in_design4_main": True,
        "protein_spec": ("RXRA", 241, "A", "S"),
    },
    {
        "variant_label": "TGFBR2_N515fs",
        "gene": "TGFBR2",
        "ensembl_id": "ENSG00000163513",
        "gene_strand": "+",
        "report_label": "p.N515Kfs*26 / p.N490fs (NM_001024847.2:c.1544dupA; NM_003242.5:c.1469dupA)",
        "chrom": "chr3",
        "pos_1based": 30688454,
        "ref": "C",
        "alt": "CA",
        "indel_span_1based": (30688454, 30688457, "CAAC", "CAAAC"),
        "variant_kind": "frameshift_ins",
        "mar_2025_cfdna_pct": "0.7%",
        "oct_2025_cfdna_pct": "11.7%",
        "encoded_in_design4_main": True,
        "protein_spec": ("FRAMESHIFT_WINDOW", "EPPFGSKVREHPCVESMKDKRVERS", 19, 25),
    },
    {
        "variant_label": "KMT2D_R2789Q",
        "gene": "KMT2D",
        "ensembl_id": "ENSG00000167548",
        "gene_strand": "-",
        "report_label": "p.R2789Q (NM_003482.4:c.8366G>A)",
        "chrom": "chr12",
        "pos_1based": 49039222,
        "ref": "C",
        "alt": "T",
        "variant_kind": "missense_snv",
        "mar_2025_cfdna_pct": "0.2%",
        "oct_2025_cfdna_pct": "3.2%",
        "encoded_in_design4_main": False,
        "protein_spec": ("KMT2D", 2789, "R", "Q"),
    },
    {
        "variant_label": "AR_P555S",
        "gene": "AR",
        "ensembl_id": "ENSG00000169083",
        "gene_strand": "+",
        "report_label": "p.P555S (NM_000044.6:c.1663C>T)",
        "chrom": "chrX",
        "pos_1based": 67643302,
        "ref": "C",
        "alt": "T",
        "variant_kind": "missense_snv",
        "mar_2025_cfdna_pct": "ND",
        "oct_2025_cfdna_pct": "0.4%",
        "encoded_in_design4_main": False,
        "protein_spec": ("AR", 555, "P", "S"),
    },
    {
        "variant_label": "DHX36_L147L",
        "gene": "DHX36",
        "ensembl_id": "ENSG00000174953",
        "gene_strand": "-",
        "report_label": "p.L147L (NM_020865.3:c.439T>C)",
        "chrom": "chr3",
        "pos_1based": 154315210,
        "ref": "A",
        "alt": "G",
        "variant_kind": "synonymous_snv",
        "mar_2025_cfdna_pct": "ND",
        "oct_2025_cfdna_pct": "0.3%",
        "encoded_in_design4_main": False,
        "protein_spec": None,
    },
    {
        "variant_label": "ERBB4_D113V",
        "gene": "ERBB4",
        "ensembl_id": "ENSG00000178568",
        "gene_strand": "-",
        "report_label": "p.D113V (NM_005235.3:c.338A>T)",
        "chrom": "chr2",
        "pos_1based": 211947513,
        "ref": "T",
        "alt": "A",
        "variant_kind": "missense_snv",
        "mar_2025_cfdna_pct": "ND",
        "oct_2025_cfdna_pct": "0.3%",
        "encoded_in_design4_main": False,
        "protein_spec": ("ERBB4", 113, "D", "V"),
    },
    {
        "variant_label": "DHX36_R160K",
        "gene": "DHX36",
        "ensembl_id": "ENSG00000174953",
        "gene_strand": "-",
        "report_label": "p.R160K (NM_020865.3:c.479G>A)",
        "chrom": "chr3",
        "pos_1based": 154315170,
        "ref": "C",
        "alt": "T",
        "variant_kind": "missense_snv",
        "mar_2025_cfdna_pct": "ND",
        "oct_2025_cfdna_pct": "0.2%",
        "encoded_in_design4_main": False,
        "protein_spec": ("DHX36", 160, "R", "K"),
    },
    {
        "variant_label": "RAF1_G535S",
        "gene": "RAF1",
        "ensembl_id": "ENSG00000132155",
        "gene_strand": "-",
        "report_label": "p.G535S (NM_002880.4:c.1603G>A)",
        "chrom": "chr3",
        "pos_1based": 12585187,
        "ref": "C",
        "alt": "T",
        "variant_kind": "missense_snv",
        "mar_2025_cfdna_pct": "ND",
        "oct_2025_cfdna_pct": "0.2%",
        "encoded_in_design4_main": False,
        "protein_spec": ("RAF1", 535, "G", "S"),
    },
    {
        "variant_label": "DNMT3A_C351fs",
        "gene": "DNMT3A",
        "ensembl_id": "ENSG00000119772",
        "gene_strand": "-",
        "report_label": "p.C351Afs*56 (Guardant report); NM_022552.5:c.1048_1051del (p.F350Afs*56)",
        "chrom": "chr2",
        "pos_1based": 25247121,
        "ref": "CAAAA",
        "alt": "C",
        "indel_span_1based": (25247121, 25247126, "CAAAAC", "CC"),
        "variant_kind": "frameshift_del",
        "mar_2025_cfdna_pct": "0.2%",
        "oct_2025_cfdna_pct": "0.1%",
        "encoded_in_design4_main": False,
        "protein_spec": None,
    },
)


def _classify_read_allele(read, pos0, ref, alt, indel_span_1based=None):
    """Classify a single aligned read at 0-based VCF anchor pos0 for SNVs and homopolymer-flanked indels."""
    qseq = read.query_sequence
    quals = read.query_qualities
    if qseq is None or quals is None:
        return None
    if len(ref) == 1 and len(alt) == 1 and indel_span_1based is None:
        for qpos, rpos in read.get_aligned_pairs(matches_only=True):
            if rpos == pos0:
                if quals[qpos] < 20:
                    return None
                base = qseq[qpos]
                if base == ref:
                    return "REF"
                if base == alt:
                    return "ALT"
                return "OTHER"
        return None

    left_1, right_1, exp_ref, exp_alt = indel_span_1based
    l0, r0 = left_1 - 1, right_1 - 1
    if read.reference_start > l0 or (read.reference_end or 0) <= r0:
        return None
    # Reject RNA splice skips (CIGAR op 3 == BAM_CREF_SKIP, 'N') overlapping [l0, r0] so an intron skip is never mistaken for a deletion
    ref_cursor = read.reference_start
    for cig_op, cig_len in (read.cigartuples or ()):
        if cig_op == 3:  # BAM_CREF_SKIP ('N')
            if ref_cursor <= r0 and ref_cursor + cig_len > l0:
                return None
            ref_cursor += cig_len
        elif cig_op in (0, 2, 7, 8):  # M, D, =, X consume reference
            ref_cursor += cig_len
    q_left, q_right = None, None
    for qpos, rpos in read.get_aligned_pairs(matches_only=True):
        if rpos == l0:
            q_left = qpos
        elif rpos == r0:
            q_right = qpos
            break
    if q_left is None or q_right is None or q_right < q_left:
        return None
    if any(quals[i] < 20 for i in range(q_left, q_right + 1)):
        return None
    observed_segment = qseq[q_left : q_right + 1]
    if observed_segment == exp_alt:
        return "ALT"
    if observed_segment == exp_ref:
        return "REF"
    return "OTHER"


def _pileup_locus(bam, chrom, pos_1based, ref, alt, *, indel_span_1based=None, stranded_gene_strand=None):
    pos0 = pos_1based - 1
    frag_calls = collections.defaultdict(set)
    strand_calls = collections.defaultdict(set)
    for read in bam.fetch(chrom, max(0, pos0 - 2), pos0 + len(ref) + 2):
        if (read.is_unmapped or read.is_duplicate or read.is_secondary
                or read.is_supplementary or read.is_qcfail or read.mapping_quality < 20):
            continue
        call = _classify_read_allele(read, pos0, ref, alt, indel_span_1based=indel_span_1based)
        if call is None:
            continue
        frag_calls[read.query_name].add(call)
        if stranded_gene_strand is not None:
            # fr-firststrand: Read 1 antisense to transcript, Read 2 sense to transcript
            tx_strand = ("-" if not read.is_reverse else "+") if read.is_read1 else ("+" if not read.is_reverse else "-")
            strand_calls[read.query_name].add((call, "coding" if tx_strand == stranded_gene_strand else "opposite"))

    counts = collections.Counter()
    for qname, cset in frag_calls.items():
        counts[next(iter(cset)) if len(cset) == 1 else "CONFLICT"] += 1

    out = {
        "ref": counts["REF"],
        "alt": counts["ALT"],
        "other": counts["OTHER"],
        "conflict": counts["CONFLICT"],
    }
    if stranded_gene_strand is not None:
        sc = collections.Counter()
        for qname, sset in strand_calls.items():
            if len(frag_calls[qname]) == 1:
                call = next(iter(frag_calls[qname]))
                strands = {st for _, st in sset}
                st_label = next(iter(strands)) if len(strands) == 1 else "mixed"
                sc[(call, st_label)] += 1
        out["coding_ref"] = sc[("REF", "coding")]
        out["coding_alt"] = sc[("ALT", "coding")]
        out["opposite_ref"] = sc[("REF", "opposite")]
        out["opposite_alt"] = sc[("ALT", "opposite")]
    return out


def execute_audit(normal_bam_src, normal_bai, dna_bam_src, dna_bai, rna_bam_src, rna_bai, paired_vcf, models_dir):
    import pysam
    import torch
    from mhcflurry import Class1AffinityPredictor

    os.environ["OMP_NUM_THREADS"] = "2"
    os.environ.setdefault("HTS_CURL_TIMEOUT", "60")
    torch.set_num_threads(2)

    if not paired_vcf or not Path(paired_vcf).is_file():
        raise FileNotFoundError("Provide an existing paired Mutect2 VCF via --paired-vcf when running --execute")
    vcf_lookup = {}
    opener = gzip.open if str(paired_vcf).endswith(".gz") else open
    with opener(paired_vcf, "rt") as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            key = (cols[0], int(cols[1]), cols[3], cols[4])
            fmt = cols[8].split(":")
            tum = dict(zip(fmt, cols[9].split(":")))
            nrm = dict(zip(fmt, cols[10].split(":")))
            vcf_lookup[key] = {
                "filter": cols[6],
                "tumour_ad": tum.get("AD", ""),
                "tumour_af": tum.get("AF", ""),
                "normal_ad": nrm.get("AD", ""),
            }

    # Query one BAM handle at a time so remote HTTPS connections never sit idle and time out
    n_piles, d_piles, r_piles = {}, {}, {}
    with pysam.AlignmentFile(normal_bam_src, "rb", index_filename=str(normal_bai)) as nbam:
        for locus in CTDNA_LOCI:
            n_piles[locus["variant_label"]] = _pileup_locus(
                nbam, locus["chrom"], locus["pos_1based"], locus["ref"], locus["alt"],
                indel_span_1based=locus.get("indel_span_1based"),
            )
    with pysam.AlignmentFile(dna_bam_src, "rb", index_filename=str(dna_bai)) as dbam:
        for locus in CTDNA_LOCI:
            d_piles[locus["variant_label"]] = _pileup_locus(
                dbam, locus["chrom"], locus["pos_1based"], locus["ref"], locus["alt"],
                indel_span_1based=locus.get("indel_span_1based"),
            )
    with pysam.AlignmentFile(rna_bam_src, "rb", index_filename=str(rna_bai)) as rbam:
        for locus in CTDNA_LOCI:
            r_piles[locus["variant_label"]] = _pileup_locus(
                rbam, locus["chrom"], locus["pos_1based"], locus["ref"], locus["alt"],
                indel_span_1based=locus.get("indel_span_1based"),
                stranded_gene_strand=locus["gene_strand"],
            )

    pileup_rows = []
    for locus in CTDNA_LOCI:
        vid = locus["variant_label"]
        ch, pos1, ref, alt = locus["chrom"], locus["pos_1based"], locus["ref"], locus["alt"]
        npile, dpile, rpile = n_piles[vid], d_piles[vid], r_piles[vid]
        vrec = vcf_lookup.get((ch, pos1, ref, alt), {})
        derived_status = _derive_observation_status(locus, npile, dpile, rpile, vrec)
        pileup_rows.append({
            "variant_label": vid,
            "gene": locus["gene"],
            "ensembl_id": locus["ensembl_id"],
            "gene_strand": locus["gene_strand"],
            "report_label": locus["report_label"],
            "grch38_variant": f"{ch}:{pos1}:{ref}>{alt}",
            "variant_kind": locus["variant_kind"],
            "mar_2025_cfdna": locus["mar_2025_cfdna_pct"],
            "oct_2025_cfdna": locus["oct_2025_cfdna_pct"],
            "encoded_in_design4_main": locus["encoded_in_design4_main"],
            "normal_wbc_ref_fragments": npile["ref"],
            "normal_wbc_alt_fragments": npile["alt"],
            "normal_wbc_other_or_conflict": npile["other"] + npile["conflict"],
            "liver_wes_ref_fragments": dpile["ref"],
            "liver_wes_alt_fragments": dpile["alt"],
            "liver_wes_other_or_conflict": dpile["other"] + dpile["conflict"],
            "paired_mutect2_filter": vrec.get("filter", "not_called_in_supplied_vcf"),
            "paired_mutect2_tumour_ad": vrec.get("tumour_ad", ""),
            "paired_mutect2_tumour_af": vrec.get("tumour_af", ""),
            "liver_wts_ref_fragments": rpile["ref"],
            "liver_wts_alt_fragments": rpile["alt"],
            "liver_wts_coding_ref_fragments": rpile["coding_ref"],
            "liver_wts_coding_alt_fragments": rpile["coding_alt"],
            "liver_wts_opposite_ref_fragments": rpile["opposite_ref"],
            "liver_wts_opposite_alt_fragments": rpile["opposite_alt"],
            "sample_observation_status": derived_status,
        })

    # Extract 25-aa protein contexts from bundled GENCODE v46 translations and run MHCflurry 2.3.14
    gencode_fa = resolve_gencode_fasta("gencode.v46.pc_translations.fa.gz", ROOT / "references")
    gene_seqs = collections.defaultdict(list)
    with gzip.open(gencode_fa, "rt") as stream:
        hdr, parts = None, []
        for line in stream:
            if line.startswith(">"):
                if hdr is not None:
                    p = hdr[1:].strip().split("|")
                    if len(p) >= 7:
                        gene_seqs[p[6]].append((p[1], "".join(parts)))
                hdr, parts = line, []
            else:
                parts.append(line.strip())
        if hdr is not None:
            p = hdr[1:].strip().split("|")
            if len(p) >= 7:
                gene_seqs[p[6]].append((p[1], "".join(parts)))

    alleles = load(ROOT / "exports/source-case.json")["hla_class_I"]
    ap = Class1AffinityPredictor.load(str(Path(models_dir) / "affinity_predictor"))

    pred_rows = []
    for locus in CTDNA_LOCI:
        pspec = locus["protein_spec"]
        if pspec is None:
            continue
        if pspec[0] == "FRAMESHIFT_WINDOW":
            _, mut_win, t_lo, t_hi = pspec
            tx_id = "ENST00000359013"
            wt_win = None
            for _, full_s in gene_seqs["TGFBR2"]:
                if mut_win[:t_lo] in full_s:
                    idx = full_s.index(mut_win[:t_lo])
                    cand = full_s[idx : idx + len(mut_win)]
                    if len(cand) == len(mut_win):
                        wt_win = cand
                        break
        else:
            gsym, pos_aa, wt_aa, mut_aa = pspec
            tx_id, wt_win, mut_win, t_lo, t_hi = None, None, None, None, None
            for tid, full_s in gene_seqs[gsym]:
                if len(full_s) >= pos_aa and full_s[pos_aa - 1] == wt_aa:
                    lo = max(0, min(len(full_s) - 25, pos_aa - 13))
                    hi = min(len(full_s), lo + 25)
                    w = full_s[lo:hi]
                    idx = pos_aa - 1 - lo
                    tx_id = tid
                    wt_win = w
                    mut_win = w[:idx] + mut_aa + w[idx + 1:]
                    t_lo, t_hi = idx, idx + 1
                    break
            if mut_win is None:
                raise ValueError("Failed to locate reference protein window for " + locus["variant_label"])

        pep_pairs = []
        for k in (8, 9, 10, 11):
            for st in range(len(mut_win) - k + 1):
                en = st + k
                if st < t_hi and en > t_lo:
                    m_pep = mut_win[st:en]
                    w_pep = wt_win[st:en] if wt_win and len(wt_win) >= en else ""
                    if m_pep != w_pep:
                        pep_pairs.append((m_pep, w_pep))

        unique_peps = sorted({p for pair in pep_pairs for p in pair if p})
        aff_lookup = {}
        for al in alleles:
            df = ap.predict_to_dataframe(peptides=unique_peps, allele=al)
            for row in df.itertuples():
                aff_lookup[(row.peptide, al)] = (float(row.prediction), float(row.prediction_percentile))

        best_per_allele = []
        for al in alleles:
            best_pair = min(pep_pairs, key=lambda pw: aff_lookup[(pw[0], al)][0])
            m_nm, m_rk = aff_lookup[(best_pair[0], al)]
            w_nm, w_rk = aff_lookup[(best_pair[1], al)] if best_pair[1] else (None, None)
            rec = {
                "variant_label": locus["variant_label"],
                "gene": locus["gene"],
                "transcript_id": tx_id,
                "allele": al,
                "best_mutant_peptide": best_pair[0],
                "mutant_affinity_nm": round(m_nm, 4),
                "mutant_percentile_rank": round(m_rk, 4),
                "paired_wildtype_peptide": best_pair[1],
                "wildtype_affinity_nm": "" if w_nm is None else round(w_nm, 4),
                "wildtype_percentile_rank": "" if w_rk is None else round(w_rk, 4),
            }
            pred_rows.append(rec)
            best_per_allele.append(rec)

    def _bam_identity(src_str):
        p = Path(src_str)
        if p.is_file():
            return {"access_mode": "local_whole_file_verified", "bytes": p.stat().st_size, "sha256": sha(p)}
        return {
            "access_mode": "indexed_remote_http_range_unverified_whole_payload",
            "url": src_str,
            "whole_file_sha256_independently_verified": False,
        }

    table(ROOT / "exports/ctdna-locus-pileups.tsv", pileup_rows)
    table(ROOT / "exports/ctdna-mhcflurry-predictions.tsv", pred_rows)
    save(ROOT / "exports/ctdna-audit-receipt.json", {
        "alignment_filters": {
            "count_unit": "Read-name fragment union; conflicting overlapping mates counted separately under other_or_conflict; not UMI molecules.",
            "duplicates_excluded": True,
            "indel_rule": (
                "Flanking-anchor read-span comparison across [span_start_1based, span_end_1based] "
                "(chr3:30688454-30688457 CAAC->CAAAC for TGFBR2; chr2:25247121-25247126 CAAAAC->CC for DNMT3A) "
                "excluding reads with CIGAR reference skips ('N') overlapping [span_start, span_end] and requiring "
                "continuous alignment of both anchors and min base quality >=20 across the spanned bases."
            ),
            "min_focal_base_quality": 20,
            "min_mapping_quality": 20,
            "qcfail_secondary_supplementary_excluded": True,
            "rna_strand_orientation_rule": "fr-firststrand (Read 1 antisense to transcript, Read 2 sense to transcript).",
        },
        "environment": {
            "machine": platform.machine(),
            "mhcflurry": version("mhcflurry"),
            "pysam": pysam.__version__,
            "python": platform.python_version(),
            "system": platform.system(),
        },
        "input_identities": {
            "gencode_v46_pc_translations_sha256": sha(gencode_fa),
            "mhcflurry_affinity_manifest_sha256": sha(Path(models_dir) / "affinity_predictor/manifest.csv"),
            "normal_bam": _bam_identity(normal_bam_src),
            "normal_bai_sha256": sha(normal_bai),
            "paired_mutect2_all_unique_vcf_sha256": sha(paired_vcf),
            "source_case_sha256": sha(ROOT / "exports/source-case.json"),
            "tumour_dna_bam": _bam_identity(dna_bam_src),
            "tumour_dna_bai_sha256": sha(dna_bai),
            "tumour_rna_bam": _bam_identity(rna_bam_src),
            "tumour_rna_bai_sha256": sha(rna_bai),
        },
        "loci_audited": len(pileup_rows),
        "output_sha256": {
            "ctdna-locus-pileups.tsv": sha(ROOT / "exports/ctdna-locus-pileups.tsv"),
            "ctdna-mhcflurry-predictions.tsv": sha(ROOT / "exports/ctdna-mhcflurry-predictions.tsv"),
        },
        "peptide_hla_rows": len(pred_rows),
        "scope": "Reproducible locus-level BAM fragment pileups and MHCflurry Class-I predictions for the 13 nucleotide-resolved Guardant360 cfDNA alterations. Does not establish clonality, anatomical lesion origin, mRNA translation or clinical immunogenicity.",
        "script_sha256": sha(Path(__file__)),
    })


def verify_offline():
    receipt = load(ROOT / "exports/ctdna-audit-receipt.json")
    for fname, expected in receipt["output_sha256"].items():
        if sha(ROOT / "exports" / fname) != expected:
            raise ValueError("ctDNA locus audit checksum mismatch: " + fname)
    p_rows = rows(ROOT / "exports/ctdna-locus-pileups.tsv")
    m_rows = rows(ROOT / "exports/ctdna-mhcflurry-predictions.tsv")
    if len(p_rows) != len(CTDNA_LOCI) or len(p_rows) != receipt["loci_audited"]:
        raise ValueError("Unexpected ctDNA pileup row count")
    if len(m_rows) != receipt["peptide_hla_rows"]:
        raise ValueError("Unexpected ctDNA MHCflurry row count")
    return {
        "ctdna_loci_pileups_verified": len(p_rows),
        "ctdna_mhcflurry_rows_verified": len(m_rows),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execute", action="store_true", help="Run live indexed BAM pileups and MHCflurry inference")
    parser.add_argument("--normal-bam", default=DEFAULT_URLS["normal_bam"])
    parser.add_argument("--normal-bai", type=Path)
    parser.add_argument("--tumour-dna-bam", default=DEFAULT_URLS["tumour_dna_bam"])
    parser.add_argument("--tumour-dna-bai", type=Path)
    parser.add_argument("--tumour-rna-bam", default=DEFAULT_URLS["tumour_rna_bam"])
    parser.add_argument("--tumour-rna-bai", type=Path)
    parser.add_argument("--paired-vcf", type=Path)
    parser.add_argument("--models", type=Path)
    args = parser.parse_args()
    if args.execute:
        for req_name in ("normal_bai", "tumour_dna_bai", "tumour_rna_bai", "paired_vcf", "models"):
            if getattr(args, req_name) is None:
                parser.error(f"--{req_name.replace('_', '-')} is required with --execute")
        execute_audit(
            args.normal_bam, args.normal_bai,
            args.tumour_dna_bam, args.tumour_dna_bai,
            args.tumour_rna_bam, args.tumour_rna_bai,
            args.paired_vcf, args.models,
        )
    print(json.dumps(verify_offline(), indent=2))


if __name__ == "__main__":
    main()
