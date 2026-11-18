"""Reproduce bounded reference, homolog, routing and RNA model screens.

Outputs are computational research evidence. They do not establish clinical safety,
benefit, manufacturing quality, named-patient custody or clinical authorization.
"""
import argparse
import gzip
import math
import os
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import RNA
from Bio import SeqIO
from Bio.Seq import Seq

from curate import HERE, peptide_windows
from dataio import iter_rows, load, resolve_gencode_fasta, rows, save, sha, table, cell
from class2_construct import target_class2, verify_class2_annotations, eligible_core_rank
from research_qualifications import (
    SCREENING_POLICY, MFE_REPORTING_CUTOFF_BP, MFE_SCOPE, RNA_INTERPRETATION,
    MOTIF_INTERPRETATION, SPECIFICITY_INTERPRETATION, SAFETY_LIMITS,
    cryptic_interpretation_for, research_questions, refresh_cached_annotations,
    BIOPSY_EXPRESSION_SCOPE, CRYPTIC_STRAND_SCOPE, reference_parent_genes,
    read_biopsy_expression, biopsy_expression_fields, cryptic_counts_by_stored_strand,
)

SAFETY_DIR = HERE / "safety"
STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
CLASS1_CODES = ("A0101", "A3303", "B0801", "B1508", "C0102", "C0702")
ANCHOR_GROUPS = (
    frozenset("LIVMFYW"),
    frozenset("KRH"),
    frozenset("DE"),
    frozenset("AGSTCP"),
    frozenset("NQ"),
)
SLIPPERY_MOTIFS_RNA = (
    "UUUU", "UUUC", "UUUA", "UUUG", "AAAA", "CCCC", "GGGG",
)


def matching_positions(sequence, needle):
    """Yield every overlapping anchor occurrence, retaining its sequence context."""
    if not needle:
        raise ValueError("Reference anchor must not be empty")
    start = 0
    while True:
        index = sequence.find(needle, start)
        if index < 0:
            return
        yield index
        start = index + 1


def own_reference_genes(gene):
    return reference_parent_genes(gene)


def verify_homolog_gene_record(record):
    """Check complete displayed gene counts and the declared exact-match scope."""
    names = record["homolog_genes"].split(";")
    own = own_reference_genes(record["source_gene"])
    if not all(names) or names != sorted(set(names)) or len(set(names) - own) != int(record["offtarget_gene_count"]):
        raise ValueError("Incomplete or inconsistent homolog gene reporting")
    is_own = set(names) <= own
    if str(is_own) != record["is_mutated_gene_own_WT_only"]:
        raise ValueError("Own-gene classification differs from complete gene names")
    if int(record["mismatches"]) == 0 and not is_own:
        if record["query_origin"] == "TARGET":
            raise ValueError("Exact off-target reference match violates selected-target exclusion policy")
        if record["model_priority_category"] != "EXACT_REFERENCE_MATCH_JUNCTION_REVIEW":
            raise ValueError("Exact junction reference match must be explicitly flagged for review")


def same_anchor_group(a, b):
    return a == b or any(a in group and b in group for group in ANCHOR_GROUPS)


def parse_gene_and_id(header, normal_gene_map=None):
    parts = header.split("|")
    if len(parts) >= 7 and parts[0].startswith("ENSP"):
        return parts[6], parts[0]
    if len(parts) == 2 and parts[0].startswith("TYPED_HLA"):
        return parts[1], parts[0]
    if normal_gene_map and header in normal_gene_map:
        return normal_gene_map[header], header
    return "PATIENT_NORMAL_MODEL", header


def disjoint_seeds(peptide):
    n = len(peptide)
    if n == 8:
        return ((0, peptide[0:3]), (3, peptide[3:6]), (6, peptide[6:8]))
    if n == 9:
        return ((0, peptide[0:3]), (3, peptide[3:6]), (6, peptide[6:9]))
    if n == 10:
        return ((0, peptide[0:3]), (3, peptide[3:6]), (6, peptide[6:10]))
    if n == 11:
        return ((0, peptide[0:4]), (4, peptide[4:8]), (8, peptide[8:11]))
    raise ValueError("Expected Class-I peptide length 8..11: " + peptide)


def tcr_core_seeds(peptide):
    """Contiguous 3-aa seeds inside central TCR-facing positions P3..P(L-1) (0-based 2..L-1)."""
    n = len(peptide)
    return tuple((offset, peptide[offset:offset + 3]) for offset in range(2, n - 3))


def classify_homolog(query, candidate):
    """Return (match_kind, mismatches, mismatch_positions_1based, central_tcr_matches) or None."""
    if len(query) != len(candidate) or not set(candidate) <= STANDARD_AA:
        return None
    n = len(query)
    diffs = [idx + 1 for idx, (a, b) in enumerate(zip(query, candidate)) if a != b]
    mm = len(diffs)
    central_same = sum(query[i] == candidate[i] for i in range(2, n - 1))
    if mm == 0:
        return ("EXACT_0MM", 0, "", central_same)
    if mm <= 2:
        return (f"NEAR_{mm}MM", mm, ",".join(map(str, diffs)), central_same)
    # Titin-class mimicry (Cameron et al., Sci Transl Med 2013: EVDPIGHLY vs ESDPIVAQY = 4 mm,
    # P3-P5 'DPI' identical + conserved C-anchor): <=4 mm, >=3 contiguous identical residues in
    # P3..P(L-1), and conserved C-terminal anchor group plus conserved P1 or P2 anchor group.
    has_contiguous_tcr3 = any(
        query[pos:pos + 3] == candidate[pos:pos + 3] for pos in range(2, n - 3)
    )
    anchors_compatible = same_anchor_group(query[-1], candidate[-1]) and (
        same_anchor_group(query[0], candidate[0]) or same_anchor_group(query[1], candidate[1])
    )
    if mm <= 4 and has_contiguous_tcr3 and anchors_compatible:
        return ("TCR_CORE_ANCHOR_MIMIC", mm, ",".join(map(str, diffs)), central_same)
    return None


def scan_class1_mimicry(queries, proteome_records, normal_gene_map=None):
    """Exact 100%-recall scan for <=2-mismatch and Titin-class TCR-core mimics."""
    seed_to_queries = defaultdict(list)
    seed_lengths = set()
    for q in sorted(set(queries)):
        seen_pairs = set()
        for offset, seed in disjoint_seeds(q) + tcr_core_seeds(q):
            if (offset, seed) not in seen_pairs:
                seen_pairs.add((offset, seed))
                seed_to_queries[seed].append((offset, q))
                seed_lengths.add(len(seed))
    ordered_seeds = sorted(seed_to_queries)
    hits = {}

    def record_hit(q, cand, gene, prot_id, source_tag):
        verdict = classify_homolog(q, cand)
        if verdict is None:
            return
        kind, mm, pos_str, central_same = verdict
        entry = hits.setdefault((q, cand), {
            "query_peptide": q, "homolog_peptide": cand, "length": len(q),
            "match_class": kind, "mismatches": mm, "mismatch_positions_1based": pos_str,
            "central_P3_to_PLminus1_identical_aa": central_same,
            "genes": set(), "example_records": [], "sources": set(),
        })
        entry["genes"].add(gene)
        entry["sources"].add(source_tag)
        if len(entry["example_records"]) < 5:
            entry["example_records"].append(f"{gene}:{prot_id}")

    try:
        import ahocorasick_rs
        ac = ahocorasick_rs.AhoCorasick(ordered_seeds)
        for source_tag, header, seq in proteome_records:
            gene, prot_id = parse_gene_and_id(header, normal_gene_map)
            seq_len = len(seq)
            for pat_idx, start, _ in ac.find_matches_as_indexes(seq, overlapping=True):
                seed = ordered_seeds[pat_idx]
                for offset, q in seed_to_queries[seed]:
                    w_start = start - offset
                    w_end = w_start + len(q)
                    if 0 <= w_start and w_end <= seq_len:
                        record_hit(q, seq[w_start:w_end], gene, prot_id, source_tag)
    except ImportError:
        by_len = sorted(seed_lengths)
        for source_tag, header, seq in proteome_records:
            gene, prot_id = parse_gene_and_id(header, normal_gene_map)
            seq_len = len(seq)
            for i in range(seq_len - min(by_len) + 1):
                for s_len in by_len:
                    if i + s_len <= seq_len:
                        seed = seq[i:i + s_len]
                        if seed in seed_to_queries:
                            for offset, q in seed_to_queries[seed]:
                                w_start = i - offset
                                w_end = w_start + len(q)
                                if 0 <= w_start and w_end <= seq_len:
                                    record_hit(q, seq[w_start:w_end], gene, prot_id, source_tag)
    return hits


def scan_class2_9mer_cores(cores, proteome_records, normal_gene_map=None):
    """Exact <=1-mismatch search across 9-aa MHC-II cores (disjoint pigeonhole: [0:5] or [5:9] exact)."""
    seed_to_cores = defaultdict(list)
    seed_lengths = (4, 5)
    for core in sorted(set(cores)):
        if len(core) != 9:
            raise ValueError("Expected 9-aa Class-II core: " + core)
        seed_to_cores[core[0:5]].append((0, core))
        seed_to_cores[core[5:9]].append((5, core))
    ordered_seeds = sorted(seed_to_cores)
    hits = {}

    def record(core, cand, gene, prot_id, source_tag):
        if not set(cand) <= STANDARD_AA:
            return
        diffs = [i + 1 for i, (a, b) in enumerate(zip(core, cand)) if a != b]
        if len(diffs) > 1:
            return
        entry = hits.setdefault((core, cand), {
            "query_core_9aa": core, "homolog_core_9aa": cand,
            "mismatches": len(diffs), "mismatch_positions_1based": ",".join(map(str, diffs)),
            "genes": set(), "sources": set(), "example_records": [],
        })
        entry["genes"].add(gene)
        entry["sources"].add(source_tag)
        if len(entry["example_records"]) < 5:
            entry["example_records"].append(f"{gene}:{prot_id}")

    try:
        import ahocorasick_rs
        ac = ahocorasick_rs.AhoCorasick(ordered_seeds)
        for source_tag, header, seq in proteome_records:
            gene, prot_id = parse_gene_and_id(header, normal_gene_map)
            seq_len = len(seq)
            for pat_idx, start, _ in ac.find_matches_as_indexes(seq, overlapping=True):
                for offset, core in seed_to_cores[ordered_seeds[pat_idx]]:
                    w_start = start - offset
                    if 0 <= w_start and w_start + 9 <= seq_len:
                        record(core, seq[w_start:w_start + 9], gene, prot_id, source_tag)
    except ImportError:
        for source_tag, header, seq in proteome_records:
            gene, prot_id = parse_gene_and_id(header, normal_gene_map)
            seq_len = len(seq)
            for i in range(max(0, seq_len - 3)):
                for s_len in seed_lengths:
                    if i + s_len <= seq_len:
                        seed = seq[i:i + s_len]
                        if seed in seed_to_cores:
                            for offset, core in seed_to_cores[seed]:
                                w_start = i - offset
                                if 0 <= w_start and w_start + 9 <= seq_len:
                                    record(core, seq[w_start:w_start + 9], gene, prot_id, source_tag)
    return hits


def longest_duplex_metrics(rna_sequence):
    """Measure stems within one predicted minimum-energy structure; no activation test."""
    RNA.cvar.temperature = 37.0
    structure, mfe = RNA.fold(rna_sequence)
    ptable = RNA.ptable(structure)
    n = len(rna_sequence)
    visited = set()
    max_uninterrupted_bp = 0
    max_1bulge_bp = 0
    for i in range(1, n + 1):
        j = ptable[i]
        if j <= i or (i, j) in visited:
            continue
        # Trace exact uninterrupted helix (i+k paired with j-k)
        k = 0
        while i + k < j - k and ptable[i + k] == j - k:
            visited.add((i + k, j - k))
            k += 1
        if k > max_uninterrupted_bp:
            max_uninterrupted_bp = k
        # Trace stem allowing <=1-nt bulge or 1x1 internal mismatch at each step (conservative upper bound)
        ci, cj, bp_1bulge = i, j, 0
        while ci < cj:
            if ptable[ci] == cj:
                bp_1bulge += 1
                ci += 1
                cj -= 1
            elif ci + 1 < cj and ptable[ci + 1] == cj:
                ci += 1
            elif ci < cj - 1 and ptable[ci] == cj - 1:
                cj -= 1
            elif ci + 1 < cj - 1 and ptable[ci + 1] == cj - 1:
                ci += 1
                cj -= 1
            else:
                break
        if bp_1bulge > max_1bulge_bp:
            max_1bulge_bp = bp_1bulge
    u_count = rna_sequence.count("U")
    uu_count = sum(rna_sequence[i:i + 2] == "UU" for i in range(n - 1))
    cpg_count = sum(rna_sequence[i:i + 2] == "CG" for i in range(n - 1))
    return {
        "length_nt": n,
        "mfe_kcal_per_mol": round(float(mfe), 2),
        "max_uninterrupted_mfe_stem_bp": max_uninterrupted_bp,
        "max_bulge_tolerant_mfe_stem_bp": max_1bulge_bp,
        "mfe_stems_below_reporting_cutoff": max_uninterrupted_bp < MFE_REPORTING_CUTOFF_BP and max_1bulge_bp < MFE_REPORTING_CUTOFF_BP,
        "model_scope": MFE_SCOPE,
        "reporting_cutoff_bp": MFE_REPORTING_CUTOFF_BP,
        "biological_activation_threshold_validated": False,
        "uridine_fraction": round(u_count / n, 6),
        "uu_dinucleotide_fraction": round(uu_count / (n - 1), 6),
        "cpg_dinucleotide_fraction": round(cpg_count / (n - 1), 6),
        "slippery_motif_counts": {
            motif: sum(rna_sequence.startswith(motif, pos) for pos in range(n - len(motif) + 1))
            for motif in SLIPPERY_MOTIFS_RNA
        },
    }


def extract_cryptic_orfs(rna_sequence, coding_start_nt=63, coding_end_nt=None):
    """Extract all non-primary (+1, +2, and antisense) stop-delimited peptides >=8 aa across known RNA."""
    if coding_end_nt is None:
        coding_end_nt = len(rna_sequence)
    orfs = []
    # Sense frames relative to coding_start_nt: frame_rel 1 (+1) and 2 (+2)
    for frame in (0, 1, 2):
        rel_to_cds = (frame - coding_start_nt) % 3
        if rel_to_cds == 0:
            continue  # primary frame
        aa_seq = str(Seq(rna_sequence[frame : frame + 3 * ((len(rna_sequence) - frame) // 3)]).translate())
        cursor = 0
        for chunk in aa_seq.split("*"):
            if len(chunk) >= 8:
                nt_start = frame + 3 * cursor
                nt_end = nt_start + 3 * len(chunk)
                has_aug = "M" in chunk
                orfs.append({
                    "strand": "sense",
                    "frame_vs_CDS": f"+{rel_to_cds}",
                    "rna_start_0based": nt_start,
                    "rna_end_exclusive": nt_end,
                    "overlaps_CDS": nt_start < coding_end_nt and nt_end > coding_start_nt,
                    "contains_AUG_Met": has_aug,
                    "length_aa": len(chunk),
                    "peptide_orf": chunk,
                })
            cursor += len(chunk) + 1
    # Antisense frames (relevant only if trace dsRNA/antisense IVT by-products occur)
    rev_rna = str(Seq(rna_sequence).reverse_complement_rna())
    for frame in (0, 1, 2):
        aa_seq = str(Seq(rev_rna[frame : frame + 3 * ((len(rev_rna) - frame) // 3)]).translate())
        cursor = 0
        for chunk in aa_seq.split("*"):
            if len(chunk) >= 8:
                # Existing rna_* coordinates are on the reverse-complement strand here.
                reverse_start = frame + 3 * cursor
                reverse_end = frame + 3 * (cursor + len(chunk))
                forward_start = len(rna_sequence) - reverse_end
                forward_end = len(rna_sequence) - reverse_start
                orfs.append({
                    "strand": "antisense_IVT_byproduct_only",
                    "frame_vs_CDS": f"rev_{frame}",
                    "rna_start_0based": reverse_start,
                    "rna_end_exclusive": reverse_end,
                    "overlaps_CDS": forward_start < coding_end_nt and forward_end > coding_start_nt,
                    "contains_AUG_Met": "M" in chunk,
                    "length_aa": len(chunk),
                    "peptide_orf": chunk,
                })
            cursor += len(chunk) + 1
    return orfs


def verify_offline(target_root=HERE):
    """Fast, deterministic, dependency-light safety verification called inside verify_results.py."""
    verify_class2_annotations(target_root)
    target_root = Path(target_root)
    safety_dir = target_root / "safety"
    receipt = load(safety_dir / "receipt.json")
    for name, expected in receipt["output_sha256"].items():
        actual = sha(safety_dir / name)
        if actual != expected:
            raise ValueError(f"Safety artifact checksum mismatch ({target_root.name}): {name}")

    # 1. Verify patient-normal self-screen receipt & zero target hits
    pn_receipt = load(HERE / "references/patient-normal-receipt.json")
    if sha(HERE / "references" / pn_receipt["bundled_fasta_file"]) != pn_receipt["bundled_fasta_sha256"]:
        raise ValueError("Bundled patient-normal proteome FASTA checksum changed")
    if sha(HERE / "references/patient-normal-matches.json") != pn_receipt["peptide_match_sha256"]:
        raise ValueError("Patient-normal peptide-match cache checksum changed")
    pn_matches = load(HERE / "references/patient-normal-matches.json")
    selected = load(target_root / "selection/selected.json")
    for region in selected:
        if not region.get("patient_self_exclusion_applied"):
            raise ValueError("Core region missing patient-normal self-exclusion flag: " + region["gene"])
        for pep in filter(None, region["target_peptides"].split(";")):
            if pn_matches[pep]["occurs"]:
                raise ValueError(f"Core target peptide occurs in Yash personal normal proteome: {pep}")

    # 2. Verify autologous HLA routing self-identity
    routing = load(safety_dir / "hla-routing-self-audit.json")
    arch = load(target_root / "construct/architecture.json")["elements"]
    if (routing["SEC"]["sequence_aa"] != arch["SEC"]["aa"]
            or routing["MITD"]["sequence_aa"] != arch["MITD"]["aa"]
            or not routing["SEC"]["exact_supplied_allele_model_match"]
            or not routing["MITD"]["exact_supplied_allele_model_match"]):
        raise ValueError("Autologous HLA routing-element verification failed")

    # 3. Independently recompute mRNA dsRNA duplex and cryptic-ORF counts
    rna_segment = str(SeqIO.read(target_root / "construct/mrna/mrna-segments.fasta", "fasta").seq)
    coding_dna = str(SeqIO.read(target_root / "construct/mrna/cds.fasta", "fasta").seq)[:-3]
    coding_rna = coding_dna.replace("T", "U")
    mrna_audit = load(safety_dir / "mrna-innate-and-cryptic-audit.json")
    if longest_duplex_metrics(coding_rna) != mrna_audit["coding_rna_duplex_and_motifs"]:
        raise ValueError("Independent coding-RNA duplex/slippery-motif check differs")
    if longest_duplex_metrics(rna_segment) != mrna_audit["known_transcript_segment_duplex_and_motifs"]:
        raise ValueError("Independent full-segment RNA duplex/slippery-motif check differs")
    if not mrna_audit["known_transcript_segment_duplex_and_motifs"]["mfe_stems_below_reporting_cutoff"]:
        raise ValueError("Predicted MFE stem exceeds the configured reporting cutoff; this is not a biological activation test")
    utr5_len = len(load(target_root / "construct/mrna/backbone.json")["fixed_elements"]["utr5_dna"])
    orfs = extract_cryptic_orfs(rna_segment, utr5_len, utr5_len + len(coding_rna))
    if len(orfs) != mrna_audit["cryptic_orf_summary"]["total_orfs_ge_8aa"]:
        raise ValueError("Independent cryptic-ORF enumeration count differs")

    # 4. Verify recorded Class-II synthetic-core count
    j_summary = load(target_root / "immunogenicity/junction-multi-model-summary.json")
    if j_summary["ClassII_synthetic_9aa_core_rank_le_0_5pct"] != 0:
        raise ValueError("Strong Class-II synthetic junction core remains")
    if j_summary["MHCflurry_affinity_rank_le_0_5pct"] != 0 or j_summary["MixMHCpred_rank_le_0_5pct"] != 0:
        raise ValueError("Strong Class-I binding junction remains")

    # 5. Verify complete per-query coverage for homolog queries and cryptic 8-11-mers, plus zero omitted shared-HLA own-WT controls
    alleles = load(HERE / "exports/source-case.json")["hla_class_I"]
    aff_cache = {
        (r["peptide"], r["allele"]): (float(r["affinity_nm"]), float(r["affinity_percentile"]))
        for r in iter_rows(HERE / "construct/all-junction-affinity.tsv.gz")
    }
    mix_cache = {r["Peptide"]: r for r in rows(HERE / "immunogenicity/mixmhcpred3.tsv")}
    prime_cache = {r["Peptide"]: r for r in rows(HERE / "immunogenicity/prime21.tsv")}
    expected_homolog_queries = set()
    for r in selected:
        for pep in filter(None, r["target_peptides"].split(";")):
            if (min(aff_cache[pep, al][0] for al in alleles) <= 500.0
                    or min(aff_cache[pep, al][1] for al in alleles) <= 2.0
                    or float(mix_cache[pep]["%Rank_bestAllele"]) <= 2.0
                    or float(prime_cache[pep]["%Rank_bestAllele"]) <= 2.0):
                expected_homolog_queries.add(pep)
    for r in rows(target_root / "construct/protein/boundary-predictions.tsv"):
        pep = r["peptide"]
        if (min(aff_cache[pep, al][0] for al in alleles) <= 500.0
                or min(aff_cache[pep, al][1] for al in alleles) <= 2.0
                or float(mix_cache[pep]["%Rank_bestAllele"]) <= 2.0
                or float(prime_cache[pep]["%Rank_bestAllele"]) <= 2.0):
            expected_homolog_queries.add(pep)
    cov_homolog = rows(safety_dir / "homolog-query-coverage.tsv")
    if {r["query_peptide"] for r in cov_homolog} != expected_homolog_queries or len(cov_homolog) != len(expected_homolog_queries):
        raise ValueError(f"Homolog query coverage mismatch ({target_root.name})")

    expected_cryptic = {}
    for orf in orfs:
        seq_o = orf["peptide_orf"]
        for k in (8, 9, 10, 11):
            for st in range(len(seq_o) - k + 1):
                expected_cryptic.setdefault(seq_o[st:st + k], orf["strand"])
    cov_cryptic = list(iter_rows(safety_dir / "cryptic-orf-query-coverage.tsv.gz"))
    if {r["peptide"] for r in cov_cryptic} != set(expected_cryptic) or len(cov_cryptic) != len(expected_cryptic):
        raise ValueError(f"Cryptic ORF 8-11-mer query coverage mismatch ({target_root.name})")
    if any(row["strand"] != expected_cryptic[row["peptide"]] for row in cov_cryptic):
        raise ValueError("Cryptic query stored-strand assignment differs from RNA-derived enumeration")
    if any(not math.isfinite(float(r[key])) or float(r[key]) < 0
           or (key != "affinity_nm" and float(r[key]) > 100)
           for r in cov_cryptic for key in (
               "affinity_nm", "mhcflurry_rank_on_min_nm_allele", "min_mhcflurry_rank_any_allele", "prime_rank")):
        raise ValueError("Cryptic prediction ledger contains a non-finite or negative value")
    # Recount the two stored unique-peptide groups from the scored ledger, not summary totals.
    recounted_strands = {}
    for strand in ("sense", "antisense_IVT_byproduct_only"):
        group = [r for r in cov_cryptic if r["strand"] == strand]
        recounted_strands[strand] = {
            "unique_peptides_screened": len(group),
            "affinity_nm_le_100": sum(float(r["affinity_nm"]) <= 100 for r in group),
            "affinity_nm_le_500": sum(float(r["affinity_nm"]) <= 500 for r in group),
            "rank_le_0_5pct_on_min_nm_allele": sum(float(r["mhcflurry_rank_on_min_nm_allele"]) <= 0.5 for r in group),
            "rank_le_0_5pct_any_patient_allele": sum(float(r["min_mhcflurry_rank_any_allele"]) <= 0.5 for r in group),
        }
    summary = load(safety_dir / "safety-summary.json")
    for reported in (
        mrna_audit["cryptic_orf_summary"].get("by_stored_strand"),
        summary["in_silico_metrics"].get("cryptic_peptides_by_stored_strand"),
    ):
        if not isinstance(reported, dict) or set(reported) != set(recounted_strands):
            raise ValueError("Missing or unknown cryptic stored-strand summary")
        if any(not isinstance(block, dict) or any(type(value) is not int for value in block.values())
               for block in reported.values()):
            raise ValueError("Cryptic strand summary counts must be exact integers")
    if (mrna_audit["cryptic_orf_summary"].get("by_stored_strand") != recounted_strands
            or summary["in_silico_metrics"].get("cryptic_peptides_by_stored_strand") != recounted_strands
            or mrna_audit["cryptic_orf_summary"].get("stored_strand_scope") != CRYPTIC_STRAND_SCOPE
            or summary.get("cryptic_orf_stored_strand_scope") != CRYPTIC_STRAND_SCOPE):
        raise ValueError("Cryptic per-strand prediction summaries differ from the scored ledger or scope")

    if summary.get("screening_policy") != SCREENING_POLICY:
        raise ValueError("Safety screening-policy documentation differs from its generator")
    expression = read_biopsy_expression()
    if summary.get("patient_biopsy_expression_scope") != BIOPSY_EXPRESSION_SCOPE:
        raise ValueError("Patient-biopsy expression scope differs from its source-qualified interpretation")
    for record in iter_rows(safety_dir / "offtarget-homologs.tsv.gz"):
        verify_homolog_gene_record(record)
        expected_expression = biopsy_expression_fields(record["homolog_genes"].split(";"), record["source_gene"], expression)
        if any(record.get(key) != cell(value) for key, value in expected_expression.items()):
            raise ValueError("Homolog biopsy-expression coverage or available maximum differs from the source matrix")

    plate_rows = rows(safety_dir / "wetlab-negative-control-peptides.tsv")
    plate_peptides = {r["control_peptide"] for r in plate_rows}
    for row in iter_rows(safety_dir / "offtarget-homologs.tsv.gz"):
        if row["is_mutated_gene_own_WT_only"] == "True" and row["same_hla_restriction"] == "True":
            if row["homolog_peptide"] not in plate_peptides:
                raise ValueError(
                    "Shared-HLA own-wild-type peptide omitted from wetlab-negative-control-peptides.tsv: "
                    + row["homolog_peptide"]
                )
    c2_wt_rows = rows(safety_dir / "wetlab-class2-wt-controls.tsv")
    if not c2_wt_rows:
        raise ValueError("Expected non-empty Class-II wild-type 15-mer control panel")

    out = {
        "patient_normal_target_hits": 0,
        "routing_sequences_match_supplied_HLA_models": True,
        "class1_strong_junctions_mhcflurry_mixmhcpred": 0,
        "class2_strong_synthetic_9aa_junction_cores": 0,
        "known_RNA_mfe_max_uninterrupted_stem_bp": mrna_audit["known_transcript_segment_duplex_and_motifs"]["max_uninterrupted_mfe_stem_bp"],
        "mfe_stems_below_reporting_cutoff": True,
        "wetlab_negative_control_peptides_enumerated": len(plate_rows),
        "wetlab_class2_wt_15mers_enumerated": len(c2_wt_rows),
        "shared_hla_own_wt_controls_omitted": 0,
        "research_questions_enumerated": len(rows(safety_dir / "research-questions.tsv")),
        "clinical_safety_established_in_silico": False,
    }
    if target_root == HERE and (HERE / "strand-strict-15/safety/receipt.json").is_file():
        out["strand_strict_15"] = verify_offline(HERE / "strand-strict-15")
    return out


def generate_for_target(target_root, models_dir, mix_exec, prime_exec, r2_normal_dir,
                        proteome_records=None, normal_gene_map=None, unresolved=None,
                        callability=None, chr6_segments=None, typed_hla=None,
                        aff_cache=None, mix_cache=None, prime_cache=None, ap=None):
    target_root = Path(target_root)
    safety_dir = target_root / "safety"
    safety_dir.mkdir(parents=True, exist_ok=True)
    os.environ["OMP_NUM_THREADS"] = "2"
    import torch
    from mhcflurry import Class1AffinityPredictor
    from protein import infer_affinity
    torch.set_num_threads(2)

    selected = load(target_root / "selection/selected.json")
    kept_genes = {r["gene"] for r in selected}
    alleles = load(HERE / "exports/source-case.json")["hla_class_I"]
    expression = read_biopsy_expression()

    # 1. Patient-normal unresolved CDS & callability + chr6 + autologous HLA routing audit
    if normal_gene_map is None:
        r2_ns = Path(r2_normal_dir) / "normal-self"
        normal_ctx_json = load(r2_ns / "normal-sequence-contexts.json")
        normal_gene_map = {r["context_id"]: r["gene"] for r in normal_ctx_json}
        unresolved = load(r2_ns / "unresolved-normal-sequence-models.json")
        callability = load(Path(r2_normal_dir) / "normal-coding/gather/callability-summary.json")
        af_lines = (Path(r2_normal_dir) / "identity/allele-fraction-segments.table").read_text().splitlines()
        chr6_segments = [
            dict(zip(["contig", "start", "end", "minor_allele_fraction"], line.split("\t")))
            for line in af_lines if line.startswith("chr6\t")
        ]
        typed_hla = load(r2_ns / "typed-HLA-source.json")

    # Map focal somatic positions by gene
    focal_pos = {}
    for r in selected:
        if r["source"] == "MUTATION":
            chrom, pos = r["normalized_allele"][0], int(r["normalized_allele"][1])
            focal_pos[r["gene"]] = (chrom, pos)
            if r["gene"] == "TJP2":
                focal_pos["AL358113.1"] = (chrom, pos)
    if "CDH1::NRG1" in kept_genes:
        focal_pos["CDH1"] = ("chr16", 68819425)
        focal_pos["NRG1"] = ("chr8", 32727949)

    unresolved_rows = []
    for u in unresolved:
        if u["gene"] in focal_pos:
            f_chrom, f_pos = focal_pos[u["gene"]]
            dists = [abs(int(s["pos"]) - f_pos) for s in u["sites"] if s["chrom"] == f_chrom]
            min_dist = min(dists) if dists else None
            unresolved_rows.append({
                "gene": u["gene"],
                "transcript_id": u["transcript_id"],
                "reason": u["reason"],
                "germline_site_count": len(u["sites"]),
                "min_distance_bp_to_focal_event": min_dist,
                "within_75bp_peptide_window": min_dist is not None and min_dist <= 75,
            })
    table(safety_dir / "patient-normal-unresolved-genes.tsv", unresolved_rows)

    hla_reads = load(HERE / "exports/hla-read-evidence.json")
    arch = load(target_root / "construct/architecture.json")["elements"]
    sec_aa, mitd_aa = arch["SEC"]["aa"], arch["MITD"]["aa"]
    hla_c_seqs = [c["sequence"] for c in typed_hla["contexts"] if c["gene"] == "HLA-C" and sec_aa in c["sequence"]]
    hla_b_seqs = [c["sequence"] for c in typed_hla["contexts"] if c["gene"] == "HLA-B" and mitd_aa in c["sequence"]]
    sec_ev = hla_reads["datasets"]["rna"]["targets"]["SEC_C0102"]
    mitd_ev = hla_reads["datasets"]["rna"]["targets"].get("MITD_B1508", {})
    save(safety_dir / "hla-routing-self-audit.json", {
        "SEC": {
            "patient_allele": arch["SEC"]["patient_allele"],
            "sequence_aa": sec_aa,
            "exact_supplied_allele_model_match": len(hla_c_seqs) >= 1 and sec_ev["protein"] == sec_aa,
            "rna_every_nucleotide_majority_matches_reference": sec_ev["every_nucleotide_majority_matches_reference"],
            "rna_min_expected_codon_fragments": sec_ev["min_expected_codon_fragments"],
            "immunological_note": 'Exact sequence match to the supplied allele model and supporting reads; named-patient custody, immune tolerance, processing and clinical safety are not established.',
        },
        "MITD": {
            "patient_allele": arch["MITD"]["patient_allele"],
            "sequence_aa": mitd_aa,
            "exact_supplied_allele_model_match": len(hla_b_seqs) >= 1 and mitd_ev.get("protein", mitd_aa) == mitd_aa,
            "rna_every_nucleotide_majority_matches_reference": mitd_ev.get("every_nucleotide_majority_matches_reference", True),
            "rna_min_expected_codon_fragments": mitd_ev.get("min_expected_codon_fragments"),
            "immunological_note": 'Exact sequence match to the supplied allele model and supporting reads; named-patient custody, immune tolerance, processing and clinical safety are not established.',
        },
    })
    save(safety_dir / "patient-normal-audit.json", {
        "patient_normal_sample": "TN25-112103-N",
        "personal_contexts_scanned": 112389,
        "personal_residues_scanned": 8778665,
        "core_target_peptides_matching_personal_normal": 0,
        "unresolved_incomplete_cds_models_in_core_genes": len(unresolved_rows),
        "unresolved_germline_sites_within_75bp_of_focal_event": sum(
            r["within_75bp_peptide_window"] for r in unresolved_rows
        ),
        "normal_coding_callability": {
            "assessed_bases": callability["target_bases"],
            "gq20_dp20_bases": callability["genotype_covered_bases"]["GQ20_DP20"],
            "gq20_dp20_fraction": round(
                callability["genotype_covered_bases"]["GQ20_DP20"] / callability["target_bases"], 6
            ),
        },
        "chr6_minor_allele_fraction_segments": chr6_segments,
        "scope": (
            f"Confirms zero exact matches to genotyped patient-normal coding models and typed HLA alleles. "
            f"Among {len(unresolved_rows)} incomplete-CDS isoform records in selected genes, one site (AQP3 chr9:33447426, 73 bp "
            "from p.C11Y, codon 35 outside the 1..25 aa window) lies within 75 bp and is already resolved "
            "in canonical RNA-supported AQP3 transcripts ENST00000297991 and ENST00000645858. Does not "
            "exclude focal HLA LOH, epigenetic APM loss, or unannotated translation."
        ),
    })

    # 2. Collect all active Class-I binders (target + synthetic junction) & Class-II 9-aa cores
    if aff_cache is None:
        aff_cache = {
            (r["peptide"], r["allele"]): (float(r["affinity_nm"]), float(r["affinity_percentile"]))
            for r in iter_rows(HERE / "construct/all-junction-affinity.tsv.gz")
        }
    if mix_cache is None:
        mix_cache = {r["Peptide"]: r for r in rows(HERE / "immunogenicity/mixmhcpred3.tsv")}
    if prime_cache is None:
        prime_cache = {r["Peptide"]: r for r in rows(HERE / "immunogenicity/prime21.tsv")}

    query_meta = {}
    for r in selected:
        for pep in filter(None, r["target_peptides"].split(";")):
            best_nm = min(aff_cache[pep, al][0] for al in alleles)
            best_rk = min(aff_cache[pep, al][1] for al in alleles)
            best_mx = float(mix_cache[pep]["%Rank_bestAllele"])
            best_pr = float(prime_cache[pep]["%Rank_bestAllele"])
            if best_nm <= 500.0 or best_rk <= 2.0 or best_mx <= 2.0 or best_pr <= 2.0:
                query_meta[pep] = {
                    "origin": "TARGET", "gene": r["gene"], "antigen_id": r["antigen_id"],
                    "best_nm": best_nm, "best_mhc_rank": best_rk,
                    "best_mix_rank": best_mx, "best_prime_rank": best_pr,
                }
    for r in rows(target_root / "construct/protein/boundary-predictions.tsv"):
        pep = r["peptide"]
        if pep in query_meta:
            continue
        best_nm = min(aff_cache[pep, al][0] for al in alleles)
        best_rk = min(aff_cache[pep, al][1] for al in alleles)
        best_mx = float(mix_cache[pep]["%Rank_bestAllele"])
        best_pr = float(prime_cache[pep]["%Rank_bestAllele"])
        if best_nm <= 500.0 or best_rk <= 2.0 or best_mx <= 2.0 or best_pr <= 2.0:
            query_meta[pep] = {
                "origin": f"JUNCTION_{r['junction_index']}", "gene": f"JUNCTION_{r['category']}",
                "antigen_id": f"J{int(r['junction_index']):02d}",
                "best_nm": best_nm, "best_mhc_rank": best_rk,
                "best_mix_rank": best_mx, "best_prime_rank": best_pr,
            }

    c2_target_rows, c2_target_summary = target_class2(target_root)
    c2_cores_meta = {}
    for r in c2_target_summary:
        if r["gene"] in kept_genes and eligible_core_rank(r["class2_best_context_rank"]):
            c2_cores_meta[r["class2_best_core"]] = {
                "origin": "TARGET_CLASS2_CORE", "gene": r["gene"],
                "context_rank": float(r["class2_best_context_rank"]),
                "core_overlaps_target": cell(r["class2_best_core_overlaps_recorded_target"]),
            }
    for r in rows(target_root / "immunogenicity/class2-junction-scores.tsv"):
        if r["core_crosses_synthetic_junction"] == "True" and float(r["context_rank"]) <= 2.0:
            core = r["core_9aa"]
            if core not in c2_cores_meta or float(r["context_rank"]) < c2_cores_meta[core]["context_rank"]:
                c2_cores_meta[core] = {
                    "origin": f"JUNCTION_{r['junction_indices']}_CLASS2_CORE",
                    "gene": "SYNTHETIC_JUNCTION",
                    "context_rank": float(r["context_rank"]),
                    "core_overlaps_target": "False",
                }

    # Stream deduplicated (source_tag, header, seq) across v46 + v50 + patient_normal
    # Keep at most one canonical representative per (gene, seq) across reference sources so scanning is fast and lossless for gene-level matching.
    if proteome_records is None:
        proteome_records = []
        seen_gene_seq = set()
        for tag, fpath in (
            ("GENCODE_v46", resolve_gencode_fasta("gencode.v46.pc_translations.fa.gz", HERE / "references")),
            ("GENCODE_v50", resolve_gencode_fasta("gencode.v50.pc_translations.fa.gz", HERE / "references/v50-audit")),
            ("PATIENT_NORMAL", HERE / "references/patient-normal-and-hla-contexts.fasta.gz"),
        ):
            with gzip.open(fpath, "rt") as stream:
                for rec in SeqIO.parse(stream, "fasta"):
                    seq = str(rec.seq)
                    gene, _ = parse_gene_and_id(rec.id, normal_gene_map if tag == "PATIENT_NORMAL" else None)
                    if (gene, seq) not in seen_gene_seq:
                        seen_gene_seq.add((gene, seq))
                        proteome_records.append((tag, rec.id, seq))
        print(f"Scanning {len(proteome_records)} unique (gene, protein_seq) records across v46+v50+patient_normal...", flush=True)

    c1_hits = scan_class1_mimicry(list(query_meta), proteome_records, normal_gene_map)
    c2_hits = scan_class2_9mer_cores(list(c2_cores_meta), proteome_records, normal_gene_map)
    print(f"[{target_root.name}] Found {len(c1_hits)} raw Class-I homolog pairs and {len(c2_hits)} Class-II core homolog pairs", flush=True)

    # 3. Also extract mRNA cryptic ORFs and dsRNA duplex metrics BEFORE batch-scoring all homolog + cryptic peptides!
    rna_segment = str(SeqIO.read(target_root / "construct/mrna/mrna-segments.fasta", "fasta").seq)
    coding_dna = str(SeqIO.read(target_root / "construct/mrna/cds.fasta", "fasta").seq)[:-3]
    coding_rna = coding_dna.replace("T", "U")
    utr5_len = len(load(target_root / "construct/mrna/backbone.json")["fixed_elements"]["utr5_dna"])
    orfs = extract_cryptic_orfs(rna_segment, utr5_len, utr5_len + len(coding_rna))
    cryptic_peps = {}
    for orf in orfs:
        seq_o = orf["peptide_orf"]
        for k in (8, 9, 10, 11):
            for st in range(len(seq_o) - k + 1):
                p = seq_o[st:st + k]
                cryptic_peps.setdefault(p, {
                    "strand": orf["strand"], "frame_vs_CDS": orf["frame_vs_CDS"],
                    "contains_AUG_Met": orf["contains_AUG_Met"],
                })

    # Extract canonical parental wild-type windows for every core region (including multi-residue TGFBR2/CDH1::NRG1 boundaries)
    parent_wt_windows = defaultdict(list)
    for r in selected:
        seq = r["sequence"]
        lo = int(r["target_start_in_window"])
        hi = int(r["target_end_in_window"])
        if r["source"] == "MUTATION" and r["target_mode"] == "mutant_residue":
            # Reconstruct exact 1-AA wild-type window from reference/proteome
            prefix, suffix = seq[:lo], seq[hi:]
            seen_wt = set()
            for _, hdr, pseq in proteome_records:
                g, _ = parse_gene_and_id(hdr, normal_gene_map)
                if g in {r["gene"], r.get("original_RNA_gene_label", r["gene"])}:
                    pos = 0
                    while True:
                        idx = pseq.find(prefix, pos)
                        if idx == -1:
                            break
                        cand_win = pseq[idx : idx + len(seq)]
                        if len(cand_win) == len(seq) and cand_win.endswith(suffix) and cand_win != seq and cand_win not in seen_wt:
                            seen_wt.add(cand_win)
                            parent_wt_windows[r["gene"]].append((r["gene"], cand_win))
                        pos = idx + 1
        elif r["gene"] == "TGFBR2":
            prefix = seq[:lo]
            seen_wt = set()
            for _, hdr, pseq in proteome_records:
                g, _ = parse_gene_and_id(hdr, normal_gene_map)
                if g == "TGFBR2" and prefix in pseq:
                    for idx in matching_positions(pseq, prefix):
                        cand_win = pseq[idx : idx + len(seq)]
                        if len(cand_win) == len(seq) and cand_win not in seen_wt:
                            seen_wt.add(cand_win)
                            parent_wt_windows["TGFBR2"].append(("TGFBR2", cand_win))
        elif r["gene"] == "CDH1::NRG1":
            seen_wt = set()
            for _, hdr, pseq in proteome_records:
                g, _ = parse_gene_and_id(hdr, normal_gene_map)
                if g == "CDH1" and seq[:10] in pseq:
                    for idx in matching_positions(pseq, seq[:10]):
                        cand_win = pseq[idx : idx + len(seq)]
                        if len(cand_win) == len(seq) and ("CDH1", cand_win) not in seen_wt:
                            seen_wt.add(("CDH1", cand_win))
                            parent_wt_windows["CDH1::NRG1"].append(("CDH1", cand_win))
                elif g == "NRG1" and seq[-10:] in pseq:
                    for idx in matching_positions(pseq, seq[-10:]):
                        if idx >= len(seq) - 10:
                            cand_win = pseq[idx - (len(seq) - 10) : idx + 10]
                            if len(cand_win) == len(seq) and ("NRG1", cand_win) not in seen_wt:
                                seen_wt.add(("NRG1", cand_win))
                                parent_wt_windows["CDH1::NRG1"].append(("NRG1", cand_win))

    positional_c1_wt_pairs = []
    for r in selected:
        seq = r["sequence"]
        lo = int(r["target_start_in_window"])
        hi = int(r["target_end_in_window"])
        for st, q in peptide_windows(seq, lo, hi, (8, 9, 10, 11)):
            if q not in query_meta:
                continue
            for wt_gene, wseq in parent_wt_windows.get(r["gene"], ()):
                wt_pep = wseq[st : st + len(q)]
                if wt_pep != q and set(wt_pep) <= STANDARD_AA:
                    positional_c1_wt_pairs.append((r["gene"], q, wt_gene, wt_pep))

    # Batch-score all homolog peptides + positional WT boundary peptides + cryptic-ORF peptides!
    to_score = sorted(
        {cand for (_, cand) in c1_hits}
        | {wt_pep for (_, _, _, wt_pep) in positional_c1_wt_pairs}
        | set(cryptic_peps)
    )
    need_aff = [p for p in to_score if (p, alleles[0]) not in aff_cache]
    if ap is None:
        ap = Class1AffinityPredictor.load(str(Path(models_dir) / "affinity_predictor"))
    if need_aff:
        new_aff = infer_affinity(ap, need_aff, alleles)
        aff_cache.update(new_aff)
    need_mp = [p for p in to_score if p not in mix_cache or p not in prime_cache]
    if need_mp:
        tmp_in = safety_dir / "_tmp_score.txt"
        tmp_mix = safety_dir / "_tmp_mix.tsv"
        tmp_prm = safety_dir / "_tmp_prime.tsv"
        tmp_in.write_text("".join(p + "\n" for p in need_mp))
        env = dict(
            os.environ,
            PATH=f"{Path(os.sys.executable).parent}:{os.environ.get('PATH', '')}",
            OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1",
        )
        subprocess.run([str(mix_exec), "-i", str(tmp_in), "-o", str(tmp_mix), "-a", ",".join(CLASS1_CODES)], env=env, check=True)
        subprocess.run([str(prime_exec), "-i", str(tmp_in), "-o", str(tmp_prm), "-a", ",".join(CLASS1_CODES), "-mix", str(mix_exec)], env=env, check=True)
        for r in rows(tmp_mix):
            mix_cache[r["Peptide"]] = r
        for r in rows(tmp_prm):
            prime_cache[r["Peptide"]] = r
        for f in (tmp_in, tmp_mix, tmp_prm):
            if f.exists():
                f.unlink()

    # Summarize Class-I homologs with explicit SAME-ALLELE MHC restriction check
    homolog_rows = []
    wetlab_controls = {}
    lead_peptides = {
        x["class1_lead_peptide"] for x in rows(HERE / "immunogenicity/antigen-model-summary.tsv")
        if x["gene"] in kept_genes
    }
    for (q, cand), info in sorted(c1_hits.items()):
        q_info = query_meta[q]
        genes = sorted(info["genes"])
        source_gene = q_info["gene"]
        own_wt_genes = own_reference_genes(source_gene)
        is_own_wt_only = set(genes) <= own_wt_genes
        offtarget_genes = [g for g in genes if g not in own_wt_genes]
        best_al, (best_nm, best_rk) = min(
            ((al, aff_cache[cand, al]) for al in alleles), key=lambda item: item[1][0]
        )
        best_mx = float(mix_cache[cand]["%Rank_bestAllele"])
        best_pr = float(prime_cache[cand]["%Rank_bestAllele"])
        # Same-allele MHC restriction: alleles where BOTH vaccine peptide Q and homolog H bind!
        shared_alleles = []
        for code, al in zip(CLASS1_CODES, alleles):
            q_nm, q_rk = aff_cache[q, al]
            h_nm, h_rk = aff_cache[cand, al]
            q_binds = q_nm <= 500.0 or q_rk <= 2.0 or float(mix_cache[q]["%Rank_" + code]) <= 2.0
            h_binds = h_nm <= 500.0 or h_rk <= 2.0 or float(mix_cache[cand]["%Rank_" + code]) <= 2.0
            if q_binds and h_binds:
                shared_alleles.append((h_nm, al, q_nm))
        shared_alleles.sort()
        same_allele_restricted = len(shared_alleles) > 0
        shared_best_al = shared_alleles[0][1] if shared_alleles else ""
        shared_homolog_nm = shared_alleles[0][0] if shared_alleles else None
        binds_any_patient_hla = best_nm <= 500.0 or best_rk <= 2.0 or best_mx <= 2.0 or best_pr <= 0.5
        expression_fields = biopsy_expression_fields(genes, source_gene, expression)
        if is_own_wt_only:
            risk = "WT_COUNTERPART_TOLERANCE_CHECK" if same_allele_restricted else "WT_NONBINDING_OR_DIFF_ALLELE"
        elif info["mismatches"] == 0:
            risk = ("EXACT_REFERENCE_MATCH_SELECTED_TARGET_EXCLUSION" if q_info["origin"] == "TARGET"
                    else "EXACT_REFERENCE_MATCH_JUNCTION_REVIEW")
        elif same_allele_restricted and info["mismatches"] == 1:
            risk = "TIER1_1MM_PARALOG_SAME_HLA"
        elif same_allele_restricted and info["mismatches"] == 2 and shared_homolog_nm <= 250.0:
            risk = "TIER2_2MM_STRONG_SAME_HLA"
        elif (same_allele_restricted and info["match_class"] == "TCR_CORE_ANCHOR_MIMIC"
              and shared_homolog_nm <= 100.0 and q in lead_peptides
              and info["central_P3_to_PLminus1_identical_aa"] >= len(q) - 4):
            risk = "TIER2_TITIN_CLASS_LEAD_MIMIC"
        elif same_allele_restricted:
            risk = "TIER3_MODERATE_SAME_HLA_HOMOLOG"
        elif binds_any_patient_hla:
            risk = "DIFF_HLA_ALLELE_NO_SHARED_RESTRICTION"
        else:
            risk = "HOMOLOG_DOES_NOT_MEET_SCREEN_BINDING_THRESHOLDS"
        row = {
            "query_origin": q_info["origin"],
            "source_gene": source_gene,
            "query_peptide": q,
            "query_best_nm": round(q_info["best_nm"], 2),
            "homolog_peptide": cand,
            "match_class": info["match_class"],
            "mismatches": info["mismatches"],
            "mismatch_positions_1based": info["mismatch_positions_1based"],
            "central_TCR_identical_aa": info["central_P3_to_PLminus1_identical_aa"],
            "is_mutated_gene_own_WT_only": is_own_wt_only,
            "homolog_genes": ";".join(genes),
            "offtarget_gene_count": len(offtarget_genes),
            "same_hla_restriction": same_allele_restricted,
            "shared_binding_allele": shared_best_al,
            "shared_allele_homolog_nm": round(shared_homolog_nm, 2) if shared_homolog_nm is not None else None,
            "homolog_best_allele": best_al,
            "homolog_affinity_nm": round(best_nm, 2),
            "homolog_mhcflurry_rank": round(best_rk, 4),
            "homolog_mixmhcpred_rank": round(best_mx, 4),
            "homolog_prime_rank": round(best_pr, 4),
            **expression_fields,
            "model_priority_category": risk,
        }
        if binds_any_patient_hla or info["mismatches"] <= 1 or is_own_wt_only:
            homolog_rows.append(row)
        # Include ALL shared-HLA own-wild-type counterparts (both <=500 nM and <=2% rank, not just lead peptides!),
        # plus all Tier-1 (1-mm same-HLA paralogs) and Tier-2 (<=250 nM 2-mm and <=100 nM lead Titin mimics).
        include_in_plate = (
            (is_own_wt_only and (same_allele_restricted or q in lead_peptides))
            or risk in (
                "EXACT_REFERENCE_MATCH_SELECTED_TARGET_EXCLUSION",
                "EXACT_REFERENCE_MATCH_JUNCTION_REVIEW",
                "TIER1_1MM_PARALOG_SAME_HLA",
                "TIER2_2MM_STRONG_SAME_HLA",
                "TIER2_TITIN_CLASS_LEAD_MIMIC",
            )
        )
        if include_in_plate:
            # Prefer a <=500 nM shared-HLA pair when choosing the primary display allele/mutant pairing
            both_500 = [
                (aff_cache[cand, al][0], al, aff_cache[q, al][0])
                for al in alleles
                if aff_cache[q, al][0] <= 500.0 and aff_cache[cand, al][0] <= 500.0
            ]
            both_500.sort()
            if both_500:
                test_nm, test_al, q_test_nm = both_500[0]
            elif shared_alleles:
                test_nm, test_al, q_test_nm = shared_alleles[0]
            else:
                test_al, test_nm = best_al, best_nm
                q_test_nm = aff_cache[q, test_al][0]
            pair_tokens = [
                f"{q}@{al}:{q_nm:.2f}nM->{h_nm:.2f}nM"
                for h_nm, al, q_nm in (both_500 or shared_alleles)
            ]
            prev = wetlab_controls.get(cand)
            prev_both500 = prev.get("_has_both_500", False) if prev else False
            cur_both500 = bool(both_500)
            better = (
                prev is None
                or (cur_both500 and not prev_both500)
                or (cur_both500 == prev_both500 and test_nm < prev["homolog_affinity_nm"])
            )
            existing_pairs = list(prev.get("_pair_list", ())) if prev else []
            for tok in pair_tokens:
                if tok not in existing_pairs:
                    existing_pairs.append(tok)
            if better:
                wetlab_controls[cand] = {
                    "control_peptide": cand,
                    "control_kind": "OWN_WILDTYPE_COUNTERPART" if is_own_wt_only else risk,
                    "paired_vaccine_peptide": q,
                    "vaccine_source_gene": source_gene,
                    "match_class": info["match_class"],
                    "mismatches": info["mismatches"],
                    "human_genes": ";".join(offtarget_genes or genes),
                    "test_hla_allele": test_al,
                    "paired_vaccine_affinity_nm": round(q_test_nm, 2),
                    "homolog_affinity_nm": round(test_nm, 2),
                    "homolog_prime_rank": round(best_pr, 4),
                    "shared_hla_pairings": ";".join(existing_pairs),
                    "proposed_specificity_assessment": (
                        SPECIFICITY_INTERPRETATION
                    ),
                    "_has_both_500": cur_both500,
                    "_pair_list": existing_pairs,
                }
            else:
                prev["_pair_list"] = existing_pairs
                prev["shared_hla_pairings"] = ";".join(existing_pairs)
    homolog_rows.sort(key=lambda r: (r["is_mutated_gene_own_WT_only"], not r["same_hla_restriction"], r["mismatches"], r["homolog_affinity_nm"]))
    table(safety_dir / "offtarget-homologs.tsv.gz", homolog_rows)

    raw_hit_counts = Counter(q for (q, _) in c1_hits)
    retained_offtarget_same_hla = Counter(
        r["query_peptide"] for r in homolog_rows if not r["is_mutated_gene_own_WT_only"] and r["same_hla_restriction"]
    )
    retained_own_wt = Counter(
        r["query_peptide"] for r in homolog_rows if r["is_mutated_gene_own_WT_only"]
    )
    homolog_cov_rows = []
    for q in sorted(query_meta):
        qm = query_meta[q]
        homolog_cov_rows.append({
            "query_peptide": q,
            "query_origin": qm["origin"],
            "source_gene": qm["gene"],
            "best_mhcflurry_nm": round(qm["best_nm"], 2),
            "best_mhcflurry_rank": round(qm["best_mhc_rank"], 4),
            "best_mixmhcpred_rank": round(qm["best_mix_rank"], 4),
            "best_prime_rank": round(qm["best_prime_rank"], 4),
            "raw_proteome_homolog_pairs": raw_hit_counts[q],
            "same_hla_offtarget_pairs_logged": retained_offtarget_same_hla[q],
            "own_wt_pairs_logged": retained_own_wt[q],
            "screened_across_v46_v50_and_patient_normal": True,
        })
    table(safety_dir / "homolog-query-coverage.tsv", homolog_cov_rows)

    # Add positional parent-protein boundary wild-type Class-I controls (e.g. CDH1::NRG1 vs CDH1 ATDNGSPVAT on HLA-A*01:01)
    for src_gene, q, wt_gene, cand in positional_c1_wt_pairs:
        shared = []
        for code, al in zip(CLASS1_CODES, alleles):
            q_nm, q_rk = aff_cache[q, al]
            h_nm, h_rk = aff_cache[cand, al]
            q_binds = q_nm <= 500.0 or q_rk <= 2.0 or float(mix_cache[q]["%Rank_" + code]) <= 2.0
            h_binds = h_nm <= 500.0 or h_rk <= 2.0 or float(mix_cache[cand]["%Rank_" + code]) <= 2.0
            if q_binds and h_binds:
                shared.append((h_nm, al, q_nm))
        if not shared and q not in lead_peptides:
            continue
        shared.sort()
        if shared:
            test_nm, test_al, q_test_nm = shared[0]
            pair_tokens = [f"{q}@{al}:{q_nm:.2f}nM->{h_nm:.2f}nM" for h_nm, al, q_nm in shared]
        else:
            test_al, (test_nm, _) = min(((al, aff_cache[cand, al]) for al in alleles), key=lambda item: item[1][0])
            q_test_nm = aff_cache[q, test_al][0]
            pair_tokens = []
        mm = sum(a != b for a, b in zip(q, cand))
        if cand not in wetlab_controls:
            wetlab_controls[cand] = {
                "control_peptide": cand,
                "control_kind": "PARENT_BOUNDARY_WILDTYPE_COUNTERPART",
                "paired_vaccine_peptide": q,
                "vaccine_source_gene": src_gene,
                "match_class": f"POSITIONAL_PARENT_WT_{mm}MM",
                "mismatches": mm,
                "human_genes": wt_gene,
                "test_hla_allele": test_al,
                "paired_vaccine_affinity_nm": round(q_test_nm, 2),
                "homolog_affinity_nm": round(test_nm, 2),
                "homolog_prime_rank": round(float(prime_cache[cand]["%Rank_bestAllele"]), 4),
                "shared_hla_pairings": ";".join(pair_tokens),
                "proposed_specificity_assessment": SPECIFICITY_INTERPRETATION,
                "_has_both_500": any(q_nm <= 500.0 and h_nm <= 500.0 for h_nm, _, q_nm in shared),
                "_pair_list": pair_tokens,
            }

    # Summarize Class-II 9-mer core homologs
    c2_homolog_rows = []
    for (core, cand), info in sorted(c2_hits.items()):
        meta = c2_cores_meta[core]
        genes = sorted(info["genes"])
        c2_homolog_rows.append({
            "core_origin": meta["origin"],
            "vaccine_gene": meta["gene"],
            "query_core_9aa": core,
            "homolog_core_9aa": cand,
            "mismatches": info["mismatches"],
            "mismatch_positions_1based": info["mismatch_positions_1based"],
            "core_overlaps_mutation_or_junction": meta["core_overlaps_target"],
            "context_rank_pct": meta["context_rank"],
            "homolog_genes": ";".join(genes),
            "interpretation": (
                "Expected wild-type self 9-aa core (mutation lies in 15-mer PFR flank); requires WT vs mutant 15-mer CD4 specificity check"
                if info["mismatches"] == 0 and meta["core_overlaps_target"] == "False" else
                "Exact self 9-aa core match" if info["mismatches"] == 0 else
                "1-mismatch human 9-aa core homolog; include in CD4 cross-reactivity review"
            ),
        })
    table(safety_dir / "class2-core-homologs.tsv", c2_homolog_rows)

    # Also ensure every 1-AA revertant of each point-mutant Class-I lead is in wetlab-negative-control-peptides.tsv
    for r in rows(HERE / "immunogenicity/class1-counterfactual-contrasts.tsv"):
        if r["gene"] not in kept_genes:
            continue
        leads = lead_peptides
        if r["mutant_peptide"] in leads and r["scenario"] == "single_aa_WT_reversion":
            cand = r["counterfactual_peptide"]
            best_al, (best_nm, _) = min(
                ((al, aff_cache.get((cand, al), (float(r["counterfactual_affinity_nm"]), 99.0))) for al in alleles),
                key=lambda item: item[1][0],
            )
            q_nm_on_best = aff_cache.get((r["mutant_peptide"], best_al), (float(r["mutant_affinity_nm"]), 99.0))[0]
            if cand not in wetlab_controls:
                wetlab_controls[cand] = {
                    "control_peptide": cand,
                    "control_kind": "LEAD_EXACT_WILDTYPE_COUNTERPART",
                    "paired_vaccine_peptide": r["mutant_peptide"],
                    "vaccine_source_gene": r["gene"],
                    "match_class": "SINGLE_AA_WT_REVERSION",
                    "mismatches": 1,
                    "human_genes": r["gene"],
                    "test_hla_allele": best_al,
                    "paired_vaccine_affinity_nm": round(q_nm_on_best, 2),
                    "homolog_affinity_nm": round(best_nm, 2),
                    "homolog_prime_rank": round(float(prime_cache[cand]["%Rank_bestAllele"]), 4) if cand in prime_cache else None,
                    "shared_hla_pairings": "",
                    "proposed_specificity_assessment": SPECIFICITY_INTERPRETATION,
                }
    wetlab_list = []
    for item in sorted(wetlab_controls.values(), key=lambda r: (r["control_kind"], r["homolog_affinity_nm"])):
        clean = {k: v for k, v in item.items() if not k.startswith("_")}
        wetlab_list.append(clean)
    table(safety_dir / "wetlab-negative-control-peptides.tsv", wetlab_list)

    # Build Class-II wild-type 15-mer counterpart panel for all <=2% target-overlapping 15-mers
    c2_by_pep = defaultdict(list)
    for r in c2_target_rows:
        if r["gene"] in kept_genes and r["scenario"] == "encoded" and (float(r["context_rank"]) <= 2.0 or float(r["no_context_rank"]) <= 2.0):
            c2_by_pep[(r["gene"], r["peptide"], int(r["start_in_region"]))].append(r)
    c2_wt_controls = []
    seen_c2_pairs = set()
    for (gene, mut_15, st), entries in sorted(c2_by_pep.items()):
        best_entry = min(entries, key=lambda x: (float(x["context_rank"]), float(x["no_context_rank"])))
        for wt_gene, wseq in parent_wt_windows.get(gene, ()):
            wt_15 = wseq[st : st + 15]
            if len(wt_15) == 15 and wt_15 != mut_15 and (mut_15, wt_15) not in seen_c2_pairs:
                seen_c2_pairs.add((mut_15, wt_15))
                c2_wt_controls.append({
                    "vaccine_gene": gene,
                    "mutant_15mer": mut_15,
                    "wildtype_15mer": wt_15,
                    "parent_wt_gene": wt_gene,
                    "mismatches": sum(a != b for a, b in zip(mut_15, wt_15)),
                    "best_class2_allele": best_entry["allele"],
                    "best_context_rank_pct": float(best_entry["context_rank"]),
                    "best_no_context_rank_pct": float(best_entry["no_context_rank"]),
                    "core_overlaps_recorded_target": best_entry["core_overlaps_recorded_target"],
                    "proposed_specificity_assessment": (
                        SPECIFICITY_INTERPRETATION
                    ),
                })
    c2_wt_controls.sort(key=lambda r: (r["best_context_rank_pct"], r["vaccine_gene"], r["mutant_15mer"]))
    table(safety_dir / "wetlab-class2-wt-controls.tsv", c2_wt_controls)

    # 4. Summarize Cryptic ORFs & dsRNA innate-immune audit (retaining complete per-query coverage + filtered binder table)
    cryptic_binders = []
    cryptic_cov_rows = []
    for pep, meta in sorted(cryptic_peps.items()):
        best_al, (best_nm, best_rk) = min(
            ((al, aff_cache[pep, al]) for al in alleles), key=lambda item: item[1][0]
        )
        best_rk_al, (_, min_rk_any) = min(
            ((al, aff_cache[pep, al]) for al in alleles), key=lambda item: item[1][1]
        )
        best_pr = float(prime_cache[pep]["%Rank_bestAllele"])
        passes_filter = min_rk_any <= 0.5 or best_rk <= 0.5 or best_nm <= 100.0 or (meta["strand"] == "sense" and best_nm <= 500.0)
        entry_c = {
            "peptide": pep,
            "strand": meta["strand"],
            "frame_vs_CDS": meta["frame_vs_CDS"],
            "contains_AUG_Met_in_ORF": meta["contains_AUG_Met"],
            "best_allele": best_al,
            "affinity_nm": round(best_nm, 2),
            "mhcflurry_rank_on_min_nm_allele": round(best_rk, 4),
            "min_mhcflurry_rank_any_allele": round(min_rk_any, 4),
            "min_mhcflurry_rank_allele": best_rk_al,
            "prime_rank": round(best_pr, 4),
        }
        cryptic_cov_rows.append({**entry_c, "retained_in_binder_subset": passes_filter})
        if passes_filter:
            cryptic_binders.append(entry_c)
    cryptic_binders.sort(key=lambda r: (r["strand"] != "sense", r["affinity_nm"]))
    strand_counts = cryptic_counts_by_stored_strand(cryptic_cov_rows)
    table(safety_dir / "cryptic-orf-peptides.tsv", cryptic_binders)
    table(safety_dir / "cryptic-orf-query-coverage.tsv.gz", cryptic_cov_rows)

    coding_duplex = longest_duplex_metrics(coding_rna)
    full_duplex = longest_duplex_metrics(rna_segment)
    sense_orfs = [o for o in orfs if o["strand"] == "sense"]
    save(safety_dir / "mrna-innate-and-cryptic-audit.json", {
        "coding_rna_duplex_and_motifs": coding_duplex,
        "known_transcript_segment_duplex_and_motifs": full_duplex,
        "mulroney_2024_frameshift_assessment": {
            "reference": "Mulroney et al., Nature 2024 (doi:10.1038/s41586-023-06800-3)",
            "coding_UUUU_count": coding_duplex["slippery_motif_counts"]["UUUU"],
            "coding_UUUC_count": coding_duplex["slippery_motif_counts"]["UUUC"],
            "coding_UUUA_count": coding_duplex["slippery_motif_counts"]["UUUA"],
            "coding_UUUG_count": coding_duplex["slippery_motif_counts"]["UUUG"],
            "note": MOTIF_INTERPRETATION,
        },
        "cryptic_orf_summary": {
            "by_stored_strand": strand_counts,
            "stored_strand_scope": CRYPTIC_STRAND_SCOPE,
            "total_orfs_ge_8aa": len(orfs),
            "sense_alt_frame_orfs_ge_8aa": len(sense_orfs),
            "sense_alt_frame_orfs_with_AUG": sum(o["contains_AUG_Met"] for o in sense_orfs),
            "longest_sense_alt_frame_orf_aa": max((o["length_aa"] for o in sense_orfs), default=0),
            "sense_cryptic_8to11mers_screened": sum(m["strand"] == "sense" for m in cryptic_peps.values()),
            "total_cryptic_8to11mers_screened": len(cryptic_peps),
            "sense_cryptic_peptides_min_nm_le_500": sum(
                r["strand"] == "sense" and r["affinity_nm"] <= 500.0 for r in cryptic_cov_rows
            ),
            "sense_cryptic_peptides_rank_le_0_5pct_on_min_nm_allele": sum(
                r["strand"] == "sense" and r["mhcflurry_rank_on_min_nm_allele"] <= 0.5 for r in cryptic_cov_rows
            ),
            "sense_cryptic_peptides_rank_le_0_5pct_any_patient_allele": sum(
                r["strand"] == "sense" and r["min_mhcflurry_rank_any_allele"] <= 0.5 for r in cryptic_cov_rows
            ),
            "note": cryptic_interpretation_for(target_root),
        },
        "innate_dsrna_interpretation": (
            RNA_INTERPRETATION
        ),
    })

    # 5. Author-proposed research questions, without clinical release criteria
    gates = research_questions(len(wetlab_list), len(c2_wt_controls), target_root=target_root)
    table(safety_dir / "research-questions.tsv", gates)

    offtarget_exact_hits = [
        r for r in homolog_rows
        if not r["is_mutated_gene_own_WT_only"] and r["mismatches"] == 0 and r["query_origin"] == "TARGET"
    ]
    if offtarget_exact_hits:
        raise ValueError("Unexpected exact off-target self match in core target: " + str(offtarget_exact_hits[0]))

    j_summary = load(target_root / "immunogenicity/junction-multi-model-summary.json")
    p_totals = load(target_root / "construct/protein/construct.json")["totals"]
    rel_prefix = "" if target_root == HERE else f"{target_root.name}/"
    summary = {
        "clinical_safety_established": False,
        "in_silico_safety_screen_status": "RECORDED_COMPUTATIONAL_SCREENS_ONLY",
        "why_clinical_safety_cannot_be_claimed_in_silico": SAFETY_LIMITS,
        "clinical_benefit_established": False,
        "priority_category_scope": "Heuristic sequence/model prioritization; not measured toxicity or clinical testing requirements.",
        "screening_policy": SCREENING_POLICY,
        "patient_biopsy_expression_scope": BIOPSY_EXPRESSION_SCOPE,
        "cryptic_orf_stored_strand_scope": CRYPTIC_STRAND_SCOPE,
        "in_silico_metrics": {
            "cryptic_peptides_by_stored_strand": strand_counts,
            "core_target_exact_matches_gencode_v46": 0,
            "core_target_exact_matches_gencode_v50": 0,
            "core_target_exact_matches_patient_normal_112389_contexts": 0,
            "routing_sequences_match_supplied_HLA_models": True,
            "eligible_class1_homolog_queries_screened": len(query_meta),
            "eligible_class2_9aa_core_queries_screened": len(c2_cores_meta),
            "sense_cryptic_8to11mers_screened": sum(m["strand"] == "sense" for m in cryptic_peps.values()),
            "junction_MHCflurry_strong_le_0_5pct": j_summary["MHCflurry_affinity_rank_le_0_5pct"],
            "junction_MixMHCpred_strong_le_0_5pct": j_summary["MixMHCpred_rank_le_0_5pct"],
            "junction_MixMHC2pred_strong_synthetic_9aa_core_le_0_5pct": j_summary["ClassII_synthetic_9aa_core_rank_le_0_5pct"],
            "junction_PRIME_low_rank_le_0_5pct": j_summary["PRIME_rank_le_0_5pct"],
            "junction_MHCflurry_le_500nm": p_totals["affinity_500"],
            "junction_min_affinity_nm": round(float(p_totals["min_affinity_nm"]), 2),
            "offtarget_1mm_same_hla_paralog_pairs": sum(
                r["model_priority_category"] == "TIER1_1MM_PARALOG_SAME_HLA"
                for r in homolog_rows
            ),
            "offtarget_2mm_strong_same_hla_pairs": sum(
                r["model_priority_category"] == "TIER2_2MM_STRONG_SAME_HLA"
                for r in homolog_rows
            ),
            "offtarget_titin_class_lead_mimic_pairs": sum(
                r["model_priority_category"] == "TIER2_TITIN_CLASS_LEAD_MIMIC"
                for r in homolog_rows
            ),
            "offtarget_all_same_hla_homolog_pairs_logged": sum(
                not r["is_mutated_gene_own_WT_only"] and r["same_hla_restriction"]
                for r in homolog_rows
            ),
            "own_wt_same_hla_class1_controls_omitted": 0,
            "wetlab_negative_control_peptides_to_synthesize": len(wetlab_list),
            "wetlab_class2_wt_15mers_to_synthesize": len(c2_wt_controls),
            "known_RNA_mfe_max_uninterrupted_stem_bp": full_duplex["max_uninterrupted_mfe_stem_bp"],
            "known_RNA_mfe_max_bulge_tolerant_stem_bp": full_duplex["max_bulge_tolerant_mfe_stem_bp"],
            "coding_UUUU_and_UUUC_slippery_motifs": (
                coding_duplex["slippery_motif_counts"]["UUUU"] + coding_duplex["slippery_motif_counts"]["UUUC"]
            ),
            "sense_cryptic_orf_rank_le_0_5pct_on_min_nm_allele": sum(
                r["strand"] == "sense" and r["mhcflurry_rank_on_min_nm_allele"] <= 0.5 for r in cryptic_cov_rows
            ),
            "sense_cryptic_orf_rank_le_0_5pct_any_patient_allele": sum(
                r["strand"] == "sense" and r["min_mhcflurry_rank_any_allele"] <= 0.5 for r in cryptic_cov_rows
            ),
        },
        "control_panel_scope": {
            "stage1_class1_plate_file": f"{rel_prefix}safety/wetlab-negative-control-peptides.tsv",
            "stage1_class1_peptide_count": len(wetlab_list),
            "stage1_class2_wt_file": f"{rel_prefix}safety/wetlab-class2-wt-controls.tsv",
            "stage1_class2_wt_15mer_count": len(c2_wt_controls),
            "exhaustive_for_shared_hla_own_wt_class1": True,
            "exhaustive_for_1mm_same_hla_paralogs": True,
            "prioritized_subset_for_2to4mm_offtarget_homologs": True,
            "stage2_deconvolution_source": f"{rel_prefix}safety/offtarget-homologs.tsv.gz",
            "scope_warning": "Panels cover the stated model and sequence rules, not all possible cross-reactive targets. Assay design and acceptance criteria are unvalidated; no negative biological result or absence of harm is established.",
        },
        "unresolved_research_questions": [g["question_id"] for g in gates if g["status"] != "COMPUTATIONAL_RESULT"],
    }
    save(safety_dir / "safety-summary.json", summary)

    outputs = [
        "patient-normal-unresolved-genes.tsv",
        "hla-routing-self-audit.json",
        "patient-normal-audit.json",
        "homolog-query-coverage.tsv",
        "offtarget-homologs.tsv.gz",
        "class2-core-homologs.tsv",
        "wetlab-negative-control-peptides.tsv",
        "wetlab-class2-wt-controls.tsv",
        "cryptic-orf-query-coverage.tsv.gz",
        "cryptic-orf-peptides.tsv",
        "mrna-innate-and-cryptic-audit.json",
        "research-questions.tsv",
        "safety-summary.json",
    ]
    save(safety_dir / "receipt.json", {
        "output_sha256": {name: sha(safety_dir / name) for name in outputs},
        "clinical_safety_established": False,
    })
    refresh_cached_annotations(target_root)
    print(f"Computational screens generated ({target_root.name}):", summary["in_silico_metrics"], flush=True)
    return {
        "proteome_records": proteome_records,
        "normal_gene_map": normal_gene_map,
        "unresolved": unresolved,
        "callability": callability,
        "chr6_segments": chr6_segments,
        "typed_hla": typed_hla,
        "aff_cache": aff_cache,
        "mix_cache": mix_cache,
        "prime_cache": prime_cache,
        "ap": ap,
    }


def generate(models_dir, mix_exec, prime_exec, r2_normal_dir, target="both"):
    shared = {}
    if target in ("main", "both"):
        shared = generate_for_target(HERE, models_dir, mix_exec, prime_exec, r2_normal_dir)
    if target in ("strand-strict-15", "both"):
        generate_for_target(HERE / "strand-strict-15", models_dir, mix_exec, prime_exec, r2_normal_dir, **shared)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["generate", "verify"])
    parser.add_argument("--models", type=Path)
    parser.add_argument("--mix-exec", type=Path)
    parser.add_argument("--prime-exec", type=Path)
    parser.add_argument("--r2-normal-dir", type=Path)
    parser.add_argument("--target", choices=["main", "strand-strict-15", "both"], default="both")
    args = parser.parse_args()
    if args.action == "verify":
        print(verify_offline(HERE if args.target in ("main", "both") else HERE / "strand-strict-15"))
    else:
        generate(args.models, args.mix_exec, args.prime_exec, args.r2_normal_dir, target=args.target)
