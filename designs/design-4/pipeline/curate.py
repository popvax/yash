"""Exploratory antigen curation from frozen local evidence, not clinical release criteria."""
import argparse
import hashlib
import math
from collections import Counter
from pathlib import Path

from dataio import load, read_text_auto, rows, save, table
from hydrophobicity import summarize

HERE = Path(__file__).resolve().parents[1]
LENGTHS = (8, 9, 10, 11)


def overlaps_target(start, end, lo, hi):
    return start < lo < end if lo == hi else start < hi and end > lo


def peptide_windows(sequence, lo, hi, lengths=LENGTHS):
    return [
        (start, sequence[start:start + length])
        for length in lengths for start in range(len(sequence) - length + 1)
        if overlaps_target(start, start + length, lo, hi)
    ]


def class2_best(sequence, lo, hi, predictions):
    candidates = [
        predictions[peptide]
        for _, peptide in peptide_windows(sequence, lo, hi, (15,))
        if peptide in predictions
    ]
    return min(candidates, key=lambda row: (float(row["%Rank_best"]), row["Peptide"])) if candidates else None


def cta_gate(row):
    values = {key: float(row[key]) for key in (
        "codons_with_ge2_fragments", "mtec_mean", "ag2_vital_max",
        "ag2_hepatocyte_sc", "ag2_testis", "ag2_testis_over_vital_max",
        "prime_best_rank", "mhc2_best_rank", "full_25aa_merged_fragments",
    )}
    return {
        "two_fragments_per_codon": values["codons_with_ge2_fragments"] == 25,
        "low_recorded_mtec_mean": values["mtec_mean"] <= 0.05,
        "low_bulk_vital_model_score": values["ag2_vital_max"] <= 0.02,
        "low_hepatocyte_model_raw_score": values["ag2_hepatocyte_sc"] <= -1.0,
        "testis_model_score": values["ag2_testis"] >= 0.30,
        "testis_vital_score_ratio": values["ag2_testis_over_vital_max"] >= 25.0,
        "prime_rank": values["prime_best_rank"] <= 0.05,
        "class2_or_full_window_support": (
            values["mhc2_best_rank"] <= 1.5 or values["full_25aa_merged_fragments"] >= 20
        ),
    }


def derive():
    pool = rows(HERE / "exports/candidate-pool.tsv")
    baseline = rows(HERE / "exports/design2-selected.tsv")
    ranked = {row["variant_id"]: row for row in rows(HERE / "exports/ranked-mutations.tsv")}
    proteins = load(HERE / "exports/rna-proteins.json")
    dna = {row["variant_id"]: row for row in load(HERE / "exports/variant-evidence.json")}
    rna = {row["variant_id"]: row for row in load(HERE / "exports/rna-evidence.json")}
    strand_records = load(HERE / "exports/rna-strand-audit.json")["mutation_strand_records"]
    reference = load(HERE / "references/peptide-matches.json")
    patient_normal = load(HERE / "references/patient-normal-matches.json")
    native_c1 = rows(HERE / "immunogenicity/native-class1-presentation.tsv")
    native_c2 = {row["Peptide"]: row for row in rows(HERE / "immunogenicity/native-class2-context.tsv")}
    binding_c2 = {row["Peptide"]: row for row in rows(HERE / "immunogenicity/native-class2-no-context.tsv")}
    expression = {}
    for line in read_text_auto(HERE / "exports/expression.tsv").splitlines():
        if line and not line.startswith("#"):
            gene, value = line.split("\t")[:2]
            expression[gene] = float(value)

    qualifying, mutation_audit = {}, []
    for protein in proteins:
        event = protein["variant_id"]
        seq = protein["sequence"]
        lo, hi = protein["mutation_start_0based"], protein["mutation_end_0based_exclusive"]
        allowed_peptides = {
            peptide for _, peptide in peptide_windows(seq, lo, hi)
            if peptide in reference and not reference[peptide]["occurs"]
            and not patient_normal[peptide]["occurs"]
        }
        c1_candidates = [
            row for row in native_c1
            if row["sequence_name"] == protein["gene"] and row["peptide"] in allowed_peptides
            and overlaps_target(int(row["pos"]), int(row["pos"]) + len(row["peptide"]), lo, hi)
        ]
        best_c1 = min((float(row["affinity"]) for row in c1_candidates), default=math.inf)
        c2 = class2_best(seq, lo, hi, native_c2)
        best_c2 = float(c2["%Rank_best"]) if c2 else math.inf
        strand_info = strand_records[event]
        coding_alt = int(strand_info["coding_strand_alt_fragments"])
        opposite_alt = int(strand_info["opposite_strand_alt_fragments"])
        gates = {
            "paired_normal_supported": (
                dna[event]["primary_DNA_eligible"]
                and dna[event]["somatic_status"] == "SUPPORTED_PAIRED_SOMATIC"
            ),
            "unstranded_isovar_context": (
                protein["usable_mutant_context"] and protein["supporting_fragments"] >= 2
            ),
            "model_nomination": best_c1 <= 500.0 or best_c2 <= 2.0,
        }
        accepted = all(gates.values())
        mutation_audit.append({
            "event_id": event, "gene": protein["gene"],
            "native_best_target_class1_nm": best_c1 if math.isfinite(best_c1) else None,
            "native_best_target_class2_context_rank": best_c2 if c2 else None,
            "RNA_ALT_fragments": rna[event]["rna_alt_fragments"],
            "coding_strand_RNA_ALT_fragments": coding_alt,
            "opposite_strand_RNA_ALT_fragments": opposite_alt,
            "rna_strand_tier": strand_info["rna_strand_tier"],
            "RNA_context_fragments": protein["supporting_fragments"],
            "gates": gates,
            "coding_strand_supported": coding_alt >= 1,
            "selected": accepted,
        })
        if accepted:
            qualifying[event] = protein

    order = [row["event_id"] for row in baseline if row["source"] == "MUTATION"]
    order += [
        "mn_chr6_27146459_C_T", "mn_chr9_69246779_T_C",
        "mn_chr14_94383173_G_T", "mn_chr15_90895499_C_T",
    ]
    if set(order) != set(qualifying) or len(order) != len(set(order)):
        raise ValueError("Frozen evidence does not produce the declared exploratory mutation panel")
    baseline_by_id = {row["event_id"]: row for row in baseline if row["source"] == "MUTATION"}
    selected = []
    for event in order:
        protein, evidence = qualifying[event], dna[event]
        row = dict(baseline_by_id.get(event, ranked.get(event, {})))
        full = protein["sequence"]
        if protein["gene"] == "KCTD19":
            start, end = 8, 31
        elif protein["gene"] == "SERPINA1":
            start, end = 16, 41
        elif protein["gene"] == "FES":
            start, end = 0, len(full)
        else:
            start, end = int(row["window_start_0based"]), int(row["window_end_exclusive"])
        sequence = full[start:end]
        lo = protein["mutation_start_0based"] - start
        hi = min(len(sequence), protein["mutation_end_0based_exclusive"] - start)
        if not (0 <= lo < hi <= len(sequence)) or "*" in sequence:
            raise ValueError("Mutant window or translated context is inconsistent")
        windows = peptide_windows(sequence, lo, hi)
        if any(peptide not in reference or peptide not in patient_normal for _, peptide in windows):
            raise ValueError("Missing frozen exact-reference or patient-normal check for a mutant peptide")
        targets = sorted({
            peptide for _, peptide in windows
            if not reference[peptide]["occurs"] and not patient_normal[peptide]["occurs"]
        })
        c1_rows = [
            value for value in native_c1
            if value["sequence_name"] == protein["gene"]
            and value["peptide"] in targets and start <= int(value["pos"])
            and int(value["pos"]) + len(value["peptide"]) <= end
        ]
        best = min(c1_rows, key=lambda value: (float(value["affinity"]), value["peptide"]))
        gene = "TJP2" if protein["gene"] == "AL358113.1" else protein["gene"]
        change = "p.S886P" if gene == "TJP2" else protein["DNA_annotation_protein_change"]
        c2 = class2_best(sequence, lo, hi, native_c2)
        b2 = class2_best(sequence, lo, hi, binding_c2)
        normal = evidence["normal_evidence"]["read_evidence"]
        strand_info = strand_records[event]
        coding_alt = int(strand_info["coding_strand_alt_fragments"])
        opposite_alt = int(strand_info["opposite_strand_alt_fragments"])
        if gene == "SERPINA1":
            mut_tier = "signal-peptide-validation-hold"
            mut_basis = (
                "Paired-normal DNA + 470 coding-strand mutant RNA fragments + native-context Class II rank <=2%; "
                "held in reserve because S22Y lies in cleaved native signal peptide 1-24 and Class I exceeds 500 nM"
            )
        elif coding_alt == 0:
            mut_tier = "somatic-dna-only-opposite-strand-rna-watch"
            mut_basis = (
                "Somatic DNA-supported (paired-normal verified, AF 0.318; clonality unestablished) + Class I <=500 nM; "
                "0 coding-strand and 2 opposite-strand RNA ALT fragments (unstranded Isovar artifact; 0/4 coding-strand "
                "depth does not prove absence of expression; coding-RNA status is unconfirmed. SG-04 proposes strand-specific RT-ddPCR as one confirmation approach)"
            )
        elif coding_alt == 1:
            mut_tier = "somatic-single-coding-fragment-low-depth"
            mut_basis = (
                "Paired-normal DNA + 1 coding-strand and 1 opposite-strand RNA ALT fragment (pooled to 2 in "
                "unstranded Isovar) + spliced host-gene transcription + Class I <=500 nM"
            )
        else:
            mut_tier = "somatic-research"
            mut_basis = (
                "Paired-normal DNA + >=2 coding-strand mutant RNA fragments + translated RNA context + Class I <=500 nM"
            )
        row.update({
            "source": "MUTATION", "gene": gene, "event_id": event, "variant_id": event,
            "sequence": sequence, "protein_change": change, "protein_id": protein["protein_id"],
            "DNA_annotation_protein_change": protein["DNA_annotation_protein_change"],
            "DNA_annotation_transcript_id": protein["DNA_annotation_transcript_id"],
            "original_RNA_gene_label": protein["gene"],
            "RNA_supported_transcript_ids": protein["transcript_ids"],
            "gene_assignment_unambiguous": protein["gene_assignment_unambiguous"],
            "annotation_resolution": (
                "Display canonical TJP2; retain original readthrough gene ambiguity and DNA isoform annotation."
                if gene == "TJP2" else "Source annotation retained"
            ),
            "region_id": "REG_" + hashlib.sha256(sequence.encode()).hexdigest()[:20],
            "window_start_0based": start, "window_end_exclusive": end,
            "target_start_in_window": lo, "target_end_in_window": hi,
            "target_mode": "frameshift_tail" if protein["frameshift"] else "mutant_residue",
            "target_peptides": ";".join(targets), "distinct_target_peptides": len(targets),
            "exact_self_target_occurrences": sum(reference[p]["occurs"] for _, p in windows),
            "exact_patient_normal_target_occurrences": sum(patient_normal[p]["occurs"] for _, p in windows),
            "nonself_target_occurrences": sum(
                not reference[p]["occurs"] and not patient_normal[p]["occurs"] for _, p in windows
            ),
            "normalized_allele": evidence["normalized_allele"],
            "primary_DNA_eligible": evidence["primary_DNA_eligible"],
            "DNA_evidence_tier": evidence["evidence_tier"], "somatic_status": evidence["somatic_status"],
            "DNA_AD_ref": evidence["DNA_AD"][0], "DNA_AD_alt": evidence["DNA_AD"][1],
            "DNA_model_AF": evidence["DNA_model_AF"],
            "normal_REF_fragments": normal["fragment_counts"].get("REF", 0),
            "normal_ALT_fragments": evidence["normal_evidence"]["target_alternate_normal_fragments"],
            "gene_coding_strand": strand_info["gene_coding_strand"],
            "RNA_ALT_fragments": rna[event]["rna_alt_fragments"],
            "coding_strand_RNA_ALT_fragments": coding_alt,
            "opposite_strand_RNA_ALT_fragments": opposite_alt,
            "RNA_total_fragments": rna[event]["rna_total_fragments"],
            "coding_strand_RNA_REF_fragments": int(strand_info["coding_strand_ref_fragments"]),
            "opposite_strand_RNA_REF_fragments": int(strand_info["opposite_strand_ref_fragments"]),
            "coding_strand_spliced_reads_pm20kb": int(strand_info["coding_strand_spliced_reads_pm20kb"]),
            "rna_strand_tier": strand_info["rna_strand_tier"],
            "RNA_context_fragments": protein["supporting_fragments"],
            "rna_alternate_fraction": rna[event]["rna_alt_fraction"],
            "supplied_gene_value": expression.get(gene),
            "best_target_peptide": best["peptide"], "best_target_allele": best["best_allele"],
            "best_target_affinity_nm": float(best["affinity"]),
            "passes_500nm_gate": float(best["affinity"]) <= 500.0,
            "native_class2_context_peptide": c2["Peptide"] if c2 else "",
            "native_class2_context_rank": float(c2["%Rank_best"]) if c2 else None,
            "class2_no_context_peptide": b2["Peptide"] if b2 else "",
            "class2_no_context_rank": float(b2["%Rank_best"]) if b2 else None,
            "hydrophobicity": summarize(sequence),
            "primary_window_eligible": True,
            "primary_reference_exclusion_applied": True,
            "patient_self_exclusion_applied": True,
            "self_screen_scope": "Exact GENCODE v46 + v50 + patient matched-normal WBC coding/typed-HLA search (112,389 personal contexts; excludes 8,987 incomplete-CDS transcript models and unmeasured noncanonical translation)",
            "design_tier": mut_tier,
            "research_only": True,
            "clinical_release_ready": False,
            "eligibility_basis": mut_basis,
        })
        if gene == "KCTD19":
            original = baseline_by_id[event]
            old_targets = {p for _, p in peptide_windows(
                original["sequence"], int(original["target_start_in_window"]),
                int(original["target_end_in_window"]), (8, 9, 10, 11, 15),
            )}
            new_targets = {p for _, p in peptide_windows(sequence, lo, hi, (8, 9, 10, 11, 15))}
            if old_targets != new_targets:
                raise ValueError("KCTD19 trimming lost a target-overlapping peptide")
        selected.append(row)

    splice_audit = []
    for row in sorted((r for r in pool if r["source"] == "SPLICE"),
                      key=lambda r: (-float(r["weighted_score_left"] or 0), r["event_id"])):
        support = int(row["fragment_count"] or 0)
        agree = int(row["strand_agreeing_fragments"] or 0)
        gates = {
            "no_recorded_exact_self_target": int(row["exact_self_target_occurrences"] or -1) == 0,
            "single_gene_locus": row["event_id"].count("ENSG") == 1,
            "full_window_fragments": support >= 10,
            "strand_agreement": support > 0 and agree / support >= 0.70,
            "complete_rna_window": row["patient_RNA_sequence_confirmed"] == "True",
        }
        accepted = all(gates.values())
        splice_audit.append({
            "gene": row["gene"], "event_id": row["event_id"], "gates": gates,
            "selected": accepted, "fragment_count": support,
            "strand_agreement_fraction": agree / support if support else None,
        })
        if accepted:
            selected.append({
                **row, "hydrophobicity": summarize(row["sequence"]),
                "design_tier": "splice-validation-hold", "research_only": True,
                "clinical_release_ready": False,
                "eligibility_basis": "Exploratory splice evidence filter; noncanonical motif still requires orthogonal validation",
            })

    cta_audit = []
    for audit in sorted(rows(HERE / "exports/cta-audit.tsv"), key=lambda r: int(r["pool_rank"])):
        gates = cta_gate(audit)
        accepted = all(gates.values())
        cta_audit.append({**audit, "gates": gates, "selected": accepted})
        if accepted:
            matches = [r for r in pool if r["source"] == "CTA/SELF" and r["region_id"] == audit["region_id"]]
            if len(matches) != 1 or matches[0]["sequence"] != audit["sequence"]:
                raise ValueError("CTA aggregate does not identify a unique source sequence")
            selected.append({
                **matches[0], "hydrophobicity": summarize(audit["sequence"]),
                "design_tier": "intentional-self-validation-hold", "research_only": True,
                "clinical_release_ready": False,
                "eligibility_basis": "Retrospective exploratory eight-criterion filter, NOT a validated safety gate",
                "RNA_full_window_fragments": int(audit["full_25aa_merged_fragments"]),
            })

    request = load(HERE / "exports/inframe-request.json")
    inframe = load(HERE / "exports/inframe-selection.json")
    for source in (r for r in baseline if r["source"] == "FUSION"):
        row = dict(source)
        if source["gene"] == "CDH1::NRG1":
            lo = inframe["junction_codon_0based"] - inframe["source_protein_window_0based"][0]
            fusion_targets = sorted(
                p for p in set(request["peptides"]) - set(request["reference_matching_peptides"])
                if not patient_normal[p]["occurs"]
            )
            row.update(
                target_start_in_window=lo, target_end_in_window=lo + 1,
                target_mode="recorded_junction_codon",
                target_peptides=";".join(fusion_targets),
                exact_patient_normal_target_occurrences=sum(
                    patient_normal[p]["occurs"] for _, p in peptide_windows(source["sequence"], lo, lo + 1)
                ),
                patient_self_exclusion_applied=True,
                design_tier="inframe-fusion-research",
                RNA_junction_reads=inframe["original_junction_read_count"],
                RNA_spanning_fragments=inframe["original_spanning_fragment_count"],
            )
        else:
            row["design_tier"] = "alternative-fusion-hypothesis-hold"
        row.update(research_only=True, clinical_release_ready=False)
        selected.append(row)

    for index, row in enumerate(selected, 1):
        row["antigen_id"] = f"AG{index:02d}"
        row.setdefault("protein_change", "")
        row.setdefault("context_id", "")
    counts = Counter(row["source"] for row in selected)
    if counts != {"MUTATION": 16, "SPLICE": 1, "CTA/SELF": 5, "FUSION": 2}:
        raise ValueError("Exploratory panel counts changed: " + str(counts))
    summary = {
        "selected_regions": len(selected), "source_counts": dict(counts),
        "antigen_residues_aa": sum(len(row["sequence"]) for row in selected),
        "mutation_coding_strand_RNA_ALT_fragment_count_sum": sum(
            int(row["coding_strand_RNA_ALT_fragments"]) for row in selected if row["source"] == "MUTATION"
        ),
        "mutation_opposite_strand_RNA_ALT_fragment_count_sum": sum(
            int(row["opposite_strand_RNA_ALT_fragments"]) for row in selected if row["source"] == "MUTATION"
        ),
        "mutation_RNA_ALT_fragment_count_sum": sum(
            int(row["RNA_ALT_fragments"]) for row in selected if row["source"] == "MUTATION"
        ),
        "mutation_context_fragment_count_sum": sum(
            int(row["RNA_context_fragments"]) for row in selected if row["source"] == "MUTATION"
        ),
        "cta_audit_scope": "40 tabulated candidates: 36 codon-supported candidates among ranks 1-84 plus four rejected baseline candidates; not all possible antigens.",
        "cta_score_ratio_formula": "ag2_testis / (ag2_vital_max + 0.01); arbitrary pseudocount, not a tissue-expression or safety ratio",
        "clinical_release_ready": False,
    }
    return selected, mutation_audit, splice_audit, cta_audit, summary


def core_regions(selected):
    return [
        row for row in selected
        if (row["source"] == "MUTATION" and row["passes_500nm_gate"])
        or (row["source"] == "FUSION" and row["gene"] == "CDH1::NRG1")
    ]


def strand_strict_core_regions(selected):
    """Contingency core requiring >=1 coding-strand mutant RNA fragment (excises KCTD19)."""
    return [
        row for row in core_regions(selected)
        if row["source"] != "MUTATION" or int(row.get("coding_strand_RNA_ALT_fragments", 0)) >= 1
    ]


def write(output, panel="core"):
    selected, mutations, splices, ctas, summary = derive()
    output = Path(output)
    core = core_regions(selected)
    if len(core) != 16 or len(selected) - len(core) != 8:
        raise ValueError("Core/reserve partition changed")
    if panel == "core":
        active = core
    elif panel == "strand-strict-core":
        active = strand_strict_core_regions(selected)
        if len(active) != 15:
            raise ValueError("Strand-strict core expected 15 regions")
    else:
        active = selected
    active_ids = {row["region_id"] for row in active}
    reserve = [row for row in selected if row["region_id"] not in active_ids]
    table(output / "selection/selected.tsv", active)
    save(output / "selection/selected.json", active)
    (output / "selection/selected.fasta").write_text("".join(
        f">{r['antigen_id']}|{r['source']}|{r['gene']}\n{r['sequence']}\n" for r in active
    ))
    table(output / "selection/reserve.tsv", reserve)
    save(output / "selection/reserve.json", reserve)
    table(output / "selection/exploratory-selected.tsv", selected)
    save(output / "selection/exploratory-selected.json", selected)
    # "selected" is the exploratory-gate outcome (paired-normal + RNA context + model nomination),
    # not an encoding decision: a gate-passing mutation can still be held in reserve. "encoded"
    # records whether this event's window is in the panel written to selected.tsv for this run.
    encoded_events = {row["event_id"] for row in active if row["source"] == "MUTATION"}
    mutations = [dict(row, encoded=row["event_id"] in encoded_events) for row in mutations]
    table(output / "selection/mutation-audit.tsv", mutations)
    table(output / "selection/splice-audit.tsv", splices)
    table(output / "selection/cta-audit.tsv", ctas)
    save(output / "selection/exploratory-summary.json", summary)
    active_muts = [row for row in active if row["source"] == "MUTATION"]
    fragment_note = (
        f"In this {len(active_muts)}-mutation selection, "
        f"{sum(int(r['RNA_ALT_fragments']) for r in active_muts)} pooled ALT fragments decompose into "
        f"{sum(int(r['coding_strand_RNA_ALT_fragments']) for r in active_muts)} coding-strand + "
        f"{sum(int(r['opposite_strand_RNA_ALT_fragments']) for r in active_muts)} opposite-strand fragments. "
        "Pooled Isovar counts require strand interpretation in this inferred fr-firststrand library. "
        "Counts are not UMI molecules, tumour-cell coverage, clonality or vaccine efficacy."
    )

    current = {
        **summary, "fragment_sum_interpretation": fragment_note, "panel": panel, "selected_regions": len(active),
        "source_counts": dict(Counter(row["source"] for row in active)),
        "antigen_residues_aa": sum(len(row["sequence"]) for row in active),
        "mutation_coding_strand_RNA_ALT_fragment_count_sum": sum(
            int(row["coding_strand_RNA_ALT_fragments"]) for row in active_muts
        ),
        "mutation_opposite_strand_RNA_ALT_fragment_count_sum": sum(
            int(row["opposite_strand_RNA_ALT_fragments"]) for row in active_muts
        ),
        "mutation_RNA_ALT_fragment_count_sum": sum(
            int(row["RNA_ALT_fragments"]) for row in active_muts
        ),
        "mutation_context_fragment_count_sum": sum(
            int(row["RNA_context_fragments"]) for row in active_muts
        ),
        "core_mutations_by_coding_strand_tier": {
            "coding_strand_ge2_alt": sum(int(r["coding_strand_RNA_ALT_fragments"]) >= 2 for r in active_muts),
            "single_coding_strand_alt": sum(int(r["coding_strand_RNA_ALT_fragments"]) == 1 for r in active_muts),
            "zero_coding_strand_alt_dna_only_watch": sum(int(r["coding_strand_RNA_ALT_fragments"]) == 0 for r in active_muts),
        },
        "reserve_regions": len(reserve),
        "core_policy": (
            "Main research selection: 12 paired-normal somatic mutation candidates with >=2 coding-strand mutant "
            "RNA fragments, FES and SYT12 with one each, and KCTD19 with somatic DNA support but unconfirmed "
            "coding-RNA expression and clonality. KCTD19 is retained as an unconfirmed research slot; "
            "shallow 0/4 coding-strand coverage does not establish biological silence. The in-frame CDH1::NRG1 "
            "fusion is an explicit inclusion exception. No RNA-only splice or intentional-self/CTA region is encoded. "
            "A separately optimized 15-region protein alternative excludes KCTD19."
            if panel == "core" else
            "15-region research alternative: excludes KCTD19, retains 14 mutation candidates "
            "with at least one coding-strand mutant RNA fragment plus the in-frame CDH1::NRG1 fusion, and "
            "reoptimizes the remaining order, linkers and known-segment PopRNA/GenBank nucleotide encoding."
            if panel == "strand-strict-core" else
            "Exploratory 24-candidate research pool; includes reserves that are not encoded in the main protein."
        ),
        "SERPINA1_policy": (
            "Retain as a reserve: S22Y is in the normally cleaved signal peptide and "
            "Class I exceeds 500 nM. Conditional Class II scores do not establish "
            "presentation of this signal-derived tumor antigen."
        ),
        "KCTD19_reassessment": (
            "KCTD19 p.G920R is somatic DNA-supported (134/419 tumour DNA reads, AF 0.318; 49/0 WBC normal REF; "
            "clonality unestablished without purity, copy-number and mutant-copy-multiplicity deconvolution) and predicts "
            "Class I <=500 nM (RSVSYSILR 36.9 nM HLA-A*33:03), but has 0 coding-strand ('-') and 2 opposite-strand ('+') "
            "mutant RNA fragments (0 '-' spliced reads in +-20 kb; AlphaGenome 2 tumour RNA score 0.0024). It entered "
            "frozen Isovar solely via unstranded read pooling, whereas 7 other paired-normal somatic DNA missense "
            "mutations (AF ~0.31-0.38; clonality unestablished) with 0 coding-strand ALT reads at 1-4x local depth were "
            "excluded. Shallow 4-fragment coding-strand REF depth does not prove KCTD19 is unexpressed across all "
            "lesions, so KCTD19 is downgraded to somatic DNA-only / unconfirmed coding-RNA watch with an author-proposed SG-04 "
            "strand-specific RT-ddPCR and a pre-solved 15-region KCTD19-excised protein alternative with complete cached Class-II context coverage in strand-strict-15/."
        ),
    }
    save(output / "selection/summary.json", current)
    return current


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE)
    parser.add_argument("--panel", choices=["core", "strand-strict-core", "exploratory"], default="core")
    args = parser.parse_args()
    print(write(args.output, args.panel))
