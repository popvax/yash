"""Generate bounded public evidence summaries from the supplied Design 4 artifacts."""
from pathlib import Path

from curate import HERE, peptide_windows
from dataio import load, rows, save, sha, table
from research_qualifications import refresh_cached_annotations


NOTES = {
    "KLF5": "Normal KLF-family homologs warrant specificity testing; cross-recognition is unmeasured.",
    "SF1": "The best recorded Class-II core is unchanged; mutation-specific CD4 recognition is unproven.",
    "AQP3": "Two coding-strand mutant fragments provide limited RNA sampling.",
    "RXRA": "Wild-type and RXR-family homolog comparisons remain necessary; clinical isoform labels differ.",
    "HSPG2": "DNA/RNA support does not establish malignant-cell presentation.",
    "TGFBR2": "Frameshift numbering differs by isoform; clonality is unestablished.",
    "PPIC": "Class-II support is weak and the best recorded core is unchanged.",
    "KCTD19": "Somatic DNA-supported; clonality and coding-RNA expression remain unestablished. The separate protein alternative excludes it.",
    "VPS11": "Predicted presentation and mutation-specific T-cell recognition remain untested.",
    "ARHGAP40": "Model affinity does not establish presentation or recognition.",
    "PLEC": "Mutant-versus-wild-type recognition remains untested.",
    "SYT12": "One coding-strand and one opposite-strand mutant fragment; RNA support is shallow.",
    "H2BC12": "Supplied gene-level expression value is missing; histone homologs warrant specificity testing.",
    "TJP2": "Readthrough/isoform annotation remains ambiguous; the observed peptide window is the stable identity.",
    "FES": "Only the observed 15-aa window is encoded; one coding-strand mutant fragment. Regional spliced reads are not mutant support.",
    "CDH1::NRG1": "Approximately 509 nM Class-I prediction is an explicit inclusion exception. The source records 200 junction reads and one spanning fragment in distinct, non-UMI units; peptide presentation is untested.",
    "SERPINA1": "S22Y lies in the annotated signal peptide; Class-I nomination exceeds 500 nM. Signal-derived presentation remains unproven.",
    "POLE": "RNA-only splice hypothesis with a noncanonical CT-GC motif; translation and tumour specificity remain unconfirmed.",
    "PRAME": "Unmutated self/CTA candidate; no full-window merged fragment in the recorded audit. Tumour specificity and normal-tissue safety are unestablished.",
    "NUTM1": "Unmutated self/CTA candidate; full-window RNA support does not establish malignant-cell origin.",
    "MAGEA6": "Unmutated self/CTA candidate; homologous mapping and T-cell specificity require validation.",
    "TSPYL6": "Unmutated self/CTA candidate; normal-tissue presentation and useful T-cell recognition are unestablished.",
    "MAGEA1": "Unmutated self/CTA candidate; normal-tissue presentation and useful T-cell recognition are unestablished.",
    "CDH1--NRG1": "Alternative out-of-frame fusion hypothesis, not a second independently established driver.",
}


def markdown(headers, records):
    def safe(value):
        return str(value if value is not None else "n/a").replace("|", "/").replace("\n", " ")
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join("---" for _ in headers) + " |"]
    lines.extend("| " + " | ".join(safe(row.get(key)) for key in headers) + " |" for row in records)
    return "\n".join(lines) + "\n"


def results_summary_snapshot(target_root=HERE):
    """Read the deterministic source snapshots embedded in results-summary.json.

    This function does not refresh annotations, infer models or write artifacts.
    Verification uses the same source mapping as generation so hashing a stale
    review file cannot substitute for checking its dependencies.
    """
    target_root = Path(target_root)

    def mrna_snapshot(directory):
        encoding = load(directory / "construct.json")
        audit = load(directory / "optimizer-audit.json")
        verification = load(directory / "verification.json")
        return {
            "metrics": encoding["optimization_metrics"], "lengths": encoding["lengths"],
            "known_segment_annotation": verification["genbank_annotation"],
            "original_alt_frame_stops_preserved": len(audit["alternative_frame_stop_intervals"]),
            "accessibility": audit["initiation_accessibility"],
            "structure_rescue": audit["bounded_structure_rescue"],
        }

    return {
        "core": load(target_root / "selection/summary.json"),
        "protein": load(target_root / "construct/protein/construct.json"),
        "junction_other_models": load(target_root / "immunogenicity/junction-multi-model-summary.json"),
        "safety": load(target_root / "safety/safety-summary.json"),
        "strand_strict_15_safety": load(target_root / "strand-strict-15/safety/safety-summary.json"),
        "alphagenome2_additional": load(target_root / "alphagenome2/additional-evidence-receipt.json"),
        "mrna": mrna_snapshot(target_root / "construct/mrna"),
        "strand_strict_15_mrna": mrna_snapshot(target_root / "strand-strict-15/construct/mrna"),
        "clinical_release_ready": False,
        "clinical_efficacy_or_superiority_established": False,
        "interpretation": "Recorded prediction and sequence checks do not establish clinical performance. The computational screens and unresolved research questions do not supply clinical release criteria.",
    }


def main():
    refresh_cached_annotations(HERE)
    if (HERE / "strand-strict-15/safety/receipt.json").is_file():
        refresh_cached_annotations(HERE / "strand-strict-15")
    selected = load(HERE / "selection/selected.json")
    reserve = load(HERE / "selection/reserve.json")
    models = {r["region_id"]: r for r in rows(HERE / "immunogenicity/antigen-model-summary.tsv")}
    cta = {r["gene"]: r for r in rows(HERE / "exports/cta-audit.tsv")}
    core_ids = {r["region_id"] for r in selected}
    audit, core_table, reserve_table = [], [], []
    for region in selected + reserve:
        core = region["region_id"] in core_ids
        model = models.get(region["region_id"], {})
        full = {
            "panel": "CORE_RESEARCH_NOT_CLINICAL" if core else "RESERVE_NOT_ENCODED",
            "antigen_id": region["antigen_id"], "gene": region["gene"], "source": region["source"],
            "region_id": region["region_id"], "event_id": region["event_id"],
            "protein_change": region.get("protein_change", ""),
            "sequence_length_aa": len(region["sequence"]),
            **{key: region.get(key) for key in (
                "DNA_AD_ref", "DNA_AD_alt", "DNA_model_AF", "normal_REF_fragments", "normal_ALT_fragments",
                "RNA_ALT_fragments", "coding_strand_RNA_ALT_fragments", "opposite_strand_RNA_ALT_fragments",
                "rna_strand_tier", "RNA_total_fragments", "RNA_context_fragments",
            )},
            "RNA_full_window_fragments": region.get("RNA_full_window_fragments", cta.get(region["gene"], {}).get("full_25aa_merged_fragments")),
            "ClassI_lead": model.get("class1_lead_peptide"),
            "ClassI_lead_allele": model.get("class1_lead_allele"),
            "ClassI_lead_model_affinity_nm": model.get("class1_lead_affinity_nm"),
            "ClassI_lead_presentation_rank": model.get("class1_lead_presentation_rank"),
            "best_target_PRIME_rank": model.get("best_target_PRIME_rank"),
            "ClassII_context_rank": model.get("class2_best_context_rank"),
            "ClassII_best_core_overlaps_target": model.get("class2_best_core_overlaps_recorded_target"),
            "ClassII_allele_qualification": model.get("class2_allele_qualification"),
            "interpretation": NOTES[region["gene"]], "clinical_release_ready": False,
        }
        audit.append(full)
        if core:
            lead = model["class1_lead_peptide"]
            allele = model["class1_lead_allele"]
            affinity = float(model["class1_lead_affinity_nm"])
            core_table.append({
                "Gene / change": f"{region['gene']} {region.get('protein_change') or ''}".strip(),
                "AA": len(region["sequence"]),
                "Coding / opposite mutant fragments": (
                    f"{region['coding_strand_RNA_ALT_fragments']} / {region['opposite_strand_RNA_ALT_fragments']}"
                    if region.get("coding_strand_RNA_ALT_fragments") is not None else "n/a (fusion)"
                ),
                "Predicted lead / HLA / nM": f"{lead} / {allele} / {affinity:.2f}",
                "Limitation": NOTES[region["gene"]],
            })
        else:
            reserve_table.append({"Gene": region["gene"], "Source": region["source"], "AA": len(region["sequence"]), "Reason for reserve": NOTES[region["gene"]]})
    table(HERE / "review/antigen-review.tsv", audit)
    (HERE / "review/antigen-review.md").write_text(
        "# Candidate assessment\n\n"
        "The main protein encodes 15 mutation windows and one in-frame fusion window. Eight candidates remain in reserve. "
        "All are research candidates; presentation, specificity and benefit have not been established.\n\n"
        "## Encoded regions\n\n"
        + markdown(["Gene / change", "AA", "Coding / opposite mutant fragments", "Predicted lead / HLA / nM", "Limitation"], core_table)
        + "\nAffinities and HLA assignments in this table come from `immunogenicity/antigen-model-summary.tsv`. "
        "They describe the recorded lead in the engineered context, not measured binding. Per-region target and context scores are retained in "
        "[the TSV](antigen-review.tsv) and [prediction files](../immunogenicity/). Fragment counts are read-name-based, not UMI molecules. "
        "Fusion support is reported separately: 200 recorded junction reads and one spanning fragment in distinct, non-UMI units; these are not 201 independently counted molecules.\n\n"
        "The 12 mutation windows with at least two coding-strand mutant fragments sum to 234; FES and SYT12 bring the coding total to 236. "
        "The nine opposite-strand fragments include KCTD19's two. These sums do not measure tumour coverage or clonality. "
        "The KCTD19-excluded [481-aa protein](../strand-strict-15/construct/protein/cassette.fasta) and its "
        "[1,586-nt known RNA segment](../strand-strict-15/construct/mrna/mrna-segments.fasta) remain a computational alternative.\n\n"
        "## Reserves\n\n"
        + markdown(["Gene", "Source", "AA", "Reason for reserve"], reserve_table)
        + "\nReserves are not encoded. Their earlier expression and model scores do not establish malignant-cell origin or normal-organ safety. "
        "See [the comparison](comparison.md), [methods](../methods.md), [safety notes](../SAFETY.md) and [disclaimer](../DISCLAIMER.md).\n"
    )

    d2_selected = rows(HERE / "exports/design2-selected.tsv")
    d2 = load(HERE / "exports/design2-protein.json")
    d2_mrna = load(HERE / "exports/design2-mrna-manifest.json")
    d3_selected = rows(HERE / "exports/design3-selected.tsv")
    d3 = load(HERE / "exports/design3-protein.json")
    d3_mrna = load(HERE / "exports/design3-mrna-manifest.json")
    current = load(HERE / "construct/protein/construct.json")
    current_mrna = load(HERE / "construct/mrna/construct.json")
    strict_selected = load(HERE / "strand-strict-15/selection/selected.json")
    strict = load(HERE / "strand-strict-15/construct/protein/construct.json")
    strict_mrna = load(HERE / "strand-strict-15/construct/mrna/construct.json")
    strict_receipt = load(HERE / "strand-strict-15/contingency-receipt.json")
    strand = load(HERE / "exports/rna-strand-audit.json")["mutation_strand_records"]
    groups = [d2_selected, d3_selected, selected, strict_selected]
    paired_mutation_ids = [
        [r["event_id"] for r in group if r["source"] == "MUTATION" and r["event_id"] in strand]
        for group in groups
    ]
    guardant_vus_counts = [
        sum(r["source"] == "MUTATION" and r["event_id"].startswith("guardant_") for r in group)
        for group in groups
    ]
    proteins = [d2, d3, current, strict]
    mrnas = [d2_mrna, d3_mrna, current_mrna, strict_mrna]
    added_linkers = [
        d2["solver"]["objective_tuple"][3],
        d3["solver"]["objective_tuple"][3],
        current["solver"]["objective_tuple"][3],
        strict_receipt["solver_objectives"]["added_residues"],
    ]
    comparison = [
        ("Encoded regions", *[len(group) for group in groups]),
        ("Paired-DNA/RNA mutation windows", *[len(ids) for ids in paired_mutation_ids]),
        ("Reference-derived Guardant cfDNA VUS windows", *guardant_vus_counts),
        ("Paired-DNA mutations with >=2 coding-strand mutant fragments", *[sum(strand[x]["coding_strand_alt_fragments"] >= 2 for x in ids) for ids in paired_mutation_ids]),
        ("Paired-DNA mutations with one coding-strand mutant fragment", *[sum(strand[x]["coding_strand_alt_fragments"] == 1 for x in ids) for ids in paired_mutation_ids]),
        ("Paired-DNA mutations with no coding-strand mutant fragments", *[sum(strand[x]["coding_strand_alt_fragments"] == 0 for x in ids) for ids in paired_mutation_ids]),
        ("Coding-strand mutant fragments across paired-DNA mutations", *[sum(strand[x]["coding_strand_alt_fragments"] for x in ids) for ids in paired_mutation_ids]),
        ("RNA-only splice windows", *[sum(r["source"] == "SPLICE" for r in group) for group in groups]),
        ("Intentional-self/CTA windows", *[sum(r["source"] == "CTA/SELF" for r in group) for group in groups]),
        ("Fusion windows", *[sum(r["source"] == "FUSION" for r in group) for group in groups]),
        ("Protein length, aa (inserted linker residues)", *[f"{p['protein_length_aa']} ({lk})" for p, lk in zip(proteins, added_linkers)]),
        ("Sense coding, nt / known RNA segment, nt", *[f"{m['lengths']['sense_coding_nt']} / {m['lengths']['known_rna_segments_nt'][0]}" for m in mrnas]),
        ("Longest coding paired-bracket run, nt (target <=14 nt)", *[m["optimization_metrics"]["coding_longest_paired_bracket_run"] for m in mrnas]),
        ("Class-I junction peptide-HLA occurrences", *[r["totals"]["peptide_hla_rows"] for r in proteins]),
        ("MHCflurry affinity percentile <=0.5% rows", *[r["totals"]["strong_rank"] for r in proteins]),
        ("MHCflurry affinity percentile <=2% rows", *[r["totals"]["mild_rank"] for r in proteins]),
        ("MHCflurry predicted affinity <=500 nM rows", *[r["totals"]["affinity_500"] for r in proteins]),
        ("<=2% rows / all junction-HLA rows (%)", *[round(100 * r["totals"]["mild_rank"] / r["totals"]["peptide_hla_rows"], 4) for r in proteins]),
        ("<=500 nM rows / all junction-HLA rows (%)", *[round(100 * r["totals"]["affinity_500"] / r["totals"]["peptide_hla_rows"], 4) for r in proteins]),
        ("Lowest predicted junction affinity, nM", *[round(r["totals"]["min_affinity_nm"], 2) for r in proteins]),
    ]
    headers = ["Metric", "Design 2", "Design 3", "Design 4 main", "Design 4 without KCTD19", "How to read this"]
    values = {row[0]: row[1:] for row in comparison}
    def comparison_row(label, key, interpretation):
        return dict(zip(headers, (label, *values[key], interpretation)))
    composition = [
        comparison_row("Encoded regions", "Encoded regions", "No known clinical direction: more candidates may add breadth and uncertainty."),
        comparison_row("Paired-DNA/RNA mutation windows", "Paired-DNA/RNA mutation windows", "Somatic mutation windows from the paired liver WES/WTS and matched-normal WBC pool; not independent cancer clones."),
        comparison_row("Reference-derived Guardant cfDNA VUS windows", "Reference-derived Guardant cfDNA VUS windows", "Design 3 adds five Oct 2025 plasma cfDNA VUS hypotheses (AR P555S, DHX36 R160K, ERBB4 D113V, KMT2D R2789Q, RAF1 G535S) translated from reference sequences rather than reconstructed from liver RNA."),
        comparison_row("Paired-DNA mutations with >=2 coding-strand mutant fragments", "Paired-DNA mutations with >=2 coding-strand mutant fragments", "Evaluated in Jan 2025 liver WTS under fr-firststrand orientation; does not show expression in every cancer cell or unbiopsied lesion."),
        comparison_row("Paired-DNA mutations with one coding-strand mutant fragment", "Paired-DNA mutations with one coding-strand mutant fragment", "SYT12 in Designs 2/3/4 plus the 15-aa FES window in Design 4; shallow RNA sampling, not a safety score."),
        comparison_row("Paired-DNA mutations with no coding-strand mutant fragments", "Paired-DNA mutations with no coding-strand mutant fragments", "KCTD19 (0 coding-strand and 2 opposite-strand mutant fragments in Jan 2025 liver WTS). Design 4 without KCTD19 removes this unconfirmed expression assumption; clinical impact is unknown."),
        comparison_row("Coding-strand mutant fragments across paired-DNA mutations", "Coding-strand mutant fragments across paired-DNA mutations", "Read-name fragment sum across included paired-DNA mutation events in Jan 2025 liver WTS (Design 3 adds H2BC12 [+10]; Design 4 also adds TJP2 [+29] and FES [+1]); not UMI molecules or clonality."),
        comparison_row("RNA-only splice windows", "RNA-only splice windows", "Design 2 encodes F5, MUC5B, MUC6, MUC3A and HLA-B; Design 3 uses junction-only windows for F5, HLA-B, HSF1, NBPF14 and KMT2D; Design 4 omits RNA-only splice windows."),
        comparison_row("Intentional-self/CTA windows", "Intentional-self/CTA windows", "Designs 2 and 3 include PRAME, SPATA31A6, LRRC36, YBX2 and TEKT5; Design 4 omits unmutated CTA/self windows."),
        comparison_row("Fusion windows", "Fusion windows", "Designs 2 and 3 encode both CDH1::NRG1 hypotheses; Design 4 encodes only the in-frame chr16:68819425:+ -> chr8:32727949:+ join (Caris prints exon 10:exon 6 on NM_001317184.2/NM_001159999.3; the NRG1 exon ordinal is transcript-dependent, so the breakpoint rather than an exon number is the stable identifier; see review/clinical-reconciliation.json)."),
        comparison_row("Protein length, aa (inserted linker residues)", "Protein length, aa (inserted linker residues)", "Includes fixed 26-aa SEC and 55-aa MITD routing domains; shorter length is an engineering property, not established safety or benefit."),
        comparison_row("Sense coding, nt / known RNA segment, nt", "Sense coding, nt / known RNA segment, nt", "Known RNA segment before unspecified PolyA (63-nt 5' UTR + sense coding + 9-nt triple stop + 71-nt 3' UTR)."),
        comparison_row("Longest coding paired-bracket run, nt (target <=14 nt)", "Longest coding paired-bracket run, nt (target <=14 nt)", "ViennaRNA 37 C MFE coding structure after repeat removal. Design 4 applies post-repeat synonymous partner-codon rescue (structure_rescue.py) while preserving all protected alternative-frame stop nucleotides; not a measured stability or PKR assay."),
    ]
    engineering = [comparison_row("Class-I peptide-HLA occurrences assessed", "Class-I junction peptide-HLA occurrences", "The denominator of this screen; not useful tumour targets or measured immune responses.")]
    engineering.append(comparison_row("MHCflurry affinity percentile <=0.5% flags", "MHCflurry affinity percentile <=0.5% rows", "Fewer is preferred by this computational screen only."))
    for count_key, label in (("mild_rank", "MHCflurry affinity percentile <=2% flags"), ("affinity_500", "MHCflurry predicted affinity <=500 nM flags")):
        cells = [f"{r['totals'][count_key]} ({100*r['totals'][count_key]/r['totals']['peptide_hla_rows']:.4f}%)" for r in proteins]
        engineering.append(dict(zip(headers, (label, *cells, "Fewer flags is a computational goal. Fractions describe different denominators; they are not the solver objective or a toxicity rate."))))
    engineering.append(comparison_row("Lowest predicted junction affinity, nM", "Lowest predicted junction affinity, nM", "Higher means weaker predicted binding of an unintended joining peptide. This does not establish safety."))
    engineering.append(dict(zip(headers, ("Lowest MHCflurry junction affinity percentile", *[f"{r['totals']['min_rank']:.4f}%" for r in proteins], "Higher is preferred by this screen. Affinity percentiles and affinity in nM can have different orderings across alleles."))))
    main_models = load(HERE / "immunogenicity/junction-multi-model-summary.json")
    alt_models = load(HERE / "strand-strict-15/immunogenicity/junction-multi-model-summary.json")
    main_safety = load(HERE / "safety/safety-summary.json")["in_silico_metrics"]
    alt_safety = load(HERE / "strand-strict-15/safety/safety-summary.json")["in_silico_metrics"]
    extra_headers = ["Metric", "Design 4 main", "Design 4 without KCTD19", "How to read this"]
    model_rows = [
        ("MixMHCpred rank <=0.5% flags", main_models["MixMHCpred_rank_le_0_5pct"], alt_models["MixMHCpred_rank_le_0_5pct"], "Fewer is the model-screen preference; no relative difference here."),
        ("PRIME rank <=0.5% flags", main_models["PRIME_rank_le_0_5pct"], alt_models["PRIME_rank_le_0_5pct"], "One flagged junction remains in each option; it is not a measured immune response."),
        ("Class-II strong synthetic-core rows", main_models["ClassII_synthetic_9aa_core_rank_le_0_5pct"], alt_models["ClassII_synthetic_9aa_core_rank_le_0_5pct"], "Strong means rank <=0.5%; zero does not exclude unpredicted presentation."),
        ("Class-II mild synthetic-core rows", main_models["ClassII_synthetic_9aa_core_rank_le_2pct"], alt_models["ClassII_synthetic_9aa_core_rank_le_2pct"], "Mild means rank <=2%; fewer is a finite optimization preference."),
        ("Class-II mild input-core overhang rows", main_models["ClassII_native_core_overhang_rank_le_2pct"], alt_models["ClassII_native_core_overhang_rank_le_2pct"], "The core stays within an input antigen region while the 15-mer touches an engineered join. It is not necessarily a healthy-tissue self peptide; clinical direction is unknown."),
        ("All mild junction-touching Class-II rows", main_models["ClassII_15mer_junction_context_rank_le_2pct"], alt_models["ClassII_15mer_junction_context_rank_le_2pct"], "Synthetic-core plus input-core-overhang rows. This descriptive total is not a calibrated risk score."),
        ("Eligible Class-I homolog queries screened", main_safety["eligible_class1_homolog_queries_screened"], alt_safety["eligible_class1_homolog_queries_screened"], "Target + junction binders screened across GENCODE v46, v50 and 112,389 patient-normal contexts (including all 24 alternative-specific junction queries)."),
        ("Same-HLA 1-mismatch paralog / <=250 nM 2-mismatch / <=100 nM lead Titin-class pairs", f"{main_safety['offtarget_1mm_same_hla_paralog_pairs']} / {main_safety['offtarget_2mm_strong_same_hla_pairs']} / {main_safety['offtarget_titin_class_lead_mimic_pairs']}", f"{alt_safety['offtarget_1mm_same_hla_paralog_pairs']} / {alt_safety['offtarget_2mm_strong_same_hla_pairs']} / {alt_safety['offtarget_titin_class_lead_mimic_pairs']}", "Heuristic model-prioritized homolog pairs; not measured T-cell cross-reactivity."),
        ("Stage-1 Class-I / Class-II specificity control peptides", f"{main_safety['wetlab_negative_control_peptides_to_synthesize']} / {main_safety['wetlab_class2_wt_15mers_to_synthesize']}", f"{alt_safety['wetlab_negative_control_peptides_to_synthesize']} / {alt_safety['wetlab_class2_wt_15mers_to_synthesize']}", "Includes 100% of shared-HLA own-wild-type Class-I counterparts in each construct; does not exhaust distant cross-recognition."),
        ("Sense alternative-frame 8-11-mers screened / <=0.5% rank (min-nM allele / any patient allele)", f"{main_safety['sense_cryptic_8to11mers_screened']} / {main_safety['sense_cryptic_orf_rank_le_0_5pct_on_min_nm_allele']} / {main_safety['sense_cryptic_orf_rank_le_0_5pct_any_patient_allele']}", f"{alt_safety['sense_cryptic_8to11mers_screened']} / {alt_safety['sense_cryptic_orf_rank_le_0_5pct_on_min_nm_allele']} / {alt_safety['sense_cryptic_orf_rank_le_0_5pct_any_patient_allele']}", "Evaluates every stop-delimited sense cryptic peptide across all six patient Class-I alleles, separating rank <=0.5% on each peptide's minimum-nM allele from rank <=0.5% on any patient allele; does not measure leaky scanning or frameshifting."),
        ("Hypothetical antisense-IVT-byproduct 8-11-mers screened / <=0.5% rank (min-nM allele / any patient allele)",
         " / ".join(str(main_safety["cryptic_peptides_by_stored_strand"]["antisense_IVT_byproduct_only"][key]) for key in ("unique_peptides_screened", "rank_le_0_5pct_on_min_nm_allele", "rank_le_0_5pct_any_patient_allele")),
         " / ".join(str(alt_safety["cryptic_peptides_by_stored_strand"]["antisense_IVT_byproduct_only"][key]) for key in ("unique_peptides_screened", "rank_le_0_5pct_on_min_nm_allele", "rank_le_0_5pct_any_patient_allele")),
         "Unique peptides assigned to the stored antisense_IVT_byproduct_only group (sense takes precedence for duplicates). These binding predictions are separate from admitted TARGET/JUNCTION homolog scans and do not establish IVT impurities, translation or presentation; counts use the reported rounded ledger."),
        ("Known RNA MFE max uninterrupted / bulge-tolerant stem, bp", f"{main_safety['known_RNA_mfe_max_uninterrupted_stem_bp']} / {main_safety['known_RNA_mfe_max_bulge_tolerant_stem_bp']}", f"{alt_safety['known_RNA_mfe_max_uninterrupted_stem_bp']} / {alt_safety['known_RNA_mfe_max_bulge_tolerant_stem_bp']}", "Single 37 C ViennaRNA MFE structure; not a validated PKR/MDA5 biological activation test."),
    ]
    extras = [dict(zip(extra_headers, (label, mval, aval, note))) for label, mval, aval, note in model_rows]
    records = composition + engineering
    table(HERE / "review/comparison.tsv", records)
    (HERE / "review/comparison.md").write_text(
        "# How Design 4 compares with Design 3 and Design 2\n\n"
        "[Design 3](../../design-3/README.md) (`30` regions, `885 aa`, `2,655 nt` sense coding) is PopVax's current lead candidate, and "
        "[Design 2](../../design-2/README.md) (`24` regions, `724 aa`, `2,172 nt` sense coding) is the earlier matched-normal baseline. "
        "Design 4 main (`16` regions, `507 aa`, `1,521 nt` sense coding) shares **13 exact `region_id`s** (byte-identical 25-aa sequences) with Design 3: "
        "11 paired-DNA/RNA mutations from Design 2 (`KLF5`, `SF1`, `AQP3`, `RXRA`, `HSPG2`, `TGFBR2`, `PPIC`, `VPS11`, `ARHGAP40`, `PLEC`, `SYT12`), "
        "`H2BC12` (`REG_e3aa30c984295df98692`, independently retained with missing vendor expression), and the in-frame `CDH1::NRG1` fusion (`REG_1f7babf47bbafeebeedc`). "
        "It also encodes the same `KCTD19 p.G920R` somatic event (`mn_chr16_67289592_C_T`), shortened from 25 aa (`REG_9dde223e9bf58444173b`, `MDLLIQRGLSRSVSYSILRKYLQED` in Designs 2 and 3) to 23 aa (`REG_6c716df1b65d721ee904`, `LLIQRGLSRSVSYSILRKYLQED` in Design 4 main, trimming two N-terminal unmutated residues `MD` while preserving all recorded target-overlapping peptides), "
        "adds `TJP2 p.S886P` (`REG_ef71940c452b9b68a401`, `25 aa`) and the unpadded C-terminal `FES p.R804W` window (`REG_0d4d002ad8a14e5a2f75`, `15 aa`), "
        "and omits Design 3's five reference-derived Guardant360 plasma VUS windows (`AR p.P555S`, `DHX36 p.R160K`, `ERBB4 p.D113V`, `KMT2D p.R2789Q`, `RAF1 p.G535S`), "
        "five junction-only RNA splice windows (`F5`, `HLA-B`, `HSF1`, `NBPF14`, `KMT2D`; versus Design 2's full-window splice set `F5`, `MUC5B`, `MUC6`, `MUC3A`, `HLA-B`), five intentional-self/CTA windows (`PRAME`, `SPATA31A6`, `LRRC36`, `YBX2`, `TEKT5`), and the alternative `CDH1--NRG1` fusion hypothesis (`REG_98875bfa1d09b7fab791`).\n\n"
        "**Design 4 without KCTD19 is a separate 15-region construct alternative.** It removes that DNA-supported candidate because coding-strand mutant RNA is unconfirmed (`0` coding-strand and `2` opposite-strand mutant fragments in Jan 2025 liver WTS), "
        "then reoptimizes the order, linkers and known-segment PopRNA/GenBank nucleotide encoding (`481 aa`, `1,443 nt` sense coding, `1,586 nt` known RNA segment before unspecified PolyA). "
        "The column is not a list of exclusions, and neither option is established as clinically safer or ready for clinical use.\n\n"
        "## Composition, evidence and mRNA structure\n\n"
        + markdown(headers, composition)
        + "\nThe eight Design 4 reserves (`SERPINA1`, `POLE`, `PRAME`, `NUTM1`, `MAGEA6`, `TSPYL6`, `MAGEA1`, `CDH1--NRG1`) come from its exploratory candidate pool; they are not simply the candidates removed from Design 2 or Design 3. "
        "Design 2 excluded `H2BC12` after unsuccessful expression-symbol checks, `TJP2`/readthrough for missing expression and ambiguous gene attribution, and `FES` for its shorter-than-25-aa window; Design 3 independently added `H2BC12` with missing expression while retaining the 25-aa window rule and adding five reference-derived Guardant360 VUS calls. "
        "In the three-BAM locus audit (`../exports/ctdna-locus-pileups.tsv`), none of those five Guardant360 VUS calls were called by paired Mutect2 in the Jan 2025 liver biopsy (`3/2,560` raw liver WES ALT/REF fragments for `KMT2D R2789Q`, `1/870` for `ERBB4 D113V`, and `0` for `AR P555S`, `DHX36 R160K` and `RAF1 G535S`; `1/13` and `1/168` coding-strand WTS ALT/REF fragments for `AR P555S` and `RAF1 G535S`, and `0` for the other three), and their lesion origin and transcription (`0.2%-0.4%` Oct 2025 cfDNA VAF for `DHX36 R160K`, `RAF1 G535S`, `ERBB4 D113V` and `AR P555S`; `3.2%` Oct 2025 and `0.2%` Mar 2025 cfDNA VAF for `KMT2D R2789Q`) remain unestablished.\n\n"
        "## Class-I junction screen\n\n"
        "A junction is an engineered join between encoded regions or linkers. A peptide-HLA occurrence is one positioned peptide evaluated for one allele; "
        "repeated sequences may count more than once. These rows concern joining peptides being screened, rather than intended tumour targets. "
        "For these junctions, a larger nM value means weaker predicted binding. Lower nM can be desirable when selecting an intended target; "
        "neither direction proves useful T-cell recognition.\n\n"
        + markdown(headers, engineering)
        + "\n## Additional predictors for the two Design 4 options\n\n"
        "Designs 2 and 3 optimize junction order and linkers using a 4-phase MHCflurry Class-I affinity objective (in Design 3, MHCflurry Class-I presentation and CapHLA Class-II predictions are recorded as candidate annotations on source regions rather than junction optimization terms). "
        "Corresponding matched multi-model junction analyses were not freshly generated for every Design 2 or Design 3 junction query, so Designs 2 and 3 are not assigned zeros or ranked in this table. "
        "Class-II results in Design 4 use nine conditional HLA definitions with MixMHC2pred 2.1-beta1. Model percentile scales are not interchangeable or probabilities of response or harm.\n\n"
        + markdown(extra_headers, extras)
        + "\nThe 15-region alternative is not uniformly preferable: its Class-I <=2% fraction and its total mild Class-II count are higher, "
        "although some other screen counts are lower. Removing KCTD19 removes an unconfirmed expression assumption and a potential target. "
        "No option has demonstrated better patient safety, cancer control or survival. Selection for further research requires scientific judgment, not choosing a column by its largest or smallest numbers.\n\n"
        "Sources: [candidate review](antigen-review.md), [Design 3 reference receipt](../exports/design3-reference-receipt.json), [main construct](../construct/protein/construct.json.gz), "
        "[alternative construct](../strand-strict-15/construct/protein/construct.json.gz), [methods](../methods.md), [safety questions](../SAFETY.md).\n"
    )

    old_matches = load(HERE / "references/peptide-matches.json")
    new_matches = load(HERE / "references/v50-audit/peptide-matches.json")
    checked, all_target_queries, new_core_hits = set(), set(), []
    for region in selected:
        for _, peptide in peptide_windows(region["sequence"], int(region["target_start_in_window"]), int(region["target_end_in_window"]), (8, 9, 10, 11, 15)):
            all_target_queries.add(peptide)
            if not old_matches[peptide]["occurs"]:
                checked.add(peptide)
                if new_matches[peptide]["occurs"]:
                    new_core_hits.append({"gene": region["gene"], "peptide": peptide})
    if new_core_hits:
        raise ValueError("Current-reference matches require curation review: " + str(new_core_hits))
    save(HERE / "references/v50-audit/core-corroboration.json", {
        "source_v46_matches_sha256": sha(HERE / "references/peptide-matches.json"),
        "source_v50_matches_sha256": sha(HERE / "references/v50-audit/peptide-matches.json"),
        "core_novel_target_sequences_checked": len(checked), "new_core_self_hits": new_core_hits,
        "all_target_overlap_unique_queries_lengths_8_9_10_11_15": len(all_target_queries),
        "queries_excluded_for_existing_v46_match": len(all_target_queries - checked),
        "scope": "Reference-absent target-overlap queries after filtering; not every encoded peptide or an independent clinical safety endpoint.",
        "additional_self_hits_in_all_queries": sum(not old_matches[p]["occurs"] and new_matches[p]["occurs"] for p in old_matches),
        "normal_tissue_presentation_or_TCR_safety_proven": False,
    })
    save(HERE / "review/results-summary.json", results_summary_snapshot())
    ctdna_by_vid = {r["variant_label"]: r for r in rows(HERE / "alphagenome2/ctdna-progression-audit.tsv")}
    kmt2d = ctdna_by_vid["KMT2D_R2789Q"]
    save(HERE / "review/clinical-reconciliation.json", {
        "clinical_pdf_sources": [
            "stanford-pathology-ihc-2025-01.pdf",
            "stanford-stamp-2025-01-31.pdf",
            "caris-tissue-2025-02-22.pdf",
            "guardant360-2025-10-21.pdf",
        ],
        "project_and_narrative_sources": [
            "https://github.com/popvax/yash/blob/7ece8f7a77c8db935422f6502d2c4db2352d7bf0/data/clinical-timeline.md",
            "https://github.com/popvax/yash/blob/7ece8f7a77c8db935422f6502d2c4db2352d7bf0/data/samples.tsv",
            "https://rnafoundry.com/yash/",
        ],
        "reconciliation_scope": (
            "Structured clinical values reconciled against the four repository PDF reports (covering Jan-Oct 2025 "
            "pathology, tissue sequencing and plasma cfDNA testing), project narrative/sample metadata "
            "(data/clinical-timeline.md, data/samples.tsv, rnafoundry.com/yash/), DNA_TN25-112103.vcf, "
            "all-unique.vcf.gz, normal.coding.g.vcf.gz and reproducible locus pileups in exports/ctdna-locus-pileups.tsv."
        ),
        "specimen": "CT-guided liver FNA and core biopsy Block A2 (case TN25-112103), collected 2025-01-15; primary site reported as Pancreas, NOS / stage IV pancreatic ductal adenocarcinoma.",
        "histology_tumour_cellularity": (
            "Approximately 60-70% histological tumour cellularity in Stanford pathology (~60% on STAMP); "
            "remaining tissue includes background liver parenchyma/stroma. Histological cellularity does not "
            "measure the malignant vs non-malignant bulk RNA transcript mixture fraction."
        ),
        "clinically_reported_fusion": {
            "caris_report_annotation": (
                "Caris RNA-Tumor Seq (pp. 2, 17) reports pathogenic in-frame CDH1:NRG1 fusion with printed splice "
                "site exon 10:exon 6 and printed transcript IDs NM_001317184.2/NM_001159999.3."
            ),
            "canonical_refseq_transcript_equivalence": (
                "On canonical RefSeq CDH1 transcript NM_004360.5 (not printed on the Caris report), the same CDH1 "
                "3' boundary (chr16:68819425:+) corresponds to exon 11."
            ),
            "nrg1_refseq_exon_ordinal_equivalence": (
                "The NRG1 5' boundary chr8:32727949:+ is the first base of an exon in all 23 curated RefSeq NRG1 "
                "transcripts (UCSC ncbiRefSeqCurated, hg38), but the exon ordinal is transcript-dependent: exon 6 "
                "on the Caris-printed NM_001159999.3 and on NM_001160004.3, NM_001160008.2, NM_004495.4, "
                "NM_013956.5, NM_013957.5, NM_013958.4, NM_013960.5 and NM_013964.5; exon 5 on NM_001159995.3, "
                "NM_001160002.2 and NM_001322201.2; exon 4 on NM_001160001.3, NM_001160005.2, NM_001160007.2, "
                "NM_001322202.2 and NM_013962.3; and exon 2 on NM_001159996.3, NM_001322197.2, NM_001322205.2, "
                "NM_001322206.2, NM_001322207.2 and NM_013959.4. The genomic breakpoint, not an exon ordinal, is "
                "the stable identifier for this join. These ordinals are external RefSeq annotation geometry, not "
                "a measurement re-derived from this project's sequencing evidence."
            ),
            "guardant360_cfdna_detection": "Detected in plasma cfDNA at 0.4% on 2025-03-29 and 9.3% on 2025-10-14.",
            "project_genomic_breakpoint": "chr16:68819425:+ -> chr8:32727949:+ in project RNA alignment exports (1 spanning fragment in RNA_TN25-112103.bam).",
            "stanford_stamp_status": "The Stanford STAMP report does not report this fusion.",
            "immunological_caveat": "Fusion transcript detection does not establish translation or MHC presentation of the junction peptide.",
        },
        "immunohistochemistry_and_wts_context": {
            "stanford_ihc": "Stanford pathology/IHC physical pp. 1-2: CK7 positive, CK20 focal positive, MLH1/MSH2/MSH6/PMS2 intact, HER2 (4B5) negative (0), SATB2 weak subset positive, beta-catenin negative (membranous staining pattern); 2025-01-21 addendum: Claudin-18.2 (43-14A) limited expression (2-3+, 50%; below >=75% companion-diagnostic threshold).",
            "caris_ihc_and_cna": "ERBB2/HER2 (4B5) negative (Score 0) and CNA-Seq Amplification Not Detected; PD-L1 (SP142) negative (0%); MSI Stable; TMB Low (3 mut/Mb); Genomic LOH Low (8% across 552 autosomal segments). The Caris p. 12 karyotype plots complex segmental gains and losses. Representative gains include chr1q, chr13q and segments of chr7, chr20q and chr22; representative losses include 1p, chr3, 4q, 5p, distal 9p and chr21. Chr6 appears at plotted copy number 2. These are nonexhaustive visual observations, not validated gene-level or allele-specific copy-number calls.",
            "caris_wts_cohort_percentiles": "Selected cohort percentiles (p. 11): VEGFA 99, SMAD4 99, ERBB3 94, NECTIN4 94, PRAME 92, LAG3 91, ADORA2A 90, CLDN6 90, NRG1 86, PDCD1 (PD-1) 86, PDCD1LG2 78, CTLA4 75, CD274 (PD-L1) 72, CLDN18 57, ERBB2 22, MSLN 8, TACSTD2 (TROP2) 4. Bulk WTS mixes tumour, immune, stromal and liver transcripts.",
        },
        "assay_specific_TMB": "Stanford STAMP (physical pp. 1-2) reports research-use tissue TMB 7.0 mut/Mb; Caris (pp. 1-2) reports tissue TMB 3 mut/Mb; Guardant360 (p. 1) reports blood TMB 7.59 mut/Mb. Different assays, genomic scopes and specimen types prevent treating these as a calibrated longitudinal TMB series. None establishes vaccine response.",
        "pharmacogenomics": "Guardant360 (p. 2) reports UGT1A1 *1/*28 Intermediate Metabolizer (increased irinotecan toxicity risk) and CYP2D6 *4/*41 Intermediate Metabolizer (activity score 0.25); DPYD, HLA-B*57:01 and TPMT had no abnormal reportable alleles detected.",
        "longitudinal_cfDNA_and_narrative_context": {
            "historical_BCOR_call_outside_executed_locus_audit": "Guardant360 p. 3 reports BCOR L1632P at 0.2% on 2025-03-29 and Not Detected on 2025-10-14. This historical protein-only call was omitted from the executed 13-locus nucleotide audit and is disclosed here separately; no BCOR BAM or MHC prediction result is claimed. BCOR reference transcripts are available and a reference-based SNV candidate can be inferred, but the report row does not provide a transcript accession or nucleotide allele to authenticate that mapping.",
            "draws_compared": "Guardant360 plasma draws on 2025-03-29 (epigenomic Tumour Fraction 1.4%) and 2025-10-14 (reported 2025-10-21; epigenomic Tumour Fraction 31.9%, blood TMB 7.59 mut/Mb, MSI-High Not detected, HRD Not detected).",
            "project_narrative_timeline": (
                "Per data/clinical-timeline.md and rnafoundry.com/yash/ (outside the four clinical PDFs): initial "
                "response to zenocutuzumab after five FOLFIRINOX cycles (CA19-9 ~45,000 -> 18), followed in Aug-Oct "
                "2025 by lung and peritoneal progression while liver and pancreatic lesions were described as "
                "controlled; an Oct 2025 peritoneal biopsy yielded insufficient tissue; treatment switched in Nov "
                "2025 in India to gemcitabine/nab-paclitaxel plus rash-limited afatinib."
            ),
            "liver_biopsy_variants_also_detected_in_oct_2025_cfDNA": (
                "Seven Guardant360-covered somatic alterations present in the 2025-01-15 liver biopsy were also "
                "detected in 2025-10-14 cfDNA (clonality and anatomical lesion attribution remain unestablished "
                "from pooled plasma cfDNA): ROS1 F1699F (1.6% -> 21.4%, synonymous SNV), TERT c.-124C>T (0.6% -> "
                "20.5%, noncoding promoter SNV), MERTK E192* (1.0% -> 19.4%, stop-gain with 0 ALT / 25 REF liver "
                "WTS fragments [23 coding-strand REF], compatible with NMD or absent mutant transcription in the "
                "sampled biopsy rather than proof of NMD), JUN R276* (1.2% -> 14.8%, stop-gain creating no novel "
                "amino-acid sequence), RXRA A214S/A241S (0.5% -> 12.1%, encoded AG04), TGFBR2 N515Kfs*26/N490fs "
                "(0.7% -> 11.7%, encoded AG06) and CDH1::NRG1 (0.4% -> 9.3%, encoded AG24)."
            ),
            "other_oct_2025_cfDNA_and_potential_CHIP_variants": (
                f"KMT2D p.R2789Q (c.8366G>A, chr12:49039222:C>T) rose from 0.2% to 3.2% cfDNA VAF (16.0x across a "
                f"22.8x epigenomic tumour-fraction increase; cancer-cell fraction and lesion origin unestablished; "
                f"{kmt2d['jan_2025_liver_wes']} liver WES fragments; {kmt2d['jan_2025_liver_rna_all']} liver WTS "
                f"fragments; reference-model predictions show AlphaMissense={kmt2d['atlas_alphamissense']} and "
                f"unverified splice scores [Atlas MERGED_SPLICING={kmt2d['atlas_merged_splicing']}; observed lung "
                f"SPLICE_SITES max prob diff {kmt2d['lung_obs_splice_sites_max_prob_diff']}, SPLICE_SITE_USAGE max "
                f"frac diff {kmt2d['lung_obs_splice_usage_max_frac_diff']}, SPLICE_JUNCTIONS max log-count diff "
                f"{kmt2d['lung_obs_splice_junctions_max_log_diff']}] that do not establish splicing alteration "
                f"without progressing-lesion RNA). Five low-VAF calls appeared only on 2025-10-14 "
                "(AR P555S 0.4%, DHX36 L147L 0.3%, ERBB4 D113V 0.3%, DHX36 R160K 0.2%, RAF1 G535S 0.2%; lesion "
                "origin and transcription unestablished). Guardant360 annotates DNMT3A p.C351Afs*56 "
                "(chr2:25247121:CAAAA>C; 0.2% -> 0.1% cfDNA VAF; 0/58 normal WBC, 0/4,108 liver WES, 0/234 liver "
                "WTS fragments) as potential clonal hematopoiesis (potential CHIP); it is excluded from targeting."
            ),
        },
        "guardant360_panel_coverage_and_B2M_limit": (
            "Thirteen of Design 4's fifteen somatic mutation genes (KLF5, SF1, AQP3, HSPG2, PPIC, KCTD19, VPS11, ARHGAP40, PLEC, SYT12, H2BC12, TJP2, FES) are not on Guardant360's 744-gene panel (methods p. 8; gene list p. 9, Table 1), so their status in Oct 2025 plasma cfDNA and progressing lesions is unmeasured rather than excluded. Guardant360 lists B2M^ (covering SNVs, indels and promoter methylation '^') and reports no B2M alteration in this specimen; per Table 1 footnotes (p. 9), B2M is not marked for copy-number amplification ('\u2020') or deletion ('*'), and HLA-A/B/C carry '\u00a7' (SNVs/indels and somatic LOH not reported)."
        ),
        "isoform_qualification": "Protein numbering differs across reports; normalized nucleotide alleles are required to establish identity.",
        "HLA_ClassI_status": (
            "Caris tumour WES and Guardant360 blood cfDNA both report HLA-A*01:01, A*33:03, B*08:01, B*15:08, "
            "C*01:02 and C*07:02; Caris notes that tumour sequencing does not replace blood-based histocompatibility "
            "testing. No clinical blood-based histocompatibility typing report or phased Class-II typing report is "
            "supplied in the repository."
        ),
        "clonal_phylogeny_or_HLA_LOH_proven": False,
        "current_November_2026_tumour_expression_or_treatment_status_established": False,
        "matched_normal_blood_collection_date_reported_by_project": "2026-10-30",
        "metadata_source": "https://github.com/popvax/yash/blob/7ece8f7a77c8db935422f6502d2c4db2352d7bf0/data/samples.tsv",
        "raw_proprietary_reports_redistributed": False,
    })
    missing_data_rows = [
        {
            "priority_id": "MD-01",
            "open_question_or_evidence_gap": "Whether the 13 off-panel Design 4 somatic mutations persist and are transcribed in progressing lesions",
            "supplied_evidence_context": "Thirteen of Design 4's 15 somatic mutation genes are absent from Guardant360's 744-gene panel (report methods p. 8), and the Oct 2025 peritoneal biopsy narrative notes insufficient tissue while plasma epigenomic Tumour Fraction was 31.9%.",
            "bounded_scope_and_limitations": "Broader plasma cfDNA sequencing could test circulating DNA detection of off-panel mutations, but pooled cfDNA cannot assign variants to specific anatomical lesions (lung vs peritoneum vs liver/pancreas), cfDNA nondetection does not prove absence from low-shedding lesions, and DNA assays cannot measure mRNA/protein expression; no assay depth or biopsy procedure is prescribed here.",
        },
        {
            "priority_id": "MD-02",
            "open_question_or_evidence_gap": "Whether malignant cells retain intact B2M (15q21.1), HLA-A/B/C (6p21.3) copy number and cell-surface MHC presentation",
            "supplied_evidence_context": "Guardant360 reports no B2M SNV/indel/methylation alteration (B2M^) but does not mark B2M for copy-number amplification (\u2020) or deletion (*) or report somatic HLA-A/B/C LOH, and bulk biopsy DNA/RNA does not measure malignant cell-surface MHC density.",
            "bounded_scope_and_limitations": "Malignant-cell MHC-I surface expression is required for direct CD8+ T-cell recognition of tumour cells, whereas CD4+ T-cell responses to MHC-II epitopes are primarily primed and recalled via professional antigen-presenting cells rather than requiring tumour-intrinsic HLA-DR expression.",
        },
        {
            "priority_id": "MD-03",
            "open_question_or_evidence_gap": "Current (post-November 2025 through November 2026) lesion distribution, organ function and molecular evolution",
            "supplied_evidence_context": "Clinical PDF reports end on 2025-10-21 and project timeline notes end in Nov 2025, whereas matched-normal WBC DNA carries a 2026-10-30 collection date; Guardant360 also records UGT1A1 *1/*28 and CYP2D6 *4/*41.",
            "bounded_scope_and_limitations": "Identifies the 12-month temporal gap between the supplied tumour/cfDNA evidence and the matched-normal blood date; this computational package makes no clinical dosing, scheduling or patient-management recommendations.",
        },
        {
            "priority_id": "MD-04",
            "open_question_or_evidence_gap": "Confirmed 4-field germline HLA Class I/II typing and untransformed gene/isoform RNA quantification",
            "supplied_evidence_context": "Class-II alleles in Design 4 are unphased research inferences from tumour DNA/RNA (with DRB3 *02:02 vs *02:256 discordance), and vendor-transformed-expression.tsv omits H2BC12 and lacks raw count/TPM units.",
            "bounded_scope_and_limitations": "Computational HLA calling on research normal.recal.bam can cross-check germline alleles in silico but does not substitute for accredited clinical histocompatibility typing; whether raw Caris WTS count/isoform exports are obtainable outside the repository is unknown.",
        },
        {
            "priority_id": "MD-05",
            "open_question_or_evidence_gap": "Sense-strand mutant transcript presence for low-count loci (KCTD19, FES, SYT12) and T-cell mutant-vs-self discrimination",
            "supplied_evidence_context": "In Jan 2025 liver WTS, KCTD19 has 0 coding-strand and 2 opposite-strand mutant fragments while FES and SYT12 each have 1 coding-strand fragment; several encoded regions have 1-2 mismatch human homologs or shared-HLA wild-type counterparts.",
            "bounded_scope_and_limitations": "Clarifies what experimental evidence would distinguish coding-strand expression from antisense/artefactual reads and evaluate T-cell cross-reactivity against the 197/184 Class-I and 39/37 Class-II control peptides, within assay detection limits.",
        },
        {
            "priority_id": "MD-06",
            "open_question_or_evidence_gap": "Cellular origin of bulk immune/signalling transcripts in the tumour microenvironment",
            "supplied_evidence_context": "Caris IHC reports PD-L1 (SP142) 0% and HER2 Score 0, whereas bulk FFPE WTS cohort percentiles include VEGFA (99), ERBB3 (94), LAG3 (91), ADORA2A (90), NRG1 (86) and PDCD1 (86) in a liver biopsy reported at 60-70% histological tumour cellularity.",
            "bounded_scope_and_limitations": "Histological cellularity (60-70%) does not measure the malignant-versus-liver/stromal RNA transcript mixture fraction; bulk FFPE WTS mixes malignant, immune, stromal and background liver transcripts and cannot establish receptor co-expression on tumour-reactive T cells or clinical activity of any systemic or combination regimen.",
        },
    ]
    table(
        HERE / "review/missing-data-priorities.tsv",
        missing_data_rows,
        ["priority_id", "open_question_or_evidence_gap", "supplied_evidence_context", "bounded_scope_and_limitations"],
    )
    print("Generated candidate assessment, recorded comparison and bounded evidence summaries.")


if __name__ == "__main__":
    main()
