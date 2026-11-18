"""Bounded interpretation of cached research screens; no inference or clinical authorization."""
import csv
import gzip
import io
import json
import math
from pathlib import Path

from dataio import load, read_text_auto, resolve_path, sha

HERE = Path(__file__).resolve().parents[1]
MFE_REPORTING_CUTOFF_BP = 30
MFE_SCOPE = "Single ViennaRNA minimum-energy structure at 37 C with standard unmodified-RNA parameters; not all conformations or a biological activation assay."
RNA_INTERPRETATION = "Stem lengths describe one predicted minimum-energy structure. The 30-bp reporting cutoff is an author-selected computational comparison, not a validated PKR or MDA5 activation boundary. Innate immune activation and manufacturing impurities have not been measured."
CRYPTIC_INTERPRETATION = "The final CDS contains 98 alternative-frame stops and predicted stop-free runs up to 54 codons. These sequence counts do not measure leaky scanning or frameshifting. Alternative translation requires assessment with defined assay sensitivity; no experiment proves absolute absence."
MOTIF_INTERPRETATION = "The configured UUUU and UUUC motifs are absent from the CDS. This does not establish translation fidelity or safety for any nucleotide chemistry; chemistry and translation products remain untested."
SPECIFICITY_INTERPRETATION = "Candidate control for a qualified mutant/wild-type or homolog specificity study. No validated acceptance cutoff, negative biological result or absence of off-target harm is supplied."
SAFETY_LIMITS = "The reported exact-reference matches, model rankings, homolog rules and folding metrics are bounded computational results. They do not eliminate autoimmune targeting, aggregation, innate activation or other harm. Presentation, recognition, delivery, manufacturing quality, clinical safety and benefit remain unestablished."
BIOPSY_EXPRESSION_SCOPE = "Supplied vendor-transformed values from a mixed tumour-biopsy RNA matrix, not calibrated TPM or healthy-tissue expression. Missing gene symbols are unavailable, not recorded zeros. The maximum uses available genes in the stated scope and retains the existing three-decimal reporting precision; coverage and recorded-zero/positive counts refer to the source values."
CRYPTIC_STRAND_SCOPE = "Unique 8-11-mer sequences in the supplied stop-delimited ORF screen, assigned to the first enumerated strand (sense takes precedence if a sequence occurs on both). Antisense_IVT_byproduct_only queries are hypothetical reverse-complement IVT-byproduct sequences, not observed product impurities, translated antigens or measured presentation. These binding-prediction counts are separate from the homolog scan of admitted TARGET/JUNCTION queries. Threshold counts use the displayed rounded prediction ledger (nM: two decimals; ranks: four decimals), not new model inference."
SCREENING_POLICY = {
    "exact_reference_match_scope": {
        "TARGET": "An exact off-target reference match excludes a nominated TARGET query under this selection policy.",
        "JUNCTION": "Exact reference matches of junction queries are recorded for comparison; they are not automatically excluded by the nominated-target rule. Neither classification establishes clinical safety or harm.",
    },
    "prime_rank_scope": {
        "query_admission_max_percentile": 2.0,
        "homolog_classification_max_percentile": 0.5,
        "shared_hla_predicate_uses_prime": False,
        "interpretation": "Broad query admission and downstream homolog classification use different declared PRIME cutoffs. A query admitted by PRIME alone need not satisfy a homolog-priority rule.",
    },
    "shared_hla_predicate": "Each peptide on the same allele satisfies MHCflurry nM <=500 OR MHCflurry percentile <=2 OR MixMHCpred rank <=2. PRIME is not part of this predicate.",
    "gene_field_scope": {
        "homolog_genes": "All matching reference genes, including the query's own parent gene where matched; no display truncation.",
        "offtarget_gene_count": "Count of matching genes after excluding the query's own parent-gene/counterpart set; it need not equal the number of names in homolog_genes.",
        "human_genes": "In the specificity-comparison plate, off-target genes when present, otherwise the own-WT parent/counterpart genes. This is reference annotation, not measured tissue expression or cross-recognition.",
        "patient_biopsy_expression": BIOPSY_EXPRESSION_SCOPE,
    },
}


def reference_parent_genes(gene):
    """The existing own-parent mapping, shared by expression reporting and classification."""
    own = {gene}
    if gene == "TJP2":
        own.update(("AL358113.1", "FPGT-TJP2", "ENSG00000285130"))
    if gene == "CDH1::NRG1":
        own.update(("CDH1", "NRG1"))
    return own


def read_biopsy_expression(path=HERE / "exports/expression.tsv"):
    values = {}
    for line in read_text_auto(path).splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 2:
            raise ValueError("Expected exactly two source biopsy-expression columns")
        gene, raw = fields
        value = float(raw)
        if not gene or gene in values or not math.isfinite(value) or value < 0:
            raise ValueError("Invalid or duplicate source biopsy-expression record: " + gene)
        values[gene] = value
    return values


def biopsy_expression_fields(homolog_genes, source_gene, expression):
    genes = sorted(set(homolog_genes))
    off = sorted(set(genes) - reference_parent_genes(source_gene))
    considered = off or genes
    available = [g for g in considered if g in expression]
    missing = [g for g in considered if g not in expression]
    return {
        "max_patient_biopsy_wts_expr": round(max(expression[g] for g in available), 3) if available else None,
        "expression_gene_scope": "OFFTARGET_GENES" if off else "ALL_MATCHED_REFERENCE_GENES",
        "expression_genes_considered": ";".join(considered),
        "expression_genes_available": ";".join(available),
        "expression_genes_missing": ";".join(missing),
        "expression_genes_considered_count": len(considered),
        "expression_genes_available_count": len(available),
        "expression_genes_missing_count": len(missing),
        "expression_recorded_zero_gene_count": sum(expression[g] == 0 for g in available),
        "expression_recorded_positive_gene_count": sum(expression[g] > 0 for g in available),
        "expression_coverage_status": "UNAVAILABLE" if not available else ("PARTIAL" if missing else "COMPLETE"),
    }


def cryptic_counts_by_stored_strand(records):
    counts = {strand: {
        "unique_peptides_screened": 0,
        "affinity_nm_le_100": 0,
        "affinity_nm_le_500": 0,
        "rank_le_0_5pct_on_min_nm_allele": 0,
        "rank_le_0_5pct_any_patient_allele": 0,
    } for strand in ("sense", "antisense_IVT_byproduct_only")}
    for row in records:
        if row["strand"] not in counts:
            raise ValueError("Unknown cryptic stored-strand label: " + row["strand"])
        block = counts[row["strand"]]
        block["unique_peptides_screened"] += 1
        for source, key, cutoff in (
            ("affinity_nm", "affinity_nm_le_100", 100.0),
            ("affinity_nm", "affinity_nm_le_500", 500.0),
            ("mhcflurry_rank_on_min_nm_allele", "rank_le_0_5pct_on_min_nm_allele", 0.5),
            ("min_mhcflurry_rank_any_allele", "rank_le_0_5pct_any_patient_allele", 0.5),
        ):
            value = float(row[source])
            if not math.isfinite(value) or value < 0 or (source != "affinity_nm" and value > 100):
                raise ValueError("Invalid cryptic prediction value: " + source)
            block[key] += int(value <= cutoff)
        prime_rank = float(row["prime_rank"])
        if not math.isfinite(prime_rank) or not 0 <= prime_rank <= 100:
            raise ValueError("Invalid cryptic prediction value: prime_rank")
    return counts
PRIORITY_LABELS = {
    "TIER1_1MM_PARALOG_SAME_HLA_MANDATORY_ELISPOT": "TIER1_1MM_PARALOG_SAME_HLA",
    "TIER2_2MM_STRONG_SAME_HLA_MANDATORY_ELISPOT": "TIER2_2MM_STRONG_SAME_HLA",
    "TIER2_TITIN_CLASS_LEAD_MIMIC_MANDATORY_ELISPOT": "TIER2_TITIN_CLASS_LEAD_MIMIC",
}


def research_questions(class1_count, class2_count, *, target_root=HERE):
    """Describe open questions for the regions actually encoded by this panel."""
    selected = load(Path(target_root) / "selection/selected.json")
    genes = {row["gene"] for row in selected}
    kctd19_scope = (
        "KCTD19 is encoded and its coding-RNA expression is unconfirmed. "
        if "KCTD19" in genes else
        "KCTD19 is excluded from this panel; its coding-RNA question is not active for this construct. "
    )
    low_count_genes = [
        row["gene"] for row in selected
        if row["gene"] in {"FES", "SYT12"}
        and row.get("coding_strand_RNA_ALT_fragments") == 1
    ]
    low_count_scope = (
        " and ".join(sorted(low_count_genes))
        + (" each have one coding-strand mutant fragment. " if len(low_count_genes) > 1
           else " has one coding-strand mutant fragment. ")
        if low_count_genes else ""
    )
    current_target_scope = (
        kctd19_scope + low_count_scope
        + "Independent HLA confirmation, Class-II pairing and current tumour presentation remain unresolved."
    )
    topics = [
        ("SG-01", "Exact-reference coverage", "Recorded matches are zero for nominated queries retained after exact-reference exclusion.", "COMPUTATIONAL_RESULT", "This replays a selection condition, not an independent biological safety endpoint. Query selection and reference coverage are incomplete; exact absence does not exclude cross-recognition or authenticate named-patient identity."),
        ("SG-02", "Predicted junction binding", "The fixed optimization reaches the recorded MHCflurry, MixMHCpred, MixMHC2pred and PRIME threshold counts.", "COMPUTATIONAL_RESULT", "Model-defined counts and finite-search optima do not establish unwanted presentation, safety or efficacy."),
        ("SG-03", "Mutant/wild-type and off-target recognition", f"The proposed panels contain {class1_count} Class-I and {class2_count} Class-II control peptides; more distant homologs are a prioritized subset.", "UNTESTED", "Qualified specialists must define appropriate specificity assays, controls and detection limits. No universal EC50 or zero-harm acceptance criterion is established."),
        ("SG-04", "Current targets, RNA and HLA", current_target_scope, "UNRESOLVED", "DNA or cfDNA detection alone does not measure RNA expression or cell-surface HLA. No negative-test result or choice of clinical construct is assumed."),
        ("SG-05", "Complete product and translation", "Cap chemistry, nucleotide modifications, PolyA, delivery, dose and manufacturing are unspecified for this construct.", "UNSPECIFIED_AND_UNTESTED", "Platform descriptions do not establish this product's composition, translation fidelity or innate reactogenicity."),
        ("SG-06", "Manufacturing quality and toxicology", "No manufactured lot, validated product specifications, impurity assay results or product-specific toxicology are supplied.", "UNTESTED", "Applicable quality specifications and development studies require qualified product development; numerical release limits are not set by this analysis."),
        ("SG-07", "Clinical and regulatory responsibility", "Clinical investigation requires qualified oversight and the applicable ethical and regulatory authorization.", "NOT_DETERMINED", "The treating team and sponsor must establish the relevant jurisdiction and pathway. This analysis chooses no trial country, FDA route or authorization and grants no immunity from liability."),
    ]
    return [dict(zip(("question_id", "topic", "evidence_or_question", "status", "interpretation"), row)) for row in topics]


def atomic_write(path, content):
    temp = path.with_name(path.name + ".qualification-tmp")
    temp.write_bytes(content if isinstance(content, bytes) else content.encode())
    temp.replace(path)


def save_json(path, obj):
    atomic_write(path, json.dumps(obj, indent=2, sort_keys=True) + "\n")


def update_table(path, *, compressed=False, expression=None):
    resolved = resolve_path(path)
    is_gz = compressed or str(resolved).endswith(".gz")
    raw = gzip.decompress(resolved.read_bytes()).decode() if is_gz else resolved.read_text()
    reader = csv.DictReader(io.StringIO(raw), delimiter="\t")
    old_fields = reader.fieldnames
    rename = {"risk_classification": "model_priority_category", "required_wetlab_gate": "proposed_specificity_assessment"}
    fields = [rename.get(k, k) for k in old_fields]
    records = []
    for row in reader:
        row = {rename.get(k, k): PRIORITY_LABELS.get(v, v) for k, v in row.items()}
        if "proposed_specificity_assessment" in row:
            row["proposed_specificity_assessment"] = SPECIFICITY_INTERPRETATION
        if expression is not None:
            row.update(biopsy_expression_fields(row["homolog_genes"].split(";"), row["source_gene"], expression))
        records.append(row)
    fields = list(dict.fromkeys(fields + [key for row in records for key in row]))
    if expression is not None and "model_priority_category" in fields:
        # Match the fresh generator's schema: expression fields precede the priority label.
        fields = [key for key in fields if key != "model_priority_category"] + ["model_priority_category"]
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(records)
    data = output.getvalue().encode()
    atomic_write(resolved, gzip.compress(data, mtime=0) if is_gz else data)


def cryptic_interpretation_for(target_root=HERE):
    ver_path = Path(target_root) / "construct/mrna/construct.json"
    if ver_path.is_file():
        alt = json.loads(ver_path.read_text())["optimization_metrics"]["alternative_reading_frames"]
        stops = int(alt["1"]["stop_codons"]) + int(alt["2"]["stop_codons"])
        max_run = max(int(alt["1"]["longest_stop_free_run_codons"]), int(alt["2"]["longest_stop_free_run_codons"]))
        return (
            f"The final CDS contains {stops} alternative-frame stops and predicted stop-free runs up to "
            f"{max_run} codons. These sequence counts do not measure leaky scanning or frameshifting. "
            "Alternative translation requires assessment with defined assay sensitivity; no experiment proves absolute absence."
        )
    return CRYPTIC_INTERPRETATION


def refresh_cached_annotations(target_root=HERE):
    """Update interpretive labels only, preserving cached predictions and biological sequences."""
    target_root = Path(target_root)
    directory = target_root / "safety"
    receipt = json.loads((directory / "receipt.json").read_text())
    for name, expected in receipt["output_sha256"].items():
        actual = sha(directory / name)
        if actual != expected:
            raise ValueError("Refusing annotation refresh after an unaccounted artifact change: " + name)
    rna = json.loads((directory / "mrna-innate-and-cryptic-audit.json").read_text())
    for key in ("coding_rna_duplex_and_motifs", "known_transcript_segment_duplex_and_motifs"):
        metrics = rna[key]
        for old, new in (
            ("max_uninterrupted_dsrna_helix_bp", "max_uninterrupted_mfe_stem_bp"),
            ("max_1nt_bulge_tolerant_stem_bp", "max_bulge_tolerant_mfe_stem_bp"),
            ("below_PKR_30bp_duplex_threshold", "mfe_stems_below_reporting_cutoff"),
        ):
            if old in metrics:
                metrics[new] = metrics.pop(old)
        metrics.update(model_scope=MFE_SCOPE, reporting_cutoff_bp=MFE_REPORTING_CUTOFF_BP, biological_activation_threshold_validated=False)
    rna["innate_dsrna_interpretation"] = RNA_INTERPRETATION
    rna["cryptic_orf_summary"]["note"] = cryptic_interpretation_for(target_root)
    with gzip.open(directory / "cryptic-orf-query-coverage.tsv.gz", "rt", newline="") as stream:
        strand_counts = cryptic_counts_by_stored_strand(csv.DictReader(stream, delimiter="\t"))
    rna["cryptic_orf_summary"]["by_stored_strand"] = strand_counts
    rna["cryptic_orf_summary"]["stored_strand_scope"] = CRYPTIC_STRAND_SCOPE
    rna["mulroney_2024_frameshift_assessment"]["note"] = MOTIF_INTERPRETATION
    save_json(directory / "mrna-innate-and-cryptic-audit.json", rna)
    expression = read_biopsy_expression()
    for filename, compressed in (("offtarget-homologs.tsv.gz", True), ("wetlab-negative-control-peptides.tsv", False), ("wetlab-class2-wt-controls.tsv", False)):
        update_table(directory / filename, compressed=compressed,
                     expression=expression if filename == "offtarget-homologs.tsv.gz" else None)
    routing = json.loads((directory / "hla-routing-self-audit.json").read_text())
    for region in ("SEC", "MITD"):
        if "exact_patient_allele_match" in routing[region]:
            routing[region]["exact_supplied_allele_model_match"] = routing[region].pop("exact_patient_allele_match")
        routing[region]["immunological_note"] = "Exact sequence match to the supplied allele model and supporting reads; named-patient custody, immune tolerance, processing and clinical safety are not established."
    save_json(directory / "hla-routing-self-audit.json", routing)
    normal = json.loads((directory / "patient-normal-audit.json").read_text())
    normal["query_set_scope"] = "Nominated target queries retained after exact-reference filtering; zero matches replay that selection condition rather than independently validating clinical safety."
    normal["reference_limit"] = "Constructed coding/HLA reference with incomplete-CDS models and noncanonical translation omitted; not a complete normal-tissue presentation atlas."
    normal["named_patient_identity_authenticated"] = False
    save_json(directory / "patient-normal-audit.json", normal)
    summary = json.loads((directory / "safety-summary.json").read_text())
    summary["in_silico_safety_screen_status"] = "RECORDED_COMPUTATIONAL_SCREENS_ONLY"
    summary["why_clinical_safety_cannot_be_claimed_in_silico"] = SAFETY_LIMITS
    summary["clinical_safety_established"] = False
    summary["clinical_benefit_established"] = False
    summary.pop("mandatory_open_preclinical_and_clinical_gates", None)
    summary["unresolved_research_questions"] = ["SG-03", "SG-04", "SG-05", "SG-06", "SG-07"]
    summary["priority_category_scope"] = "Heuristic sequence/model prioritization; not measured toxicity or clinical testing requirements."
    summary["screening_policy"] = SCREENING_POLICY
    summary["query_set_scope"] = normal["query_set_scope"]
    summary["patient_biopsy_expression_scope"] = BIOPSY_EXPRESSION_SCOPE
    summary["cryptic_orf_stored_strand_scope"] = CRYPTIC_STRAND_SCOPE
    numeric = summary["in_silico_metrics"]
    numeric["cryptic_peptides_by_stored_strand"] = strand_counts
    for old, new in (
        ("max_intramolecular_uninterrupted_dsrna_bp", "known_RNA_mfe_max_uninterrupted_stem_bp"),
        ("max_intramolecular_1bulge_stem_bp", "known_RNA_mfe_max_bulge_tolerant_stem_bp"),
        ("autologous_HLA_routing_domains_verified", "routing_sequences_match_supplied_HLA_models"),
    ):
        if old in numeric:
            numeric[new] = numeric.pop(old)
    summary["control_panel_scope"]["scope_warning"] = "Panels cover the stated model and sequence rules, not all possible cross-reactive targets. Assay design and acceptance criteria are unvalidated; no negative biological result or absence of harm is established."
    save_json(directory / "safety-summary.json", summary)
    c1 = summary["control_panel_scope"]["stage1_class1_peptide_count"]
    c2 = summary["control_panel_scope"]["stage1_class2_wt_15mer_count"]
    questions = research_questions(c1, c2, target_root=target_root)
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=list(questions[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(questions)
    atomic_write(directory / "research-questions.tsv", output.getvalue())
    legacy = directory / "release-gates.tsv"
    if legacy.exists():
        legacy.unlink()
    receipt["output_sha256"].pop("release-gates.tsv", None)
    names = list(receipt["output_sha256"])
    if "research-questions.tsv" not in names:
        names.append("research-questions.tsv")
    receipt["output_sha256"] = {name: sha(directory / name) for name in names}
    receipt["clinical_safety_established"] = False
    receipt["scope"] = "Integrity of cached computational screens and bounded research annotations; not a clinical certification. Peptides and numerical model scores were not regenerated by this annotation refresh."
    save_json(directory / "receipt.json", receipt)


if __name__ == "__main__":
    refresh_cached_annotations(HERE)
    if (HERE / "strand-strict-15/safety/receipt.json").is_file():
        refresh_cached_annotations(HERE / "strand-strict-15")
