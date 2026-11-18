import json
import sys
import unittest
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "pipeline"), str(ROOT / "tools")]

from alphagenome2_audit import verify_offline as verify_ag2
from audit_ctdna_loci import verify_offline as verify_ctdna_loci
from build_strand_strict_contingency import cached_junction_check_record
from dataio import iter_rows, load, rows, sha
from replay_helpers import canonical_path, read_edges, read_triplets


class AlphaGenomeAndStrictMrnaTests(unittest.TestCase):
    def test_zero_percentile_class2_core_is_not_dropped(self):
        from class2_construct import eligible_core_rank
        self.assertTrue(eligible_core_rank(0.0))
        self.assertTrue(eligible_core_rank(2.0))
        self.assertFalse(eligible_core_rank(None))
        self.assertFalse(eligible_core_rank(2.001))
        self.assertFalse(eligible_core_rank(-0.1))

    def test_class2_annotations_use_each_construct_context(self):
        from class2_construct import verify_class2_annotations
        verify_class2_annotations(ROOT)
        verify_class2_annotations(ROOT / "strand-strict-15")
        controls = rows(ROOT / "strand-strict-15/safety/wetlab-class2-wt-controls.tsv")
        got = {r["mutant_15mer"]: float(r["best_context_rank_pct"]) for r in controls}
        self.assertEqual(got["VSKGTKAVTKYTSAK"], 0.881)
        self.assertEqual(got["YGEMLHIRYRLLRQA"], 1.80)

    def test_fragment_summary_uses_active_selection(self):
        for base in (ROOT, ROOT / "strand-strict-15"):
            s = load(base / "selection/summary.json")
            self.assertIn(str(s["source_counts"]["MUTATION"])+"-mutation", s["fragment_sum_interpretation"])
            self.assertIn(str(s["mutation_RNA_ALT_fragment_count_sum"])+" pooled", s["fragment_sum_interpretation"])

    def test_recorder_does_not_claim_raw_identity_execution(self):
        from record_verification import _build_receipt_payload
        r = _build_receipt_payload([], 0, {}, {})
        self.assertFalse(r["data_access_limits"]["whole_RNA_BAM_checksum_checked_in_this_workflow"])
        self.assertFalse(r["identity_check"]["portable_checker_reexecuted_in_this_workflow"])
        self.assertIn("Historical", r["identity_check"]["evidence_status"])

    def test_alphagenome2_offline_receipt_and_observed_only_per_track_retention(self):
        summary = verify_ag2()
        self.assertEqual(summary["conditions_with_imputation_metadata"], 12)
        self.assertEqual(summary["ctdna_longitudinal_variants_audited"], 13)
        self.assertEqual(summary["metastatic_and_tme_genes_profiled"], 32)
        self.assertGreater(summary["retained_variant_track_rows"], 1000)
        self.assertGreater(summary["retained_interval_track_rows"], 3000)
        self.assertGreater(summary["embedding_cosine_similarity"], 0.99)
        conds = {r["label"]: r for r in rows(ROOT / "alphagenome2/condition-imputation-metadata.tsv")}
        for organ in ("testis", "heart", "lung", "liver", "kidney_cortex", "brain", "pancreas_tissue"):
            self.assertEqual(conds[organ]["rna_seq"], "observed")
            self.assertIn("not_patient_measurement", conds[organ]["metadata_semantics"])
        for cell in ("pancreatic_ductal_cell", "liver_hepatocyte", "peritoneal_mesothelial"):
            self.assertEqual(conds[cell]["sc_gene_expression"], "observed")
        # Ensure single-cell conditions never run imputed_exploratory splicing scorers and units stay separate
        ctdna = rows(ROOT / "alphagenome2/ctdna-progression-audit.tsv")
        self.assertIn("lung_obs_splice_sites_max_prob_diff", ctdna[0])
        self.assertIn("lung_obs_splice_usage_max_frac_diff", ctdna[0])
        self.assertIn("lung_obs_splice_junctions_max_log_diff", ctdna[0])
        self.assertNotIn("ductal_max_splice_delta", ctdna[0])
        for tr in iter_rows(ROOT / "alphagenome2/variant-track-scores.tsv.gz"):
            if tr["condition_label"] in ("ductal", "peritoneum"):
                self.assertNotIn("SPLICE", tr["output_type"])

    def test_strand_strict_15_mrna_fresh_recheck_and_honest_cached_mode(self):
        prot = str(SeqIO.read(ROOT / "strand-strict-15/construct/protein/cassette.fasta", "fasta").seq)
        cds = str(SeqIO.read(ROOT / "strand-strict-15/construct/mrna/cds.fasta", "fasta").seq)
        self.assertEqual(len(prot), 481)
        self.assertEqual(len(cds), 1446)
        self.assertEqual(str(Seq(cds[:-3]).translate()), prot)
        ver = load(ROOT / "strand-strict-15/construct/mrna/verification.json")
        self.assertEqual(ver["status"], "KNOWN_NUCLEOTIDE_SEGMENTS_VERIFIED_POLYA_UNSPECIFIED")
        self.assertEqual(ver["upstream_junction_check_status"], "FRESH_FULL_JUNCTION_CHECK_PASSED")
        self.assertFalse(ver["complete_construct_verified"])
        self.assertEqual(ver["repeated_12nt_pairs"], 0)
        self.assertLessEqual(ver["final_paired_bracket_run_nt"], 14)

        # Verify cached assembly mode never fabricates FRESH_FULL_JUNCTION_CHECK_PASSED or 0.0/True fresh flags
        selected_all = load(ROOT / "selection/selected.json")
        strict_selected = load(ROOT / "strand-strict-15/selection/selected.json")
        kctd_idx = next(i for i, r in enumerate(selected_all) if r["gene"] == "KCTD19")
        def remap(idx):
            if idx == kctd_idx:
                return None
            if idx == len(selected_all):
                return len(strict_selected)
            return idx if idx < kctd_idx else idx - 1
        edges15 = [
            dict(e, left_index=remap(e["left_index"]), right_index=remap(e["right_index"]))
            for e in read_edges(ROOT / "construct/protein/edge-costs.tsv")
            if remap(e["left_index"]) is not None and remap(e["right_index"]) is not None
        ]
        kept_ids = {e["edge_id"] for e in edges15}
        triplets15 = [
            dict(t, region_index=remap(t["region_index"]))
            for t in read_triplets(ROOT / "construct/protein/processing-triplets.tsv")
            if remap(t["region_index"]) is not None
            and t["incoming_edge_id"] in kept_ids and t["outgoing_edge_id"] in kept_ids
        ]
        rec15 = load(ROOT / "strand-strict-15/contingency-receipt.json")
        path15 = canonical_path(len(strict_selected), edges15, rec15["chosen_edge_ids"])
        alleles = load(ROOT / "exports/source-case.json")["hla_class_I"]
        values = {
            (r["peptide"], r["allele"]): (float(r["affinity_nm"]), float(r["affinity_percentile"]))
            for r in iter_rows(ROOT / "construct/all-junction-affinity.tsv.gz")
        }
        cached_rec = cached_junction_check_record(ROOT / "strand-strict-15", strict_selected, path15, triplets15, alleles, values)
        self.assertEqual(cached_rec["status"], "CACHED_ASSEMBLY_JUNCTION_CHECK_PASSED")
        self.assertFalse(cached_rec["fresh_repeat_model_inference_performed"])
        self.assertIsNone(cached_rec["max_relative_affinity_difference"])
        self.assertIsNone(cached_rec["threshold_classifications_identical"])
        self.assertIsNone(cached_rec["processing_objective_matches_full_cassette"])

        # Verify strand-strict-15 safety audit covers all 1,782 sense cryptic 8-11-mers & 174 homolog queries
        s15_mrna_audit = load(ROOT / "strand-strict-15/safety/mrna-innate-and-cryptic-audit.json")
        self.assertEqual(s15_mrna_audit["cryptic_orf_summary"]["sense_cryptic_8to11mers_screened"], 1782)
        self.assertIn("sense_cryptic_peptides_rank_le_0_5pct_on_min_nm_allele", s15_mrna_audit["cryptic_orf_summary"])
        self.assertIn("sense_cryptic_peptides_rank_le_0_5pct_any_patient_allele", s15_mrna_audit["cryptic_orf_summary"])
        s15_safety = load(ROOT / "strand-strict-15/safety/safety-summary.json")["in_silico_metrics"]
        self.assertEqual(s15_safety["eligible_class1_homolog_queries_screened"], 174)
        self.assertIn("sense_cryptic_orf_rank_le_0_5pct_on_min_nm_allele", s15_safety)
        self.assertIn("sense_cryptic_orf_rank_le_0_5pct_any_patient_allele", s15_safety)

    def test_cigar_splice_skip_is_never_mistaken_for_indel_deletion(self):
        from audit_ctdna_loci import _classify_read_allele

        class _FakeSpliceRead:
            query_sequence = "CC"
            query_qualities = [37, 37]
            reference_start = 25247120
            reference_end = 25247126
            # 1M (pos 25247120) + 4N intron skip (25247121..25247124) + 1M (pos 25247125)
            cigartuples = [(0, 1), (3, 4), (0, 1)]

            def get_aligned_pairs(self, matches_only=False):
                return [(0, 25247120), (1, 25247125)]

        class _FakeTrueDelRead(_FakeSpliceRead):
            # 1M + 4D true deletion (op 2) + 1M
            cigartuples = [(0, 1), (2, 4), (0, 1)]

        span = (25247121, 25247126, "CAAAAC", "CC")
        self.assertIsNone(_classify_read_allele(_FakeSpliceRead(), 25247120, "CAAAA", "C", indel_span_1based=span))
        self.assertEqual(_classify_read_allele(_FakeTrueDelRead(), 25247120, "CAAAA", "C", indel_span_1based=span), "ALT")

    def test_upstream_source_ctdna_pileups_and_bounded_clinical_questions(self):
        sv = load(ROOT / "review/source-verification.json")
        self.assertEqual(sv["checked_upstream_commit"], "7ece8f7a77c8db935422f6502d2c4db2352d7bf0")
        self.assertEqual(sv["source_receipt_sha256"], sha(ROOT / "exports/source-receipt.json"))
        self.assertEqual(sv["repository_files_matching"], 20)

        ctdna_check = verify_ctdna_loci()
        self.assertEqual(ctdna_check["ctdna_loci_pileups_verified"], 13)
        self.assertEqual(ctdna_check["ctdna_mhcflurry_rows_verified"], 42)

        clin = load(ROOT / "review/clinical-reconciliation.json")
        clin_blob = json.dumps(clin)
        for banned in ("truncal", "subclone", "splice-neutral", "testis-restricted", "OCR", "did not detect", "unperformed", "all negative", "15q14", "amplification (+)"):
            self.assertNotIn(banned, clin_blob)
        self.assertNotIn("NM_004360.5", clin["clinically_reported_fusion"]["caris_report_annotation"])
        self.assertIn("NM_004360.5", clin["clinically_reported_fusion"]["canonical_refseq_transcript_equivalence"])
        self.assertIn("clinical_pdf_sources", clin)
        self.assertIn("project_and_narrative_sources", clin)

        priorities = rows(ROOT / "review/missing-data-priorities.tsv")
        self.assertEqual([r["priority_id"] for r in priorities], [f"MD-0{i}" for i in range(1, 7)])
        self.assertIn("open_question_or_evidence_gap", priorities[0])
        self.assertIn("15q21.1", priorities[1]["open_question_or_evidence_gap"])
        for row in priorities:
            self.assertNotIn("200-300x", json.dumps(row))
            self.assertNotIn("15q14", json.dumps(row))

    def test_verification_receipt_separates_current_and_historical_runs(self):
        ver = load(ROOT / "review/verification.json")
        self.assertIn("current_execution", ver)
        self.assertIn("historical_upstream_execution_reference", ver)
        self.assertIn(ver["overall_status"], ("IN_PROGRESS", "COMPLETED"))
        cur = ver["current_execution"]
        for chk in cur["checks"]:
            self.assertIn(chk["status"], ("IN_PROGRESS", "COMPLETED"))
            if chk["status"] == "COMPLETED":
                self.assertEqual(chk["exit_code"], 0)
            else:
                self.assertEqual(chk["check"], "final-tests")
                self.assertIsNone(chk["exit_code"])
            for arg in chk["args"]:
                self.assertNotIn("/dev/shm", arg)
            log_path = ROOT / chk["log_file"]
            self.assertTrue(log_path.is_file())
            self.assertEqual(sha(log_path), chk["log_sha256"])
        for rel, exp_hash in ver["executed_scientific_inputs_sha256"].items():
            self.assertEqual(sha(ROOT / rel), exp_hash, f"Input hash drift: {rel}")


if __name__ == "__main__":
    unittest.main()
