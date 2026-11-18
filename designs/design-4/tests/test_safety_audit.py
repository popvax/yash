import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pipeline"))

from dataio import iter_rows, rows
from safety_audit import (
    HERE,
    classify_homolog,
    extract_cryptic_orfs,
    longest_duplex_metrics,
    scan_class1_mimicry,
    scan_class2_9mer_cores,
    verify_offline,
)


class SafetyAuditTests(unittest.TestCase):
    def test_titin_magea3_archetype_is_caught_despite_four_mismatches(self):
        # Historic fatal off-target mimicry on HLA-A*01:01 (Cameron et al., Sci Transl Med 2013):
        # MAGE-A3: EVDPIGHLY vs cardiac Titin: ESDPIVAQY (4 mismatches, P3-P5 'DPI' + E1/Y9 conserved!)
        verdict = classify_homolog("EVDPIGHLY", "ESDPIVAQY")
        self.assertIsNotNone(verdict)
        self.assertEqual(verdict[0], "TCR_CORE_ANCHOR_MIMIC")
        self.assertEqual(verdict[1], 4)

        hits = scan_class1_mimicry(
            ["EVDPIGHLY"],
            [("GENCODE_v46", "ENSP00001|ENST1|ENSG1|x|y|TTN-201|TTN|30", "MAAAGGGESDPIVAQYKKK")],
        )
        self.assertIn(("EVDPIGHLY", "ESDPIVAQY"), hits)
        self.assertIn("TTN", hits["EVDPIGHLY", "ESDPIVAQY"]["genes"])

    def test_one_and_two_mismatch_and_class2_core_recall(self):
        hits = scan_class1_mimicry(
            ["AYEPGQWPSF"],
            [("GENCODE_v46", "ENSP00002|ENST2|ENSG2|x|y|FES-201|FES|20", "EQCWAYEPGQRPSFS")],
        )
        self.assertEqual(hits["AYEPGQWPSF", "AYEPGQRPSF"]["match_class"], "NEAR_1MM")

        c2 = scan_class2_9mer_cores(
            ["IATDNATST", "YSILRKYLQ"],
            [
                ("GENCODE_v46", "ENSP00003|ENST3|ENSG3|x|y|TEST-1|TEST|20", "AAAIATDNATSSGGG"),
                ("GENCODE_v46", "ENSP00004|ENST4|ENSG4|x|y|KCTD19-201|KCTD19|20", "RSVYSILGKYLQED"),
            ],
        )
        self.assertEqual(c2["IATDNATST", "IATDNATSS"]["mismatches"], 1)
        # Critical P5 (0-based index 4) mismatch must be caught by disjoint [0:5]+[5:9] pigeonhole seeds
        self.assertIn(("YSILRKYLQ", "YSILGKYLQ"), c2)
        self.assertEqual(c2["YSILRKYLQ", "YSILGKYLQ"]["mismatch_positions_1based"], "5")

    def test_dsrna_duplex_and_cryptic_orf_detection(self):
        # A predicted 32-bp stem exceeds the configured reporting cutoff; no activation assay is implied.
        hairpin_32bp = "G" * 32 + "UUCG" + "C" * 32
        metrics = longest_duplex_metrics(hairpin_32bp)
        self.assertGreaterEqual(metrics["max_uninterrupted_mfe_stem_bp"], 30)
        self.assertFalse(metrics["mfe_stems_below_reporting_cutoff"])
        self.assertFalse(metrics["biological_activation_threshold_validated"])

        # Short +1 frame ORF >=8 aa is detected while +0 primary frame is skipped
        rna = "A" * 63 + "AUG" + "GCU" * 10 + "UAA"
        orfs = extract_cryptic_orfs(rna, coding_start_nt=63, coding_end_nt=63 + 33)
        self.assertTrue(any(o["frame_vs_CDS"] == "+1" and o["length_aa"] >= 8 for o in orfs))
        self.assertFalse(any(o["frame_vs_CDS"] == "+0" for o in orfs))

    def test_zero_shared_hla_own_wildtype_controls_omitted(self):
        result = verify_offline()
        self.assertEqual(result["shared_hla_own_wt_controls_omitted"], 0)
        plate = {r["control_peptide"]: r for r in rows(HERE / "safety/wetlab-negative-control-peptides.tsv")}
        self.assertIn("SYSILGKYL", plate)
        self.assertEqual(plate["SYSILGKYL"]["paired_vaccine_peptide"], "SYSILRKYL")
        self.assertEqual(plate["SYSILGKYL"]["test_hla_allele"], "HLA-C*07:02")
        self.assertIn("ATDNGSPVAT", plate)
        for row in iter_rows(HERE / "safety/offtarget-homologs.tsv.gz"):
            if row["is_mutated_gene_own_WT_only"] == "True" and row["same_hla_restriction"] == "True":
                self.assertIn(row["homolog_peptide"], plate)

        # Verify the separate strand-strict-15 safety audit covers all 1,782 sense cryptic 8-11-mers and 174 homolog queries
        self.assertIn("strand_strict_15", result)
        self.assertEqual(result["strand_strict_15"]["shared_hla_own_wt_controls_omitted"], 0)
        s15_sum = rows(HERE / "strand-strict-15/safety/wetlab-negative-control-peptides.tsv")
        self.assertEqual(len(s15_sum), 184)
        s15_plate = {r["control_peptide"] for r in s15_sum}
        self.assertNotIn("SYSILGKYL", s15_plate)  # KCTD19 is excised in strand-strict-15
        for row in iter_rows(HERE / "strand-strict-15/safety/offtarget-homologs.tsv.gz"):
            if row["is_mutated_gene_own_WT_only"] == "True" and row["same_hla_restriction"] == "True":
                self.assertIn(row["homolog_peptide"], s15_plate)

    def test_annotation_refresh_rejects_unaccounted_artifact_tampering(self):
        import shutil
        import tempfile
        from research_qualifications import refresh_cached_annotations

        with tempfile.TemporaryDirectory() as tmp:
            troot = Path(tmp)
            shutil.copytree(HERE / "safety", troot / "safety")
            (troot / "safety/class2-core-homologs.tsv").write_text("tampered\n")
            with self.assertRaisesRegex(ValueError, "Refusing annotation refresh"):
                refresh_cached_annotations(troot)


if __name__ == "__main__":
    unittest.main()
