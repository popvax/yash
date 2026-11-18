import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pipeline"))

from curate import core_regions, derive, strand_strict_core_regions
from strand_audit import verify_offline as verify_strand_offline


class CoreReserveTests(unittest.TestCase):
    def test_weak_signal_mutation_and_unvalidated_self_sources_are_reserves(self):
        candidates = [
            {"source": "MUTATION", "gene": "SERPINA1", "passes_500nm_gate": False, "coding_strand_RNA_ALT_fragments": 470},
            {"source": "MUTATION", "gene": "H2BC12", "passes_500nm_gate": True, "coding_strand_RNA_ALT_fragments": 10},
            {"source": "MUTATION", "gene": "KCTD19", "passes_500nm_gate": True, "coding_strand_RNA_ALT_fragments": 0},
            {"source": "CTA/SELF", "gene": "PRAME"},
            {"source": "SPLICE", "gene": "POLE"},
            {"source": "FUSION", "gene": "CDH1--NRG1"},
            {"source": "FUSION", "gene": "CDH1::NRG1"},
        ]
        self.assertEqual(
            [row["gene"] for row in core_regions(candidates)], ["H2BC12", "KCTD19", "CDH1::NRG1"]
        )
        self.assertEqual(
            [row["gene"] for row in strand_strict_core_regions(candidates)], ["H2BC12", "CDH1::NRG1"]
        )

    def test_strand_audit_reconciles_kctd19_fes_and_syt12(self):
        check = verify_strand_offline()
        self.assertGreaterEqual(check["spliced_loci_examined"], 200)
        self.assertEqual(check["kctd19_coding_vs_opposite_alt"], [0, 2])
        self.assertEqual(check["fes_coding_vs_opposite_alt"], [1, 1])
        self.assertEqual(check["syt12_coding_vs_opposite_alt"], [1, 1])
        derived = {row["gene"]: row for row in core_regions(derive()[0])}
        self.assertEqual(derived["KCTD19"]["design_tier"], "somatic-dna-only-opposite-strand-rna-watch")
        self.assertEqual(derived["FES"]["design_tier"], "somatic-single-coding-fragment-low-depth")
        self.assertEqual(derived["SYT12"]["design_tier"], "somatic-single-coding-fragment-low-depth")


if __name__ == "__main__":
    unittest.main()
