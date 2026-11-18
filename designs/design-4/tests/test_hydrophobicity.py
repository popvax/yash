import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pipeline"))

from hydrophobicity import protein_screen


class ShortAntigenHydropathyTests(unittest.TestCase):
    def test_exact_15aa_antigen_is_not_padded_or_individually_passed(self):
        antigen = "EQCWAYEPGQWPSFS"
        protein = "G" * 19 + antigen + "G" * 19
        parts = [
            {"name": "linker-1", "kind": "linker", "start_aa": 0, "end_aa": 19, "aa": "G" * 19},
            {"name": "FES", "kind": "antigen", "start_aa": 19, "end_aa": 34, "aa": antigen},
            {"name": "linker-2", "kind": "linker", "start_aa": 34, "end_aa": 53, "aa": "G" * 19},
        ]

        windows, report = protein_screen(protein, parts)

        fes = next(row for row in report["antigen_regions"] if row["name"] == "FES")
        self.assertEqual(fes["length_aa"], 15)
        self.assertEqual(fes["status"], "NOT_SCORABLE_SHORTER_THAN_19_AA")
        self.assertEqual(fes["windows_scanned"], 0)
        self.assertEqual(len(windows), len(protein) - 18)
        self.assertFalse(report["all_antigen_regions_individually_pass"])
        self.assertEqual(report["policy"]["edge_padding"], False)
        self.assertEqual(report["status"], "PASS")

    def test_direct_join_is_a_valid_zero_length_linker(self):
        protein = "G" * 40
        parts = [
            {"name": "a", "kind": "antigen", "start_aa": 0, "end_aa": 20, "aa": "G" * 20},
            {"name": "direct", "kind": "linker", "start_aa": 20, "end_aa": 20, "aa": ""},
            {"name": "b", "kind": "antigen", "start_aa": 20, "end_aa": 40, "aa": "G" * 20},
        ]
        _, report = protein_screen(protein, parts)
        self.assertEqual(report["status"], "PASS")

    def test_individual_hydrophobic_antigen_cannot_hide_in_routing_overlap(self):
        protein = "G" * 19 + "I" * 19
        parts = [
            {"name": "SEC", "kind": "signal_peptide", "start_aa": 0, "end_aa": 19, "aa": "G" * 19},
            {"name": "a", "kind": "antigen", "start_aa": 19, "end_aa": 38, "aa": "I" * 19},
        ]
        _, report = protein_screen(protein, parts)
        self.assertEqual(report["status"], "FAIL")


if __name__ == "__main__":
    unittest.main()
