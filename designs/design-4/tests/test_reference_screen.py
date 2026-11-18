import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pipeline"))

from reference_screen import exact_search


class ReferenceSearchTests(unittest.TestCase):
    def test_overlapping_hits_and_final_complete_window_are_counted(self):
        wanted = {"AAAAAAAA", "AAAAAAAAA", "CDEFGHIK"}
        found, _, _ = exact_search([
            ("repeat", "AAAAAAAAA"), ("terminal", "LCDEFGHIK"),
        ], wanted)
        self.assertEqual(found["AAAAAAAA"]["occurrence_count"], 2)
        self.assertEqual(found["AAAAAAAAA"]["occurrence_count"], 1)
        self.assertEqual(found["CDEFGHIK"]["first_hits"][0]["start_aa"], 1)

    def test_absent_peptide_is_reported_without_an_invented_hit(self):
        found, _, _ = exact_search([("reference", "ACDEFGHIKLMN")], {"YYYYYYYY"})
        self.assertFalse(found["YYYYYYYY"]["occurs"])
        self.assertEqual(found["YYYYYYYY"]["first_hits"], [])


if __name__ == "__main__":
    unittest.main()
