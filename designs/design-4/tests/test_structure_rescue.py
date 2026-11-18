import sys
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "construct/mrna"))
from structure_rescue import rescue


class FrozenStemRescueTests(unittest.TestCase):
    def setUp(self):
        self.codons = ["AUG", "GCC", "GGU", "CUG", "GUG", "AAG", "CAG",
                       "GAC", "CCU", "UUC", "AGC", "CGC", "ACC"]
        self.protein = "MAGLVKQDPFSRT"
        self.coding = "".join(self.codons)
        self.allowed = {aa: [codon] for aa, codon in zip(self.protein, self.codons)}
        self.allowed["D"] = ["GAC", "GAU"]
        self.frequency = {aa: {codon: 1.0 for codon in codons}
                          for aa, codons in self.allowed.items()}
        self.bad = "(" * 15 + "." * 6 + ")" * 8 + "." * 3 + ")" * 7
        self.good = "(" * 14 + "." * 7 + ")" * 7 + "." * 4 + ")" * 7

    def test_partner_codon_can_repair_a_fully_frozen_long_bracket_run(self):
        def fold(sequence):
            return (self.good if sequence[21:24] == "GAU" else self.bad, -1.0)
        with patch("structure_rescue.RNA.fold", side_effect=fold):
            result, audit = rescue(
                self.coding, self.protein, self.allowed, self.frequency, [], set(range(18))
            )
        self.assertEqual(result[:18], self.coding[:18])
        self.assertEqual(result[21:24], "GAU")
        self.assertTrue(audit["target_reached"])

    def test_impossible_frozen_target_is_reported_not_fabricated(self):
        with patch("structure_rescue.RNA.fold", return_value=(self.bad, -1.0)):
            result, audit = rescue(
                self.coding, self.protein, self.allowed, self.frequency, [],
                set(range(len(self.coding))),
            )
        self.assertEqual(result, self.coding)
        self.assertFalse(audit["target_reached"])


if __name__ == "__main__":
    unittest.main()
