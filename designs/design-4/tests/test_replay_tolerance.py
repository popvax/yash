import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pipeline"))

from replay_helpers import compare_json_with_tolerance, compare_optimizer_log


class ReplayToleranceTests(unittest.TestCase):
    def test_sub_1e15_float_drift_passes_while_discrete_or_sequence_drift_fails(self):
        expected = {
            "sequence_sha256": {"rna_before_PolyA": "abc123"},
            "lengths": {"sense_coding_nt": 1521},
            "all_constraints_pass": True,
            "initiation_accessibility": {
                "mean_marginal_unpaired_probability": 0.36083389524675435,
                "per_base": [0.5210000000000001, 0.6123456789012345],
            },
        }
        observed = {
            "sequence_sha256": {"rna_before_PolyA": "abc123"},
            "lengths": {"sense_coding_nt": 1521},
            "all_constraints_pass": True,
            "initiation_accessibility": {
                "mean_marginal_unpaired_probability": 0.3608338952467546,
                "per_base": [0.5210000000000004, 0.6123456789012342],
            },
        }
        diff = compare_json_with_tolerance(observed, expected, rel_tol=1e-12, abs_tol=1e-12)
        self.assertGreater(diff, 0.0)
        self.assertLess(diff, 1e-14)

        # Any sequence hash, integer, or boolean change must fail immediately
        bad_seq = {**observed, "sequence_sha256": {"rna_before_PolyA": "abc124"}}
        with self.assertRaises(ValueError):
            compare_json_with_tolerance(bad_seq, expected)

        bad_int = {**observed, "lengths": {"sense_coding_nt": 1524}}
        with self.assertRaises(ValueError):
            compare_json_with_tolerance(bad_int, expected)

        bad_bool = {**observed, "all_constraints_pass": 1}
        with self.assertRaises(ValueError):
            compare_json_with_tolerance(bad_bool, expected)

        # Reverse bool/int substitution (True in place of integer 1) must also fail
        with self.assertRaises(ValueError):
            compare_json_with_tolerance({"lengths": {"sense_coding_nt": True}}, {"lengths": {"sense_coding_nt": 1}})
        with self.assertRaises(ValueError):
            compare_json_with_tolerance({"lengths": {"count": False}}, {"lengths": {"count": 0}})

        # Fractional or exact float in place of integer (and integer in place of float) must fail
        with self.assertRaises(ValueError):
            compare_json_with_tolerance({"count": 1.0}, {"count": 1})
        with self.assertRaises(ValueError):
            compare_json_with_tolerance({"count": 1.0 + 1e-16}, {"count": 1})
        with self.assertRaises(ValueError):
            compare_json_with_tolerance({"prob": 1}, {"prob": 1.0})

        # Float drift exceeding tolerance must fail
        bad_float = {
            **observed,
            "initiation_accessibility": {
                **observed["initiation_accessibility"],
                "mean_marginal_unpaired_probability": 0.3608349,
            },
        }
        with self.assertRaises(ValueError):
            compare_json_with_tolerance(bad_float, expected, rel_tol=1e-12, abs_tol=1e-12)

    def test_optimizer_log_comparison(self):
        exp = 'Coding DNA generated: 1521 nt\n{"codon_distribution_error": 0.20109512462109444, "seed": 0}\nEND\n'
        obs = 'Coding DNA generated: 1521 nt\n{"codon_distribution_error": 0.20109512462109466, "seed": 0}\nEND\n'
        self.assertLess(compare_optimizer_log(obs, exp), 1e-14)
        with self.assertRaises(ValueError):
            compare_optimizer_log(obs.replace("1521", "1524"), exp)


if __name__ == "__main__":
    unittest.main()
