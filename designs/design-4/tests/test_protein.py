import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pipeline"))

from dataio import load, rows
from protein import (
    HERE,
    assembled_class2_context,
    assess_assembled_class2,
    eligible_options,
    solve,
    target_candidates,
)


class EventSpecificityTests(unittest.TestCase):
    def test_inframe_fusion_self_flank_is_not_counted_as_a_target(self):
        sequence = "STYTALIIATDNATSTSTTGTSHLV"
        region = {
            "source": "FUSION", "sequence": sequence,
            "target_start_in_window": 12, "target_end_in_window": 13,
        }
        reference = {
            sequence[start:start + length]: {"occurs": False}
            for length in (8, 9, 10, 11)
            for start in range(len(sequence) - length + 1)
        }
        reference["NATSTSTT"]["occurs"] = True
        targets = {peptide for _, peptide in target_candidates(region, reference)}
        self.assertNotIn("TSTTGTSHL", targets)
        self.assertNotIn("NATSTSTT", targets)
        self.assertIn("ALIIATDNA", targets)


class TwoNeighborProcessingTests(unittest.TestCase):
    def test_pruning_respects_model_agreement_before_mild_binding_count(self):
        common = {"left_index": 0, "right_index": 1, "strong_rank": 0,
                  "hydrophobicity_pass": True, "prime_low_rank": 0}
        options = [
            {**common, "edge_id": "primary-only", "mild_rank": 0, "mix_strong_rank": 1},
            {**common, "edge_id": "agreement", "mild_rank": 1, "mix_strong_rank": 0},
        ]
        self.assertEqual([row["edge_id"] for row in eligible_options(options)], ["agreement"])

    def test_two_flank_term_is_counted_once_for_the_selected_neighbor_pair(self):
        edges = [
            {
                "edge_id": f"{i}:{j}", "left_index": i, "right_index": j,
                "strong_rank": 0, "mild_rank": 0, "min_rank": 3.0,
                "added_aa": 0, "intended_proc_score_x10000": 0,
                "hydrophobicity_pass": True,
            }
            for i in range(3) for j in range(3) if i != j
        ]
        triplets = [
            {"incoming_edge_id": "0:1", "outgoing_edge_id": "1:2",
             "processing_score_x10000": 8000},
            {"incoming_edge_id": "2:1", "outgoing_edge_id": "1:0",
             "processing_score_x10000": 1000},
        ]
        path, result = solve(2, edges, triplets, seconds=5)
        self.assertEqual([edge["edge_id"] for edge in path], ["2:0", "0:1", "1:2"])
        self.assertEqual(result["phases"][-1]["value"], 8000)

    def test_fallback_15_has_complete_assembled_class2_context_coverage(self):
        strict_dir = HERE / "strand-strict-15"
        seq15 = "".join(
            line for line in (strict_dir / "construct/protein/cassette.fasta").read_text().splitlines()
            if not line.startswith(">")
        )
        joins15 = rows(strict_dir / "construct/protein/junctions.tsv")
        c2_ctx = {(r["Peptide"], r["Context"]): r for r in rows(strict_dir / "immunogenicity/mixmhc2-context.tsv")}
        n2 = {r["Peptide"]: r for r in rows(strict_dir / "immunogenicity/mixmhc2-no-context.tsv")}
        self.assertIn(("GEQCWAYEPGQWPSF", "GGGGEQPSFSGI"), c2_ctx)
        for st in range(len(seq15) - 14):
            pep = seq15[st:st + 15]
            ctx = assembled_class2_context(seq15, st, st + 15)
            self.assertIn((pep, ctx), c2_ctx)
        _, summary = assess_assembled_class2(seq15, joins15, c2_ctx, n2)
        recorded = load(strict_dir / "immunogenicity/junction-multi-model-summary.json")
        self.assertEqual(summary["ClassII_synthetic_9aa_core_rank_le_0_5pct"], 0)
        for key, val in summary.items():
            self.assertEqual(recorded[key], val)


if __name__ == "__main__":
    unittest.main()
