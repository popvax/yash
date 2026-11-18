"""Typed cached scores, canonical path checks, and numerical-tolerance JSON comparison."""
import json
import math

from dataio import rows
from protein import eligible_options


def compare_json_with_tolerance(observed, expected, rel_tol=1e-12, abs_tol=1e-12, path="$"):
    """Verify exact Python type and value equality on discrete fields and bounded tolerance on floats."""
    if type(observed) is not type(expected):
        raise ValueError(
            f"Type mismatch at {path}: {type(observed).__name__} ({observed!r}) != "
            f"{type(expected).__name__} ({expected!r})"
        )
    if type(expected) is float:
        if not (math.isfinite(observed) and math.isfinite(expected)):
            raise ValueError(f"Non-finite float at {path}: {observed!r} vs {expected!r}")
        if not math.isclose(observed, expected, rel_tol=rel_tol, abs_tol=abs_tol):
            raise ValueError(
                f"Float diagnostic exceeded tolerance (rel={rel_tol}, abs={abs_tol}) at {path}: "
                f"{observed!r} vs {expected!r}"
            )
        return abs(observed - expected)
    if type(expected) is dict:
        if set(observed) != set(expected):
            raise ValueError(f"JSON key mismatch at {path}")
        return max(
            (compare_json_with_tolerance(observed[k], expected[k], rel_tol, abs_tol, f"{path}.{k}")
             for k in expected),
            default=0.0,
        )
    if type(expected) is list:
        if len(observed) != len(expected):
            raise ValueError(f"JSON list length mismatch at {path}")
        return max(
            (compare_json_with_tolerance(o, e, rel_tol, abs_tol, f"{path}[{i}]")
             for i, (o, e) in enumerate(zip(observed, expected))),
            default=0.0,
        )
    if observed != expected:
        raise ValueError(f"Discrete value mismatch at {path}: {observed!r} != {expected!r}")
    return 0.0


def compare_optimizer_log(observed_text, expected_text, rel_tol=1e-12, abs_tol=1e-12):
    """Verify exact discrete optimizer trace lines and float tolerance on the JSON summary line."""
    obs_lines = observed_text.splitlines()
    exp_lines = expected_text.splitlines()
    if len(obs_lines) != len(exp_lines):
        raise ValueError("Optimizer log line count changed")
    max_diff = 0.0
    for idx, (o_line, e_line) in enumerate(zip(obs_lines, exp_lines)):
        if o_line.startswith("{") and e_line.startswith("{"):
            diff = compare_json_with_tolerance(
                json.loads(o_line), json.loads(e_line), rel_tol, abs_tol, f"optimizer.log[{idx}]"
            )
            max_diff = max(max_diff, diff)
        elif o_line != e_line:
            raise ValueError(f"Optimizer discrete log line {idx} changed: {o_line!r} != {e_line!r}")
    return max_diff


def read_edges(path):
    records = rows(path)
    integers = (
        "left_index", "right_index", "strong_rank", "mild_rank", "affinity_50",
        "affinity_500", "added_aa", "peptide_occurrences", "intended_proc_score_x10000",
        "mix_strong_rank", "prime_low_rank", "mixmhc2_strong_core",
        "mixmhc2_mild_core", "mixmhc2_all_mild",
    )
    for row in records:
        for key in integers:
            if key in row:
                row[key] = int(row[key])
        for key in ("min_rank", "min_affinity_nm", "max_pair_mean_kd"):
            row[key] = float(row[key])
        row["hydrophobicity_pass"] = row["hydrophobicity_pass"] == "True"
    return records


def read_triplets(path):
    records = rows(path)
    for row in records:
        row["region_index"] = int(row["region_index"])
        row["processing_score_x10000"] = int(row["processing_score_x10000"])
    return records


def canonical_path(n, edges, identifiers):
    chosen = [edge for edge in edges if edge["edge_id"] in set(identifiers)]
    if len(chosen) != n + 1:
        raise ValueError("Wrong number of canonical edges")
    outgoing = {edge["left_index"]: edge for edge in chosen}
    path, node = [], n
    for _ in range(n + 1):
        edge = outgoing[node]
        path.append(edge)
        node = edge["right_index"]
    if node != n or sorted(edge["right_index"] for edge in path[:-1]) != list(range(n)):
        raise ValueError("Canonical path does not visit every antigen exactly once")
    return path


def objective_values(edges, path, triplets):
    ranks = sorted({edge["min_rank"] for edge in eligible_options(edges)})
    chosen = {edge["edge_id"] for edge in path}
    processing = sum(edge["intended_proc_score_x10000"] for edge in path)
    processing += sum(
        row["processing_score_x10000"] for row in triplets
        if row["incoming_edge_id"] in chosen and row["outgoing_edge_id"] in chosen
    )
    values = {
        "strong_count": sum(edge["strong_rank"] for edge in path),
        "mild_count": sum(edge["mild_rank"] for edge in path),
        "worst_rank_ordinal": ranks.index(min(edge["min_rank"] for edge in path)),
        "added_residues": sum(edge["added_aa"] for edge in path),
        "intended_boundary_processing_x10000": processing,
    }
    for key, name in (("mix_strong_rank", "mixmhcpred_strong_count"),
                      ("mixmhc2_strong_core", "mixmhc2_strong_core_count"),
                      ("prime_low_rank", "prime_low_rank_count"),
                      ("affinity_500", "affinity_500_count")):
        if all(key in edge for edge in path):
            values[name] = sum(edge[key] for edge in path)
    if all("mixmhc2_mild_core" in edge and "affinity_500" in edge for edge in path):
        values["total_mild_and_500nm_burden"] = sum(
            edge["mild_rank"] + edge["affinity_500"] + edge["mixmhc2_mild_core"] for edge in path
        )
    return values
