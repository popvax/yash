"""Replay fictional claims and optionally check the existing expression source."""

import argparse
import hashlib
import json
import platform
import subprocess
import sys
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

import source_check


HERE = Path(__file__).resolve().parent
MATCHES = {"value_matches", "value_matches_after_rounding"}


def execute(command):
    exit_code = None
    stdout = ""
    stderr = ""
    try:
        completed = subprocess.run(
            command, cwd=HERE, capture_output=True, text=True, timeout=30
        )
        exit_code = completed.returncode
        stdout = completed.stdout
        stderr = completed.stderr
        report = json.loads(stdout)
        if not isinstance(report, dict):
            raise ValueError("Checker output must be a JSON object.")
    except (OSError, subprocess.TimeoutExpired, ValueError) as error:
        report = {"status": "execution_error", "reason": str(error)}
    return exit_code, report, stdout, stderr


def source_checks():
    identity = json.loads((HERE / "source-identity.json").read_text(encoding="utf-8"))
    path = source_check.ROOT / identity["source"]["path"]
    content = path.read_bytes()
    if hashlib.sha256(content).hexdigest() != identity["source"]["sha256"]:
        raise ValueError("Existing expression source bytes differ from the pinned version.")
    rows = [line for line in content.decode("utf-8").split("\n") if line and not line.startswith("#")]
    gene, raw_value = rows[0].removesuffix("\r").split("\t")
    sample = identity["source"]["sample_id"]
    rounded = str(Decimal(raw_value).quantize(Decimal("0.01"), rounding=ROUND_HALF_UP))
    cases = [
        ("existing-exact", gene, sample, raw_value, [], "value_matches"),
        ("existing-rounded", gene, sample, rounded, ["--decimal-places", "2"], "value_matches_after_rounding"),
        ("existing-wrong-number", gene, sample, str(Decimal(raw_value) + 1), [], "value_differs"),
        ("existing-wrong-sample", gene, "FICTIONAL-OTHER-SAMPLE", raw_value, [], "gene_or_sample_differs"),
        ("existing-unknown-scale", gene, sample, raw_value, ["--scale", "unknown"], "insufficient_information"),
        ("documented-missing-lookup", "H2BC12", sample, "6", [], "measurement_missing_from_source"),
    ]
    results = []
    for case_id, gene, sample, value, extra, expected in cases:
        command = [
            sys.executable, str(HERE / "source_check.py"),
            "--gene", gene, "--sample", sample, "--value", value,
        ] + extra
        exit_code, report, _, _ = execute(command)
        expected_exit = 0 if expected in MATCHES else 1
        results.append({
            "id": case_id,
            "expected_status": expected,
            "actual_status": report.get("status"),
            "exit_code": exit_code,
            "matches_expected": report.get("status") == expected and exit_code == expected_exit,
        })
    return {
        "scope": "Read-only integration. Claims are constructed from the pinned input; independent scientific labels are pending.",
        "expression_sha256": identity["source"]["sha256"],
        "reference_audit_bytes_verified": False,
        "output_scope": "Statuses only. Source measurements and sample values are omitted from this receipt.",
        "results": results,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="A fresh result file in an existing directory.")
    parser.add_argument("--with-source", action="store_true", help="Also read the existing pinned repository expression table.")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Choose a new result filename to preserve earlier evidence.")

    identity = json.loads((HERE / "fixture-identity.json").read_text(encoding="utf-8"))
    if hashlib.sha256((HERE / "source_check.py").read_bytes()).hexdigest() != identity["checker_sha256"]:
        raise ValueError("Checker fingerprint changed. Review the change before replaying.")
    for profile in identity["profiles"].values():
        if hashlib.sha256((HERE / profile["path"]).read_bytes()).hexdigest() != profile["sha256"]:
            raise ValueError("Fictional fixture bytes changed: " + profile["path"])

    cases = json.loads((HERE / "cases.json").read_text(encoding="utf-8"))
    results = []
    for case in cases:
        command = [
            sys.executable, str(HERE / "fixture_check.py"),
            "--fixture-profile", case.get("fixture_profile", "bundle"),
        ] + case.get("arguments", [])
        if case.get("claim_file"):
            command += ["--claim", str(HERE / case["claim_file"])]
        exit_code, report, stdout, stderr = execute(command)
        matches = (
            report.get("status") == case["expected_status"]
            and exit_code == case["expected_exit_code"]
            and all(report.get(key) == value for key, value in case.get("expected_fields", {}).items())
        )
        results.append({
            "id": case["id"],
            "expected_status": case["expected_status"],
            "expected_fields": case.get("expected_fields", {}),
            "exit_code": exit_code,
            "matches_expected": matches,
            "report": report,
            "stdout": stdout,
            "stderr": stderr,
        })

    integration = source_checks() if args.with_source else None
    source_results = integration["results"] if integration else []
    summary = {
        "synthetic_cases": len(results),
        "synthetic_expected_outcomes": sum(row["matches_expected"] for row in results),
        "source_integration_cases": len(source_results),
        "source_expected_outcomes": sum(row["matches_expected"] for row in source_results),
        "python_version": platform.python_version(),
        "independent_human_review": "not_performed",
        "domain_validation": "not_performed",
    }
    files = [HERE / filename for filename in (
        "source_check.py", "fixture_check.py", "run.py", "cases.json",
        "fixture-identity.json", "source-identity.json",
    )] + sorted((HERE / "claims").glob("*.json"))
    receipt = {
        "summary": summary,
        "file_sha256": {
            path.relative_to(HERE).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in files
        },
        "fixture_identity": identity,
        "synthetic_results": results,
        "source_integration": integration,
    }
    with args.output.open("x", encoding="utf-8") as handle:
        json.dump(receipt, handle, indent=2)
        handle.write("\n")
    print(json.dumps(summary, indent=2))
    return 0 if all(row["matches_expected"] for row in results + source_results) else 1


if __name__ == "__main__":
    try:
        exit_code = main()
    except (OSError, ValueError, KeyError) as error:
        print(json.dumps({"status": "replay_input_error", "reason": str(error)}, indent=2))
        exit_code = 2
    raise SystemExit(exit_code)
