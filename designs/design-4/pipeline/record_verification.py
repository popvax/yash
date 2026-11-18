"""Execute cross-platform verification, unit tests, reference rescans and full dual replay.

Uses tempfile.TemporaryDirectory() (never OS-specific /dev/shm), writes sanitized logs to
review/logs/, separates current execution evidence from historical upstream Linux provenance
in review/verification.json, and refreshes manifest.json cleanly.
"""
from importlib.metadata import version
import json
import platform
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from dataio import load, save, sha
from package_review import HERE, included_files, write_manifest


LOG_DIR = HERE / "review/logs"
VER_PATH = HERE / "review/verification.json"


def _scientific_input_hashes():
    excluded = {"review/verification.json"}
    out = {}
    for p in included_files():
        rel = str(p.relative_to(HERE))
        if rel in excluded or rel.startswith("review/logs/"):
            continue
        out[rel] = sha(p)
    return out


def _sanitize_text(text, replay_dir_str=None):
    cleaned = text.replace(str(HERE), "<DESIGN4_ROOT>")
    if replay_dir_str:
        cleaned = cleaned.replace(replay_dir_str, "<TEMP_REPLAY_DIR>")
    # Normalize unittest wall-clock line ("Ran N tests in X.XXXs") only in stdout/stderr if needed?
    # Keep actual measured output intact — only local filesystem prefixes are anonymized.
    return cleaned


def _build_receipt_payload(checks, unit_test_count, replay_summary, input_hashes, overall_status="COMPLETED"):
    return {
        "overall_status": overall_status,
        "clinical_benefit_established": False,
        "clinical_release_ready": False,
        "clinical_safety_established": False,
        "clinical_validation": False,
        "guarantee_of_zero_harm_cure_or_legal_immunity": False,
        "named_patient_chain_of_custody_authenticated": False,
        "current_check_does_not_repeat": [
            "Licensed predictor inference (unless recheck_construct.py, safety_audit.py generate, or audit_ctdna_loci.py --execute is run)",
            "Live AlphaGenome 2 API calls (unless alphagenome2_audit.py --execute is run with ALPHAGENOME_API_KEY)",
            "Raw alignment or variant calling",
            "Whole tumour/normal DNA fingerprint analysis",
            "Laboratory or clinical experiments",
        ],
        "data_access_limits": {
            "DNA_indexed_queries_distinct_from_whole_file_verification": True,
            "all_data_scientifically_reanalysed": False,
            "whole_RNA_BAM_checksum_checked_in_this_workflow": False,
            "historical_RNA_checksum_evidence": "review/identity-concordance.json",
            "whole_tumour_and_normal_DNA_BAM_payload_checksums_independently_verified": False,
        },
        "current_execution": {
            "environment": {
                "system": platform.system(),
                "machine": platform.machine(),
                "python": platform.python_version(),
                "packages": {
                    "ViennaRNA": version("ViennaRNA"),
                    "biopython": version("biopython"),
                    "numpy": version("numpy"),
                    "ortools": version("ortools"),
                    "python-codon-tables": version("python-codon-tables"),
                },
            },
            "replay_directory_policy": "Cross-platform Python tempfile.TemporaryDirectory() (no /dev/shm dependency)",
            "unit_tests_passed": unit_test_count,
            "checks": checks,
            "replay_summary": replay_summary,
        },
        "historical_upstream_execution_reference": {
            "note": (
                "Retained strictly as historical provenance from the earlier 224-file Linux revision; "
                "these timings and /dev/shm paths are NOT from the current execution recorded above."
            ),
            "historical_platform": "Linux",
            "historical_checks": [
                {"check": "verify_results.py --rescan-reference", "elapsed_seconds": 81.262, "exit_code": 0},
                {"check": "unittest discover -s tests -p test_*.py", "elapsed_seconds": 5.372, "exit_code": 0},
                {"check": "run.py --output <linux-tmpfs>", "elapsed_seconds": 123.752, "exit_code": 0},
            ],
        },
        "identity_check": {
            "exact_threshold_matches": 93,
            "markers": 110,
            "normal_heterozygote_RNA_homozygote_calls": 7,
            "opposite_homozygote_conflicts": 0,
            "opposite_homozygote_eligible_markers": 8,
            "portable_checker_reexecuted_in_this_workflow": False,
            "evidence_status": "Historical result from the referenced identity receipt; not reexecuted by this workflow.",
            "receipt": "identity-concordance.json",
            "scope": "Bounded donor concordance, not named-patient authentication or an identity probability.",
            "uncertain_RNA_calls": 10,
        },
        "executed_scientific_inputs_sha256": input_hashes,
    }


def main():
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    input_hashes = _scientific_input_hashes()

    # Seed placeholder files + manifest so verify_results.py sees a complete manifest during execution
    placeholder_logs = [
        ("final-reference-rescan", ["verify_results.py", "--rescan-reference"], "review/logs/01-verify-rescan.log"),
        ("final-tests", ["-m", "unittest", "discover", "-s", "tests", "-p", "test_*.py"], "review/logs/02-unittest.log"),
        ("final-full-replay", ["run.py", "--output", "<TEMP_REPLAY_DIR>"], "review/logs/03-full-replay.log"),
    ]
    seed_checks = []
    for cname, cargs, lrel in placeholder_logs:
        lp = HERE / lrel
        if not lp.is_file():
            lp.write_text("pending\n")
        seed_checks.append({
            "check": cname, "args": cargs, "status": "PENDING_NOT_EXECUTED", "exit_code": None,
            "elapsed_seconds": None, "log_file": lrel, "log_sha256": sha(lp),
        })
    save(VER_PATH, _build_receipt_payload(seed_checks, None, {}, input_hashes, overall_status="IN_PROGRESS"))
    write_manifest()

    recorded_checks = []

    # 1. Run verify_results.py --rescan-reference
    t0 = time.perf_counter()
    p1 = subprocess.run(
        [sys.executable, str(HERE / "verify_results.py"), "--rescan-reference"],
        cwd=str(HERE), capture_output=True, text=True, check=True,
    )
    dt1 = round(time.perf_counter() - t0, 3)
    log1 = HERE / "review/logs/01-verify-rescan.log"
    log1.write_text(_sanitize_text(p1.stdout + p1.stderr))
    write_manifest()
    recorded_checks.append({
        "check": "final-reference-rescan",
        "args": ["verify_results.py", "--rescan-reference"],
        "status": "COMPLETED",
        "exit_code": p1.returncode,
        "elapsed_seconds": dt1,
        "log_file": "review/logs/01-verify-rescan.log",
        "log_sha256": sha(log1),
    })

    # 2. Run full dual-construct replay into a cross-platform TemporaryDirectory
    with tempfile.TemporaryDirectory(prefix="yash-design4-replay-") as tmp_parent:
        replay_target = Path(tmp_parent) / "replay-out"
        t0 = time.perf_counter()
        p3 = subprocess.run(
            [sys.executable, str(HERE / "run.py"), "--output", str(replay_target)],
            cwd=str(HERE), capture_output=True, text=True, check=True,
        )
        dt3 = round(time.perf_counter() - t0, 3)
        replay_res = load(replay_target / "replay-result.json")
        log3 = HERE / "review/logs/03-full-replay.log"
        log3.write_text(_sanitize_text(p3.stdout + p3.stderr, str(replay_target)))
        recorded_checks.append({
            "check": "final-full-replay",
            "args": ["run.py", "--output", "<TEMP_REPLAY_DIR>"],
            "status": "COMPLETED",
            "exit_code": p3.returncode,
            "elapsed_seconds": dt3,
            "log_file": "review/logs/03-full-replay.log",
            "log_sha256": sha(log3),
        })

    # Write provisional receipt marking final-tests IN_PROGRESS (exit_code: None) while unittest executes
    log2 = HERE / "review/logs/02-unittest.log"
    provisional_checks = [
        recorded_checks[0],
        {
            "check": "final-tests",
            "args": ["-m", "unittest", "discover", "-s", "tests", "-p", "test_*.py"],
            "status": "IN_PROGRESS",
            "exit_code": None,
            "elapsed_seconds": None,
            "log_file": "review/logs/02-unittest.log",
            "log_sha256": sha(log2),
        },
        recorded_checks[1],
    ]
    save(VER_PATH, _build_receipt_payload(provisional_checks, None, replay_res, input_hashes, overall_status="IN_PROGRESS"))
    write_manifest()

    # 3. Run unittest suite
    t0 = time.perf_counter()
    p2 = subprocess.run(
        [sys.executable, "-m", "unittest", "discover", "-s", "tests", "-p", "test_*.py"],
        cwd=str(HERE), capture_output=True, text=True, check=True,
    )
    dt2 = round(time.perf_counter() - t0, 3)
    log2.write_text(_sanitize_text(p2.stdout + p2.stderr))
    m_tests = re.search(r"Ran (\d+) tests", p2.stderr + p2.stdout)
    n_tests = int(m_tests.group(1)) if m_tests else 0

    final_checks = [
        recorded_checks[0],
        {
            "check": "final-tests",
            "args": ["-m", "unittest", "discover", "-s", "tests", "-p", "test_*.py"],
            "status": "COMPLETED",
            "exit_code": p2.returncode,
            "elapsed_seconds": dt2,
            "log_file": "review/logs/02-unittest.log",
            "log_sha256": sha(log2),
        },
        recorded_checks[1],
    ]
    save(VER_PATH, _build_receipt_payload(final_checks, n_tests, replay_res, input_hashes, overall_status="COMPLETED"))
    write_manifest()

    # 4. Final read-only invariant check confirming COMPLETED status + manifest + verification.json fixed point
    final_ver = load(VER_PATH)
    if final_ver["overall_status"] != "COMPLETED" or any(c["status"] != "COMPLETED" or c["exit_code"] != 0 for c in final_ver["current_execution"]["checks"]):
        raise RuntimeError("Final verification receipt is incomplete")
    subprocess.run([sys.executable, str(HERE / "verify_results.py")], cwd=str(HERE), check=True, capture_output=True, text=True)
    subprocess.run([sys.executable, "-m", "unittest", "discover", "-s", "tests", "-p", "test_*.py"], cwd=str(HERE), check=True, capture_output=True, text=True)
    print(json.dumps({
        "unit_tests_passed": n_tests,
        "rescan_seconds": dt1,
        "unittest_seconds": dt2,
        "replay_seconds": dt3,
        "replay_max_float_diff": replay_res["diagnostic_max_abs_float_diff"],
    }, indent=2))


if __name__ == "__main__":
    main()
