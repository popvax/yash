"""Replay the specificity-first design offline from included evidence.

Prerequisites (see README.md): MHCflurry 2.3.14 pan-allele weights, the
Ensembl 103 proteome and GTF under /workspace/ref, the recorded normal coding
GVCF, and a pinned RNA venv (.venv-mrna) for the PopRNA stage. No cloud jobs.
"""
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
STAGES = [
    "splice_specificity_screen.py",
    "build_selected.py",
    "build_regions.py",
    "build_prediction_cache.py",
    "build_terminal_cache.py",
    "fresh_junction_check.py",
    "hla_robustness.py",
    "prepare_mrna_stage.py",
    "compare_to_design2.py",
]
RNA_STAGE = HERE / "construct" / "mrna" / "run.py"
RNA_VENV = REPO / ".venv-mrna" / "bin" / "python"


def main():
    for script in STAGES:
        print("== " + script, flush=True)
        subprocess.run([sys.executable, HERE / script], check=True)
    if RNA_VENV.exists():
        print("== construct/mrna/run.py (pinned RNA venv)", flush=True)
        subprocess.run([str(RNA_VENV), RNA_STAGE], check=True, cwd=RNA_STAGE.parent)
    else:
        print("== RNA stage skipped: create the pinned venv first "
              "(numpy==1.26.4 ViennaRNA==2.7.2 python-codon-tables==0.1.18 "
              "biopython==1.86) and rerun.", flush=True)


if __name__ == "__main__":
    main()
