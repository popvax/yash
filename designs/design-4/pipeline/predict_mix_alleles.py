#!/usr/bin/env python3
"""Allele-only MixMHCpred CLI adapter using literal argv, never an eval-based alignment path."""
import argparse
import os
import subprocess
import sys
from pathlib import Path

root = Path(os.environ["MIXMHCPRED_ROOT"]).resolve()
if "MixMHCpred3.0" not in (root / "MixMHCpred").read_text():
    raise ValueError("This adapter requires the documented MixMHCpred 3.0 layout")
if "--help" in sys.argv or "-h" in sys.argv:
    print("MixMHCpred3.0")
    print("Allele-scoring adapter: -i INPUT -o OUTPUT -a A0101,A3303,...")
    raise SystemExit(0)
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-i", "--input", required=True)
parser.add_argument("-o", "--output", required=True)
parser.add_argument("-a", "--alleles", required=True)
parser.add_argument("-m", default="0")
args = parser.parse_args()
raise SystemExit(subprocess.run([
    sys.executable, str(root / "code/main.py"), "-i", args.input, "-o", args.output,
    "-a", args.alleles, "-m", args.m, "-l", str(root / "lib"),
]).returncode)
