"""Independently audit the included patient-RNA splice evidence offline."""
import argparse
import json
from pathlib import Path
from yash_vaccine.splice_evidence import verify

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(verify(args.results, args.output), indent=2))
