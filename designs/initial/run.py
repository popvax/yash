"""Reproduce Design #1, from selection through RNA encoding."""
import argparse
import json
from pathlib import Path
from yash_vaccine.initial import reproduce


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="New directory for reproduced files")
    args = parser.parse_args()
    print(json.dumps(reproduce(Path(__file__).resolve().parents[2], args.output), indent=2))
