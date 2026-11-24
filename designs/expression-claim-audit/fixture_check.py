"""Run the unchanged checker against one declared fictional review source."""
import json
import sys
from pathlib import Path
import source_check as checker

HERE = Path(__file__).resolve().parent

def main():
    try:
        identity = json.loads((HERE / "fixture-identity.json").read_text(encoding="utf-8"))
        parser = checker.ClaimArgumentParser(add_help=False, allow_abbrev=False)
        parser.add_argument("--fixture-profile", choices=sorted(identity["profiles"]), action=checker.SingleValue)
        options, remaining = parser.parse_known_args()
        profile = identity["profiles"][options.fixture_profile or "bundle"]
        checker.ROOT = HERE
        checker.SOURCE = HERE / profile["path"]
        checker.EXPECTED_HASH = profile["sha256"]
        checker.PRACTICE_SOURCE = checker.SOURCE
        checker.PRACTICE_HASH = checker.EXPECTED_HASH
        sys.argv = [sys.argv[0]] + remaining
        return checker.main()
    except (OSError, ValueError, KeyError) as error:
        print(json.dumps({"status": "input_error", "reason": str(error)}, indent=2))
        return 2

if __name__ == "__main__":
    raise SystemExit(main())
