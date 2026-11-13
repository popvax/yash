"""Pin the deterministic PopRNA coding-sequence output into backbone.json.

Reads the CDS produced by run.py (coding sequence = CDS minus the first stop
triplet), records its SHA-256 as the design's recorded coding-sequence pin, and
prints it. A rerun of run.py then verifies bit-identical reproduction.
"""
import csv
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
MRNA = HERE / "construct" / "mrna"


def main():
    cds = "".join(
        line.strip() for line in (MRNA / "cds.fasta").read_text().splitlines()[1:])
    assert len(cds) % 3 == 0 and cds[-3:] in ("TAA", "TAG", "TGA")
    coding = cds[:-3]
    digest = hashlib.sha256(coding.encode("ascii")).hexdigest()
    backbone_path = MRNA / "backbone.json"
    backbone = json.loads(backbone_path.read_text())
    backbone["optimization"]["recorded_coding_sequence_sha256"] = digest
    backbone_path.write_text(json.dumps(backbone, indent=2) + "\n")
    print(json.dumps({"recorded_coding_sequence_sha256": digest,
                      "coding_nt": len(coding)}, indent=2))


if __name__ == "__main__":
    main()
