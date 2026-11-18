"""Record new local RNA parameters and coding hash, then execute the independent builder."""
import hashlib
import subprocess
import sys
import tempfile
from pathlib import Path

from dataio import load, save, sha
from construct_identifiers import configure_identifiers

HERE = Path(__file__).resolve().parents[1]
MRNA = HERE / "construct/mrna"


def main():
    cfg = load(MRNA / "backbone.json")
    check = load(HERE / "construct/final-check.json")
    cfg.update(
        protein_length_aa=check["protein_length_aa"],
        protein_file_sha256=sha(HERE / "construct/protein/cassette.fasta"),
        protein_components_sha256=sha(HERE / "construct/protein/components.tsv"),
    )
    configure_identifiers(cfg, "MAIN16")
    settings = cfg["optimization"]
    settings.update(
        accessibility_trials=12,
        optimizer_script_sha256=sha(MRNA / "optimizer.py"),
        method="Source PopRNA heuristics, complete alternative-frame-stop nucleotide protection, then bounded full-known-segment ensemble initiation-accessibility refinement.",
    )
    settings["implementation_changes"] = list(dict.fromkeys(settings["implementation_changes"] + [
        "Independently enumerate and freeze every alternative-frame stop nucleotide after stop insertion.",
        "Record every bounded initiation-accessibility trial and verify all protected nucleotides after refinement.",
    ]))
    sources = load(MRNA / "annotation-sources.json")
    sources["record_date"] = "12-NOV-2026"
    for key in ("selected_regions", "routing_architecture"):
        sources[key]["sha256"] = sha(MRNA / sources[key]["path"])
    save(MRNA / "annotation-sources.json", sources)
    protein = "".join(line for line in (HERE / "construct/protein/cassette.fasta").read_text().splitlines()
                      if not line.startswith(">"))
    fixed = cfg["fixed_elements"]
    with tempfile.TemporaryDirectory(prefix="design4-mrna-record-") as directory:
        path = Path(directory) / "coding.txt"
        command = [
            sys.executable, str(MRNA / "optimizer.py"), "--protein", "-i", protein,
            "--output", str(path), "--organism", settings["species_table"],
            "--cutoff", str(settings["codon_cutoff"]), "--seed", str(settings["seed"]),
            "--fold-iterations", str(settings["structure_iterations_per_pass"]),
            "--utr5-rna", fixed["utr5_dna"].replace("T", "U"),
            "--fixed-suffix-rna", (fixed["stop_block_dna"] + fixed["utr3_body_dna"]).replace("T", "U"),
            "--accessibility-trials", str(settings["accessibility_trials"]),
            "--sequences", *settings["avoid_sequences_rna"],
        ]
        subprocess.run(command, check=True)
        coding = path.read_text().strip()
    settings["recorded_coding_sequence_sha256"] = hashlib.sha256(coding.encode()).hexdigest()
    save(MRNA / "backbone.json", cfg)
    subprocess.run([sys.executable, str(MRNA / "run.py")], check=True)


if __name__ == "__main__":
    main()
