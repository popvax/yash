"""Unconstrained LinearDesign references for the same fixed protein and codon table."""
import argparse
import csv
import hashlib
import importlib.util
import re
import subprocess
import sys
from pathlib import Path

import RNA
from Bio.Seq import Seq

from dataio import fasta, load, save, sha, table

HERE = Path(__file__).resolve().parents[1]
MRNA = HERE / "construct/mrna"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--linear-design-root", type=Path, required=True)
    args = parser.parse_args()
    tool = args.linear_design_root.resolve()
    output = MRNA / "benchmarks"
    output.mkdir(parents=True, exist_ok=True)
    codons = load(MRNA / "human-codon-usage.json")["codon_usage"]
    codon_file = output / "human-codon-usage.csv"
    with codon_file.open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["#", "", ""])
        for aa, entries in sorted(codons.items()):
            for codon, frequency in sorted(entries.items()):
                writer.writerow([codon.replace("T", "U"), aa, frequency])
    protein = "".join(line for line in (HERE / "construct/protein/cassette.fasta").read_text().splitlines()
                      if not line.startswith(">"))
    cfg = load(MRNA / "backbone.json")
    audit = load(MRNA / "optimizer-audit.json")
    snapshot = audit["coding_rna_after_stop_insertion"]
    sys.path.insert(0, str(MRNA))
    spec = importlib.util.spec_from_file_location("mrna_metrics", MRNA / "run.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    RNA.cvar.temperature = 37.0
    benchmarks = []
    for value in (0, 0.5, 1, 2, 4):
        label = "lambda-" + str(value).replace(".", "_")
        command = [str(tool / "bin/LinearDesign_2D"), str(value), "1", str(codon_file)]
        result = subprocess.run(command, input=protein + "\n", cwd=tool,
                                text=True, capture_output=True, timeout=600)
        (output / (label + ".log")).write_text(result.stdout + result.stderr)
        match = re.search(r"mRNA sequence:\s+([ACGU]+)", result.stdout)
        if result.returncode or match is None:
            raise RuntimeError("LinearDesign reference failed: " + label + " " + result.stderr)
        coding = match.group(1)
        if str(Seq(coding).translate()) != protein or len(coding) != 3 * len(protein):
            raise ValueError("LinearDesign changed the protein")
        fasta(output / (label + ".fasta"), "LinearDesign_reference_" + label, coding)
        metrics = module.metrics(coding.replace("U", "T"), protein, codons, cfg["optimization"])
        energy = re.search(r"folding free energy:\s*([-0-9.]+)", result.stdout)
        cai = re.search(r"mRNA CAI:\s*([0-9.]+)", result.stdout)
        row = {
            "lambda": value, "protein_length_aa": len(protein),
            "coding_rna_sha256": hashlib.sha256(coding.encode()).hexdigest(),
            "LinearDesign_reported_energy": float(energy.group(1)) if energy else None,
            "LinearDesign_reported_CAI": float(cai.group(1)) if cai else None,
            **metrics,
            "source_stop_and_initial_nt_preserved": all(
                coding[position] == snapshot[position]
                for position in audit["protected_nt_positions"]
            ),
            "clinical_release_ready": False,
        }
        benchmarks.append(row)
    table(output / "lineardesign-benchmark.tsv", benchmarks)
    save(output / "receipt.json", {
        "protein_sha256": sha(HERE / "construct/protein/cassette.fasta"),
        "executable_sha256": sha(tool / "bin/LinearDesign_2D"),
        "codon_table_sha256": sha(codon_file),
        "export_identifier_revision": {
            "construct_variant": "MAIN16",
            "current_protein_file_sha256": "6759bb7225980e5aad288bbd8f63a1c300002d807adc8ea3c5051c24e3d3a7d6",
            "previous_protein_file_sha256": "3c8026244fdbc8020f721cbcb0d8ac1eef5e8a66d6aa42387e4e278ada2e4734",
            "protein_sequence_sha256": "48eb33bf91f9e8ed852b120a534d1dcd765e20ea2755afed87116a576fad1960",
            "scope": "Benchmark outputs and before-rescue snapshots were recorded prior to the MAIN16 FASTA-header revision; all 507 protein residues and 1,524 CDS bases are unchanged.",
        },
        "lambda_values": [0, 0.5, 1, 2, 4],
        "same_fixed_protein_and_same_human_codon_frequencies": True,
        "scope": "Five unconstrained algorithm reference points; not a proven global Pareto frontier, and not alternatives satisfying the source motif/repeat/stop constraints.",
        "clinical_validation": False,
        "output_sha256": {path.name: sha(path) for path in sorted(output.iterdir())
                          if path.is_file() and path.name != "receipt.json"},
    })
    print("Saved five same-protein, same-table LinearDesign reference points")


if __name__ == "__main__":
    main()
