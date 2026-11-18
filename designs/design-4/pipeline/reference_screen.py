"""Exact encoded-peptide search in a checksum-recorded public reference proteome."""
import argparse
import gzip
import re
import os
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

from curate import HERE, peptide_windows
from dataio import load, resolve_gencode_fasta, rows, save, sha


def wanted_peptides():
    wanted = set()
    for protein in load(HERE / "exports/rna-proteins.json"):
        seq = protein["sequence"]
        lo, hi = protein["mutation_start_0based"], protein["mutation_end_0based_exclusive"]
        wanted.update(p for _, p in peptide_windows(seq, lo, hi, (8, 9, 10, 11, 15)))
        match = re.fullmatch(r"p\.([A-Z])\d+([A-Z])", protein["DNA_annotation_protein_change"])
        if match and hi == lo + 1:
            if seq[lo] != match.group(2):
                raise ValueError("Counterfactual reversion disagrees with the source amino-acid annotation")
            wild = seq[:lo] + match.group(1) + seq[hi:]
            wanted.update(p for _, p in peptide_windows(wild, lo, hi, (8, 9, 10, 11, 15)))
    for row in rows(HERE / "exports/candidate-pool.tsv"):
        if row["source"] != "MUTATION":
            seq = row["sequence"]
            wanted.update(
                seq[start:start + length] for length in (8, 9, 10, 11, 15)
                for start in range(len(seq) - length + 1)
            )
    request = load(HERE / "exports/inframe-request.json")
    wanted.update(request["peptides"])
    return wanted


def exact_search(records, wanted):
    prefixes = defaultdict(set)
    for peptide in wanted:
        prefixes[peptide[:3]].add(len(peptide))
    counts = defaultdict(int)
    first_hits = defaultdict(list)
    proteins, residues = 0, 0
    for identifier, sequence in records:
        proteins += 1
        residues += len(sequence)
        for start in range(max(0, len(sequence) - 7)):
            for length in prefixes.get(sequence[start:start + 3], ()):
                if start + length > len(sequence):
                    continue
                peptide = sequence[start:start + length]
                if peptide in wanted:
                    counts[peptide] += 1
                    if len(first_hits[peptide]) < 10:
                        first_hits[peptide].append({
                            "protein_record": identifier, "start_aa": start,
                            "end_aa": start + length,
                        })
    return {
        peptide: {
            "occurs": counts[peptide] > 0,
            "occurrence_count": counts[peptide],
            "first_hits": first_hits[peptide],
            "source_provenance_complete": True,
        }
        for peptide in sorted(wanted)
    }, proteins, residues


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, default=resolve_gencode_fasta("gencode.v46.pc_translations.fa.gz", HERE / "references"))
    parser.add_argument("--output", type=Path, default=HERE / "references")
    args = parser.parse_args()
    with gzip.open(args.fasta, "rt") as stream:
        matches, proteins, residues = exact_search(
            ((record.id, str(record.seq)) for record in SeqIO.parse(stream, "fasta")),
            wanted_peptides(),
        )
    checksum = sha(args.fasta)
    version = re.search(r"\.v(\d+)\.", args.fasta.name)
    if version is None:
        raise ValueError("Reference release must be explicit in the FASTA filename")
    release = version.group(1)
    save(args.output / "peptide-matches.json", matches)
    save(args.output / "reference-receipt.json", {
        "reference_url": f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{release}/gencode.v{release}.pc_translations.fa.gz",
        "reference_release": f"GENCODE human v{release} (GRCh38), protein-coding translations",
        "fasta_file": os.path.relpath(args.fasta, args.output), "fasta_sha256": checksum,
        "proteins_scanned": proteins, "residues_scanned": residues,
        "peptides_queried": len(matches),
        "peptide_match_sha256": sha(args.output / "peptide-matches.json"),
        "scope": (
            "Exact encoded sequence in this reference only. Not exhaustive patient "
            "normal isoforms, noncanonical translation, post-translational modifications, "
            "normal tissue expression, HLA presentation, or TCR cross-reactivity."
        ),
    })
    print("Reference screen:", proteins, "proteins;", len(matches), "peptides;", residues, "residues")


if __name__ == "__main__":
    main()
