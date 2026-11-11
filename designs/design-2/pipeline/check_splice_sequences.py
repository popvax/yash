"""Check splice-region reference matches and optionally the focal normal-DNA base.

Uses the provider Ensembl 103 protein FASTA. No selection or predictions change.
The optional normal check requires pysam and uses indexed public data.
"""
from collections import Counter, defaultdict
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import json
import re

ROOT = Path(__file__).resolve().parents[3]
PROTEIN_SHA256 = "f8df3d7cfc85718e9161579bd9168a3598547ee50022b5e1c272ed66293f0505"
CHANGES = [("F5", 8, "F", "S"), ("MUC5B", 12, "I", "S"),
           ("MUC3A", 14, "F", "I"), ("MUC6", 6, "M", "T")]


def translations(path):
    header, parts = None, []
    with gzip.open(path, "rt") as stream:
        for line in stream:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header, parts = line[1:].strip(), []
            else:
                parts.append(line.strip())
    if header is not None:
        yield header, "".join(parts)


def primary_coding(header):
    match = re.search(r"chromosome:GRCh38:([^:]+):", header)
    chromosomes = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}
    return bool(match and match[1] in chromosomes
                and "gene_biotype:protein_coding " in header
                and "transcript_biotype:protein_coding " in header)


def overlapping_peptides(sequence, lo, hi):
    return {sequence[start:start + length] for length in range(8, 12)
            for start in range(len(sequence) - length + 1)
            if start < hi and start + length > lo}


def normal_base():
    import pysam
    access = json.loads((ROOT / "data/catalog/access.json").read_text())
    key = "raw/TN25-112103/normal/normal.recal.bam"
    url = access["public_base_url"].rstrip("/") + "/" + key
    observations = defaultdict(set)
    with pysam.AlignmentFile(url, "rb", index_filename=url + ".bai") as bam:
        chrom = "chr11" if "chr11" in bam.references else "11"
        for read in bam.fetch(chrom, 1017575, 1017576):
            if (read.is_unmapped or read.is_duplicate or read.is_secondary
                    or read.is_supplementary or read.is_qcfail
                    or read.mapping_quality < 20 or read.query_qualities is None):
                continue
            for query, reference in read.get_aligned_pairs(matches_only=True):
                if reference == 1017575 and read.query_qualities[query] >= 20:
                    observations[read.query_name].add(read.query_sequence[query])
    counts = Counter(next(iter(bases)) if len(bases) == 1 else "CONFLICT"
                     for bases in observations.values())
    return dict(bam_url=url, chromosome="chr11", position_1based=1017576,
                reference_base="G", methionine_associated_alternate_base="A",
                minimum_mapping_quality=20, minimum_base_quality=20,
                duplicate_secondary_supplementary_qcfail_excluded=True,
                count_unit="fragment names; overlapping mates counted once",
                base_counts={base: counts[base] for base in ["A", "C", "G", "T", "N", "CONFLICT"]})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-proteins", type=Path, required=True)
    parser.add_argument("--check-normal", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Choose a new output path")
    if hashlib.sha256(args.reference_proteins.read_bytes()).hexdigest() != PROTEIN_SHA256:
        raise ValueError("Reference FASTA checksum differs from Ensembl 103")
    records = list(translations(args.reference_proteins))
    primary = "!".join(sequence for header, sequence in records if primary_coding(header))
    extended = "!".join(sequence for _, sequence in records)
    with (ROOT / "designs/design-2/selection/selected.tsv").open() as stream:
        selected = {row["gene"]: row for row in csv.DictReader(stream, delimiter="\t")}
    comparisons = []
    for gene, position, observed, reference in CHANGES:
        row = selected[gene]
        sequence = row["sequence"]
        if sequence[position - 1] != observed:
            raise ValueError("Selected residue differs: " + gene)
        reverted = sequence[:position - 1] + reference + sequence[position:]
        lo, hi = int(row["target_start_in_window"]), int(row["target_end_in_window"])
        peptides = overlapping_peptides(sequence, lo, hi)
        reverted_peptides = overlapping_peptides(reverted, lo, hi)
        comparisons.append(dict(gene=gene, region_id=row["region_id"], sequence=sequence,
            residue_position_1based=position, observed_residue=observed, reference_residue=reference,
            reverted_sequence=reverted, target_interval_0based=[lo, hi],
            original_primary_nonself_peptides=sum(p not in primary for p in peptides),
            reverted_primary_nonself_peptides=sum(p not in primary for p in reverted_peptides),
            original_extended_nonself_peptides=sum(p not in extended for p in peptides),
            full_window_extended_matches=[header for header, protein in records if sequence in protein]))
    result = dict(reference_protein_fasta_sha256=PROTEIN_SHA256,
                  reference_translations=len(records),
                  primary_coding_translations=sum(primary_coding(h) for h, _ in records),
                  comparisons=comparisons)
    if args.check_normal:
        result["MUC6_normal_base"] = normal_base()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print("Checked four splice sequences" + (" and the focal MUC6 normal-DNA base" if args.check_normal else ""))


if __name__ == "__main__":
    main()
