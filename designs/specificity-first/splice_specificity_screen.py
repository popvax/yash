"""Specificity-attribution screen over the full recorded splice pool.

For every splice window in the Design #2 candidate pool:
  1. enumerate junction-overlapping 8-11mers and count exact matches against
     the Ensembl 103 primary protein-coding and extended (all-translation)
     reference scopes, using the lead design's own matching definitions;
  2. fetch the reference genomic sequence spanning each window's recorded
     junction, rebuild the spliced reference coding sequence, and identify
     every extra substitution (codon differences outside the junction codon);
  3. build counterfactual windows with the extra substitutions reverted and
     recount nonself peptides, attributing peptide novelty to the coordinate
     join versus the extra substitutions.

Computational screen only; no claim about RNA-sequence artifact status.
"""
import csv
import gzip
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PROTEIN_SHA256 = "f8df3d7cfc85718e9161579bd9168a3598547ee50022b5e1c272ed66293f0505"
OUT = Path(__file__).resolve().parent / "splice-specificity-attribution.tsv"

CODON = {}
bases = "TCAG"
aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
i = 0
for b1 in bases:
    for b2 in bases:
        for b3 in bases:
            CODON[b1 + b2 + b3] = aas[i]
            i += 1


def translate(nt):
    aa = []
    for i in range(0, len(nt) - len(nt) % 3, 3):
        aa.append(CODON.get(nt[i:i + 3], "X"))
    return "".join(aa)


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
    return sorted({sequence[start:start + length] for length in range(8, 12)
                   for start in range(len(sequence) - length + 1)
                   if start < hi and start + length > lo})


def revcomp(seq):
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def rest_sequence(region):
    url = f"https://rest.ensembl.org/sequence/region/human/{region}"
    request = urllib.request.Request(url, headers={"Content-Type": "text/plain"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(request, timeout=30) as response:
                return response.read().decode().strip()
        except Exception as error:
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)


def window_reference_cds(chrom, positions_0based, strand_minus):
    """Fetch the reference genome across the window footprint and rebuild the
    spliced coding sequence (gap removed) in transcript orientation."""
    ordered = sorted(set(positions_0based))
    lo, hi = ordered[0], ordered[-1]
    span = rest_sequence(f"{chrom}:{lo + 1}-{hi + 1}:1")
    blocks, current = [], [ordered[0]]
    for previous, position in zip(ordered, ordered[1:]):
        if position - previous == 1:
            current.append(position)
        else:
            blocks.append(current)
            current = [position]
    blocks.append(current)
    parts = []
    for block in blocks:
        blo, bhi = min(block), max(block)
        segment = span[blo - lo:bhi - lo + 1]
        if len(segment) != bhi - blo + 1:
            raise ValueError("reference span extraction failed")
        parts.append(segment)
    cds = "".join(parts)
    return revcomp(cds) if strand_minus else cds


def main():
    proteome = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT.parent / "ref" / "pep.all.fa.gz"
    digest = __import__("hashlib").sha256(proteome.read_bytes()).hexdigest()
    if digest != PROTEIN_SHA256:
        raise ValueError("Reference FASTA checksum differs from Ensembl 103")
    records = list(translations(proteome))
    primary = "!".join(seq for header, seq in records if primary_coding(header))
    extended = "!".join(seq for _, seq in records)
    n_primary = sum(primary_coding(h) for h, _ in records)
    print(f"translations={len(records)} primary_coding={n_primary}", flush=True)

    with gzip.open(ROOT / "data/shared-rna/splicing/reconstruction/patient-windows.jsonl.gz", "rt") as stream:
        windows_by_id = {}
        for line in stream:
            record = json.loads(line)
            windows_by_id[record["window_id"]] = record

    pool_path = ROOT / "designs/design-2/inputs/integrated/candidate-pool.tsv"
    with pool_path.open() as stream:
        pool = [row for row in csv.DictReader(stream, delimiter="\t") if row["source"] == "SPLICE"]
    print(f"splice windows in pool: {len(pool)}", flush=True)

    rows_out = []
    for row in pool:
        sequence = row["sequence"]
        window_id = row["window_id"]
        record = windows_by_id.get(window_id)
        if record is None:
            raise ValueError(f"window {window_id} not found in reconstruction")
        lo = int(row["target_start_in_window"])
        hi = int(row["target_end_in_window"])
        peptides = overlapping_peptides(sequence, lo, hi)
        original_primary = sum(p not in primary for p in peptides)
        original_extended = sum(p not in extended for p in peptides)

        coordinates = record["coordinates"]
        match = re.match(r"^(chr[^:]+):(\d+)-(\d+)\(([+-])\)$", coordinates)
        if not match:
            raise ValueError("unparsed coordinates: " + coordinates)
        chrom, _, _, strand = match.groups()
        strand_minus = strand == "-"
        positions = record["genomic_positions_0based"]
        cds_observed = record["coding_sequence"]
        junction_nt = int(record["junction_nt_offset"])
        start_relative = int(record["start_nt_relative_to_join"])

        reference_cds = window_reference_cds(chrom, positions, strand_minus)
        # The observed coding sequence starts `start_relative` nt relative to the
        # join; align the reference spliced CDS to the observed window by the join.
        join_in_observed = junction_nt - start_relative if start_relative < 0 else junction_nt
        # Simpler: both sequences span the same genomic footprint minus the gap,
        # so compare directly when lengths agree.
        extra_base_changes = []
        junction_codon_index = None
        if len(reference_cds) == len(cds_observed):
            # locate the join: the genomic discontinuity, in transcript orientation
            ordered = sorted(set(positions))
            blocks, current = [], [ordered[0]]
            for previous, position in zip(ordered, ordered[1:]):
                if position - previous == 1:
                    current.append(position)
                else:
                    blocks.append(current)
                    current = [position]
            blocks.append(current)
            if len(blocks) == 2:
                join_nt = len(blocks[1]) if strand_minus else len(blocks[0])
                junction_codon_index = join_nt // 3
                for index, (observed_base, reference_base) in enumerate(zip(cds_observed, reference_cds)):
                    if observed_base != reference_base:
                        codon = index // 3
                        extra_base_changes.append(dict(
                            nt_index=index, codon_index=codon,
                            observed_base=observed_base, reference_base=reference_base,
                            in_junction_codon=(codon == junction_codon_index)))
            else:
                join_nt = None
        else:
            join_nt = None

        observed_aa = translate(cds_observed) if len(cds_observed) % 3 == 0 else ""
        reference_aa = translate(reference_cds) if len(reference_cds) == len(cds_observed) and len(reference_cds) % 3 == 0 else ""
        extra_aa = []
        if observed_aa and reference_aa and len(observed_aa) == len(reference_aa):
            for index, (obs, ref) in enumerate(zip(observed_aa, reference_aa)):
                if obs != ref and (junction_codon_index is None or index != junction_codon_index):
                    extra_aa.append(dict(position_0based=index, observed=obs, reference=ref))

        # counterfactuals: revert all extra aa differences, and each singly
        def reverted_sequence(changes):
            aa = list(observed_aa)
            for change in changes:
                aa[change["position_0based"]] = change["reference"]
            return "".join(aa)

        counterfactuals = []
        if extra_aa:
            all_reverted = reverted_sequence(extra_aa)
            peptides_reverted = overlapping_peptides(all_reverted, lo, hi)
            counterfactuals.append(dict(
                scope="all_extra_reverted",
                reverted=[f"{c['reference']}{c['position_0based'] + 1}{c['observed']}" for c in extra_aa],
                reverted_sequence=all_reverted,
                primary_nonself=sum(p not in primary for p in peptides_reverted),
                extended_nonself=sum(p not in extended for p in peptides_reverted)))
            for change in extra_aa:
                single = reverted_sequence([change])
                peptides_single = overlapping_peptides(single, lo, hi)
                counterfactuals.append(dict(
                    scope="single_reverted",
                    reverted=[f"{change['reference']}{change['position_0based'] + 1}{change['observed']}"],
                    reverted_sequence=single,
                    primary_nonself=sum(p not in primary for p in peptides_single),
                    extended_nonself=sum(p not in extended for p in peptides_single)))

        if original_primary == 0:
            classification = "self_no_novelty"
        elif not extra_aa and len(reference_cds) == len(cds_observed) and len(cds_observed) % 3 == 0:
            classification = "junction_only_novelty"
        elif original_primary > 0 and counterfactuals and counterfactuals[0]["primary_nonself"] == 0:
            classification = "substitution_dependent"
        elif original_primary > 0 and counterfactuals:
            classification = "junction_only_novelty"
        else:
            classification = "novelty_attribution_incomplete"

        rows_out.append(dict(
            window_id=window_id, gene=row["gene"], gene_id=row["gene_id"],
            coordinates=coordinates, sequence=sequence,
            target_start_in_window=lo, target_end_in_window=hi,
            junction_nt_offset=junction_nt, start_nt_relative_to_join=start_relative,
            fragment_count=int(row["fragment_count"]),
            supplied_transformed_gene_tpm=row["supplied_transformed_gene_tpm"],
            best_target_affinity_nm=row["best_target_affinity_nm"],
            native_affinity_sum=row["native_affinity_sum"],
            junction_is_annotated=row["junction_is_annotated"],
            coordinate_topology=row["coordinate_topology"],
            native_event_review_tier=row["native_event_review_tier"],
            region_id=row["region_id"],
            original_primary_nonself=original_primary,
            original_extended_nonself=original_extended,
            distinct_junction_peptides=len(peptides),
            extra_aa_changes=";".join(f"{c['reference']}{c['position_0based'] + 1}{c['observed']}" for c in extra_aa),
            extra_base_change_count=len(extra_base_changes),
            cds_length_match=len(reference_cds) == len(cds_observed),
            classification=classification,
            all_reverted_primary_nonself=counterfactuals[0]["primary_nonself"] if counterfactuals else "",
            all_reverted_extended_nonself=counterfactuals[0]["extended_nonself"] if counterfactuals else "",
            counterfactuals_json=json.dumps(counterfactuals),
            observed_aa=observed_aa, reference_aa=reference_aa))
        print(f"{row['gene']} {window_id}: primary_nonself={original_primary} "
              f"extra_aa={len(extra_aa)} -> {classification}", flush=True)

    fields = list(rows_out[0].keys())
    with OUT.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows_out)
    summary = dict(
        windows=len(rows_out),
        classification={c: sum(r["classification"] == c for r in rows_out) for c in {r["classification"] for r in rows_out}},
        reference_translations=len(records), primary_coding_translations=n_primary,
        reference_protein_fasta_sha256=PROTEIN_SHA256)
    (OUT.parent / "splice-specificity-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
