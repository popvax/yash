"""Independent local patient-normal self check for selected splice windows.

Model semantics follow the lead design's normal-self scope: annotated coding
transcripts (Ensembl 103 GTF here; the lead used GENCODE 37/Ensembl 103)
combined with the patient's genotyped normal coding variants (public normal
coding GVCF), enumerated as bounded unphased local possibilities. For each
selected splice window we collect annotated CDS segments overlapping the
window's exonic blocks, translate them in frame, apply local genotype
possibilities, and test whether any junction-overlapping 8-11mer of the
window occurs. A hit is conservatively treated as potential patient self.
"""
import csv
import gzip
import json
import re
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GTF = Path("/workspace/Homo_sapiens.GRCh38.103.gtf.gz")
GVCF = Path("/workspace/normal.coding.g.vcf.gz")
MAX_COMBINATIONS = 4096
CLUSTER_SEPARATION = 33
MARGIN = 50

CODON = {}
_bases = "TCAG"
_aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_i = 0
for _b1 in _bases:
    for _b2 in _bases:
        for _b3 in _bases:
            CODON[_b1 + _b2 + _b3] = _aas[_i]
            _i += 1


def translate(nt):
    return "".join(CODON.get(nt[i:i + 3], "X") for i in range(0, len(nt) - len(nt) % 3, 3))


def revcomp(seq):
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def overlapping_peptides(sequence, lo, hi):
    return sorted({sequence[start:start + length] for length in range(8, 12)
                   for start in range(len(sequence) - length + 1)
                   if start < hi and start + length > lo})


def fetch_region(region):
    url = f"https://rest.ensembl.org/sequence/region/human/{region}"
    request = urllib.request.Request(url, headers={"Content-Type": "text/plain"})
    with urllib.request.urlopen(request, timeout=60) as response:
        return response.read().decode().strip().upper()


def gtf_cds_segments(gene_name, chrom):
    """CDS segments for every transcript of a gene: transcript -> list of
    (start_0based, end_exclusive, phase) in transcript order."""
    transcripts = {}
    with gzip.open(GTF, "rt") as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if fields[2] != "CDS":
                continue
            attributes = dict(re.findall(r'(\S+) "([^"]*)";', fields[8]))
            if attributes.get("gene_name") != gene_name:
                continue
            transcript = attributes["transcript_id"]
            start, end, strand, frame = int(fields[3]), int(fields[4]), fields[6], int(fields[7])
            transcripts.setdefault(transcript, dict(strand=strand, segments=[]))
            transcripts[transcript]["segments"].append((start - 1, end, frame))
    for transcript in transcripts.values():
        reverse = transcript["strand"] == "-"
        transcript["segments"].sort(key=lambda s: s[0], reverse=reverse)
    return transcripts


def gvcf_variants(chrom, lo_0based, hi_0based):
    import pysam
    vcf = pysam.VariantFile(str(GVCF))
    sample = vcf.header.samples[0]
    variants = []
    for record in vcf.fetch(chrom, lo_0based, hi_0based):
        alts = [a for a in (record.alts or []) if a != "<NON_REF>" and "*" not in a]
        if not alts:
            continue
        gt = record.samples[sample].get("GT") or ()
        alleles = record.samples[sample].alleles or ()
        choices = {}
        for allele in gt:
            if allele is None or allele == 0 or allele >= len(alleles):
                continue
            alt = alleles[allele]
            if alt in ("<NON_REF>",) or "*" in alt:
                continue
            choices[allele - 1] = alt
        if choices:
            variants.append(dict(pos=record.pos - 1, ref=record.ref, choices=choices))
    return variants


def local_possibilities(sequence, variants):
    if not variants:
        return [sequence]
    variants = sorted(variants, key=lambda v: v["pos"])
    clusters, current = [], [variants[0]]
    for variant in variants[1:]:
        if variant["pos"] - current[-1]["pos"] <= CLUSTER_SEPARATION:
            current.append(variant)
        else:
            clusters.append(current)
            current = [variant]
    clusters.append(current)
    sequences = [sequence]
    for cluster in clusters:
        combinations = [[]]
        for variant in cluster:
            combinations = [combo + [choice] for combo in combinations
                            for choice in variant["choices"].items()]
            combinations = combinations[:MAX_COMBINATIONS]
        next_sequences = []
        for sequence in sequences:
            for combo in combinations:
                mutated = sequence
                try:
                    for (allele_index, alt), variant in zip(combo, cluster):
                        offset = variant["pos"] - variant["origin"]
                        if offset < 0 or mutated[offset:offset + len(variant["ref"])] != variant["ref"]:
                            raise ValueError
                        mutated = mutated[:offset] + alt + mutated[offset + len(variant["ref"]):]
                except ValueError:
                    continue
                next_sequences.append(mutated)
        sequences = next_sequences or [sequence]
    return sequences


def main():
    selected = list(csv.DictReader(
        (Path(__file__).resolve().parent / "selected-splice-windows.tsv").open(), delimiter="\t"))
    # primary reference scope, same definition as the specificity screen
    import hashlib
    proteome = ROOT.parent / "ref" / "pep.all.fa.gz"
    if hashlib.sha256(proteome.read_bytes()).hexdigest() != \
            "f8df3d7cfc85718e9161579bd9168a3598547ee50022b5e1c272ed66293f0505":
        raise ValueError("Reference FASTA checksum differs from Ensembl 103")

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

    primary = "!".join(seq for header, seq in translations(proteome)
                       if re.search(r"chromosome:GRCh38:([^:]+):", header)
                       and re.search(r"chromosome:GRCh38:([^:]+):", header).group(1)
                       in {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}
                       and "gene_biotype:protein_coding " in header
                       and "transcript_biotype:protein_coding " in header)
    for row in selected:
        row["_nonself_exempt"] = primary

    with gzip.open(ROOT / "data/shared-rna/splicing/reconstruction/patient-windows.jsonl.gz", "rt") as stream:
        windows_by_id = {}
        for line in stream:
            record = json.loads(line)
            windows_by_id[record["window_id"]] = record

    results = []
    # precompute primary-reference nonself peptides per window (screen results)
    attribution = {r["window_id"]: r for r in csv.DictReader(
        (Path(__file__).resolve().parent / "splice-specificity-attribution.tsv").open(),
        delimiter="\t")}
    for row in selected:
        record = windows_by_id[row["window_id"]]
        match = re.match(r"^(chr[^:]+):(\d+)-(\d+)\(([+-])\)$", record["coordinates"])
        chrom = match.group(1)
        positions = sorted(set(record["genomic_positions_0based"]))
        blocks, current = [], [positions[0]]
        for previous, position in zip(positions, positions[1:]):
            if position - previous == 1:
                current.append(position)
            else:
                blocks.append(current)
                current = [position]
        blocks.append(current)
        block_lo = min(min(b) for b in blocks)
        block_hi = max(max(b) for b in blocks)

        # local genomic sequence covering both blocks plus margin (for variant application)
        span = fetch_region(f"{chrom}:{block_lo + 1 - MARGIN}-{block_hi + MARGIN}:1")
        span_offset = block_lo - MARGIN  # genomic 0-based offset of span[0]

        transcripts = gtf_cds_segments(row["gene"], chrom)
        all_targets = overlapping_peptides(row["sequence"], int(row["target_start_in_window"]),
                                           int(row["target_end_in_window"]))
        # only nonself peptides matter: a reference-self peptide occurring in an
        # annotated transcript is expected and not a patient-normal concern
        nonself_targets = [p for p in all_targets if p not in row["_nonself_exempt"]]
        hits = {}
        n_possibilities = 0
        n_variants = 0
        transcripts_checked = 0
        for transcript_id, info in transcripts.items():
            reverse = info["strand"] == "-"
            segments = info["segments"]
            # cumulative CDS offsets to derive the frame at each segment
            cumulative = 0
            segment_phase = []
            for start, end, frame in segments:
                segment_phase.append(cumulative % 3)
                cumulative += end - start
            # CDS segments overlapping the local region
            overlapping = [(start, end, phase)
                           for (start, end, _), phase in zip(segments, segment_phase)
                           if end > block_lo - MARGIN and start < block_hi + MARGIN]
            if not overlapping:
                continue
            transcripts_checked += 1
            # spliced local CDS in transcript orientation
            pieces, piece_positions = [], []
            for start, end, _phase in overlapping:
                piece = span[start - span_offset:end - span_offset]
                pieces.append(piece)
                piece_positions.extend(range(start, end))
            local_cds = "".join(pieces)
            if reverse:
                local_cds = revcomp(local_cds)
                piece_positions = list(reversed(piece_positions))
            # variants overlapping the local CDS bases
            local_variants = []
            for variant in gvcf_variants(chrom, block_lo - MARGIN, block_hi + MARGIN):
                for offset, genomic in enumerate(piece_positions):
                    if genomic == variant["pos"]:
                        local_variants.append(dict(variant, origin=offset - variant["pos"] + variant["pos"], pos=offset))
                        break
            n_variants = max(n_variants, len(local_variants))
            for possibility in local_possibilities(local_cds, local_variants):
                n_possibilities += 1
                aa = translate(possibility)
                for peptide in nonself_targets:
                    if peptide in aa:
                        hits.setdefault(peptide, []).append(transcript_id)

        results.append(dict(
            gene=row["gene"], window_id=row["window_id"], coordinates=record["coordinates"],
            annotated_transcripts_checked=transcripts_checked,
            local_variant_sites=n_variants,
            possibilities_evaluated=n_possibilities,
            junction_peptides_checked=len(all_targets),
            nonself_peptides_checked=len(nonself_targets),
            peptides_potential_normal_self=sorted(hits),
            hit_transcripts=sorted({t for ts in hits.values() for t in ts}),
            interpretation="POTENTIAL_PATIENT_SELF" if hits else "NO_HIT_IN_LOCAL_ANNOTATED_NORMAL_MODEL",
            proves_normal_tissue_absence=False))
        print(f"{row['gene']}: transcripts={transcripts_checked} possibilities={n_possibilities} "
              f"self-hits={len(hits)} -> {results[-1]['interpretation']}", flush=True)

    out = Path(__file__).resolve().parent / "normal-model-check.json"
    out.write_text(json.dumps(results, indent=2) + "\n")
    print("wrote", out)


if __name__ == "__main__":
    main()
