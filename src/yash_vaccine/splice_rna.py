"""Patient-sequence splice windows, without filling gaps from a reference.

Coordinates within a fragment are query-base offsets from the RNA join. Mates
may extend one another only when their aligned paths overlap consistently.
Different fragments are never combined into a synthetic consensus haplotype.
"""
from collections import defaultdict
import hashlib

from Bio.Seq import Seq

AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def eligible(read, min_mapq=20):
    return (not (read.is_unmapped or read.is_secondary or read.is_supplementary
                 or read.is_duplicate or read.is_qcfail)
            and read.mapping_quality >= min_mapq
            and (not read.has_tag("NH") or read.get_tag("NH") == 1))


def aligned_path(read):
    """Query bases and reference positions, including query insertions.

    Soft-clipped bases are excluded. An N advances the reference, not the query.
    Returning actual query bases avoids the reference-codon error being fixed.
    """
    query, reference = 0, read.reference_start
    path, joins = {}, []
    for op, length in read.cigartuples or ():
        if op in (0, 7, 8):
            for i in range(length):
                path[query + i] = (reference + i, read.query_sequence[query + i],
                                   int(read.query_qualities[query + i])
                                   if read.query_qualities is not None else 0)
            query += length
            reference += length
        elif op == 1:
            for i in range(length):
                path[query + i] = (None, read.query_sequence[query + i],
                                   int(read.query_qualities[query + i])
                                   if read.query_qualities is not None else 0)
            query += length
        elif op == 2:
            reference += length
        elif op == 3:
            joins.append((reference, reference + length, query))
            reference += length
        elif op == 4:
            query += length
        elif op not in (5, 6):
            raise ValueError(f"Unsupported CIGAR operation {op}")
    return path, joins


def orient_path(path, query_join, strand):
    if strand == "+":
        return {q - query_join: value for q, value in path.items()}
    return {query_join - 1 - q: (p, base.translate(COMPLEMENT), quality)
            for q, (p, base, quality) in path.items()}


def cds_index(parts, position, strand):
    offset = -parts[0][2]
    for start, end, _ in parts:
        if start <= position <= end:
            return offset + (position - start if strand == "+" else end - position)
        offset += end - start + 1
    return None


def anchored_frames(path, transcript_models, strand, min_anchor=10, index_cache=None):
    """Anchor translation in a read-observed upstream annotated CDS block.

    Use the earliest upstream CDS block within 75 nt of the join, rather than
    restarting the frame after a local read insertion/deletion. All supported
    transcript phases are retained. This is a local frame, not a full-ORF claim.
    """
    result = defaultdict(list)
    for transcript, parts in transcript_models.items():
        blocks, block = [], []
        for relative in range(-75, 0):
            value = path.get(relative)
            index = None
            if value and value[0] is not None and value[2] >= 20:
                key = (transcript, value[0])
                if index_cache is not None and key in index_cache:
                    index = index_cache[key]
                else:
                    index = cds_index(parts, value[0] + 1, strand)
                    if index_cache is not None:
                        index_cache[key] = index
            if index is None or index < 0:
                if block:
                    blocks.append(block)
                    block = []
                continue
            if block and (relative != block[-1][0] + 1 or index != block[-1][1] + 1):
                blocks.append(block)
                block = []
            block.append((relative, index))
        if block:
            blocks.append(block)
        blocks = [b for b in blocks if len(b) >= min_anchor]
        if not blocks:
            continue
        anchor = blocks[0][0]
        first_codon = anchor[0] + (-anchor[1]) % 3
        end_codon = first_codon + 3 * ((-first_codon) // 3)
        prefix_positions = range(first_codon, end_codon)
        if not all(q in path and path[q][2] >= 20 and path[q][1] in "ACGT"
                   for q in prefix_positions):
            continue
        prefix = "".join(path[q][1] for q in prefix_positions)
        if any(prefix[i:i+3] in {"TAA", "TAG", "TGA"} for i in range(0, len(prefix), 3)):
            continue
        # CDS index at the first base downstream of the join, modulo three.
        phase = (anchor[1] - anchor[0]) % 3
        result[phase].append(dict(transcript_id=transcript,
                                 upstream_query_offset=anchor[0],
                                 upstream_CDS_index_0based=anchor[1],
                                 anchor_bases=len(blocks[0]),
                                 first_complete_codon_query_offset=first_codon,
                                 upstream_coding_sequence=prefix,
                                 upstream_coding_positions_0based=[path[q][0] for q in prefix_positions]))
    return dict(result)


def extend_with_mate(join_path, mate_path, strand, min_overlap=10):
    """Project an overlapping mate onto the join path, or leave it unextended."""
    reference_to_relative = {p: relative for relative, (p, _, _) in join_path.items()
                             if p is not None}
    direction = 1 if strand == "+" else -1
    offsets = [reference_to_relative[p] - direction * q
               for q, (p, _, _) in mate_path.items() if p in reference_to_relative]
    if len(offsets) < min_overlap or len(set(offsets)) != 1:
        return dict(join_path), "NO_CONSISTENT_MATE_OVERLAP"
    offset = offsets[0]
    projected = {direction * q + offset:
                 (p, base if strand == "+" else base.translate(COMPLEMENT), quality)
                 for q, (p, base, quality) in mate_path.items()}
    merged = dict(join_path)
    for relative, value in projected.items():
        if relative in merged:
            old = merged[relative]
            if old[0] != value[0]:
                return dict(join_path), "MATE_PATH_CONFLICT"
            if old[2] >= 20 and value[2] >= 20 and old[1] != value[1]:
                # A conflict becomes a hole, never an arbitrarily chosen base.
                merged[relative] = (old[0], "N", 0)
                continue
            if value[2] <= old[2]:
                continue
        merged[relative] = value
    return merged, "CONSISTENT_OVERLAPPING_MATE"


def patient_windows(path, phase, length=25, min_baseq=20):
    """All fully observed, phase-anchored windows crossing this RNA join."""
    for start in range(-3 * length + 1, 0):
        end = start + 3 * length
        if (start + phase) % 3 or end <= 0:
            continue
        positions = range(start, end)
        if not all(p in path and path[p][2] >= min_baseq and path[p][1] in "ACGT"
                   for p in positions):
            continue
        dna = "".join(path[p][1] for p in positions)
        sequence = str(Seq(dna).translate())
        if len(sequence) != length or not set(sequence) <= AMINO_ACIDS:
            continue
        boundary = -start
        target_start = boundary // 3
        target_end = target_start + bool(boundary % 3)
        yield dict(sequence=sequence, coding_sequence=dna,
                   start_nt_relative_to_join=start, junction_nt_offset=boundary,
                   target_start_in_window=target_start,
                   target_end_in_window=target_end,
                   target_mode="junction_codon" if boundary % 3 else "junction_boundary",
                   genomic_positions_0based=[path[p][0] for p in positions])


def fragment_agrees_with_window(read_paths, window, strand, min_baseq=20):
    """Reject contradictory high-quality mate bases within this local path.

    A mate can disagree in alignment outside the window, preventing whole-read
    extension, while still exposing a conflicting base inside the window. Such
    a fragment must not be counted by choosing just the matching read.
    """
    positions = window["genomic_positions_0based"]
    indexes = {p: i for i, p in enumerate(positions) if p is not None}
    direction = 1 if strand == "+" else -1
    for path in read_paths:
        offsets = {indexes[p] - direction * q for q, (p, _, _) in path.items() if p in indexes}
        placements = []
        for offset in offsets:
            covered = [(q, direction * q + offset) for q in path
                       if 0 <= direction * q + offset < len(positions)]
            if covered and all(path[q][0] == positions[i] for q, i in covered):
                placements.append(covered)
        if len(placements) != 1:
            continue
        for q, i in placements[0]:
            _, base, quality = path[q]
            if strand == "-":
                base = base.translate(COMPLEMENT)
            if quality >= min_baseq and base != window["coding_sequence"][i]:
                return False
    return True


def target_peptides(window):
    sequence = window["sequence"]
    left, right = window["target_start_in_window"], window["target_end_in_window"]
    for length in (8, 9, 10, 11):
        for start in range(len(sequence) - length + 1):
            if start < right and start + length > left:
                yield dict(peptide=sequence[start:start + length], start=start,
                           end=start + length)


def stable_id(event, phase, start, coding_sequence):
    key = f"{event}|{phase}|{start}|{coding_sequence}"
    return "RNA_SPLICE_" + hashlib.sha256(key.encode()).hexdigest()[:20]
