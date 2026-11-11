"""Reconstruct the complete normal-filtered LENS splice pool from patient RNA.

Run only in the recorded mounted environment through run.py --execute. Native
LENS event review is retained; patient sequences replace its reference-only
peptide consistency check. No old HLA prediction is used to nominate events.
"""
import ast
from collections import Counter, defaultdict
import csv
import datetime
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import sys

from Bio import SeqIO
from Bio.Seq import Seq
import pysam

from splice_rna import (aligned_path, anchored_frames, eligible,
                        extend_with_mate, fragment_agrees_with_window,
                        orient_path, patient_windows, stable_id)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def table(path, records):
    records = list(records)
    fields = list(dict.fromkeys(k for row in records for k in row))
    with Path(path).open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)


def save(path, data):
    Path(path).write_text(json.dumps(data, indent=2) + "\n")


def main(root, shard=0, shards=1):
    policy = json.loads((root / "policy.json").read_text())
    lens = Path(policy["lens_root"])
    bam_path, gtf_path = Path(policy["rna_bam"]), Path(policy["gtf"])
    genome_path = Path(policy["genome"])
    out = root / "reconstruction" if shards == 1 else root / "reconstruction-shards" / str(shard)
    out.mkdir(parents=True, exist_ok=True)
    # Use the two original MSNAF coordinate functions without loading its model
    # dependencies. Their source file and all inputs are hashed in the receipt.
    source = root / "msnaf-core.8-11.py"
    syntax = ast.parse(source.read_text())
    namespace = {}
    functions = [node for node in syntax.body if isinstance(node, ast.FunctionDef)
                 and node.name in {"exon_coords_to_dict", "uid_to_coord"}]
    assert len(functions) == 2
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(source), "exec"), namespace)
    exon_table = Path(policy["msnaf_exons"])
    exon_coordinates = namespace["exon_coords_to_dict"](str(exon_table))
    with (lens / "splice/NeoJunction_statistics_maxmin.txt").open() as handle:
        normal_rows = [row for row in csv.DictReader(handle, delimiter="\t")
                       if row["cond"] == "True" and row["cond_add_tcga_control"] == "True"]
    events = []
    for row in normal_rows:
        event = row[""]
        events.append(dict(event_id=event, gene_id=event.split(":")[0],
                           coordinates=namespace["uid_to_coord"](event, exon_coordinates),
                           **{k: v for k, v in row.items() if k}))
    assert len(events) == len({x["event_id"] for x in events})
    # Independently check the coordinate adapter against every existing MSNAF
    # translated event, including events that never passed an HLA screen.
    recorded = {r["splice_description"]: r["splice_coords"] for r in
                csv.DictReader((lens / "splice/junction-review-candidates.tsv").open(), delimiter="\t")}
    assert all(r["coordinates"] == recorded[r["event_id"]]
               for r in events if r["event_id"] in recorded)
    assert set(recorded) <= {r["event_id"] for r in events}
    total_events = len(events)
    events = [event for i, event in enumerate(events) if i % shards == shard]
    table(out / "input-events.tsv", events)
    wanted = {r["gene_id"] for r in events}
    tx_parts, tx_exons, by_gene, tx_meta, gene_names = (defaultdict(list), defaultdict(list),
                                                      defaultdict(dict), {}, {})
    with gtf_path.open() as handle:
        for line in handle:
            columns = line.rstrip().split("\t")
            if len(columns) != 9 or columns[2] not in ("CDS", "exon", "gene"):
                continue
            attrs = dict(re.findall(r'(\w+) "([^\"]*)"', columns[8]))
            gene = attrs.get("gene_id", "").split(".")[0]
            if gene not in wanted:
                continue
            gene_names[gene] = attrs.get("gene_name", gene)
            tx = attrs.get("transcript_id")
            if not tx:
                continue
            tx_meta[tx] = (gene, columns[0], columns[6])
            if columns[2] == "CDS":
                tx_parts[tx].append((int(columns[3]), int(columns[4]), int(columns[7])))
            elif columns[2] == "exon":
                tx_exons[tx].append((int(columns[3]), int(columns[4])))
    proteins = {}
    with gzip.open(root / "gencode37.proteins.fasta.gz", "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            parts = record.id.split("|")
            if parts[2].split(".")[0] in wanted:
                proteins[parts[1]] = str(record.seq)
    # Match the native LENS junction dictionary, stripping only gene versions.
    junctions = {"gene_junctions": defaultdict(set), "tx_junctions": defaultdict(set)}
    for tx, parts in tx_exons.items():
        gene, chrom, strand = tx_meta[tx]
        parts = sorted(parts)
        for left, right in zip(parts, parts[1:]):
            if left[1] + 1 > right[0] - 1:
                continue
            key = (chrom, left[1] + 1, right[0] - 1, strand)
            junctions["gene_junctions"][gene].add(key)
            junctions["tx_junctions"][tx].add(key)
    genome = pysam.FastaFile(str(genome_path))
    transcript_audit = []
    for tx, parts in tx_parts.items():
        if tx not in proteins:
            continue
        gene, chrom, strand = tx_meta[tx]
        parts = sorted(parts, reverse=strand == "-")
        pieces = [genome.fetch(chrom, a - 1, b).upper() for a, b, _ in parts]
        if strand == "-":
            pieces = [str(Seq(p).reverse_complement()) for p in pieces]
        dna = "".join(pieces)[parts[0][2]:]
        protein = str(Seq(dna[:len(dna) // 3 * 3]).translate()).rstrip("*")
        match = protein == proteins[tx].rstrip("*")
        transcript_audit.append(dict(transcript_id=tx, gene_id=gene,
                                     reference_translation_exact=match))
        if match:
            by_gene[(gene, chrom, strand)][tx] = parts
    table(out / "reference-transcripts.tsv", transcript_audit)
    native_source = lens / "metadata.splice-review.py"
    sys.path.insert(0, str(lens / "native-source/lenstools"))
    spec = importlib.util.spec_from_file_location("lens_splice_review", native_source)
    native = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(native)
    os.environ["LENS_MSNAF_EXON_BOUNDS"] = "1"
    cfg_path = lens / "review_scoring_defaults.expression-aware.yaml"
    cfg = native._load_review_config(str(cfg_path))
    rna = pysam.AlignmentFile(str(bam_path), "rb")
    expression = {}
    for line in Path(policy["expression"]).read_text().splitlines():
        if line and not line.startswith("#"):
            gene, value = line.split("\t")[:2]
            expression[gene] = float(value)
    audit, windows, evidence, counts = [], [], {}, Counter()
    for number, event in enumerate(events):
        if number % 50 == 0:
            print("PROGRESS", shard, number, "/", len(events), "windows", len(windows), flush=True)
        uid, gene = event["event_id"], event["gene_id"]
        match = re.fullmatch(r"([^:]+):(\d+)-(\d+)\(([+-])\)", event["coordinates"])
        base = dict(event, gene=gene_names.get(gene, gene))
        if not match:
            audit.append(dict(base, status="COORDINATES_UNRESOLVED"))
            continue
        chrom, low, high, strand = match[1], int(match[2]), int(match[3]), match[4]
        models = by_gene.get((gene, chrom, strand), {})
        if not models:
            audit.append(dict(base, status="NO_VERIFIED_ANNOTATED_CDS"))
            continue
        # Fetch only boundary neighborhoods. Reads crossing a long intron are
        # still found without scanning every intervening genomic base.
        alignments = {}
        for a, b in [(max(0, low - 201), low + 200), (max(0, high - 201), high + 200)]:
            for read in rna.fetch(chrom, a, b):
                if not eligible(read, policy["min_mapq"]):
                    continue
                key = (read.query_name, read.flag, read.reference_start, read.cigarstring)
                alignments[key] = read
        fragments, paths, seeds = defaultdict(list), {}, []
        for key, read in alignments.items():
            path, joins = aligned_path(read)
            paths[key] = path
            fragments[read.query_name].append(key)
            if high == low + 1:
                pos = {p: q for q, (p, _, _) in path.items() if p is not None}
                if low - 1 in pos and high - 1 in pos and pos[high - 1] == pos[low - 1] + 1:
                    joins.append((low, high - 1, pos[high - 1]))
            for first, last, query_join in joins:
                if (first, last) == (low, high - 1):
                    seeds.append((key, query_join))
        sequence_support, representatives = defaultdict(set), {}
        path_support = defaultdict(lambda: defaultdict(set))
        path_anchors = defaultdict(lambda: defaultdict(dict))
        frame_support, mate_outcomes, index_cache = defaultdict(dict), Counter(), {}
        conflicting_fragment_windows = set()
        for seed_number, (key, query_join) in enumerate(seeds):
            read = alignments[key]
            path = orient_path(paths[key], query_join, strand)
            for other in fragments[read.query_name]:
                if other != key:
                    path, outcome = extend_with_mate(path, paths[other], strand,
                                                    policy["min_mate_overlap_nt"])
                    mate_outcomes[outcome] += 1
            frames = anchored_frames(path, models, strand, policy["min_upstream_CDS_anchor_nt"], index_cache)
            if seed_number == 0:
                assert frames == anchored_frames(path, models, strand, policy["min_upstream_CDS_anchor_nt"])
            anchor_agreement = {}
            consistent_frames = {}
            for phase, anchors in frames.items():
                kept = []
                for anchor in anchors:
                    key = (anchor["upstream_coding_sequence"], tuple(anchor["upstream_coding_positions_0based"]))
                    if key not in anchor_agreement:
                        anchor_agreement[key] = fragment_agrees_with_window(
                            [paths[k] for k in fragments[read.query_name]],
                            dict(coding_sequence=key[0], genomic_positions_0based=key[1]),
                            strand, policy["min_baseq"])
                    if anchor_agreement[key]:
                        kept.append(anchor)
                if kept:
                    consistent_frames[phase] = kept
            frames = consistent_frames
            for phase, anchors in frames.items():
                for window in patient_windows(path, phase, policy["region_length"], policy["min_baseq"]):
                    wid = stable_id(uid, phase, window["start_nt_relative_to_join"], window["coding_sequence"])
                    if not fragment_agrees_with_window(
                            [paths[k] for k in fragments[read.query_name]], window, strand, policy["min_baseq"]):
                        conflicting_fragment_windows.add((wid, read.query_name))
                        continue
                    sequence_support[wid].add(read.query_name)
                    representatives[wid] = dict(window, window_id=wid, context_phase=phase)
                    signature = tuple(window["genomic_positions_0based"])
                    path_support[wid][signature].add(read.query_name)
                    path_anchors[wid][signature][read.query_name] = anchors
                    for anchor in anchors:
                        frame_support[wid][anchor["transcript_id"]] = anchor
        supported = {wid: names for wid, names in sequence_support.items()
                     if len(names) >= policy["min_full_window_fragments"]}
        rec = dict(base, junction_fragment_names=len({alignments[k].query_name for k, _ in seeds}),
                   candidate_haplotype_windows=len(sequence_support),
                   supported_haplotype_windows=len(supported),
                   ambiguous_fragment_window_pairs=len(conflicting_fragment_windows),
                   mate_extension_counts=json.dumps(dict(mate_outcomes), sort_keys=True))
        if not supported:
            audit.append(dict(rec, status="NO_FULL_PATIENT_WINDOW_WITH_MINIMUM_SUPPORT"))
            continue
        # Preserve native event review. Do not compare patient DNA bases to a
        # reference-only peptide model: that was the circular check being fixed.
        sample, event_evidence = native._extract_splice_evidence(
            rna, chrom, low, high, strand, junctions, uid, cfg["extract"]["min_mapq"],
            cfg["windows"]["splice_window_bp"], cfg["extract"]["count_mode"], "")
        normalized = event_evidence["normalized_junction"]
        motif, valid = native._junction_motif(genome, chrom, normalized["start"], normalized["end"], strand)
        event_evidence["junction"] = dict(motif=motif, motif_valid=valid, strand=strand)
        ev = dict(tumor_rna=sample, event_specific=event_evidence)
        flags, hard_fail = native._flags_from_evidence(dict(source_type="SPLICE", peptide=""), ev, cfg)
        score = max(0, 100 - sum(cfg["scoring"]["penalties"].get(f, 0) for f in flags))
        tier = ("FAIL" if hard_fail or score < cfg["scoring"]["tier_thresholds"]["review_min"] else
                "REVIEW" if score < cfg["scoring"]["tier_thresholds"]["pass_min"] else "PASS")
        evidence[uid] = dict(ev, flags=flags, score=score, tier=tier)
        rec.update(native_review_tier=tier, native_review_score=score,
                   native_review_flags=";".join(flags), motif=motif)
        for wid, names in sorted(supported.items()):
            window = representatives[wid]
            witnesses = []
            for signature, fragment_names in sorted(path_support[wid].items(), key=lambda item: str(item[0])):
                anchor_groups = {}
                for name in sorted(fragment_names):
                    for anchor in path_anchors[wid][signature][name]:
                        key = json.dumps(anchor, sort_keys=True)
                        if key not in anchor_groups:
                            anchor_groups[key] = dict(anchor, fragment_names=[])
                        anchor_groups[key]["fragment_names"].append(name)
                witnesses.append(dict(genomic_positions_0based=list(signature),
                                      fragment_names=sorted(fragment_names),
                                      upstream_frame_anchors=list(anchor_groups.values())))
            assert set().union(*(set(w["fragment_names"]) for w in witnesses)) == names
            window["genomic_positions_0based"] = witnesses[0]["genomic_positions_0based"]
            windows.append(dict(window, source="SPLICE", event_id=uid, gene_id=gene,
                                gene=base["gene"], coordinates=event["coordinates"],
                                fragment_count=len(names), fragment_names=sorted(names),
                                genomic_path_support=witnesses, genomic_path_count=len(witnesses),
                                upstream_frame_anchors=list(frame_support[wid].values()),
                                supporting_transcripts=";".join(sorted(frame_support[wid])),
                                supplied_transformed_gene_tpm=expression.get(base["gene"]),
                                native_event_review_tier=tier, native_event_review_score=score,
                                native_event_review_flags=flags, motif=motif))
        audit.append(dict(rec, status="PATIENT_WINDOWS_RECONSTRUCTED"))
        if number % 50 == 0 or len(windows) < 100:
            print("EVENT", number, "/", len(events), base["gene"], "windows", len(windows), flush=True)
            save(out / "status.json", dict(stage="reconstructing", processed_events=number + 1,
                                            input_events=len(events), windows=len(windows)))
    for record in audit:
        counts[record["status"]] += 1
    with gzip.open(out / "patient-windows.jsonl.gz", "wt") as handle:
        for window in windows:
            handle.write(json.dumps(window) + "\n")
    with gzip.open(out / "native-event-evidence.json.gz", "wt") as handle:
        json.dump(evidence, handle)
    table(out / "event-audit.tsv", audit)
    receipt = dict(status="COMPLETED", completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                   policy=policy, shard_index=shard, shards=shards, total_pool_events=total_events,
                   event_counts=dict(counts), input_events=len(events),
                   patient_windows=len(windows),
                   patient_window_events=len({w["event_id"] for w in windows}),
                   native_PASS_window_events=len({w["event_id"] for w in windows if w["native_event_review_tier"] == "PASS"}),
                   old_MSNAF_coordinates_checked=len(recorded),
                   verified_transcripts=sum(x["reference_translation_exact"] for x in transcript_audit),
                   input_sha256={str(p): digest(p) for p in [bam_path, gtf_path, genome_path, exon_table,
                       cfg_path, native_source, source, root / "splice_rna.py", root / "reconstruct.py",
                       root / "policy.json", root / "gencode37.proteins.fasta.gz",
                       lens / "splice/NeoJunction_statistics_maxmin.txt", Path(policy["expression"])]})
    save(out / "receipt.json", receipt)
    save(out / "status.json", dict(stage="reconstruction_complete", **receipt))
    print(json.dumps({k: v for k, v in receipt.items() if k not in {"input_sha256", "policy"}}, indent=2), flush=True)


if __name__ == "__main__":
    main(Path(sys.argv[1]), int(sys.argv[2]) if len(sys.argv) > 2 else 0,
         int(sys.argv[3]) if len(sys.argv) > 3 else 1)
