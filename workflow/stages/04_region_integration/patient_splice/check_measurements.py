"""Read-only coordinate-control and limited orientation measurements on Modal.

Executed explicitly through checks.py. No DNA calling, extra-substitution
assessment, deduplication or sequence replacement is performed here.
"""
import argparse
import ast
from collections import Counter, defaultdict
import datetime
import gzip
import hashlib
import json
from pathlib import Path
import re


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + "\n")


def write_gzip(path, value):
    # Stable compressed bytes permit an exact offline replay.
    with Path(path).open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as compressed:
            compressed.write((json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode())


def coordinate(value):
    match = re.fullmatch(r"([^:]+):(\d+)-(\d+)\(([+-])\)", value)
    if not match:
        return None
    chrom, a, b, strand = match.groups()
    return (chrom.removeprefix("chr"), int(a), int(b)), strand


def inputs(root, source):
    policy = read_json(source / "policy.json")
    request = read_json(root / "request.json")
    assert request["source_reconstruction_sha256"] == digest(source / "reconstruction/patient-windows.jsonl.gz")
    assert request["source_events_sha256"] == digest(source / "reconstruction/input-events.tsv")
    import csv
    events = list(csv.DictReader((source / "reconstruction/input-events.tsv").open(), delimiter="\t"))
    assert len(events) == len({r["event_id"] for r in events})
    return policy, request, events


def receipt(root, mode, source, paths, outputs, details):
    write_json(root / (mode + "-receipt.json"), dict(
        status="COMPLETED", mode=mode, completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        source_root=str(source), source_reconstruction_sha256=digest(source / "reconstruction/patient-windows.jsonl.gz"),
        source_events_sha256=digest(source / "reconstruction/input-events.tsv"),
        worker_sha256=digest(__file__), checks_policy_sha256=digest(root / "checks-policy.json"),
        request_sha256=digest(root / "request.json"),
        input_sha256={str(p): digest(p) for p in paths},
        output_sha256={name: digest(root / name) for name in outputs}, **details))


def h5_strings(dataset):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in dataset[:]]


def tissue_values(var, name):
    import h5py
    column = var[name]
    if isinstance(column, h5py.Group):
        categories, codes = h5_strings(column["categories"]), column["codes"][:]
    elif "__categories" in var and name in var["__categories"]:
        categories, codes = h5_strings(var["__categories"][name]), column[:]
    else:
        return h5_strings(column)
    return [categories[int(i)] if int(i) >= 0 else "UNKNOWN" for i in codes]


def normal(root, source):
    import h5py
    import numpy as np
    source_policy, request, events = inputs(root, source)
    checks_policy = read_json(root / "checks-policy.json")
    # Use exactly the source coordinate function and complete Alt91 exon table.
    # All genes are loaded: limiting this to the selected genes misses aliases.
    msnaf = source / "msnaf-core.8-11.py"
    syntax = ast.parse(msnaf.read_text())
    functions = [n for n in syntax.body if isinstance(n, ast.FunctionDef)
                 and n.name in {"exon_coords_to_dict", "uid_to_coord"}]
    assert len(functions) == 2
    namespace = {}
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(msnaf), "exec"), namespace)
    exon_path = Path(source_policy["msnaf_exons"])
    exons = namespace["exon_coords_to_dict"](str(exon_path))
    targets = defaultdict(list)
    for event in events:
        parsed = coordinate(event["coordinates"])
        if parsed:
            targets[parsed[0]].append(event["event_id"])
    base = exon_path.parent.parent / "controls"
    result = dict(schema_version=1, source_events=len(events),
        coordinate_scope="Every row in each supplied catalog is resolved across all genes with the pinned source adapter. Unresolvable rows are counted; raw normal reads are not reprocessed.",
        source_strand="Catalog annotation; not independent molecular strand evidence", cohorts={})
    paths = [source / "policy.json", source / "reconstruction/input-events.tsv", msnaf, exon_path]
    for config in checks_policy["normal"]["cohorts"]:
        path = base / config["filename"]
        paths.append(path)
        with h5py.File(path, "r") as handle:
            obs, var = handle["obs"], handle["var"]
            ik = obs.attrs.get("_index", "_index")
            vk = var.attrs.get("_index", "_index")
            ik = ik.decode() if isinstance(ik, bytes) else ik
            vk = vk.decode() if isinstance(vk, bytes) else vk
            n = len(var[vk])
            tissues = tissue_values(var, config["tissue_column"])
            assert len(tissues) == n
            sample_names = h5_strings(var[vk])
            matrix = handle["X"]
            encoding = matrix.attrs.get("encoding-type", "")
            if isinstance(encoding, bytes):
                encoding = encoding.decode()
            assert encoding == "csr_matrix", encoding
            totals, matched_rows, event_matches = Counter(), {}, defaultdict(list)
            index = obs[ik]
            for offset in range(0, len(index), 100000):
                for i, uid in enumerate(index.asstr()[offset:offset + 100000], offset):
                    totals["catalog_rows"] += 1
                    try:
                        resolved = namespace["uid_to_coord"](uid, exons)
                        parsed = coordinate(resolved)
                    except (KeyError, ValueError, IndexError, StopIteration):
                        parsed = None
                    if parsed is None:
                        totals["unresolved_coordinate_rows"] += 1
                        continue
                    totals["resolved_coordinate_rows"] += 1
                    key, strand = parsed
                    if key not in targets:
                        continue
                    start, end = (int(x) for x in matrix["indptr"][i:i + 2])
                    columns = matrix["indices"][start:end]
                    values = matrix["data"][start:end]
                    assert np.isfinite(values).all() and (values >= 0).all()
                    assert len(set(int(x) for x in columns)) == len(columns)
                    assert uid not in matched_rows
                    matched_rows[uid] = dict(coordinates=resolved, strand=strand,
                        source_row=i, nonzero_counts=[[int(j), float(v)] for j, v in zip(columns, values) if v > 0])
                    for event in targets[key]:
                        event_matches[event].append(uid)
                    totals["matched_catalog_rows"] += 1
                print("NORMAL_SCAN", config["name"], offset + len(index[offset:offset + 100000]), "/", len(index), flush=True)
            result["cohorts"][config["name"]] = dict(
                source=str(path), tissue_column=config["tissue_column"], samples=n,
                sample_names=sample_names, tissues=tissues, catalog_counts=dict(totals),
                rows=matched_rows, event_matches=dict(event_matches))
        print("NORMAL_COMPLETE", config["name"], dict(totals), flush=True)
    write_gzip(root / "normal-counts.json.gz", result)
    receipt(root, "normal", source, paths, ["normal-counts.json.gz"],
            dict(events_measured=len(events), versions=dict(h5py=h5py.__version__, numpy=np.__version__)))


def eligible(read, min_mapq):
    return (not any((read.is_unmapped, read.is_secondary, read.is_supplementary,
                     read.is_duplicate, read.is_qcfail)) and read.mapping_quality >= min_mapq
            and (not read.has_tag("NH") or read.get_tag("NH") == 1))


def read1_orientation(read, strand):
    if (not read.is_paired or read.mate_is_unmapped or read.is_read1 == read.is_read2
            or read.reference_id != read.next_reference_id):
        return "unevaluable"
    # FR paired reads should imply the same read-1 direction independently.
    if read.is_reverse == read.mate_is_reverse:
        return "unevaluable"
    reverse = read.is_reverse if read.is_read1 else not read.is_reverse
    return "sense" if reverse == (strand == "-") else "antisense"


def has_boundary(read, low, high):
    if high == low + 1:
        positions = {p: q for q, p in read.get_aligned_pairs(matches_only=True) if p is not None and q is not None}
        return low - 1 in positions and high - 1 in positions and positions[high - 1] == positions[low - 1] + 1
    pos = read.reference_start
    for operation, length in read.cigartuples or ():
        if operation == 3 and pos == low and pos + length == high - 1:
            return True
        if operation in (0, 2, 3, 7, 8):
            pos += length
    return False


def strand(root, source):
    import pysam
    source_policy, request, events = inputs(root, source)
    checks = read_json(root / "checks-policy.json")["strand"]
    bam_path, gtf_path, genome_path = [Path(source_policy[k]) for k in ("rna_bam", "gtf", "genome")]
    rna, genome = pysam.AlignmentFile(str(bam_path), "rb"), pysam.FastaFile(str(genome_path))
    transcripts = defaultdict(list)
    for line in gtf_path.open():
        c = line.rstrip().split("\t")
        if len(c) != 9 or c[2] != "CDS" or c[0] not in rna.references:
            continue
        attrs = dict(re.findall(r'(\w+) "([^\"]*)"', c[8]))
        if "transcript_id" in attrs and "gene_id" in attrs:
            transcripts[attrs["transcript_id"]].append((c[0], int(c[3]) - 1, int(c[4]), c[6], attrs["gene_id"]))
    junctions = defaultdict(set)
    for parts in transcripts.values():
        parts.sort(key=lambda x: x[1])
        for a, b in zip(parts, parts[1:]):
            if a[0] == b[0] and b[1] - a[2] >= 80:
                junctions[(a[0], a[2], b[1])].add((a[3], a[4]))
    counts, names, genes, intervals, library_witnesses = Counter(), set(), set(), 0, []
    for (chrom, low, high), sources in sorted(junctions.items()):
        directions, ids = {s for s, _ in sources}, {g for _, g in sources}
        if len(directions) != 1 or ids & genes:
            continue
        intervals += 1
        direction, added = next(iter(directions)), 0
        for read in rna.fetch(chrom, low - 1, low + 1):
            if not eligible(read, source_policy["min_mapq"]) or read.query_name in names:
                continue
            if not has_boundary(read, low, high + 1):
                continue
            orientation = read1_orientation(read, direction)
            if orientation == "unevaluable":
                continue
            counts[orientation] += 1
            names.add(read.query_name)
            library_witnesses.append(dict(fragment=read.query_name, genes=sorted(ids),
                coordinates=[chrom, low, high, direction], read1_orientation=orientation))
            added += 1
            if added >= checks["library_sample_max_fragments_per_gene"]:
                break
        if added:
            genes.update(ids)
        if len(genes) >= checks["library_sample_max_genes"] or intervals >= checks["library_sample_max_intervals"]:
            break
    windows, by_event = {}, defaultdict(list)
    with gzip.open(source / "reconstruction/patient-windows.jsonl.gz", "rt") as handle:
        for line in handle:
            window = json.loads(line)
            windows[window["window_id"]] = window
            by_event[window["event_id"]].append(window)
    event_by_id = {x["event_id"]: x for x in events}
    result = dict(schema_version=1, library=dict(counts=dict(counts), fragments=len(names),
        genes=len(genes), intervals_examined=intervals, witnesses=library_witnesses), events={}, windows={})
    canonical = {"+": {"GT-AG", "GC-AG", "AT-AC"}, "-": {"CT-AC", "CT-GC", "GT-AT"}}
    for number, (event, event_windows) in enumerate(sorted(by_event.items()), 1):
        record = event_by_id[event]
        parsed = coordinate(record["coordinates"])
        assert parsed
        (chrom, low, high), direction = parsed
        chrom = record["coordinates"].split(":")[0]
        observations, xs = defaultdict(set), Counter()
        for read in rna.fetch(chrom, max(0, low - 201), low + 200):
            if not eligible(read, source_policy["min_mapq"]) or not has_boundary(read, low, high):
                continue
            observations[read.query_name].add(read1_orientation(read, direction))
            if read.has_tag("XS"):
                xs[str(read.get_tag("XS"))] += 1
        fragments = {name: next(iter(values)) if len(values) == 1 else "conflicting"
                     for name, values in observations.items()}
        motif = None
        if high > low + 1:
            motif = genome.fetch(chrom, low, low + 2).upper() + "-" + genome.fetch(chrom, high - 3, high - 1).upper()
        result["events"][event] = dict(coordinates=record["coordinates"], gene_strand=direction,
            junction_fragment_orientations=dict(Counter(fragments.values())),
            genomic_plus_motif=motif, motif_agrees_with_gene=motif in canonical[direction] if motif else None,
            XS_tags=dict(xs), motif_and_XS_are_report_only=True)
        for window in event_windows:
            names = sorted(window["fragment_names"])
            witness_observations = {name: fragments.get(name, "not_recovered") for name in names}
            result["windows"][window["window_id"]] = dict(event_id=event,
                counts=dict(Counter(witness_observations.values())), fragment_orientations=witness_observations)
        if number % 100 == 0:
            print("STRAND_PROGRESS", number, "/", len(by_event), flush=True)
    assert set(result["windows"]) == set(windows)
    write_gzip(root / "strand-counts.json.gz", result)
    receipt(root, "strand", source, [source / "policy.json", bam_path, gtf_path, genome_path],
            ["strand-counts.json.gz"], dict(events_measured=len(by_event), windows_measured=len(windows),
                                         versions=dict(pysam=pysam.__version__)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["normal", "strand"])
    parser.add_argument("root", type=Path)
    parser.add_argument("source", type=Path)
    args = parser.parse_args()
    assert args.root.resolve() != args.source.resolve(), "Original source evidence must remain unchanged"
    {"normal": normal, "strand": strand}[args.mode](args.root, args.source)


if __name__ == "__main__":
    main()
