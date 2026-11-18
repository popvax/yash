"""Read-only artifact verification. Passing is NOT clinical or biological validation."""
import argparse
import gzip
import hashlib
import importlib.util
import json
import math
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "pipeline"))

import RNA
from Bio import SeqIO
from Bio.Seq import Seq

from curate import core_regions, derive, strand_strict_core_regions
from dataio import cell, iter_rows, load, rows, sha
from hydrophobicity import protein_screen, summarize
from protein import (
    CLASS2_ALLELES,
    HERE,
    assemble,
    assembled_class2_context,
    assess,
    assess_assembled_class2,
    junction_class2_counts,
    junction_peptides,
    target_candidates,
)
from replay_helpers import canonical_path, objective_values, read_edges, read_triplets
from alphagenome2_audit import verify_offline as verify_alphagenome2_offline
from safety_audit import verify_offline as verify_safety_offline
from strand_audit import verify_offline as verify_strand_offline
from build_review import results_summary_snapshot


def require(condition, message):
    if not condition:
        raise ValueError(message)


def verify_review_summary(target_root=HERE):
    """Reject stale review snapshots even when their manifest hash is current."""
    target_root = Path(target_root)
    recorded = load(target_root / "review/results-summary.json")
    expected = results_summary_snapshot(target_root)
    require(set(recorded) == set(expected), "Review results-summary fields differ from the current generator")

    def canonical(value):
        # Preserve JSON types: True must not equal 1, nor an integer a float.
        return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)

    for key, source_snapshot in expected.items():
        require(
            canonical(recorded[key]) == canonical(source_snapshot),
            "Review results-summary is stale for " + key
            + "; regenerate it from the current source artifacts with pipeline/build_review.py",
        )
    return {"current_source_snapshots_match": True, "top_level_fields_checked": len(expected)}


def mrna_module():
    directory = HERE / "construct/mrna"
    sys.path.insert(0, str(directory))
    spec = importlib.util.spec_from_file_location("verified_mrna_builder", directory / "run.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def verify(rescan_reference=False):
    require((HERE / "manifest.json").is_file(),
            "Missing bundle manifest; do not silently skip narrative and artifact integrity")
    manifest = load(HERE / "manifest.json")
    from package_review import included_files
    require(set(manifest["files_sha256"]) == {
        str(path.relative_to(HERE)) for path in included_files()
    }, "Manifest does not cover the complete published bundle")
    for relative, expected in manifest["files_sha256"].items():
        path = (HERE / relative).resolve()
        require(path.is_relative_to(HERE), "Manifest path leaves the bundle")
        require(sha(path) == expected, "Bundle checksum changed: " + relative)
    source = load(HERE / "exports/source-receipt.json")
    for relative, entry in source["files"].items():
        if not entry["adapted_after_snapshot"]:
            require(sha(HERE / relative) == entry["source_sha256"], "Frozen source input changed: " + relative)
    d3_ref = load(HERE / "exports/design3-reference-receipt.json")
    for relative, entry in d3_ref["files"].items():
        require(sha(HERE / relative) == entry["sha256"], "Frozen Design 3 reference input changed: " + relative)
        upstream_path = HERE.parents[1] / entry["upstream_path"]
        if upstream_path.is_file():
            require(sha(upstream_path) == entry["sha256"],
                    "Frozen Design 3 reference input diverges from repository designs/design-3/: " + relative)
    selected = load(HERE / "selection/selected.json")
    derived = core_regions(derive()[0])
    require(selected == derived, "The core no longer derives from the frozen source evidence")
    selected_table = rows(HERE / "selection/selected.tsv")
    for before, after in zip(selected, selected_table):
        require(all(cell(before.get(key)) == value for key, value in after.items()),
                "Selected TSV differs from the structured source records")
    require(len(selected_table) == len(selected) == 16, "Expected the declared 16-region core")
    require(Counter(row["source"] for row in selected) == {"MUTATION": 15, "FUSION": 1}, "Core source counts changed")
    require(all(row["clinical_release_ready"] is False for row in selected), "An unvalidated antigen was marked clinically ready")
    require(len(load(HERE / "selection/reserve.json")) == 8, "Reserve panel changed")

    from dataio import resolve_gencode_fasta
    receipt = load(HERE / "references/reference-receipt.json")
    require(sha(HERE / "references/peptide-matches.json") == receipt["peptide_match_sha256"], "Exact reference matches changed")
    latest_dir = HERE / "references/v50-audit"
    latest = load(latest_dir / "reference-receipt.json")
    require(sha(latest_dir / "peptide-matches.json") == latest["peptide_match_sha256"], "Latest reference matches changed")
    provider = load(latest_dir / "upstream-checksum-receipt.json")
    checksum_file = latest_dir / "upstream-MD5SUMS.txt"
    require(sha(checksum_file) == provider["provider_checksum_file_sha256"],
            "Provider checksum snapshot changed")
    line = next(
        line for line in checksum_file.read_text().splitlines()
        if line.split()[-1].lstrip("*") == provider["filename"]
    )
    require(line.split()[0] == provider["provider_published_md5"] == provider["local_md5"]
            and provider["local_sha256"] == latest["fasta_sha256"],
            "Provider checksum receipt disagrees with the recorded MD5SUMS or v50 receipt")
    v46_path = resolve_gencode_fasta(receipt["fasta_file"], HERE / "references")
    v50_path = resolve_gencode_fasta(latest["fasta_file"], latest_dir)
    if v46_path.is_file():
        require(sha(v46_path) == receipt["fasta_sha256"], "Reference proteome changed")
    if v50_path.is_file():
        require(sha(v50_path) == latest["fasta_sha256"], "Latest reference proteome changed")
        require(hashlib.md5(v50_path.read_bytes()).hexdigest() == provider["provider_published_md5"],
                "Latest reference does not match the recorded provider transfer checksum")
    old_matches = load(HERE / "references/peptide-matches.json")
    latest_matches = load(latest_dir / "peptide-matches.json")
    from curate import peptide_windows
    for region in selected:
        for _, peptide in peptide_windows(
            region["sequence"], int(region["target_start_in_window"]),
            int(region["target_end_in_window"]), (8, 9, 10, 11, 15),
        ):
            require(old_matches[peptide]["occurs"] or not latest_matches[peptide]["occurs"],
                    "The latest reference introduces an intended-target self match")
    if rescan_reference:
        require(v46_path.is_file() and v50_path.is_file(),
                "External GENCODE v46/v50 translation FASTAs are required for --rescan-reference; "
                "obtain them from the EBI URLs in references/reference-receipt.json and "
                "references/v50-audit/reference-receipt.json into work/gencode-references/ or YASH_GENCODE_DIR")
        from reference_screen import exact_search, wanted_peptides
        wanted = wanted_peptides()
        for fasta_path, directory, expected in ((v46_path, HERE / "references", receipt), (v50_path, latest_dir, latest)):
            with gzip.open(fasta_path, "rt") as stream:
                matches, protein_count, residue_count = exact_search(
                    ((record.id, str(record.seq)) for record in SeqIO.parse(stream, "fasta")), wanted
                )
            require(matches == load(directory / "peptide-matches.json"),
                    "Fresh exact proteome rescan differs: " + expected["reference_release"])
            require(protein_count == expected["proteins_scanned"]
                    and residue_count == expected["residues_scanned"],
                    "Reference record/residue count differs: " + expected["reference_release"])
        pn_receipt = load(HERE / "references/patient-normal-receipt.json")
        with gzip.open(HERE / "references" / pn_receipt["bundled_fasta_file"], "rt") as stream:
            pn_matches, pn_count, pn_residues = exact_search(
                ((record.id, str(record.seq)) for record in SeqIO.parse(stream, "fasta")), wanted
            )
        require(pn_matches == load(HERE / "references/patient-normal-matches.json"),
                "Fresh exact patient-normal proteome rescan differs")
        require(pn_count == pn_receipt["proteins_scanned"]
                and pn_residues == pn_receipt["residues_scanned"],
                "Patient-normal proteome record/residue count differs")
    ag = load(HERE / "alphagenome2/receipt.json")
    for name, expected in ag["output_sha256"].items():
        require(sha(HERE / "alphagenome2" / name) == expected, "Original AlphaGenome table changed")
    condition = ag["expression_conditioning"]
    require(sha(HERE / "alphagenome2" / condition["condition_embedding_file"]) == condition["condition_embedding_sha256"],
            "Original expression-conditioning embedding changed")
    ag2_offline = verify_alphagenome2_offline()
    ver_receipt = load(HERE / "review/verification.json")
    for rel, exp_hash in ver_receipt.get("executed_scientific_inputs_sha256", {}).items():
        require(sha(HERE / rel) == exp_hash, "Scientific input hash disagrees with verification receipt: " + rel)
    sv = load(HERE / "review/source-verification.json")
    require(sv["source_receipt_sha256"] == sha(HERE / "exports/source-receipt.json"),
            "Source verification receipt is out of sync with exports/source-receipt.json")
    # This receipt uses two hashing conventions, so every record must declare which one it means:
    # file_bytes hashes bundle_file as stored, decompressed_bytes hashes the gzipped payload that is
    # actually present. Reading them with one convention silently accepts stale values.
    for sv_group in ("adapted_from_checked_upstream_records", "previous_local_research_records", "repository_records"):
        for item in sv.get(sv_group, []):
            rel = item["bundle_file"]
            plain = HERE / rel
            target = item.get("hash_target")
            if target == "file_bytes":
                require(plain.is_file(),
                        "Record declares hash_target file_bytes but the file is absent: " + rel)
                observed = hashlib.sha256(plain.read_bytes()).hexdigest()
            elif target == "decompressed_bytes":
                packed = Path(str(plain) + ".gz")
                require(packed.is_file(),
                        "Record declares hash_target decompressed_bytes but no .gz payload is present: " + rel)
                observed = hashlib.sha256(gzip.decompress(packed.read_bytes())).hexdigest()
            else:
                require(False, "Source verification record needs an explicit hash_target: " + rel)
            require(observed == item.get("bundled_sha256", item.get("sha256")),
                    "Source verification recorded hash disagrees with the named file: " + rel)
    md_rows = rows(HERE / "review/missing-data-priorities.tsv")
    require([r["priority_id"] for r in md_rows] == [f"MD-0{i}" for i in range(1, 7)],
            "Expected six prioritized missing-data items MD-01..MD-06")

    directory = HERE / "construct/protein"
    construct = load(directory / "construct.json")
    edges = read_edges(directory / "edge-costs.tsv")
    triplets = read_triplets(directory / "processing-triplets.tsv")
    leads = load(directory / "lead-targets.json")
    boundary_query_scores = {}
    for row in rows(directory / "boundary-processing-scores.tsv"):
        key = (row["peptide"], row["n_flank"], row["c_flank"])
        require(key not in boundary_query_scores,
                "Duplicate boundary processing query in recorded table: " + row["peptide"])
        boundary_query_scores[key] = int(row["processing_score_x10000"])
    canonical = load(HERE / "construct/consensus-receipt.json")
    path = canonical_path(len(selected), edges, canonical["chosen_edge_ids"])
    sequence, parts, joins = assemble(path, selected, construct["SEC"], construct["MITD"])
    recorded_sequence = str(SeqIO.read(directory / "cassette.fasta", "fasta").seq)
    require(sequence == recorded_sequence, "Recorded protein differs from the selected path")
    for expected, observed in zip(parts, rows(directory / "components.tsv")):
        require(all(cell(expected[key]) == observed[key] for key in expected), "Protein component changed")
    expected_order = [
        {"position": str(i + 1), **{key: selected[edge["right_index"]][key]
                                    for key in ("antigen_id", "region_id", "gene", "sequence")}}
        for i, edge in enumerate(path[:-1])
    ]
    require(rows(directory / "antigen-order.tsv") == expected_order, "Antigen order table differs from selected path")
    values = {
        (row["peptide"], row["allele"]): (float(row["affinity_nm"]), float(row["affinity_percentile"]))
        for row in iter_rows(HERE / "construct/all-junction-affinity.tsv.gz")
    }
    alleles = load(HERE / "exports/source-case.json")["hla_class_I"]
    require(len(leads) == len(selected), "Lead target count differs from selected regions")
    for idx, (region, lead) in enumerate(zip(selected, leads)):
        cands = [
            (values[pep, al][0], pep, al, st)
            for st, pep in target_candidates(region, old_matches) for al in alleles
        ]
        best_nm, best_pep, best_al, best_st = min(cands)
        best_en = best_st + len(best_pep)
        require(
            lead["region_index"] == idx and lead["region_id"] == region["region_id"]
            and lead["peptide"] == best_pep and lead["allele"] == best_al
            and math.isclose(float(lead["affinity_nm"]), best_nm, rel_tol=1e-12, abs_tol=1e-12)
            and lead["start_in_region"] == best_st and lead["end_in_region"] == best_en
            and lead["needs_N_neighbor"] == (best_st < 5)
            and lead["needs_C_neighbor"] == (len(region["sequence"]) - best_en < 5),
            "Lead target record mismatch: " + region["gene"],
        )
    mix = {row["Peptide"]: row for row in rows(HERE / "immunogenicity/mixmhcpred3.tsv")}
    prime = {row["Peptide"]: row for row in rows(HERE / "immunogenicity/prime21.tsv")}
    c2_edge_map = {
        (row["Peptide"], row["Context"]): row
        for row in iter_rows(HERE / "construct/all-junction-class2.tsv.gz")
    }
    linkers = {row["name"]: row["sequence"] for row in load(HERE / "construct/linkers.json")["choices"]}

    # Boundary processing scores: every edge's intended_proc_score_x10000 is the sum of
    # round(processing_score * 10000) over all boundary-dependent lead-target (peptide, n_flank, c_flank)
    # queries on that edge. We reconstruct the exact query tuple(s) each edge carries from the selected
    # regions, SEC/MITD and lead-targets.json, and re-derive every edge's sum offline from the 128 recorded
    # per-query predictions in construct/protein/boundary-processing-scores.tsv.gz (or verify 0 when an
    # edge carries no boundary-dependent lead query).
    def edge_queries(edge):
        i, j = edge["left_index"], edge["right_index"]
        n = len(selected)
        left = construct["SEC"] if i == n else selected[i]["sequence"]
        right = construct["MITD"] if j == n else selected[j]["sequence"]
        linker = linkers[edge["linker"]]
        out = []
        if i < n:
            lead = leads[i]
            if lead["needs_C_neighbor"] and not lead["needs_N_neighbor"]:
                a, b = lead["start_in_region"], lead["end_in_region"]
                out.append((lead["peptide"], left[a - 5:a], (left[b:] + linker + right)[:5]))
        if j < n:
            lead = leads[j]
            if lead["needs_N_neighbor"] and not lead["needs_C_neighbor"]:
                a, b = lead["start_in_region"], lead["end_in_region"]
                out.append((lead["peptide"], (left + linker + right[:a])[-5:], right[b:b + 5]))
        return out

    for edge in edges:
        require(edge["linker"] in linkers, "Edge names an unknown linker: " + edge["edge_id"])
    # Every boundary processing query the edges raise must appear exactly once in the recorded
    # per-query table, and the table must carry no row that no edge uses - so a tampered or stale
    # row cannot hide in an unsummed entry.
    used_queries = set()
    for edge in edges:
        used_queries.update(edge_queries(edge))
    require(set(boundary_query_scores) == used_queries,
            "Recorded per-query processing table does not match the edges' boundary queries exactly")
    proc_counts = Counter()

    for edge in edges:
        require(linkers[edge["linker"]] == edge["linker_sequence"], "Edge linker changed")
        left = construct["SEC"] if edge["left_index"] == len(selected) else selected[edge["left_index"]]["sequence"]
        right = construct["MITD"] if edge["right_index"] == len(selected) else selected[edge["right_index"]]["sequence"]
        require(edge["left"] == left and edge["right"] == right, "Edge source sequence changed")
        peptides = junction_peptides(left, edge["linker_sequence"], right)
        scores = [values[peptide, allele] for peptide in peptides for allele in alleles]
        c2_sc, c2_mc, c2_all = junction_class2_counts(
            edge["left_index"], edge["right_index"], left, edge["linker_sequence"],
            right, len(selected), c2_edge_map,
        )
        expected = {
            "strong_rank": sum(rank <= 0.5 for _, rank in scores),
            "mild_rank": sum(rank <= 2.0 for _, rank in scores),
            "affinity_50": sum(nm <= 50 for nm, _ in scores),
            "affinity_500": sum(nm <= 500 for nm, _ in scores),
            "added_aa": len(edge["linker_sequence"]),
            "peptide_occurrences": len(peptides),
            "mix_strong_rank": sum(
                float(mix[peptide]["%Rank_" + allele[4:].replace("*", "").replace(":", "")]) <= 0.5
                for peptide in peptides for allele in alleles
            ),
            "prime_low_rank": sum(
                float(prime[peptide]["%Rank_" + allele[4:].replace("*", "").replace(":", "")]) <= 0.5
                for peptide in peptides for allele in alleles
            ),
            "mixmhc2_strong_core": c2_sc,
            "mixmhc2_mild_core": c2_mc,
            "mixmhc2_all_mild": c2_all,
        }
        require(all(edge[key] == value for key, value in expected.items()), "Cached edge cost changed")
        require(edge["min_rank"] == min(rank for _, rank in scores), "Edge rank minimum changed")
        require(math.isclose(edge["min_affinity_nm"], min(nm for nm, _ in scores), rel_tol=1e-12, abs_tol=1e-12),
                "Edge affinity minimum changed")
        hydro = summarize(left + edge["linker_sequence"] + right)
        require(math.isclose(edge["max_pair_mean_kd"], hydro["max_mean_kd"], rel_tol=1e-12, abs_tol=1e-12),
                "Edge max pair mean KD changed")
        require(edge["hydrophobicity_pass"] == (edge["category"] != "internal" or hydro["status"] == "PASS"),
                "Edge hydropathy eligibility changed")
        queries = edge_queries(edge)
        recorded_proc = edge["intended_proc_score_x10000"]
        if not queries:
            require(recorded_proc == 0,
                    "Edge with no qualifying boundary query must score 0: " + edge["edge_id"])
            proc_counts["zero_by_construction"] += 1
        else:
            # Re-derive the edge's processing score by summing the recorded per-query table over the
            # edge's reconstructed queries. Each per-query score is one rounded probability in
            # [0, 10000]. This covers every edge exactly, chosen or not, with no self-referential gap.
            for query in queries:
                require(0 <= boundary_query_scores[query] <= 10000,
                        "Recorded per-query processing score out of valid range: " + edge["edge_id"])
            require(recorded_proc == sum(boundary_query_scores[query] for query in queries),
                    "Edge boundary processing score disagrees with recorded per-query table: "
                    + edge["edge_id"])
            proc_counts["re_derived_from_recorded_queries"] += 1
    edge_by_id = {edge["edge_id"]: edge for edge in edges}
    triplet_query_scores = {}
    for trip in triplets:
        t_lead = leads[trip["region_index"]]
        region_seq = selected[trip["region_index"]]["sequence"]
        a, b = t_lead["start_in_region"], t_lead["end_in_region"]
        before, after = edge_by_id[trip["incoming_edge_id"]], edge_by_id[trip["outgoing_edge_id"]]
        require(
            t_lead["needs_N_neighbor"] and t_lead["needs_C_neighbor"]
            and before["right_index"] == trip["region_index"]
            and after["left_index"] == trip["region_index"]
            and trip["peptide"] == t_lead["peptide"]
            and trip["n_flank"] == (before["left"] + before["linker_sequence"] + region_seq[:a])[-5:]
            and trip["c_flank"] == (region_seq[b:] + after["linker_sequence"] + after["right"])[:5]
            and 0 <= trip["processing_score_x10000"] <= 10000,
            "Triplet processing query does not re-derive from the recorded construct",
        )
        t_key = (trip["peptide"], trip["n_flank"], trip["c_flank"])
        require(triplet_query_scores.setdefault(t_key, trip["processing_score_x10000"])
                == trip["processing_score_x10000"],
                "Identical triplet processing queries carry different recorded scores")
    predictions, totals = assess(sequence, joins, path, alleles, values)
    objectives = objective_values(edges, path, triplets)
    for phase in construct["solver"]["phases"]:
        require(phase["status"] == "OPTIMAL" and phase["value"] == phase["bound"], "A solver phase lacks its recorded optimal bound")
        require(objectives[phase["objective"]] == phase["value"], "Canonical path objective changed")
    fresh = rows(HERE / "construct/selected-inference/predictions.tsv")
    require(len(fresh) == len(predictions), "Fresh repeat model inference row count changed")
    for expected, observed in zip(predictions, fresh):
        require(expected["peptide"] == observed["peptide"] and expected["allele"] == observed["allele"], "Fresh repeat identity changed")
        require(abs(float(observed["affinity_nm"]) / expected["affinity_nm"] - 1) < 1e-4, "Fresh repeat numerical tolerance failed")
        require(
            [float(observed["affinity_nm"]) <= 50, float(observed["affinity_nm"]) <= 500,
             float(observed["affinity_percentile"]) <= 0.5, float(observed["affinity_percentile"]) <= 2.0]
            == [expected["affinity_nm"] <= 50, expected["affinity_nm"] <= 500,
                expected["affinity_percentile"] <= 0.5, expected["affinity_percentile"] <= 2.0],
            "Main fresh repeat threshold classification differs",
        )
    _, hydro = protein_screen(sequence, parts)
    require(hydro == load(directory / "hydrophobicity.json") and hydro["status"] == "PASS", "Full-chain hydropathy record differs")
    short = [row for row in hydro["antigen_regions"] if row["length_aa"] < 19]
    require(len(short) == 1 and short[0]["status"] == "NOT_SCORABLE_SHORTER_THAN_19_AA", "Short FES was incorrectly padded or individually passed")
    final = load(HERE / "construct/final-check.json")
    main_leads = rows(HERE / "immunogenicity/placed-class1-leads.tsv")
    main_placed_proc = sum(
        round(float(r["processing_score"]) * 10000) for r in main_leads
        if r["needs_N_neighbor"] == "True" or r["needs_C_neighbor"] == "True"
    )
    require(final["protein_sha256"] == sha(directory / "cassette.fasta")
            and final["components_sha256"] == sha(directory / "components.tsv")
            and final["threshold_classifications_identical"] is True
            and final["processing_objective_matches_full_cassette"] is True
            and abs(main_placed_proc - construct["solver"]["phases"][-1]["value"]) <= 1,
            "Fresh recheck refers to another construct or contradicts retained lead/classification evidence")
    for name, expected in load(HERE / "immunogenicity/receipt.json")["output_sha256"].items():
        require(sha(HERE / "immunogenicity" / name) == expected, "Immunogenicity output changed: " + name)

    mrna = HERE / "construct/mrna"
    nucleotide = load(mrna / "verification.json")
    for group in ("input_sha256", "output_sha256"):
        for relative, expected in nucleotide[group].items():
            require(sha(mrna / relative) == expected, "Nucleotide " + group + " changed: " + relative)
    cfg = load(mrna / "backbone.json")
    coding = str(SeqIO.read(mrna / "cds.fasta", "fasta").seq)[:-3]
    require(str(Seq(coding).translate()) == sequence, "RNA encoding changes an amino acid")
    audit = load(mrna / "optimizer-audit.json")
    snapshot = audit["coding_rna_after_stop_insertion"]
    expected_stops = [
        [offset, start, start + 3] for offset in (1, 2)
        for start in range(offset, len(snapshot) - 2, 3)
        if snapshot[start:start + 3] in ("UAA", "UAG", "UGA")
    ]
    require(audit["alternative_frame_stop_intervals"] == expected_stops, "An original alternative-frame stop was omitted")
    protected = set(audit["protected_nt_positions"])
    require(set(range(18)) <= protected and all(
        set(range(start, end)) <= protected for _, start, end in expected_stops
    ), "Stop or initial-codon protection is incomplete")
    coding_rna = coding.replace("T", "U")
    require(all(coding_rna[i] == snapshot[i] for i in protected), "A protected nucleotide changed")
    RNA.cvar.temperature = 37.0
    module = mrna_module()
    metrics = module.metrics(coding, sequence, load(mrna / "human-codon-usage.json")["codon_usage"], cfg["optimization"])
    require(metrics == load(mrna / "construct.json")["optimization_metrics"], "Independent nucleotide metrics differ")
    require(metrics["rare_codons_below_cutoff"] == metrics["repeated_12nt_pairs"] == 0
            and not any(metrics["coding_forbidden_motif_occurrences"].values()), "A hard coding constraint failed")
    require(nucleotide["complete_construct_verified"] is False
            and nucleotide["polyA_status"] == "unspecified", "Unknown molecular architecture was fabricated")
    require(sha(mrna / "optimizer.py") == cfg["optimization"]["optimizer_script_sha256"],
            "Optimizer implementation hash in the configuration is stale")
    for name, expected in load(mrna / "benchmarks/receipt.json")["output_sha256"].items():
        require(sha(mrna / "benchmarks" / name) == expected, "RNA benchmark changed: " + name)
    require(load(mrna / "benchmarks/before-rescue-audit.json")["coding_rna_after_stop_insertion"] == snapshot,
            "Benchmark stop-protection reference no longer matches this encoding")
    raw_components = rows(mrna / "components.tsv")
    components, cursor = [], 0
    before = len(cfg["template_context"]["prefix_before_transcribed_sequence_dna"])
    for raw in raw_components:
        row = dict(raw)
        unknown = row["kind"] == "polyA_unspecified"
        for key in ("rna_start_0based", "rna_end_exclusive", "template_start_0based",
                    "template_end_exclusive", "length_nt"):
            row[key] = (0 if key == "length_nt" else
                        cursor if key.startswith("rna_") else before + cursor) if unknown else int(row[key])
        cursor = row["rna_end_exclusive"]
        components.append(row)
    from genbank import build_record, load_inputs, verify_export
    sources, annotations, architecture, _ = load_inputs(mrna)
    dna = str(SeqIO.read(mrna / "template-segments.fasta", "fasta").seq)
    rna_segment = str(SeqIO.read(mrna / "mrna-segments.fasta", "fasta").seq)
    from accessibility import ensemble_accessibility
    calculated_accessibility = ensemble_accessibility(rna_segment, len(cfg["fixed_elements"]["utr5_dna"]))
    accessibility = audit["initiation_accessibility"]
    for name in ("mean_marginal_unpaired_probability", "AUG_mean_marginal_unpaired_probability"):
        require(abs(calculated_accessibility[name] - accessibility["selected"][name]) < 1e-8,
                "Fresh ensemble accessibility differs")
        require(accessibility["selected"][name] + 1e-9 >= accessibility["baseline"][name],
                "A claimed initiation-accessibility improvement worsened its baseline")
    expected_records = build_record(dna, sequence, components, cfg, sources, annotations, architecture)
    genbank = verify_export(mrna / "template-segments.gb", expected_records, components,
                           sequence, (rna_segment,), annotations)
    require(genbank == nucleotide["genbank_annotation"], "Independent GenBank round-trip check differs")
    strand_check = verify_strand_offline()
    # Verify exact full-chain Class-II context coverage and counts on main-16 (zero flank substitution)
    c2_main_ctx = {(r["Peptide"], r["Context"]): r for r in rows(HERE / "immunogenicity/mixmhc2-context.tsv")}
    n2_main = {r["Peptide"]: r for r in rows(HERE / "immunogenicity/mixmhc2-no-context.tsv")}
    _, main_c2_summary = assess_assembled_class2(sequence, joins, c2_main_ctx, n2_main)
    recorded_junc_mm = load(HERE / "immunogenicity/junction-multi-model-summary.json")
    for key, val in main_c2_summary.items():
        require(recorded_junc_mm[key] == val, f"Main-16 assembled Class-II summary mismatch: {key}")

    # Verify across all hydropathy-passing two-neighbor FES (ein, FES, eout) combinations that exact
    # two-neighbor contexts in all-junction-class2.tsv.gz match single-edge Class-II counts identically
    fes_idx = next(i for i, r in enumerate(selected) if r["gene"] == "FES")
    fes_in = [e for e in edges if e["hydrophobicity_pass"] and e["right_index"] == fes_idx]
    fes_out = [e for e in edges if e["hydrophobicity_pass"] and e["left_index"] == fes_idx]
    for ein in fes_in:
        for eout in fes_out:
            if ein["left_index"] == eout["right_index"]:
                continue
            tri_seq = ein["left"] + ein["linker_sequence"] + ein["right"] + eout["linker_sequence"] + eout["right"]
            b1, e1 = len(ein["left"]), len(ein["left"]) + len(ein["linker_sequence"])
            b2, e2 = e1 + len(ein["right"]), e1 + len(ein["right"]) + len(eout["linker_sequence"])
            sc_exact, mc_exact, am_exact = 0, 0, 0
            for b, e in ((b1, e1), (b2, e2)):
                for st in range(len(tri_seq) - 14):
                    en = st + 15
                    touches = (st < e and en > b) if e > b else (st < b < en)
                    if not touches:
                        continue
                    pep = tri_seq[st:en]
                    ctx = assembled_class2_context(tri_seq, st, en)
                    scored = c2_edge_map[pep, ctx]
                    for allele in CLASS2_ALLELES:
                        rk = float(scored["%Rank_" + allele])
                        if rk <= 2.0:
                            am_exact += 1
                            cp1 = int(scored["CoreP1_" + allele]) - 1
                            cs, ce = st + cp1, st + cp1 + 9
                            if (cs < e and ce > b) if e > b else (cs < b < ce):
                                mc_exact += 1
                                if rk <= 0.5:
                                    sc_exact += 1
            require(
                (sc_exact, mc_exact, am_exact) == (
                    ein["mixmhc2_strong_core"] + eout["mixmhc2_strong_core"],
                    ein["mixmhc2_mild_core"] + eout["mixmhc2_mild_core"],
                    ein["mixmhc2_all_mild"] + eout["mixmhc2_all_mild"],
                ),
                f"Two-neighbor FES Class-II context differs from pairwise sum: {ein['edge_id']} -> {eout['edge_id']}",
            )

    strict_dir = HERE / "strand-strict-15"
    strict_selected = load(strict_dir / "selection/selected.json")
    require(strict_selected == strand_strict_core_regions(derive()[0]) and len(strict_selected) == 15,
            "Strand-strict 15-region contingency selection differs")
    strict_receipt = load(strict_dir / "contingency-receipt.json")
    strict_construct = load(strict_dir / "construct/protein/construct.json")
    kctd_idx = next(i for i, r in enumerate(selected) if r["gene"] == "KCTD19")
    def remap_strict(idx):
        if idx == kctd_idx:
            return None
        if idx == len(selected):
            return len(strict_selected)
        return idx if idx < kctd_idx else idx - 1
    edges15 = [
        dict(edge, left_index=remap_strict(edge["left_index"]), right_index=remap_strict(edge["right_index"]))
        for edge in edges
        if remap_strict(edge["left_index"]) is not None and remap_strict(edge["right_index"]) is not None
    ]
    kept_ids = {edge["edge_id"] for edge in edges15}
    triplets15 = [
        dict(trip, region_index=remap_strict(trip["region_index"]))
        for trip in triplets
        if remap_strict(trip["region_index"]) is not None
        and trip["incoming_edge_id"] in kept_ids and trip["outgoing_edge_id"] in kept_ids
    ]
    strict_path = canonical_path(len(strict_selected), edges15, strict_receipt["chosen_edge_ids"])
    strict_seq, strict_parts, strict_joins = assemble(strict_path, strict_selected, strict_construct["SEC"], strict_construct["MITD"])
    require(strict_seq == str(SeqIO.read(strict_dir / "construct/protein/cassette.fasta", "fasta").seq),
            "Strand-strict 15-region contingency protein sequence mismatch")
    for expected, observed in zip(strict_parts, rows(strict_dir / "construct/protein/components.tsv")):
        require(all(cell(expected[key]) == observed[key] for key in expected), "Fallback-15 protein component changed")
    expected_strict_order = [
        {"position": str(i + 1), **{key: strict_selected[edge["right_index"]][key]
                                    for key in ("antigen_id", "region_id", "gene", "sequence")}}
        for i, edge in enumerate(strict_path[:-1])
    ]
    require(rows(strict_dir / "construct/protein/antigen-order.tsv") == expected_strict_order,
            "Fallback-15 antigen order table differs from selected path")
    strict_leads_targets = load(strict_dir / "construct/protein/lead-targets.json")
    expected_strict_leads = [
        dict(lead, region_index=remap_strict(lead["region_index"]))
        for lead in leads if remap_strict(lead["region_index"]) is not None
    ]
    require(strict_leads_targets == expected_strict_leads, "Fallback-15 lead-targets.json differs from remapped leads")
    require(objective_values(edges15, strict_path, triplets15) == strict_receipt["solver_objectives"],
            "Strand-strict 15-region contingency objective values mismatch")
    strict_preds, strict_totals = assess(strict_seq, strict_joins, strict_path, alleles, values)
    strict_fc = load(strict_dir / "construct/final-check.json")
    require(strict_fc["totals"] == strict_totals
            and strict_fc["protein_sha256"] == sha(strict_dir / "construct/protein/cassette.fasta")
            and strict_fc["components_sha256"] == sha(strict_dir / "construct/protein/components.tsv"),
            "Fallback-15 junction check totals or input hashes mismatch")
    if strict_fc["status"] == "FRESH_FULL_JUNCTION_CHECK_PASSED":
        require(strict_fc.get("fresh_repeat_model_inference_performed") is True
                and strict_fc.get("threshold_classifications_identical") is True
                and strict_fc.get("processing_objective_matches_full_cassette") is True,
                "Fallback-15 fresh junction status flags incomplete")
        strict_fresh = rows(strict_dir / "construct/selected-inference/predictions.tsv")
        require(len(strict_fresh) == len(strict_preds), "Fallback-15 fresh repeat model inference row count changed")
        for expected_p, observed_p in zip(strict_preds, strict_fresh):
            require(expected_p["peptide"] == observed_p["peptide"] and expected_p["allele"] == observed_p["allele"],
                    "Fallback-15 fresh repeat identity changed")
            require(abs(float(observed_p["affinity_nm"]) / expected_p["affinity_nm"] - 1) < 1e-4,
                    "Fallback-15 fresh repeat numerical tolerance failed")
            require(
                [float(observed_p["affinity_nm"]) <= 50, float(observed_p["affinity_nm"]) <= 500,
                 float(observed_p["affinity_percentile"]) <= 0.5, float(observed_p["affinity_percentile"]) <= 2.0]
                == [expected_p["affinity_nm"] <= 50, expected_p["affinity_nm"] <= 500,
                    expected_p["affinity_percentile"] <= 0.5, expected_p["affinity_percentile"] <= 2.0],
                "Fallback-15 fresh repeat threshold classification differs",
            )
        strict_leads = rows(strict_dir / "immunogenicity/placed-class1-leads.tsv")
        strict_placed_proc = sum(
            round(float(r["processing_score"]) * 10000) for r in strict_leads
            if r["needs_N_neighbor"] == "True" or r["needs_C_neighbor"] == "True"
        )
        require(abs(strict_placed_proc - strict_construct["solver"]["phases"][-1]["value"]) <= 1,
                "Fallback-15 placed full-cassette processing score disagrees with solver objective")
    elif strict_fc["status"] == "CACHED_ASSEMBLY_JUNCTION_CHECK_PASSED":
        require(strict_fc.get("fresh_repeat_model_inference_performed") is False
                and strict_fc.get("cached_boundary_enumeration_verified") is True
                and strict_fc.get("max_relative_affinity_difference") is None,
                "Fallback-15 cached junction check metadata inconsistent")
    else:
        require(False, "Unsupported fallback-15 junction check status: " + str(strict_fc["status"]))
    # Verify complete assembled Class-II context coverage and recounted junction scores on fallback-15
    for name, expected_hash in strict_receipt["immunogenicity_sha256"].items():
        require(sha(strict_dir / "immunogenicity" / name) == expected_hash,
                f"Strand-strict 15-region immunogenicity checksum mismatch: {name}")
    c2_strict_ctx = {(r["Peptide"], r["Context"]): r for r in rows(strict_dir / "immunogenicity/mixmhc2-context.tsv")}
    n2_strict = {r["Peptide"]: r for r in rows(strict_dir / "immunogenicity/mixmhc2-no-context.tsv")}
    for st in range(len(strict_seq) - 14):
        pep = strict_seq[st:st + 15]
        ctx = assembled_class2_context(strict_seq, st, st + 15)
        require((pep, ctx) in c2_strict_ctx and (pep, ctx) in c2_edge_map and pep in n2_strict,
                f"Missing fallback-15 assembled Class-II prediction: ({pep!r}, {ctx!r})")
    strict_c2_rows, strict_c2_summary = assess_assembled_class2(strict_seq, strict_joins, c2_strict_ctx, n2_strict)
    recorded_strict_mm = load(strict_dir / "immunogenicity/junction-multi-model-summary.json")
    for key, val in strict_c2_summary.items():
        require(recorded_strict_mm[key] == val, f"Fallback-15 assembled Class-II summary mismatch: {key}")
    require(len(strict_c2_rows) == len(rows(strict_dir / "immunogenicity/class2-junction-scores.tsv")),
            "Fallback-15 Class-II junction score row count mismatch")
    mrna15 = strict_dir / "construct/mrna"
    require(strict_receipt["mrna_verification_sha256"] == sha(mrna15 / "verification.json"),
            "Fallback-15 mRNA verification receipt checksum mismatch")
    nuc15 = load(mrna15 / "verification.json")
    for group in ("input_sha256", "output_sha256"):
        for relative, expected_hash in nuc15[group].items():
            require(sha(mrna15 / relative) == expected_hash, "Fallback-15 nucleotide " + group + " changed: " + relative)
    cfg15 = load(mrna15 / "backbone.json")
    coding15 = str(SeqIO.read(mrna15 / "cds.fasta", "fasta").seq)[:-3]
    require(str(Seq(coding15).translate()) == strict_seq, "Fallback-15 RNA encoding changes an amino acid")
    metrics15 = module.metrics(coding15, strict_seq, load(mrna15 / "human-codon-usage.json")["codon_usage"], cfg15["optimization"])
    require(metrics15 == load(mrna15 / "construct.json")["optimization_metrics"], "Fallback-15 nucleotide metrics differ")
    require(metrics15["rare_codons_below_cutoff"] == metrics15["repeated_12nt_pairs"] == 0
            and not any(metrics15["coding_forbidden_motif_occurrences"].values()), "Fallback-15 hard coding constraint failed")
    require(nuc15["complete_construct_verified"] is False and nuc15["polyA_status"] == "unspecified",
            "Fallback-15 fabricated unknown PolyA architecture")
    audit15 = load(mrna15 / "optimizer-audit.json")
    snap15 = audit15["coding_rna_after_stop_insertion"]
    exp_stops15 = [
        [offset, start, start + 3] for offset in (1, 2)
        for start in range(offset, len(snap15) - 2, 3)
        if snap15[start:start + 3] in ("UAA", "UAG", "UGA")
    ]
    require(audit15["alternative_frame_stop_intervals"] == exp_stops15, "Fallback-15 original alternative-frame stop omitted")
    prot15_set = set(audit15["protected_nt_positions"])
    require(set(range(18)) <= prot15_set and all(
        set(range(start, end)) <= prot15_set for _, start, end in exp_stops15
    ), "Fallback-15 stop or initial-codon protection incomplete")
    coding15_rna = coding15.replace("T", "U")
    require(all(coding15_rna[i] == snap15[i] for i in prot15_set),
            "Fallback-15 protected nucleotide changed")
    raw_comp15 = rows(mrna15 / "components.tsv")
    comp15, cur15 = [], 0
    before15 = len(cfg15["template_context"]["prefix_before_transcribed_sequence_dna"])
    for raw in raw_comp15:
        row = dict(raw)
        unknown = row["kind"] == "polyA_unspecified"
        for key in ("rna_start_0based", "rna_end_exclusive", "template_start_0based",
                    "template_end_exclusive", "length_nt"):
            row[key] = (0 if key == "length_nt" else
                        cur15 if key.startswith("rna_") else before15 + cur15) if unknown else int(row[key])
        cur15 = row["rna_end_exclusive"]
        comp15.append(row)
    src15, ann15, arch15, _ = load_inputs(mrna15)
    dna15 = str(SeqIO.read(mrna15 / "template-segments.fasta", "fasta").seq)
    rna15_seg = str(SeqIO.read(mrna15 / "mrna-segments.fasta", "fasta").seq)
    calc_acc15 = ensemble_accessibility(rna15_seg, len(cfg15["fixed_elements"]["utr5_dna"]))
    acc15 = audit15["initiation_accessibility"]
    for name in ("mean_marginal_unpaired_probability", "AUG_mean_marginal_unpaired_probability"):
        require(abs(calc_acc15[name] - acc15["selected"][name]) < 1e-8,
                "Fallback-15 fresh ensemble accessibility differs")
        require(acc15["selected"][name] + 1e-9 >= acc15["baseline"][name],
                "Fallback-15 initiation-accessibility worsened baseline")
    exp_rec15 = build_record(dna15, strict_seq, comp15, cfg15, src15, ann15, arch15)
    gb15 = verify_export(mrna15 / "template-segments.gb", exp_rec15, comp15, strict_seq, (rna15_seg,), ann15)
    require(gb15 == nuc15["genbank_annotation"], "Fallback-15 GenBank round-trip check differs")

    safety_check = verify_safety_offline()
    review_summary_check = verify_review_summary()
    return {
        "status": "VERIFIED_RESEARCH_ARTIFACTS_NOT_CLINICAL_VALIDATION",
        "selected_regions": len(selected), "source_counts": dict(Counter(r["source"] for r in selected)),
        "protein_length_aa": len(sequence), "coding_sense_nt": len(coding),
        "junctions": totals, "all_recorded_solver_bounds_match": True,
        "final_paired_bracket_run_nt": metrics["coding_longest_paired_bracket_run"],
        "paired_bracket_soft_target_met": nucleotide["folding_heuristic_target_retained"],
        "repeated_12nt_pairs": metrics["repeated_12nt_pairs"],
        "all_original_alternative_frame_stops_preserved": len(expected_stops),
        "genbank_features": genbank["feature_count"],
        "rna_strand_audit": strand_check,
        "boundary_processing_score_derivation": {
            "edges": len(edges),
            "zero_by_construction": proc_counts["zero_by_construction"],
            "re_derived_from_recorded_queries": proc_counts["re_derived_from_recorded_queries"],
            "distinct_recorded_queries": len(boundary_query_scores),
            "recorded_query_table": "construct/protein/boundary-processing-scores.tsv.gz",
            "note": "Every edge's intended_proc_score_x10000 is re-derived offline by summing the "
                    "recorded per-query processing scores over the edge's reconstructed queries; a "
                    "zero-query edge must score 0. The per-query table was produced by re-inference "
                    "with the byte-identical recorded MHCflurry 2.3.0 bundle and reproduces every "
                    "recorded edge sum exactly (see construct/protein/boundary-processing-receipt.json).",
        },
        "strand_strict_15_contingency_verified": True,
        "strand_strict_15_mrna_sense_coding_nt": len(coding15),
        "alphagenome2_audit": ag2_offline,
        "in_silico_safety": safety_check,
        "review_summary": review_summary_check,
        "complete_molecule": False, "clinical_release_ready": False,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rescan-reference", action="store_true")
    print(json.dumps(verify(parser.parse_args().rescan_reference), indent=2))
