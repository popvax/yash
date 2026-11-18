"""Local Class-I model inference and finite, lexicographic cassette optimization."""
import argparse
import math
import os
from collections import Counter
from importlib.metadata import version
from pathlib import Path

from ortools.sat.python import cp_model

from curate import HERE, peptide_windows
from dataio import fasta, load, save, sha, table
from hydrophobicity import protein_screen, summarize
from construct_identifiers import export_identifiers, variant_for_regions

LENGTHS = (8, 9, 10, 11)
CLASS2_ALLELES = (
    "DRB1_03_01", "DRB1_08_02", "DRB3_02_02",
    "DQA1_05_01__DQB1_02_01", "DQA1_04_01__DQB1_04_02",
    "DQA1_04_01__DQB1_02_01", "DPA1_01_03__DPB1_02_01",
    "DPA1_02_01__DPB1_13_01", "DPA1_01_03__DPB1_13_01",
)


def junction_class2_counts(i, j, left, linker_seq, right, n, c2_map):
    """Single-edge Class-II junction counts; 15-aa FES two-neighbor flanks are verified full-chain via assess_assembled_class2."""
    comb = left + linker_seq + right
    begin, end_lk = len(left), len(left) + len(linker_seq)
    strong_core, mild_core, all_mild = 0, 0, 0
    for st in range(len(comb) - 14):
        en = st + 15
        touches = (st < end_lk and en > begin) if linker_seq else (st < begin < en)
        if not touches:
            continue
        n_ctx = comb[st - 3:st] if st >= 3 else ("-" * (3 - st) + comb[:st] if i == n else "G" * (3 - st) + comb[:st])
        c_ctx = comb[en:en + 3] if en + 3 <= len(comb) else (comb[en:] + "-" * (en + 3 - len(comb)) if j == n else comb[en:] + "G" * (en + 3 - len(comb)))
        scored = c2_map[comb[st:en], n_ctx + comb[st:st + 3] + comb[en - 3:en] + c_ctx]
        for allele in CLASS2_ALLELES:
            rank = float(scored["%Rank_" + allele])
            if rank <= 2.0:
                all_mild += 1
                cp1 = int(scored["CoreP1_" + allele]) - 1
                cs, ce = st + cp1, st + cp1 + 9
                if (cs < end_lk and ce > begin) if linker_seq else (cs < begin < ce):
                    mild_core += 1
                    if rank <= 0.5:
                        strong_core += 1
    return strong_core, mild_core, all_mild


def assembled_class2_context(sequence, start, end):
    peptide = sequence[start:end]
    return (
        sequence[max(0, start - 3):start].rjust(3, "-")
        + peptide[:3] + peptide[-3:]
        + sequence[end:end + 3].ljust(3, "-")
    )


def assess_assembled_class2(sequence, joins, c2_map, n2_map=None):
    """Exact full-chain Class-II junction evaluation with zero flank substitution."""
    rows_out = []
    for st in range(len(sequence) - 14):
        en = st + 15
        pep = sequence[st:en]
        ctx = assembled_class2_context(sequence, st, en)
        if (pep, ctx) not in c2_map:
            raise ValueError(f"Missing assembled Class-II context prediction: ({pep!r}, {ctx!r})")
        touched = [
            join for join in joins
            if (st < int(join["linker_end_aa"]) and en > int(join["linker_start_aa"])
                if int(join["linker_end_aa"]) > int(join["linker_start_aa"])
                else st < int(join["linker_start_aa"]) < en)
        ]
        if not touched:
            continue
        scored = c2_map[pep, ctx]
        for allele in CLASS2_ALLELES:
            cp1 = int(scored["CoreP1_" + allele]) - 1
            cs, ce = st + cp1, st + cp1 + 9
            crosses_core = any(
                (cs < int(join["linker_end_aa"]) and ce > int(join["linker_start_aa"])
                 if int(join["linker_end_aa"]) > int(join["linker_start_aa"])
                 else cs < int(join["linker_start_aa"]) < ce)
                for join in touched
            )
            row = {
                "peptide": pep,
                "context": ctx,
                "cassette_start_aa": st,
                "junction_indices": ";".join(str(join["junction_index"]) for join in touched),
                "allele": allele,
                "context_rank": float(scored["%Rank_" + allele]),
                "no_context_rank": float(n2_map[pep]["%Rank_" + allele]) if n2_map is not None else None,
                "core_9aa": pep[cp1:cp1 + 9],
                "core_crosses_synthetic_junction": crosses_core,
            }
            rows_out.append(row)
    summary = {
        "ClassII_synthetic_9aa_core_rank_le_0_5pct": sum(
            r["core_crosses_synthetic_junction"] and r["context_rank"] <= 0.5 for r in rows_out
        ),
        "ClassII_synthetic_9aa_core_rank_le_2pct": sum(
            r["core_crosses_synthetic_junction"] and r["context_rank"] <= 2.0 for r in rows_out
        ),
        "ClassII_native_core_overhang_rank_le_2pct": sum(
            not r["core_crosses_synthetic_junction"] and r["context_rank"] <= 2.0 for r in rows_out
        ),
        "ClassII_15mer_junction_context_rank_le_2pct": sum(
            r["context_rank"] <= 2.0 for r in rows_out
        ),
    }
    return rows_out, summary


def junction_peptides(left, linker, right):
    sequence = left[-10:] + linker + right[:10]
    begin, end = len(left[-10:]), len(left[-10:]) + len(linker)
    return [
        sequence[start:start + length]
        for length in LENGTHS for start in range(len(sequence) - length + 1)
        if (start < end and start + length > begin
            if linker else start < begin < start + length)
    ]


def full_scan(sequence, joins):
    occurrences = []
    for length in LENGTHS:
        for start in range(len(sequence) - length + 1):
            end = start + length
            touched = [
                join for join in joins
                if (start < join["linker_end_aa"] and end > join["linker_start_aa"]
                    if join["linker_end_aa"] > join["linker_start_aa"]
                    else start < join["linker_start_aa"] < end)
            ]
            if len(touched) > 1:
                raise ValueError("A peptide spans two joins; pairwise junction costs are invalid")
            if touched:
                occurrences.append({
                    "junction_index": touched[0]["junction_index"],
                    "category": touched[0]["category"],
                    "start_aa": start, "end_aa": end, "peptide": sequence[start:end],
                })
    return occurrences


def assemble(path, regions, sec, mitd):
    sequence, components, joins = "", [], []

    def add(kind, name, aa, source=None):
        nonlocal sequence
        source = source or {}
        components.append({
            "kind": kind, "name": name, "aa": aa,
            "start_aa": len(sequence), "end_aa": len(sequence) + len(aa),
            "region_id": source.get("region_id", ""),
            "event_id": source.get("event_id", ""), "gene": source.get("gene", ""),
        })
        sequence += aa

    add("signal_peptide", "HLA_C_01_02_SEC", sec)
    for index, edge in enumerate(path):
        start = len(sequence)
        add("linker", edge["linker"], edge["linker_sequence"])
        joins.append({
            "junction_index": index, "category": edge["category"],
            "left_id": edge["left_id"], "right_id": edge["right_id"],
            "linker": edge["linker"], "linker_sequence": edge["linker_sequence"],
            "linker_start_aa": start, "linker_end_aa": len(sequence),
            "left": edge["left"], "right": edge["right"],
        })
        if edge["right_index"] == len(regions):
            add("mitd", "HLA_B_15_08_MITD", mitd)
        else:
            region = regions[edge["right_index"]]
            add("antigen", region["antigen_id"], region["sequence"], region)
    if Counter((p["region_id"], p["aa"]) for p in components if p["kind"] == "antigen") != Counter(
        (r["region_id"], r["sequence"]) for r in regions
    ):
        raise ValueError("Assembly does not contain every exact antigen once")
    return sequence, components, joins


def summarize_predictions(predictions, allele_count=6):
    return {
        "peptide_occurrences": len(predictions) // allele_count,
        "peptide_hla_rows": len(predictions),
        "strong_rank": sum(row["affinity_percentile"] <= 0.5 for row in predictions),
        "mild_rank": sum(row["affinity_percentile"] <= 2.0 for row in predictions),
        "affinity_50": sum(row["affinity_nm"] <= 50.0 for row in predictions),
        "affinity_500": sum(row["affinity_nm"] <= 500.0 for row in predictions),
        "min_rank": min(row["affinity_percentile"] for row in predictions),
        "min_affinity_nm": min(row["affinity_nm"] for row in predictions),
    }


def assess(sequence, joins, path, alleles, values):
    occurrences = full_scan(sequence, joins)
    predictions = [{
        **row, "allele": allele,
        "affinity_nm": values[row["peptide"], allele][0],
        "affinity_percentile": values[row["peptide"], allele][1],
    } for row in occurrences for allele in alleles]
    for index, edge in enumerate(path):
        actual = Counter(row["peptide"] for row in occurrences if row["junction_index"] == index)
        if actual != Counter(junction_peptides(edge["left"], edge["linker_sequence"], edge["right"])):
            raise ValueError("Full-chain enumeration differs from the edge enumeration")
    totals = summarize_predictions(predictions, len(alleles))
    for key in ("strong_rank", "mild_rank", "affinity_50", "affinity_500", "peptide_occurrences"):
        if totals[key] != sum(edge[key] for edge in path):
            raise ValueError("Placed junction cost differs: " + key)
    return predictions, totals


def eligible_options(edges):
    passing = [edge for edge in edges if edge["hydrophobicity_pass"]]
    if passing and all("mixmhc2_strong_core" in edge for edge in passing):
        return passing
    groups = {}
    for edge in passing:
        groups.setdefault((edge["left_index"], edge["right_index"]), []).append(edge)
    allowed = []
    for _, options in sorted(groups.items()):
        def cost(edge):
            middle = tuple(edge[key] for key in ("mix_strong_rank", "prime_low_rank")
                           if key in edge)
            return (edge["strong_rank"], *middle, edge["mild_rank"])
        best = min(cost(edge) for edge in options)
        allowed.extend(edge for edge in options
                       if cost(edge) == best)
    return allowed


def solve(n, edges, triplets, seconds=600):
    allowed = eligible_options(edges)
    model = cp_model.CpModel()
    items = [(model.new_bool_var("edge_" + str(i)), edge) for i, edge in enumerate(allowed)]
    literals = {edge["edge_id"]: literal for literal, edge in items}
    model.add_circuit([
        (edge["left_index"], edge["right_index"], literal) for literal, edge in items
    ])
    ranks = sorted({edge["min_rank"] for edge in allowed})
    ordinal = {rank: index for index, rank in enumerate(ranks)}
    quality = model.new_int_var(0, len(ranks) - 1, "worst_rank_ordinal")
    for literal, edge in items:
        model.add(quality <= ordinal[edge["min_rank"]] + len(ranks) * (1 - literal))
    pair_items = []
    for index, triplet in enumerate(triplets):
        left = literals[triplet["incoming_edge_id"]]
        right = literals[triplet["outgoing_edge_id"]]
        together = model.new_bool_var("two_flanks_" + str(index))
        model.add(together <= left)
        model.add(together <= right)
        model.add(together >= left + right - 1)
        pair_items.append((together, triplet))
    processing = sum(lit * edge["intended_proc_score_x10000"] for lit, edge in items)
    processing += sum(lit * triplet["processing_score_x10000"] for lit, triplet in pair_items)
    objectives = [
        ("strong_count", sum(lit * edge["strong_rank"] for lit, edge in items), False),
    ]
    for field, name in (
        ("mix_strong_rank", "mixmhcpred_strong_count"),
        ("mixmhc2_strong_core", "mixmhc2_strong_core_count"),
        ("prime_low_rank", "prime_low_rank_count"),
    ):
        if all(field in edge for edge in allowed):
            objectives.append((name, sum(lit * edge[field] for lit, edge in items), False))
    if all("mixmhc2_mild_core" in edge and "affinity_500" in edge for edge in allowed):
        objectives.extend([
            ("total_mild_and_500nm_burden",
             sum(lit * (edge["mild_rank"] + edge["affinity_500"] + edge["mixmhc2_mild_core"])
                 for lit, edge in items), False),
            ("affinity_500_count", sum(lit * edge["affinity_500"] for lit, edge in items), False),
        ])
    objectives += [
        ("mild_count", sum(lit * edge["mild_rank"] for lit, edge in items), False),
        ("worst_rank_ordinal", quality, True),
        ("added_residues", sum(lit * edge["added_aa"] for lit, edge in items), False),
        ("intended_boundary_processing_x10000", processing, True),
    ]
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 1
    solver.parameters.random_seed = 31
    solver.parameters.max_time_in_seconds = seconds
    phases = []
    for name, expression, maximize in objectives:
        model.maximize(expression) if maximize else model.minimize(expression)
        status = solver.solve(model)
        if status != cp_model.OPTIMAL:
            raise RuntimeError("Solver phase not proven optimal: " + name + " " + solver.status_name(status))
        value = solver.value(expression)
        phase = {
            "objective": name, "status": solver.status_name(status),
            "value": value, "bound": solver.best_objective_bound,
        }
        phases.append(phase)
        print("CP-SAT:", phase, flush=True)
        model.add(expression == value)
        model.clear_hints()
        for literal, _ in items + pair_items:
            model.add_hint(literal, solver.value(literal))
    chosen = [edge for literal, edge in items if solver.value(literal)]
    outgoing = {edge["left_index"]: edge for edge in chosen}
    node, path = n, []
    for _ in range(n + 1):
        edge = outgoing[node]
        path.append(edge)
        node = edge["right_index"]
    if node != n or sorted(edge["right_index"] for edge in path[:-1]) != list(range(n)):
        raise ValueError("Solver circuit does not visit every antigen once")
    result = {
        "solver": "OR-Tools CP-SAT", "workers": 1, "seed": 31,
        "all_phases_proven_optimal": True, "phases": phases,
        "objective_tuple": [
            sum(edge["strong_rank"] for edge in path),
            sum(edge["mild_rank"] for edge in path),
            -min(edge["min_rank"] for edge in path),
            sum(edge["added_aa"] for edge in path),
            -phases[-1]["value"],
        ],
        "finite_scope": "Fixed antigen panel, recorded linker choices, model scores and integer processing-score rounding; not a biological optimum.",
        "two_neighbor_processing_terms": len(triplets),
        "consensus_objective_values": {
            phase["objective"]: phase["value"] for phase in phases
            if phase["objective"] in (
                "mixmhcpred_strong_count", "mixmhc2_strong_core_count",
                "prime_low_rank_count", "total_mild_and_500nm_burden", "affinity_500_count",
            )
        },
    }
    return path, result


def infer_affinity(predictor, peptides, alleles):
    values = {}
    for allele in alleles:
        affinity = predictor.predict(peptides=peptides, allele=allele)
        ranks = predictor.percentile_ranks(affinity, allele=allele)
        for peptide, nm, rank in zip(peptides, affinity, ranks):
            nm, rank = float(nm), float(rank)
            if not (math.isfinite(nm) and nm > 0 and math.isfinite(rank) and 0 <= rank <= 100):
                raise ValueError("Invalid predictor output")
            values[peptide, allele] = (nm, rank)
    return values


def target_candidates(region, reference):
    sequence = region["sequence"]
    lo, hi = int(region["target_start_in_window"]), int(region["target_end_in_window"])
    windows = peptide_windows(sequence, lo, hi)
    if region["source"] == "CTA/SELF":
        return windows
    return [(start, peptide) for start, peptide in windows
            if peptide in reference and not reference[peptide]["occurs"]]


def create_edges(regions, linkers, sec, mitd, alleles, values, processor, reference,
                 mix=None, prime=None, c2_map=None):
    n = len(regions)
    leads = []
    for index, region in enumerate(regions):
        candidates = [
            (values[peptide, allele][0], peptide, allele, start)
            for start, peptide in target_candidates(region, reference) for allele in alleles
        ]
        if not candidates:
            raise ValueError("No event-specific target peptide: " + region["gene"])
        affinity, peptide, allele, start = min(candidates)
        end = start + len(peptide)
        leads.append({
            "region_index": index, "region_id": region["region_id"],
            "antigen_id": region["antigen_id"], "gene": region["gene"], "source": region["source"],
            "peptide": peptide, "allele": allele, "affinity_nm": affinity,
            "start_in_region": start, "end_in_region": end,
            "needs_N_neighbor": start < 5,
            "needs_C_neighbor": len(region["sequence"]) - end < 5,
        })
    edges, queries, query_edge_ids = [], [], []
    for i in range(n + 1):
        for j in range(n + 1):
            if i == j:
                continue
            left = sec if i == n else regions[i]["sequence"]
            right = mitd if j == n else regions[j]["sequence"]
            for linker in linkers:
                edge_id = f"{i}:{j}:{linker['name']}"
                peptides = junction_peptides(left, linker["sequence"], right)
                scores = [values[peptide, allele] for peptide in peptides for allele in alleles]
                category = "SEC" if i == n else "MITD" if j == n else "internal"
                hydro = summarize(left + linker["sequence"] + right)
                edges.append({
                    "edge_id": edge_id, "left_index": i, "right_index": j,
                    "left_id": "SEC" if i == n else regions[i]["region_id"],
                    "right_id": "MITD" if j == n else regions[j]["region_id"],
                    "category": category, "left": left, "right": right,
                    "linker": linker["name"], "linker_sequence": linker["sequence"],
                    "strong_rank": sum(rank <= 0.5 for _, rank in scores),
                    "mild_rank": sum(rank <= 2.0 for _, rank in scores),
                    "affinity_50": sum(nm <= 50.0 for nm, _ in scores),
                    "affinity_500": sum(nm <= 500.0 for nm, _ in scores),
                    "min_rank": min(rank for _, rank in scores),
                    "min_affinity_nm": min(nm for nm, _ in scores),
                    "added_aa": len(linker["sequence"]), "peptide_occurrences": len(peptides),
                    "intended_proc_score_x10000": 0,
                    "max_pair_mean_kd": hydro["max_mean_kd"],
                    "hydrophobicity_pass": category != "internal" or hydro["status"] == "PASS",
                })
                if mix is not None:
                    edges[-1]["mix_strong_rank"] = sum(
                        float(mix[peptide]["%Rank_" + allele[4:].replace("*", "").replace(":", "")]) <= 0.5
                        for peptide in peptides for allele in alleles
                    )
                if prime is not None:
                    edges[-1]["prime_low_rank"] = sum(
                        float(prime[peptide]["%Rank_" + allele[4:].replace("*", "").replace(":", "")]) <= 0.5
                        for peptide in peptides for allele in alleles
                    )
                if c2_map is not None:
                    c2_sc, c2_mc, c2_all = junction_class2_counts(
                        i, j, left, linker["sequence"], right, n, c2_map
                    )
                    edges[-1].update(
                        mixmhc2_strong_core=c2_sc,
                        mixmhc2_mild_core=c2_mc,
                        mixmhc2_all_mild=c2_all,
                    )
                if i < n:
                    lead = leads[i]
                    if lead["needs_C_neighbor"] and not lead["needs_N_neighbor"]:
                        a, b = lead["start_in_region"], lead["end_in_region"]
                        queries.append((lead["peptide"], left[a - 5:a],
                                        (left[b:] + linker["sequence"] + right)[:5]))
                        query_edge_ids.append(edge_id)
                if j < n:
                    lead = leads[j]
                    if lead["needs_N_neighbor"] and not lead["needs_C_neighbor"]:
                        a, b = lead["start_in_region"], lead["end_in_region"]
                        queries.append((lead["peptide"], (left + linker["sequence"] + right[:a])[-5:],
                                        right[b:b + 5]))
                        query_edge_ids.append(edge_id)
    edge_by_id = {edge["edge_id"]: edge for edge in edges}
    if queries:
        scores = processor.predict(
            peptides=[q[0] for q in queries], n_flanks=[q[1] for q in queries],
            c_flanks=[q[2] for q in queries],
        )
        for edge_id, score in zip(query_edge_ids, scores):
            edge_by_id[edge_id]["intended_proc_score_x10000"] += round(float(score) * 10000)
    triplets, triple_queries = [], []
    allowed = eligible_options(edges)
    for index, lead in enumerate(leads):
        if not (lead["needs_N_neighbor"] and lead["needs_C_neighbor"]):
            continue
        sequence = regions[index]["sequence"]
        a, b = lead["start_in_region"], lead["end_in_region"]
        incoming = [edge for edge in allowed if edge["right_index"] == index]
        outgoing = [edge for edge in allowed if edge["left_index"] == index]
        for before in incoming:
            for after in outgoing:
                if before["left_index"] == after["right_index"]:
                    continue
                nflank = (before["left"] + before["linker_sequence"] + sequence[:a])[-5:]
                cflank = (sequence[b:] + after["linker_sequence"] + after["right"])[:5]
                triple_queries.append((lead["peptide"], nflank, cflank))
                triplets.append({
                    "region_index": index, "incoming_edge_id": before["edge_id"],
                    "outgoing_edge_id": after["edge_id"], "peptide": lead["peptide"],
                    "n_flank": nflank, "c_flank": cflank,
                })
    if triple_queries:
        scores = processor.predict(
            peptides=[q[0] for q in triple_queries], n_flanks=[q[1] for q in triple_queries],
            c_flanks=[q[2] for q in triple_queries],
        )
        for row, score in zip(triplets, scores):
            row["processing_score_x10000"] = round(float(score) * 10000)
    return edges, triplets, leads


def emit_assembly(output, regions, path, solver, sec, mitd, alleles, values):
    output = Path(output)
    directory = output / "construct/protein"
    directory.mkdir(parents=True, exist_ok=True)
    sequence, parts, joins = assemble(path, regions, sec, mitd)
    predictions, totals = assess(sequence, joins, path, alleles, values)
    windows, hydro = protein_screen(sequence, parts)
    if hydro["status"] != "PASS":
        raise ValueError("Assembled hydropathy heuristic failed")
    variant = variant_for_regions(regions)
    identifiers = export_identifiers(variant)
    fasta(directory / "cassette.fasta", identifiers["protein_fasta"], sequence)
    table(directory / "components.tsv", parts)
    table(directory / "junctions.tsv", joins)
    table(directory / "boundary-predictions.tsv", predictions)
    table(directory / "hydrophobicity-windows.tsv", windows)
    save(directory / "hydrophobicity.json", hydro)
    table(directory / "antigen-order.tsv", [
        {"position": i + 1, **{key: regions[edge["right_index"]][key]
                              for key in ("antigen_id", "region_id", "gene", "sequence")}}
        for i, edge in enumerate(path[:-1])
    ])
    save(directory / "construct.json", {
        "name": "design-4-" + variant.lower() + "-research", "regions": len(regions),
        "construct_variant": variant, "export_identifiers": identifiers,
        "protein_length_aa": len(sequence), "SEC": sec, "MITD": mitd,
        "solver": solver, "totals": totals, "hydrophobicity": hydro,
        "linker_counts": dict(Counter(edge["linker"] for edge in path)),
        "selected_regions_sha256": sha(output / "selection/selected.tsv"),
        "clinical_release_ready": False,
    })
    return sequence, parts, joins, predictions, totals


def main():
    """Standalone 5-phase MHCflurry-only solver stage; the published 9-phase multi-model consensus construct is built by pipeline/consensus.py and pipeline/recheck_construct.py."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--models", type=Path, required=True,
                        help="Local MHCflurry 2.3.0 class-I presentation model directory.")
    parser.add_argument("--seconds", type=int, default=600)
    args = parser.parse_args()
    os.environ["OMP_NUM_THREADS"] = "2"
    import torch
    from mhcflurry import Class1AffinityPredictor, Class1PresentationPredictor
    torch.set_num_threads(2)
    regions = load(HERE / "selection/selected.json")
    reference = load(HERE / "references/peptide-matches.json")
    linkers = load(HERE / "construct/linkers.json")["choices"]
    architecture = load(HERE / "construct/architecture.json")
    sec, mitd = (architecture["elements"][key]["aa"] for key in ("SEC", "MITD"))
    alleles = load(HERE / "exports/source-case.json")["hla_class_I"]
    n = len(regions)
    junctions = set()
    for i in range(n + 1):
        for j in range(n + 1):
            if i != j:
                left = sec if i == n else regions[i]["sequence"]
                right = mitd if j == n else regions[j]["sequence"]
                for linker in linkers:
                    junctions.update(junction_peptides(left, linker["sequence"], right))
    targets = {peptide for region in regions for _, peptide in target_candidates(region, reference)}
    peptides = sorted(junctions | targets)
    print("Fresh affinity inference:", len(junctions), "junction peptides;", len(peptides), "total queries", flush=True)
    affinity_predictor = Class1AffinityPredictor.load(str(args.models / "affinity_predictor"))
    presentation_predictor = Class1PresentationPredictor.load(str(args.models))
    values = infer_affinity(affinity_predictor, peptides, alleles)
    cache = HERE / "construct/all-junction-affinity.tsv.gz"
    table(cache, (
        {"peptide": peptide, "allele": allele,
         "affinity_nm": values[peptide, allele][0],
         "affinity_percentile": values[peptide, allele][1]}
        for peptide in peptides for allele in alleles
    ), fields=["peptide", "allele", "affinity_nm", "affinity_percentile"])
    edges, triplets, leads = create_edges(
        regions, linkers, sec, mitd, alleles, values,
        presentation_predictor.processing_predictor_with_flanks, reference,
    )
    table(HERE / "construct/protein/edge-costs.tsv", edges)
    table(HERE / "construct/protein/processing-triplets.tsv", triplets,
          fields=["region_index", "incoming_edge_id", "outgoing_edge_id", "peptide",
                  "n_flank", "c_flank", "processing_score_x10000"])
    save(HERE / "construct/protein/lead-targets.json", leads)
    print("Solving", len(edges), "edges and", len(triplets), "two-neighbor processing terms", flush=True)
    path, solver = solve(n, edges, triplets, args.seconds)
    sequence, parts, joins, prior, totals = emit_assembly(
        HERE, regions, path, solver, sec, mitd, alleles, values
    )
    selected_peptides = sorted({row["peptide"] for row in full_scan(sequence, joins)})
    fresh = infer_affinity(affinity_predictor, selected_peptides, alleles)
    checked, recounted = assess(sequence, joins, path, alleles, fresh)
    differences = []
    for before, after in zip(prior, checked):
        differences.append(abs(after["affinity_nm"] / before["affinity_nm"] - 1))
        if [before["affinity_nm"] <= 50, before["affinity_nm"] <= 500,
            before["affinity_percentile"] <= 0.5, before["affinity_percentile"] <= 2] != [
            after["affinity_nm"] <= 50, after["affinity_nm"] <= 500,
            after["affinity_percentile"] <= 0.5, after["affinity_percentile"] <= 2]:
            raise ValueError("Fresh boundary inference changed a threshold classification")
    if max(differences) >= 1e-4:
        raise ValueError("Fresh boundary inference differs beyond recorded tolerance")
    table(HERE / "construct/selected-inference/predictions.tsv", checked)
    full = presentation_predictor.predict_sequences(
        sequences={"Design4_v2": sequence}, alleles={"Yash": alleles}, result="all"
    )
    table(HERE / "immunogenicity/cassette-class1-presentation.tsv", full.to_dict("records"))
    by_position = {(int(row["pos"]), row["peptide"]): row for row in full.to_dict("records")}
    parts_by_id = {part["region_id"]: part for part in parts if part["kind"] == "antigen"}
    placed_leads, rounded_dynamic = [], 0
    for lead in leads:
        start = parts_by_id[lead["region_id"]]["start_aa"] + lead["start_in_region"]
        prediction = by_position[start, lead["peptide"]]
        placed_leads.append({
            **lead, "cassette_start_aa": start,
            "processing_score": float(prediction["processing_score"]),
            "presentation_score": float(prediction["presentation_score"]),
            "presentation_percentile": float(prediction["presentation_percentile"]),
            "presentation_best_allele": prediction["best_allele"],
        })
        if lead["needs_N_neighbor"] or lead["needs_C_neighbor"]:
            rounded_dynamic += round(float(prediction["processing_score"]) * 10000)
    if abs(rounded_dynamic - solver["phases"][-1]["value"]) > 1:
        raise ValueError("Processing objective does not match the placed full-cassette flanks")
    table(HERE / "immunogenicity/placed-class1-leads.tsv", placed_leads)
    protein_dir = HERE / "construct/protein"
    save(HERE / "construct/final-check.json", {
        "status": "FRESH_FULL_JUNCTION_CHECK_PASSED",
        "protein_sha256": sha(protein_dir / "cassette.fasta"),
        "components_sha256": sha(protein_dir / "components.tsv"),
        "protein_length_aa": len(sequence), "antigen_regions": n,
        "all_antigens_exactly_once": True, "totals": recounted,
        "hydrophobicity": load(protein_dir / "hydrophobicity.json"),
        "max_relative_affinity_difference": max(differences),
        "threshold_classifications_identical": True,
        "processing_objective_matches_full_cassette": True,
        "scope": "Fresh repeat model inference, not experimental processing, presentation, safety, or efficacy.",
    })
    model_files = {
        str(path.relative_to(args.models)): sha(path)
        for path in sorted(args.models.rglob("*")) if path.is_file()
    }
    save(HERE / "construct/inference-receipt.json", {
        "library": "mhcflurry", "library_version": version("mhcflurry"),
        "model_bundle_label": "2.3.0", "model_files_sha256": model_files,
        "fresh_junction_and_target_queries": len(peptides),
        "fresh_selected_junction_queries": len(selected_peptides),
        "alleles": alleles, "lengths": list(LENGTHS),
        "affinity_cache_sha256": sha(cache),
        "weights_redistributed": False, "new_cloud_jobs": False,
    })
    print("Saved model-verified cassette:", len(sequence), "aa;", totals, flush=True)


if __name__ == "__main__":
    main()
