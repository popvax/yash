"""Extend the prediction cache with terminal join peptides (SEC->antigen and
antigen->MITD) and assemble the full protein with pinned terminal elements.

Follows the lead design's architecture: patient-supported HLA-C*01:02 SEC
signal peptide fixed first, HLA-B*15:08 MITD fixed last, jointly optimized
internal order and linker choices using the same lexicographic CP-SAT
objective as the maintained solver.
"""
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
from yash_vaccine.constructs.boundaries import junction_peptides, LENGTHS  # noqa: E402
from yash_vaccine.constructs import optimize as base_optimize  # noqa: E402
from yash_vaccine.hydrophobicity import summarize, write_protein_screen  # noqa: E402
from yash_vaccine.io import write_tsv, write_json  # noqa: E402

MODELS = "/root/.local/share/mhcflurry/4/2.3.0/models_class1_pan/models.combined"
ARCHITECTURE = json.loads((ROOT / "designs/design-2/construct/architecture.json").read_text())
SEC = ARCHITECTURE["elements"]["SEC"]["aa"]
MITD = ARCHITECTURE["elements"]["MITD"]["aa"]


def extend_cache():
    regions = list(csv.DictReader((HERE / "regions.tsv").open(), delimiter="\t"))
    linkers = json.loads((ROOT / "config/linkers.json").read_text())["choices"]
    alleles = json.loads((ROOT / "config/case.json").read_text())["hla_class_I"]
    cache = np.load(HERE / "prediction-cache.npz")
    peptides = list(cache["peptides"])
    index = {p: i for i, p in enumerate(peptides)}
    terminal_queries = set()
    for row in regions:
        terminal_queries.update(junction_peptides(SEC, "", row["sequence"]))
        terminal_queries.update(junction_peptides(row["sequence"], "", MITD))
        for linker in linkers:
            terminal_queries.update(junction_peptides(SEC, linker["sequence"], row["sequence"]))
            terminal_queries.update(junction_peptides(row["sequence"], linker["sequence"], MITD))
    missing = sorted(terminal_queries - set(index))
    print(f"terminal join peptides: {len(terminal_queries)} ({len(missing)} new)", flush=True)
    if missing:
        from mhcflurry import Class1AffinityPredictor
        predictor = Class1AffinityPredictor.load(MODELS)
        affinity = cache["affinity_nm"]
        percentile = cache["affinity_percentile"]
        new_aff = np.empty((len(missing), len(alleles)))
        new_pct = np.empty((len(missing), len(alleles)))
        t0 = time.time()
        for j, allele in enumerate(alleles):
            frame = predictor.predict_to_dataframe(alleles=[allele] * len(missing), peptides=missing)
            new_aff[:, j] = frame["prediction"].to_numpy()
            new_pct[:, j] = frame["prediction_percentile"].to_numpy()
        print(f"predicted {len(missing)} x {len(alleles)} in {time.time()-t0:.0f}s", flush=True)
        affinity = np.vstack([affinity, new_aff])
        percentile = np.vstack([percentile, new_pct])
        peptides = peptides + missing
        np.savez(HERE / "prediction-cache.npz", peptides=np.array(peptides),
                 alleles=np.array(alleles), affinity_nm=affinity,
                 affinity_percentile=percentile)
        print("cache extended", flush=True)


def cost(queries, index, affinity, rank, added):
    ix = [index[p] for p in queries]
    a = affinity[ix]; r = rank[ix]
    return dict(strong_rank=int((r <= .5).sum()), mild_rank=int((r <= 2).sum()),
                affinity_500=int((a <= 500).sum()), affinity_50=int((a <= 50).sum()),
                min_rank=float(r.min()) if ix else 100.,
                min_affinity_nm=float(a.min()) if ix else '',
                added_aa=added, peptide_occurrences=len(ix))


def main():
    extend_cache()
    cache = np.load(HERE / "prediction-cache.npz")
    peptides = list(cache["peptides"])
    index = {p: i for i, p in enumerate(peptides)}
    affinity = cache["affinity_nm"]; rank = cache["affinity_percentile"]
    alleles = list(cache["alleles"])
    regions = list(csv.DictReader((HERE / "regions.tsv").open(), delimiter="\t"))
    linkers = json.loads((ROOT / "config/linkers.json").read_text())["choices"]
    n = len(regions)

    # node 0 = SEC (pinned start), nodes 1..n = antigens, node n+1 = MITD (pinned end)
    sec_starts = cost(junction_peptides(SEC, "", regions[0]["sequence"]), index, affinity, rank, 0)
    edges = []
    for i, row in enumerate(regions):
        for j, other in enumerate(regions):
            if i != j:
                for linker in linkers:
                    queries = junction_peptides(row["sequence"], linker["sequence"], other["sequence"])
                    edge = dict(left=i + 1, right=j + 1, left_id=row["antigen_id"],
                                right_id=other["antigen_id"], linker=linker["name"],
                                sequence=linker["sequence"], category="internal",
                                **cost(queries, index, affinity, rank, len(linker["sequence"])))
                    screen = summarize(row["sequence"] + edge["sequence"] + other["sequence"])
                    edge["hydrophobicity_pass"] = screen["status"] == "PASS"
                    edges.append(edge)
        for linker in linkers:
            queries = junction_peptides(SEC, linker["sequence"], row["sequence"])
            edge = dict(left=0, right=i + 1, left_id="SEC", right_id=row["antigen_id"],
                        linker=linker["name"], sequence=linker["sequence"], category="SEC",
                        **cost(queries, index, affinity, rank, len(linker["sequence"])))
            # fixed routing elements are exempt from the edge hydrophobicity gate;
            # boundary windows are separately reported by the protein-level screen
            edge["hydrophobicity_pass"] = True
            edges.append(edge)
            queries = junction_peptides(row["sequence"], linker["sequence"], MITD)
            edge = dict(left=i + 1, right=n + 1, left_id=row["antigen_id"], right_id="MITD",
                        linker=linker["name"], sequence=linker["sequence"], category="MITD",
                        **cost(queries, index, affinity, rank, len(linker["sequence"])))
            edge["hydrophobicity_pass"] = True
            edges.append(edge)
    edges = [e for e in edges if e["hydrophobicity_pass"]]
    print(f"edges after hydrophobicity filter: {len(edges)}", flush=True)

    # pinned-terminal circuit: SEC (node 0) is the only permitted start and
    # MITD (node n+1) the only permitted end.
    starts = []
    for i in range(n):
        starts.append(cost(junction_peptides(SEC, "", regions[i]["sequence"]), index, affinity, rank, 0))
    order, chosen, meta = solve_pinned(n + 2, edges, sec_starts, seconds=300)
    next_edge = {(e["left"], e["right"]): e for e in chosen}
    # export full protein: SEC + antigens with linkers + MITD
    sequence = SEC
    parts = [dict(kind="signal_peptide", name="HLA_C_01_02_SEC", aa=SEC, start_aa=0, end_aa=len(SEC))]
    for p in range(len(order) - 1):
        edge = next_edge[(order[p], order[p + 1])]
        if edge["sequence"]:
            parts.append(dict(kind="linker", name=edge["linker"], aa=edge["sequence"],
                              start_aa=len(sequence), end_aa=len(sequence) + len(edge["sequence"])))
            sequence += edge["sequence"]
        right = order[p + 1]
        if right not in (0, n + 1):
            row = regions[right - 1]
            parts.append(dict(kind="antigen", name=row["antigen_id"], aa=row["sequence"],
                              start_aa=len(sequence), end_aa=len(sequence) + len(row["sequence"])))
            sequence += row["sequence"]
    parts.append(dict(kind="mitd", name="HLA_B_1508_MITD", aa=MITD,
                      start_aa=len(sequence), end_aa=len(sequence) + len(MITD)))
    sequence += MITD
    assert all(sequence[p["start_aa"]:p["end_aa"]] == p["aa"] for p in parts)
    assert len([p for p in parts if p["kind"] == "antigen"]) == n
    out = HERE / "construct"
    out.mkdir(exist_ok=True)
    hydrophobicity = write_protein_screen(out, sequence, parts)
    pred = []
    for row in base_scan(sequence, parts):
        ix = index[row["peptide"]]
        pred.extend(dict(**row, allele=h, affinity_nm=float(affinity[ix, j]),
                         affinity_percentile=float(rank[ix, j])) for j, h in enumerate(alleles))
    (out / "cassette.fasta").write_text(
        ">specificity_first_RESEARCH_DRAFT\n" +
        "\n".join(sequence[p:p + 80] for p in range(0, len(sequence), 80)) + "\n")
    write_tsv(out / "components.tsv", parts)
    write_tsv(out / "boundary-predictions.tsv", pred)
    write_json(out / "construct.json", dict(
        name="specificity-first", regions=n, protein_length_aa=len(sequence),
        sec_aa=SEC, mitd_aa=MITD,
        order=[regions[i - 1]["antigen_id"] if 0 < i < n + 1 else ("SEC" if i == 0 else "MITD") for i in order],
        linkers=[next_edge[(order[p], order[p + 1])]["linker"] for p in range(len(order) - 1)],
        boundary_peptide_occurrences=len(pred),
        strong_rank=sum(x["affinity_percentile"] <= .5 for x in pred),
        mild_rank=sum(x["affinity_percentile"] <= 2 for x in pred),
        affinity_500=sum(x["affinity_nm"] <= 500 for x in pred),
        hydrophobicity=hydrophobicity, solver=meta))
    print("protein length:", len(sequence))
    print("order:", [regions[i - 1]["antigen_id"] if 0 < i < n + 1 else ("SEC" if i == 0 else "MITD") for i in order])


def base_scan(sequence, parts):
    from yash_vaccine.constructs.boundaries import scan_boundaries
    return list(scan_boundaries(sequence, parts))


def solve_pinned(n, edges, sec_starts, seconds=300):
    """Circuit over n nodes; node 0 (SEC) pinned as start, node n-1 (MITD) as end."""
    from ortools.sat.python import cp_model
    edge_groups = {}
    for edge in edges:
        edge_groups.setdefault((edge["left"], edge["right"]), []).append(edge)
    allowed = []
    for pair, options in sorted(edge_groups.items()):
        best = min((x["strong_rank"], x["mild_rank"]) for x in options)
        options = [x for x in options if (x["strong_rank"], x["mild_rank"]) == best]
        for x in options:
            if not any(y["min_rank"] >= x["min_rank"] and y["added_aa"] <= x["added_aa"]
                       and (y["min_rank"] > x["min_rank"] or y["added_aa"] < x["added_aa"])
                       for y in options):
                allowed.append(x)
    model = cp_model.CpModel()
    arcs = []; items = []; dummy = n
    for k, edge in enumerate(allowed):
        lit = model.new_bool_var("edge_" + str(k))
        arcs.append((edge["left"], edge["right"], lit)); items.append((lit, edge))
    # SEC (node 0) is the circuit start: dummy->0 only. MITD (node n-1) the end: ->dummy only.
    start_lit = model.new_bool_var("start_sec")
    arcs.append((dummy, 0, start_lit)); items.append((start_lit, sec_starts))
    end_lit = model.new_bool_var("end_mitd")
    arcs.append((n - 1, dummy, end_lit))
    model.add(start_lit == 1)
    model.add_circuit(arcs)
    ranks = sorted({x["min_rank"] for _, x in items}); rank_index = {v: i for i, v in enumerate(ranks)}
    quality = model.new_int_var(0, len(ranks) - 1, "worst_rank_ordinal")
    for lit, edge in items:
        model.add(quality <= rank_index[edge["min_rank"]] + len(ranks) * (1 - lit))
    expressions = [sum(lit * x["strong_rank"] for lit, x in items),
                   sum(lit * x["mild_rank"] for lit, x in items), quality,
                   sum(lit * x["added_aa"] for lit, x in items)]
    phases = []; solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 1; solver.parameters.random_seed = 31
    solver.parameters.max_time_in_seconds = seconds
    for name, expr, maximize in zip(["strong_count", "mild_count", "worst_rank_ordinal", "added_residues"],
                                    expressions, [False, False, True, False]):
        if maximize:
            model.maximize(expr)
        else:
            model.minimize(expr)
        begin = time.perf_counter(); status = solver.solve(model)
        if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            raise RuntimeError(solver.status_name(status))
        value = solver.value(expr)
        phases.append(dict(objective=name, status=solver.status_name(status), value=value,
                           bound=solver.best_objective_bound, elapsed_seconds=time.perf_counter() - begin))
        print("SOLVER", name, phases[-1], flush=True)
        if status != cp_model.OPTIMAL:
            break
        model.add(expr == value)
        model.clear_hints()
        for lit, _ in items:
            model.add_hint(lit, solver.value(lit))
    selected = [edge for lit, edge in items if solver.value(lit) and "left" in edge]
    next_edge = {x["left"]: x for x in selected}
    order = [0]; chosen = []
    while order[-1] in next_edge:
        edge = next_edge[order[-1]]
        chosen.append(edge); order.append(edge["right"])
    assert len(order) == n and sorted(order) == list(range(n))
    return order, chosen, dict(phases=phases, all_phases_proven_optimal=all(x["status"] == "OPTIMAL" for x in phases),
                               solver="OR-Tools CP-SAT pinned-terminal", workers=1, seed=31)


if __name__ == "__main__":
    main()
