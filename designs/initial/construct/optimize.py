"""Joint path/linker optimisation with scored SEC and MITD endpoint edges."""
from ortools.sat.python import cp_model


def objective(edges):
    return (sum(e["strong_rank"] for e in edges), sum(e["mild_rank"] for e in edges),
            -min(e["min_rank"] for e in edges), sum(e["added_aa"] for e in edges))


def solve(n, edges, seconds=120):
    # The anchor represents SEC on outgoing arcs and MITD on incoming arcs.
    # One circuit through it and all antigens is one SEC-to-MITD Hamiltonian path.
    groups = {}
    for edge in edges:
        groups.setdefault((edge["left_index"], edge["right_index"]), []).append(edge)
    allowed = []
    for _, options in sorted(groups.items()):
        best = min((x["strong_rank"], x["mild_rank"]) for x in options)
        options = [x for x in options if (x["strong_rank"], x["mild_rank"]) == best]
        for x in options:
            if not any(y["min_rank"] >= x["min_rank"] and y["added_aa"] <= x["added_aa"]
                       and (y["min_rank"] > x["min_rank"] or y["added_aa"] < x["added_aa"])
                       for y in options):
                allowed.append(x)
    model = cp_model.CpModel()
    items = [(model.new_bool_var("edge_" + str(i)), edge) for i, edge in enumerate(allowed)]
    model.add_circuit([(e["left_index"], e["right_index"], lit) for lit, e in items])
    ranks = sorted({e["min_rank"] for _, e in items})
    ix = {rank: i for i, rank in enumerate(ranks)}
    quality = model.new_int_var(0, len(ranks)-1, "worst_rank_ordinal")
    for lit, edge in items:
        model.add(quality <= ix[edge["min_rank"]] + len(ranks)*(1-lit))
    objectives = [sum(lit*e["strong_rank"] for lit,e in items),
                  sum(lit*e["mild_rank"] for lit,e in items), quality,
                  sum(lit*e["added_aa"] for lit,e in items)]
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 1
    solver.parameters.random_seed = 31
    solver.parameters.max_time_in_seconds = seconds
    phases = []
    for name, expr, maximize in zip(
            ["strong_count", "mild_count", "worst_rank_ordinal", "added_residues"],
            objectives, [False, False, True, False]):
        model.maximize(expr) if maximize else model.minimize(expr)
        status = solver.solve(model)
        if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            raise RuntimeError("No solution: " + solver.status_name(status))
        value = solver.value(expr)
        phase = dict(objective=name, status=solver.status_name(status), value=value,
                     bound=solver.best_objective_bound)
        phases.append(phase)
        print("SOLVER", phase, flush=True)
        if status != cp_model.OPTIMAL:
            raise RuntimeError("Not proven optimal; increase --seconds before exporting this design")
        model.add(expr == value)
        model.clear_hints()
        for lit, _ in items:
            model.add_hint(lit, solver.value(lit))
    selected = [edge for lit,edge in items if solver.value(lit)]
    outgoing = {e["left_index"]: e for e in selected}
    node = n
    path = []
    for _ in range(n+1):
        edge = outgoing[node]
        path.append(edge)
        node = edge["right_index"]
    if node != n or sorted(e["right_index"] for e in path[:-1]) != list(range(n)):
        raise ValueError("Solver path does not contain every antigen exactly once")
    return path, dict(solver="OR-Tools CP-SAT", workers=1, seed=31,
        objectives=["minimise count rank <=0.5%", "minimise count rank <=2%",
                    "maximise worst junction percentile rank", "minimise inserted linker residues"],
        phases=phases, all_phases_proven_optimal=True, considered_options=len(edges),
        nondominated_options=len(allowed), objective_tuple=objective(path))
