"""Exact lexicographic extension of native Vaxrank junction burden.

Native Vaxrank swaps linkers at a fixed junction. This module also chooses the
Hamiltonian antigen order. No RNA or antigen-rank score enters junction costs.
"""
import time
from ortools.sat.python import cp_model

def objective(order, chosen, starts):
    costs=[starts[order[0]],*chosen]
    return (sum(x['strong_rank'] for x in costs),sum(x['mild_rank'] for x in costs),
            -min(x['min_rank'] for x in costs),sum(x['added_aa'] for x in costs))

def solve(n, edges, starts, seconds=120):
    """Separate integer solves avoid arbitrary floating scalar weights.

    The percentile bottleneck uses exact ordering of observed floating values,
    so no percentile rounding affects the optimization.
    """
    edge_groups={}
    for edge in edges:edge_groups.setdefault((edge['left'],edge['right']),[]).append(edge)
    allowed=[]
    for pair, options in sorted(edge_groups.items()):
        # Exact local dominance: larger primary count cannot help any order;
        # with equal strong count, a larger mild count cannot help either.
        best=min((x['strong_rank'],x['mild_rank']) for x in options)
        options=[x for x in options if (x['strong_rank'],x['mild_rank'])==best]
        for x in options:
            if not any(y['min_rank']>=x['min_rank'] and y['added_aa']<=x['added_aa']
                       and (y['min_rank']>x['min_rank'] or y['added_aa']<x['added_aa'])
                       for y in options):allowed.append(x)
    model=cp_model.CpModel();arcs=[];items=[];dummy=n
    for k,edge in enumerate(allowed):
        lit=model.new_bool_var('edge_'+str(k))
        arcs.append((edge['left'],edge['right'],lit));items.append((lit,edge))
    start_items=[]
    for i in range(n):
        first=model.new_bool_var('start_'+str(i));last=model.new_bool_var('end_'+str(i))
        arcs.extend([(dummy,i,first),(i,dummy,last)])
        if not starts[i].get('hydrophobicity_pass',True):model.add(first==0)
        start_items.append((first,starts[i]));items.append((first,starts[i]))
    model.add_circuit(arcs)
    ranks=sorted({x['min_rank'] for _,x in items});rank_index={v:i for i,v in enumerate(ranks)}
    quality=model.new_int_var(0,len(ranks)-1,'worst_rank_ordinal')
    for lit,edge in items:
        model.add(quality<=rank_index[edge['min_rank']]+len(ranks)*(1-lit))
    expressions=[sum(lit*x['strong_rank'] for lit,x in items),
                 sum(lit*x['mild_rank'] for lit,x in items),quality,
                 sum(lit*x['added_aa'] for lit,x in items)]
    phases=[];solver=cp_model.CpSolver()
    solver.parameters.num_search_workers=1;solver.parameters.random_seed=31
    solver.parameters.max_time_in_seconds=seconds
    for name,expr,maximize in zip(['strong_count','mild_count','worst_rank_ordinal','added_residues'],expressions,[False,False,True,False]):
        if maximize:model.maximize(expr)
        else:model.minimize(expr)
        begin=time.perf_counter();status=solver.solve(model)
        if status not in (cp_model.OPTIMAL,cp_model.FEASIBLE):raise RuntimeError(solver.status_name(status))
        value=solver.value(expr)
        phases.append(dict(objective=name,status=solver.status_name(status),value=value,
                           bound=solver.best_objective_bound,elapsed_seconds=time.perf_counter()-begin))
        print('SOLVER',name,phases[-1],flush=True)
        # A feasible incumbent may be exported, but never labelled optimal.
        if status!=cp_model.OPTIMAL:break
        model.add(expr==value)
        model.clear_hints()
        for lit,_ in items:model.add_hint(lit,solver.value(lit))
    selected=[edge for lit,edge in items if solver.value(lit) and 'left' in edge]
    start=next(i for i,(lit,_) in enumerate(start_items) if solver.value(lit))
    next_edge={x['left']:x for x in selected};order=[start];chosen=[]
    while order[-1] in next_edge:
        edge=next_edge[order[-1]];chosen.append(edge);order.append(edge['right'])
    assert len(order)==n and sorted(order)==list(range(n)) and len(chosen)==n-1
    return order,chosen,dict(phases=phases,all_phases_proven_optimal=len(phases)==4 and all(x['status']=='OPTIMAL' for x in phases),
                             pair_linker_choices=len(edges),nondominated_choices=len(allowed),
                             objective_tuple=objective(order,chosen,starts),
                             solver='OR-Tools CP-SAT',workers=1,seed=31)
