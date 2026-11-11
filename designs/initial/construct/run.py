"""Optimise the fixed antigen set and routing elements; offline cached model scores."""
import argparse
from collections import Counter
import importlib.metadata
from pathlib import Path
import numpy as np
from common import (HERE, REPO, sha, read_json, write_json, rows, write_rows, fasta,
                    inputs, junction_peptides, terminal_request, load_predictions,
                    full_scan, request)
from optimize import solve
from yash_vaccine.hydrophobicity import summarize, write_protein_screen


def prediction_lookup(terminal, alleles):
    cache_path = REPO / "data/checkpoints/construct-predictions.npz"
    expected = read_json(REPO / "data/checkpoints/manifest.json")[cache_path.name]["sha256"]
    if sha(cache_path) != expected:
        raise ValueError("Frozen affinity cache changed")
    cache = np.load(cache_path)
    if cache["alleles"].tolist() != alleles:
        raise ValueError("Cache HLA mismatch")
    index = {p: i for i,p in enumerate(cache["peptides"].tolist())}
    affinity, rank = cache["affinity_nm"], cache["affinity_percentile"]
    assert len(index) == len(affinity) == len(rank)
    assert np.isfinite(affinity).all() and (affinity > 0).all()
    assert np.isfinite(rank).all() and ((rank >= 0) & (rank <= 100)).all()
    differences, rdifferences, overlap = [], [], 0
    for (peptide, allele), (a, r) in terminal.items():
        if peptide not in index:
            continue
        i, j = index[peptide], alleles.index(allele)
        ca, cr = float(affinity[i,j]), float(rank[i,j])
        differences.append(abs(a/ca-1)); rdifferences.append(abs(r-cr)); overlap += 1
        if [a<=50,a<=500,r<=.5,r<=2] != [ca<=50,ca<=500,cr<=.5,cr<=2]:
            raise ValueError("Repeated prediction changes a threshold classification")
    if max(differences, default=0) >= 1e-4 or max(rdifferences, default=0) >= 1e-4:
        raise ValueError("New terminal predictions disagree with shared cached peptides")
    def lookup(peptide, h):
        if peptide in index:
            i = index[peptide]
            return float(affinity[i,h]), float(rank[i,h])
        return terminal[peptide, alleles[h]]
    return lookup, dict(cache_sha256=expected, cache_peptides=len(index),
        fresh_terminal_pairs=len(terminal), overlap_pairs_checked=overlap,
        overlap_max_relative_affinity_difference=max(differences, default=0),
        overlap_max_absolute_rank_difference=max(rdifferences, default=0),
        overlapping_threshold_classifications_unchanged=True)


def assemble(path, regions, sec, mitd):
    sequence, components, joins = "", [], []
    def add(kind, name, aa, source=None):
        nonlocal sequence
        source = source or {}
        components.append(dict(kind=kind, name=name, aa=aa, start_aa=len(sequence),
            end_aa=len(sequence)+len(aa), region_id=source.get("region_id", ""),
            event_id=source.get("event_id", ""), gene=source.get("gene", "")))
        sequence += aa
    add("signal_peptide", "HLA_C_01_02_SEC", sec)
    for i, edge in enumerate(path):
        begin = len(sequence)
        add("linker", edge["linker"], edge["linker_sequence"])
        end = len(sequence)
        joins.append(dict(junction_index=i, category=edge["category"],
            left_id=edge["left_id"], right_id=edge["right_id"], linker=edge["linker"],
            linker_sequence=edge["linker_sequence"], linker_start_aa=begin, linker_end_aa=end,
            left=edge["left"], right=edge["right"]))
        if edge["right_index"] == len(regions):
            add("mitd", "HLA_B_15_08_MITD", mitd)
        else:
            r = regions[edge["right_index"]]
            add("antigen", r["antigen_id"], r["sequence"], r)
    assert sequence.startswith(sec) and sequence.endswith(mitd)
    assert len(sequence) == sum(len(r["sequence"]) for r in regions)+len(sec)+len(mitd)+sum(e["added_aa"] for e in path)
    assert sorted((p["region_id"],p["aa"]) for p in components if p["kind"]=="antigen") == \
           sorted((r["region_id"],r["sequence"]) for r in regions)
    assert all(sequence[p["start_aa"]:p["end_aa"]] == p["aa"] for p in components)
    return sequence, components, joins


def summarize_predictions(predictions):
    return dict(peptide_occurrences=len(predictions)//6, peptide_hla_rows=len(predictions),
        strong_rank=sum(r["affinity_percentile"]<=.5 for r in predictions),
        mild_rank=sum(r["affinity_percentile"]<=2 for r in predictions),
        affinity_50=sum(r["affinity_nm"]<=50 for r in predictions),
        affinity_500=sum(r["affinity_nm"]<=500 for r in predictions),
        min_rank=min(r["affinity_percentile"] for r in predictions),
        min_affinity_nm=min(r["affinity_nm"] for r in predictions))


def assess(sequence, joins, path, alleles, lookup):
    occurrences = full_scan(sequence, joins)
    predicted = [dict(**row, allele=allele, affinity_nm=lookup(row["peptide"],h)[0],
                      affinity_percentile=lookup(row["peptide"],h)[1])
                 for row in occurrences for h,allele in enumerate(alleles)]
    for i, edge in enumerate(path):
        actual = Counter(r["peptide"] for r in occurrences if r["junction_index"] == i)
        expected = Counter(junction_peptides(edge["left"],edge["linker_sequence"],edge["right"]))
        if actual != expected:
            raise ValueError("Full-protein scan disagrees with pairwise junction enumeration")
    totals = summarize_predictions(predicted)
    for key in ["strong_rank", "mild_rank", "affinity_50", "affinity_500", "peptide_occurrences"]:
        if totals[key] != sum(e[key] for e in path):
            raise ValueError("Full-chain cost mismatch: " + key)
    by_category = [dict(category=category, **summarize_predictions([r for r in predicted if r["category"]==category]))
                   for category in ["SEC", "internal", "MITD"]]
    return predicted, totals, by_category


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--terminal-predictions", type=Path, default=HERE / "terminal-inference")
    p.add_argument("--output", type=Path, default=HERE / "protein")
    p.add_argument("--seconds", type=int, default=120)
    args = p.parse_args()
    if args.output.exists():
        p.error("Choose a new output directory")
    versions = {x:importlib.metadata.version(x) for x in ["numpy","ortools"]}
    if versions != {"numpy":"1.26.4","ortools":"9.14.6206"}:
        raise ValueError("Use the recorded project environment: " + str(versions))
    regions, linkers, sec, mitd, alleles, hashes = inputs()
    terminal = load_predictions(args.terminal_predictions, terminal_request())
    lookup, cache_meta = prediction_lookup(terminal, alleles)
    n = len(regions)
    def edge(left_i, right_i, linker):
        left = sec if left_i==n else regions[left_i]["sequence"]
        right = mitd if right_i==n else regions[right_i]["sequence"]
        queries = junction_peptides(left, linker["sequence"], right)
        values = [lookup(q,h) for q in queries for h in range(len(alleles))]
        return dict(left_index=left_i, right_index=right_i,
            left_id="SEC" if left_i==n else regions[left_i]["region_id"],
            right_id="MITD" if right_i==n else regions[right_i]["region_id"],
            category="SEC" if left_i==n else "MITD" if right_i==n else "internal",
            left=left, right=right, linker=linker["name"], linker_sequence=linker["sequence"],
            strong_rank=sum(r<=.5 for a,r in values), mild_rank=sum(r<=2 for a,r in values),
            affinity_50=sum(a<=50 for a,r in values), affinity_500=sum(a<=500 for a,r in values),
            min_rank=min(r for a,r in values), min_affinity_nm=min(a for a,r in values),
            added_aa=len(linker["sequence"]), peptide_occurrences=len(queries))
    edges = [edge(i,j,l) for i in range(n) for j in range(n) if i!=j for l in linkers]
    edges += [edge(n,i,l) for i in range(n) for l in linkers]
    edges += [edge(i,n,l) for i in range(n) for l in linkers]
    for option in edges:
        screen=summarize(option['left']+option['linker_sequence']+option['right'])
        # Routing-domain overlaps are visible in the final unmasked full-chain
        # scan, separately from the hard gate on internal antigen/linker windows.
        option['max_pair_mean_kd']=screen['max_mean_kd']
        option['hydrophobicity_pass']=option['category']!='internal' or screen['status']=='PASS'
    eligible_edges=[option for option in edges if option['hydrophobicity_pass']]
    path, solver = solve(n,eligible_edges,args.seconds)
    solver['hydrophobicity_excluded_options']=len(edges)-len(eligible_edges)
    solver['total_evaluated_options']=len(edges)
    sequence, components, joins = assemble(path,regions,sec,mitd)
    predictions, totals, categories = assess(sequence,joins,path,alleles,lookup)
    output=args.output; output.mkdir(parents=True)
    hydrophobicity=write_protein_screen(output,sequence,components)
    fasta(output/"cassette.fasta", "Yash_joint_SEC_antigens_MITD", sequence)
    write_rows(output/"components.tsv",components)
    write_rows(output/"junctions.tsv",joins)
    write_rows(output/"boundary-predictions.tsv",predictions)
    write_rows(output/"by-category.tsv",categories)
    write_rows(output/"edge-costs.tsv",edges)
    order = [regions[e["right_index"]] for e in path[:-1]]
    write_rows(output/"antigen-order.tsv",[dict(position=i+1,antigen_id=r["antigen_id"],region_id=r["region_id"],
        gene=r["gene"],sequence=r["sequence"]) for i,r in enumerate(order)])
    required = {p for e in edges for p in junction_peptides(e["left"],e["linker_sequence"],e["right"])}
    construct = dict(name="initial", protein_length_aa=len(sequence), antigen_regions=n,
        SEC=sec, MITD=mitd, antigen_sequences_unchanged=True, exact_routing_elements_preserved=True,
        antigen_order=[r["region_id"] for r in order], first_gene=order[0]["gene"], last_gene=order[-1]["gene"],
        linker_choices=[e["linker"] for e in path], inserted_linker_residues=sum(e["added_aa"] for e in path),
        linker_counts=dict(Counter(e["linker"] for e in path)), solver=solver, prediction_cache=cache_meta,
        required_query_peptides=len(required), required_query_peptide_hla_pairs=len(required)*len(alleles),
        optimized=dict(totals=totals,by_category=categories),
        hydrophobicity=hydrophobicity,
        input_sha256={**hashes,"terminal_inference":sha(args.terminal_predictions/"inference.json")},
        code_sha256={**{name:sha(HERE/name) for name in ["run.py","common.py","optimize.py","input-locations.json","architecture.json"]},
                     "src/yash_vaccine/hydrophobicity.py":sha(REPO/'src/yash_vaccine/hydrophobicity.py')}, packages=versions,
        scope="Junction affinity objective with hydropathy feasibility gate; no cleavage, presentation or immunogenicity optimization.")
    write_json(output/"construct.json",construct)
    verify_request = request([dict(name=str(j["junction_index"]),left=j["left"],linker=j["linker_sequence"],right=j["right"])
        for j in joins], "Fresh verification of all selected junctions", {"cassette.fasta":sha(output/"cassette.fasta")})
    write_json(output/"verification-request.json",verify_request)
    names=["cassette.fasta","components.tsv","junctions.tsv","boundary-predictions.tsv","by-category.tsv",
           "edge-costs.tsv","antigen-order.tsv","construct.json","verification-request.json", "hydrophobicity-windows.tsv","hydrophobicity.json"]
    write_json(output/"verification.json",dict(status="CACHED_ASSEMBLY_VERIFIED",all_antigens_exactly_once=True,antigen_regions=n,
        all_junctions_in_objective=True,scored_junctions=n+1,independent_full_chain_scan_verified=True,
        all_solver_phases_optimal=solver["all_phases_proven_optimal"],
        hydrophobicity_pass=True,
        fresh_selected_junction_inference="separate final-check.json after verify.py",
        output_sha256={name:sha(output/name) for name in names}))
    print(__import__("json").dumps({"length":len(sequence),"optimized":categories,
        "first_gene":order[0]["gene"],"last_gene":order[-1]["gene"],"terminal_linkers":[path[0]["linker"],path[-1]["linker"]]},indent=2))


if __name__ == "__main__":
    main()
