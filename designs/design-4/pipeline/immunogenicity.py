"""Prepare, score and summarize event-specific predictions with explicit contexts."""
import argparse
import os
import re
from pathlib import Path

from curate import HERE, overlaps_target, peptide_windows
from dataio import iter_rows, load, read_text_auto, rows, save, sha, table, write_text_auto

CLASS1_ALLELES = ["A0101", "A3303", "B0801", "B1508", "C0102", "C0702"]
CLASS2_ALLELES = [
    "DRB1_03_01", "DRB1_08_02", "DRB3_02_02",
    "DQA1_05_01__DQB1_02_01", "DQA1_04_01__DQB1_04_02",
    "DQA1_04_01__DQB1_02_01", "DPA1_01_03__DPB1_02_01",
    "DPA1_02_01__DPB1_13_01", "DPA1_01_03__DPB1_13_01",
]
DIRECTORY = HERE / "immunogenicity"


def context(sequence, start, end):
    peptide = sequence[start:end]
    return (sequence[max(0, start - 3):start].rjust(3, "-") + peptide[:3]
            + peptide[-3:] + sequence[end:end + 3].ljust(3, "-"))


def prepare():
    regions = load(HERE / "selection/selected.json")
    sequence = "".join(line for line in (HERE / "construct/protein/cassette.fasta").read_text().splitlines()
                       if not line.startswith(">"))
    parts = rows(HERE / "construct/protein/components.tsv")
    parts = [row for row in parts if row["kind"] == "antigen"]
    part_by_id = {row["region_id"]: row for row in parts}
    queries = {row["peptide"] for row in iter_rows(HERE / "construct/all-junction-affinity.tsv.gz")}
    class2_requests, encoded = set(), []
    for length in (8, 9, 10, 11):
        queries.update(sequence[start:start + length] for start in range(len(sequence) - length + 1))
    for start in range(len(sequence) - 14):
        peptide, ctx = sequence[start:start + 15], context(sequence, start, start + 15)
        class2_requests.add((peptide, ctx))
        encoded.append({"peptide": peptide, "context": ctx, "cassette_start_aa": start})
    target_rows, contrasts = [], []
    for region in regions:
        left = int(part_by_id[region["region_id"]]["start_aa"])
        lo, hi = int(region["target_start_in_window"]), int(region["target_end_in_window"])
        for length in (8, 9, 10, 11, 15):
            for start, peptide in peptide_windows(region["sequence"], lo, hi, (length,)):
                target_rows.append({
                    "region_id": region["region_id"], "gene": region["gene"],
                    "scenario": "encoded", "peptide": peptide, "length": length,
                    "start_in_region": start, "cassette_start_aa": left + start,
                    "context": context(sequence, left + start, left + start + length) if length == 15 else "",
                    "target_start_in_region": lo, "target_end_in_region": hi,
                })
        change = re.fullmatch(r"p\.([A-Z])\d+([A-Z])", region["DNA_annotation_protein_change"]) if region["source"] == "MUTATION" else None
        scenarios = []
        if change and hi == lo + 1:
            if region["sequence"][lo] != change.group(2):
                raise ValueError("Mutant sequence disagrees with its source annotation")
            scenarios.append(("single_aa_WT_reversion", lo, change.group(1)))
        if region["gene"] == "CDH1::NRG1":
            sequon = region["sequence"].find("NAT")
            if sequon < 0:
                raise ValueError("Recorded fusion lacks the assessed N-glycosylation sequon")
            scenarios.append(("N_to_D_deglycosylation_unmeasured", sequon, "D"))
        for name, position, replacement in scenarios:
            counterfactual = sequence[:left + position] + replacement + sequence[left + position + 1:]
            for length in (8, 9, 10, 11, 15):
                for start, peptide in peptide_windows(region["sequence"], lo, hi, (length,)):
                    altered = counterfactual[left + start:left + start + length]
                    if length == 15:
                        ctx = context(counterfactual, left + start, left + start + length)
                        class2_requests.add((altered, ctx))
                    else:
                        ctx = ""
                        queries.add(altered)
                        contrasts.append({
                            "region_id": region["region_id"], "gene": region["gene"],
                            "scenario": name, "mutant_peptide": peptide, "counterfactual_peptide": altered,
                            "start_in_region": start,
                        })
                    target_rows.append({
                        "region_id": region["region_id"], "gene": region["gene"],
                        "scenario": name, "peptide": altered, "length": length,
                        "start_in_region": start, "cassette_start_aa": left + start,
                        "context": ctx, "target_start_in_region": lo, "target_end_in_region": hi,
                    })
    request_dir = DIRECTORY / "requests"
    request_dir.mkdir(parents=True, exist_ok=True)
    write_text_auto(request_dir / "class1.txt", "".join(peptide + "\n" for peptide in sorted(queries)))
    write_text_auto(request_dir / "class2-context.txt", "".join(
        peptide + "\t" + ctx + "\n" for peptide, ctx in sorted(class2_requests)
    ))
    write_text_auto(request_dir / "class2-no-context.txt", "".join(
        peptide + "\n" for peptide in sorted({p for p, _ in class2_requests})
    ))
    table(request_dir / "encoded-class2-occurrences.tsv", encoded)
    table(request_dir / "target-occurrences.tsv", target_rows)
    table(request_dir / "counterfactual-class1-pairs.tsv", contrasts)
    save(request_dir / "receipt.json", {
        "class1_alleles": CLASS1_ALLELES, "class2_chain_combinations": CLASS2_ALLELES,
        "class2_pair_phasing_or_surface_expression_verified": False,
        "hla_typing_source_sha256": sha(HERE / "exports/hla-typing-comparison.json"),
        "hla_typing_method": "Recorded T1K tumour DNA WES versus tumour RNA; not matched-normal or clinically validated HLA typing.",
        "DRB3_qualification": "DNA DRB3*02:02 and RNA DRB3*02:256 disagree. DRB3*02:02 scores are a conditional DNA-call scenario, not a concordant allele.",
        "class2_scope": "Nine available definitions, not all potential alpha/beta pairings. Haplotype phase and dimer formation remain unverified.",
        "class1_queries": len(queries), "class2_context_queries": len(class2_requests),
        "context_policy": "Complete engineered cassette; '-' only at actual termini. Counterfactual single-AA reversions are not measured patient-normal protein sequences.",
        "input_sha256": {
            name: sha(request_dir / name) for name in (
                "class1.txt", "class2-context.txt", "class2-no-context.txt",
                "encoded-class2-occurrences.tsv",
                "target-occurrences.tsv", "counterfactual-class1-pairs.tsv",
            )
        },
    })
    print("Prepared", len(queries), "Class-I and", len(class2_requests), "Class-II/context queries")


def wildtype_affinity(models):
    os.environ["OMP_NUM_THREADS"] = "2"
    import torch
    from mhcflurry import Class1AffinityPredictor
    from protein import infer_affinity
    torch.set_num_threads(2)
    alleles = load(HERE / "exports/source-case.json")["hla_class_I"]
    peptides = sorted({row["counterfactual_peptide"]
                       for row in rows(DIRECTORY / "requests/counterfactual-class1-pairs.tsv")})
    predictor = Class1AffinityPredictor.load(str(models / "affinity_predictor"))
    predictions = infer_affinity(predictor, peptides, alleles)
    table(DIRECTORY / "counterfactual-class1-affinity.tsv", [{
        "peptide": peptide, "allele": allele, "affinity_nm": predictions[peptide, allele][0],
        "affinity_percentile": predictions[peptide, allele][1],
    } for peptide in peptides for allele in alleles])


def summarize():
    mix = {r["Peptide"]: r for r in rows(DIRECTORY / "mixmhcpred3.tsv")}
    prime = {r["Peptide"]: r for r in rows(DIRECTORY / "prime21.tsv")}
    c2 = {(r["Peptide"], r["Context"]): r for r in rows(DIRECTORY / "mixmhc2-context.tsv")}
    n2 = {r["Peptide"]: r for r in rows(DIRECTORY / "mixmhc2-no-context.tsv")}
    expected_c1 = set(read_text_auto(DIRECTORY / "requests/class1.txt").splitlines())
    expected_c2 = {tuple(line.split("\t")) for line in
                   read_text_auto(DIRECTORY / "requests/class2-context.txt").splitlines()}
    expected_n2 = set(read_text_auto(DIRECTORY / "requests/class2-no-context.txt").splitlines())
    if set(mix) != expected_c1 or set(prime) != expected_c1 or set(c2) != expected_c2 or set(n2) != expected_n2:
        raise ValueError("Prediction queries do not exactly match the current requests")
    reference = load(HERE / "references/peptide-matches.json")
    occurrences = rows(DIRECTORY / "requests/target-occurrences.tsv")
    regions = load(HERE / "selection/selected.json")
    leads = {r["region_id"]: r for r in rows(DIRECTORY / "placed-class1-leads.tsv")}
    contrasts = rows(DIRECTORY / "requests/counterfactual-class1-pairs.tsv")
    all_affinity = {
        (r["peptide"], r["allele"]): float(r["affinity_nm"])
        for r in iter_rows(HERE / "construct/all-junction-affinity.tsv.gz")
    }
    counterfactual_affinity = {
        (r["peptide"], r["allele"]): float(r["affinity_nm"])
        for r in rows(DIRECTORY / "counterfactual-class1-affinity.tsv")
    }
    summary_rows, class2_rows, contrast_rows = [], [], []
    for region in regions:
        rid = region["region_id"]
        target = [
            r for r in occurrences if r["region_id"] == rid and r["scenario"] == "encoded"
            and (region["source"] == "CTA/SELF" or not reference[r["peptide"]]["occurs"])
        ]
        class1 = [r for r in target if int(r["length"]) != 15]
        class2 = [r for r in target if int(r["length"]) == 15]
        best_context = min(class2, key=lambda r: float(c2[r["peptide"], r["context"]]["%Rank_best"])) if class2 else None
        best_no_context = min(class2, key=lambda r: float(n2[r["peptide"]]["%Rank_best"])) if class2 else None
        lead = leads[rid]
        record = {
            "antigen_id": region["antigen_id"], "region_id": rid, "gene": region["gene"],
            "source": region["source"], "design_tier": region["design_tier"],
            "class1_lead_peptide": lead["peptide"], "class1_lead_allele": lead["allele"],
            "class1_lead_affinity_nm": lead["affinity_nm"],
            "class1_lead_processing_score": lead["processing_score"],
            "class1_lead_presentation_rank": lead["presentation_percentile"],
            "best_target_mixmhcpred_rank": min(float(mix[r["peptide"]]["%Rank_bestAllele"]) for r in class1),
            "best_target_PRIME_rank": min(float(prime[r["peptide"]]["%Rank_bestAllele"]) for r in class1),
            "class2_best_context_rank": None, "class2_best_context_peptide": "",
            "class2_best_context_allele": "", "class2_same_peptide_no_context_rank": None,
            "class2_best_no_context_rank": None, "clinical_release_ready": False,
        }
        if best_context:
            scored = c2[best_context["peptide"], best_context["context"]]
            allele = scored["BestAllele"]
            core = scored["Core_best"]
            placements = {
                start for candidate in (core, core[::-1])
                for start in range(len(best_context["peptide"]) - 8)
                if best_context["peptide"].startswith(candidate, start)
            }
            physical_overlaps = {
                overlaps_target(
                    int(best_context["start_in_region"]) + start,
                    int(best_context["start_in_region"]) + start + 9,
                    int(region["target_start_in_window"]), int(region["target_end_in_window"]),
                )
                for start in placements
            }
            record.update(
                class2_best_context_rank=float(scored["%Rank_best"]),
                class2_best_context_peptide=best_context["peptide"], class2_best_context_allele=allele,
                class2_same_peptide_no_context_rank=float(n2[best_context["peptide"]]["%Rank_" + allele]),
                class2_best_no_context_rank=float(n2[best_no_context["peptide"]]["%Rank_best"]),
                class2_best_core=scored["Core_best"],
                class2_best_core_overlaps_recorded_target=(
                    next(iter(physical_overlaps)) if len(physical_overlaps) == 1 else None
                ),
                class2_allele_qualification=(
                    "DNA_call_with_discordant_RNA_DRB3_call"
                    if allele == "DRB3_02_02" else
                    "unphased_alpha_beta_pairing_not_surface_verified"
                    if allele.startswith(("DQ", "DP")) else
                    "DNA_RNA_concordant_tool_call_not_clinically_validated"
                ),
            )
        summary_rows.append(record)
        for occurrence in (r for r in occurrences if r["region_id"] == rid and int(r["length"]) == 15):
            scored = c2[occurrence["peptide"], occurrence["context"]]
            for allele in CLASS2_ALLELES:
                core_start = int(occurrence["start_in_region"]) + int(scored["CoreP1_" + allele]) - 1
                class2_rows.append({
                    **occurrence, "allele": allele,
                    "context_rank": float(scored["%Rank_" + allele]),
                    "no_context_rank": float(n2[occurrence["peptide"]]["%Rank_" + allele]),
                    "core_overlaps_recorded_target": overlaps_target(
                        core_start, core_start + 9, int(region["target_start_in_window"]),
                        int(region["target_end_in_window"]),
                    ) if scored["SubSpec_" + allele] != "-1" else None,
                    "binding_subspecificity": scored["SubSpec_" + allele],
                    "core_coordinate_qualification": (
                        "reverse_binding_interval_not_independently_interpreted"
                        if scored["SubSpec_" + allele] == "-1" else "model_forward_CoreP1"
                    ),
                })
        for contrast in (r for r in contrasts if r["region_id"] == rid):
            mutant, other = contrast["mutant_peptide"], contrast["counterfactual_peptide"]
            if (mutant, lead["allele"]) not in all_affinity:
                continue
            for code in CLASS1_ALLELES:
                allele = "HLA-" + code[0] + "*" + code[1:3] + ":" + code[3:5]
                nm = all_affinity.get((mutant, allele))
                if nm is None:
                    continue
                wt = counterfactual_affinity[other, allele]
                contrast_rows.append({
                    **contrast, "allele": allele, "mutant_affinity_nm": nm,
                    "counterfactual_affinity_nm": wt, "predicted_affinity_ratio_WT_over_mutant": wt / nm,
                    "mutant_mix_rank": float(mix[mutant]["%Rank_" + code]),
                    "counterfactual_mix_rank": float(mix[other]["%Rank_" + code]),
                    "mutant_PRIME_rank": float(prime[mutant]["%Rank_" + code]),
                    "counterfactual_PRIME_rank": float(prime[other]["%Rank_" + code]),
                })
    table(DIRECTORY / "antigen-model-summary.tsv", summary_rows)
    table(DIRECTORY / "class2-target-scores.tsv", class2_rows)
    table(DIRECTORY / "class1-counterfactual-contrasts.tsv", contrast_rows)
    boundaries = rows(HERE / "construct/protein/boundary-predictions.tsv")
    metrics = {
        "count_unit": "Placed peptide-HLA pairs; repeated peptide occurrences count separately.",
        "MHCflurry_affinity_rank_le_0_5pct": sum(float(r["affinity_percentile"]) <= 0.5 for r in boundaries),
        "MixMHCpred_rank_le_0_5pct": sum(float(mix[r["peptide"]]["%Rank_" + r["allele"][4:].replace("*", "").replace(":", "")]) <= 0.5 for r in boundaries),
        "PRIME_rank_le_0_5pct": sum(float(prime[r["peptide"]]["%Rank_" + r["allele"][4:].replace("*", "").replace(":", "")]) <= 0.5 for r in boundaries),
        "scope": "Model-dependent rankings, not unwanted T-cell responses. PRIME and MixMHCpred are correlated; percentile scales are not interchangeable.",
    }
    save(DIRECTORY / "junction-multi-model-summary.json", metrics)
    joins = rows(HERE / "construct/protein/junctions.tsv")
    class2_junctions = []
    for row in rows(DIRECTORY / "requests/encoded-class2-occurrences.tsv"):
        start, end = int(row["cassette_start_aa"]), int(row["cassette_start_aa"]) + 15
        touched = [
            join for join in joins
            if (start < int(join["linker_end_aa"]) and end > int(join["linker_start_aa"])
                if int(join["linker_end_aa"]) > int(join["linker_start_aa"])
                else start < int(join["linker_start_aa"]) < end)
        ]
        if not touched:
            continue
        scored = c2[row["peptide"], row["context"]]
        for allele in CLASS2_ALLELES:
            cp1 = int(scored["CoreP1_" + allele]) - 1
            cs, ce = start + cp1, start + cp1 + 9
            crosses_core = any(
                (cs < int(join["linker_end_aa"]) and ce > int(join["linker_start_aa"])
                 if int(join["linker_end_aa"]) > int(join["linker_start_aa"])
                 else cs < int(join["linker_start_aa"]) < ce)
                for join in touched
            )
            class2_junctions.append({
                **row, "junction_indices": ";".join(join["junction_index"] for join in touched),
                "allele": allele, "context_rank": float(scored["%Rank_" + allele]),
                "no_context_rank": float(n2[row["peptide"]]["%Rank_" + allele]),
                "core_9aa": row["peptide"][cp1:cp1 + 9],
                "core_crosses_synthetic_junction": crosses_core,
            })
    table(DIRECTORY / "class2-junction-scores.tsv", class2_junctions)
    metrics["ClassII_synthetic_9aa_core_rank_le_0_5pct"] = sum(
        row["core_crosses_synthetic_junction"] and row["context_rank"] <= 0.5
        for row in class2_junctions
    )
    metrics["ClassII_synthetic_9aa_core_rank_le_2pct"] = sum(
        row["core_crosses_synthetic_junction"] and row["context_rank"] <= 2.0
        for row in class2_junctions
    )
    metrics["ClassII_native_core_overhang_rank_le_2pct"] = sum(
        not row["core_crosses_synthetic_junction"] and row["context_rank"] <= 2.0
        for row in class2_junctions
    )
    metrics["ClassII_15mer_junction_context_rank_le_2pct"] = sum(
        row["context_rank"] <= 2.0 for row in class2_junctions
    )
    metrics["ClassII_scope"] = (
        "15-mers across nine conditional chain combinations; synthetic junction-spanning "
        "9-aa cores are minimized in CP-SAT (0 at <=0.5%, 12 at <=2.0%), while input-region "
        "9-aa cores whose 15-mer window extends across an engineered linker/boundary "
        "(1-4 flanking linker/adjacent residues; includes mutant or wild-type input cores) "
        "are reported separately (9 at <=2.0%)."
    )
    save(DIRECTORY / "junction-multi-model-summary.json", metrics)
    save(DIRECTORY / "receipt.json", {
        "output_sha256": {name: sha(DIRECTORY / name) for name in (
            "mixmhcpred3.tsv", "prime21.tsv", "mixmhc2-context.tsv",
            "mixmhc2-no-context.tsv", "counterfactual-class1-affinity.tsv",
            "cassette-class1-presentation.tsv", "placed-class1-leads.tsv",
            "antigen-model-summary.tsv", "class2-target-scores.tsv",
            "class1-counterfactual-contrasts.tsv", "junction-multi-model-summary.json",
            "class2-junction-scores.tsv",
        )},
        "clinical_validation": False,
    })
    print("Model-specific junction metrics:", metrics)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["prepare", "wildtype-affinity", "summarize"])
    parser.add_argument("--models", type=Path)
    args = parser.parse_args()
    if args.action == "prepare":
        prepare()
    elif args.action == "wildtype-affinity":
        if args.models is None:
            parser.error("--models is required for fresh affinity inference")
        wildtype_affinity(args.models)
    else:
        summarize()
