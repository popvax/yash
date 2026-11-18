"""Build the pre-solved 15-region KCTD19-excised contingency construct (strand-strict-15/).

If strand-specific RT-ddPCR / targeted RNA-seq under Question SG-04 does not detect
coding-strand KCTD19 p.G920R expression (currently 0 coding-strand ALT / 2 opposite-strand
ALT fragments in RNA_TN25-112103.bam), this 10-phase CP-SAT-optimal 15-region construct
(14 coding-strand-supported somatic mutations + in-frame CDH1::NRG1 driver fusion) provides
both the 481-aa protein and its verified known-segment PopRNA/GenBank nucleotide encoding.
"""
import hashlib
import shutil
import subprocess
import sys

from Bio import SeqIO

from curate import HERE, write
from dataio import iter_rows, load, resolve_path, rows, save, sha, table
from protein import assemble, assess_assembled_class2, emit_assembly, solve
from replay_helpers import objective_values, read_edges, read_triplets
from safety_audit import longest_duplex_metrics
from construct_identifiers import configure_identifiers

OUT = HERE / "strand-strict-15"


def cached_junction_check_record(target_root, regions, path, triplets, alleles, values):
    """Derive a strictly cached assembly/enumeration check without claiming fresh model inference."""
    prot_dir = target_root / "construct/protein"
    prot_path = prot_dir / "cassette.fasta"
    comp_path = prot_dir / "components.tsv"
    hydro = load(prot_dir / "hydrophobicity.json")
    p_construct = load(prot_dir / "construct.json")
    prot_seq = str(SeqIO.read(prot_path, "fasta").seq)
    parts = rows(comp_path)
    actual_antigens = [(r["region_id"], r["aa"]) for r in parts if r["kind"] == "antigen"]
    expected_antigens = [(r["region_id"], r["sequence"]) for r in regions]
    if len(actual_antigens) != len(expected_antigens) or set(actual_antigens) != set(expected_antigens):
        raise ValueError("Assembled components do not contain each selected antigen once")
    joins = []
    for row in rows(prot_dir / "junctions.tsv"):
        j = dict(row)
        for k in ("junction_index", "linker_start_aa", "linker_end_aa"):
            j[k] = int(j[k])
        joins.append(j)
    from protein import assess
    recomputed_preds, recomputed_totals = assess(prot_seq, joins, path, alleles, values)
    recorded_bounds = rows(prot_dir / "boundary-predictions.tsv")
    if len(recomputed_preds) != len(recorded_bounds) or recomputed_totals != p_construct["totals"]:
        raise ValueError("Cached full-chain junction enumeration disagrees with recorded boundary predictions")
    chosen_ids = {edge["edge_id"] for edge in path}
    triplet_sum = sum(edge["intended_proc_score_x10000"] for edge in path) + sum(
        t["processing_score_x10000"] for t in triplets
        if t["incoming_edge_id"] in chosen_ids and t["outgoing_edge_id"] in chosen_ids
    )
    solver_proc = next(
        ph["value"] for ph in p_construct["solver"]["phases"]
        if ph["objective"] == "intended_boundary_processing_x10000"
    )
    if triplet_sum != solver_proc:
        raise ValueError("Selected triplet processing sum disagrees with recorded solver objective")
    return {
        "all_antigens_exactly_once": True,
        "antigen_regions": len(regions),
        "cached_boundary_enumeration_verified": True,
        "clinical_release_ready": False,
        "components_sha256": sha(comp_path),
        "fresh_repeat_model_inference_performed": False,
        "hydrophobicity": hydro,
        "max_relative_affinity_difference": None,
        "processing_objective_matches_full_cassette": None,
        "protein_length_aa": len(prot_seq),
        "protein_sha256": sha(prot_path),
        "scope": "Cached lookup-table junction enumeration and sequence checks only; fresh repeat MHCflurry inference was not executed in this offline assembly step.",
        "status": "CACHED_ASSEMBLY_JUNCTION_CHECK_PASSED",
        "threshold_classifications_identical": None,
        "totals": recomputed_totals,
        "triplet_processing_sum_matches_solver_objective": True,
    }


def encode_strict_mrna(target_root, regions, path, triplets, alleles, values, recorded_coding_sha256=None):
    prot_path = target_root / "construct/protein/cassette.fasta"
    comp_path = target_root / "construct/protein/components.tsv"
    prot_seq = str(SeqIO.read(prot_path, "fasta").seq)

    shutil.copyfile(HERE / "construct/architecture.json", target_root / "construct/architecture.json")
    existing_fc = target_root / "construct/final-check.json"
    existing_si = resolve_path(target_root / "construct/selected-inference/predictions.tsv")
    if (existing_fc.is_file() and existing_si.is_file()
            and load(existing_fc).get("status") == "FRESH_FULL_JUNCTION_CHECK_PASSED"
            and load(existing_fc).get("protein_sha256") == sha(prot_path)
            and load(existing_fc).get("components_sha256") == sha(comp_path)):
        pass
    else:
        save(
            existing_fc,
            cached_junction_check_record(target_root, regions, path, triplets, alleles, values),
        )

    mrna_dir = target_root / "construct/mrna"
    mrna_dir.mkdir(parents=True, exist_ok=True)
    for name in (
        "optimizer.py", "genbank.py", "accessibility.py", "structure_rescue.py",
        "human-codon-usage.json", "requirements.txt", "run.py",
    ):
        shutil.copyfile(HERE / "construct/mrna" / name, mrna_dir / name)

    ann = load(HERE / "construct/mrna/annotation-sources.json")
    ann["selected_regions"]["sha256"] = sha(target_root / "selection/selected.tsv")
    ann["record_date"] = "13-NOV-2026"
    save(mrna_dir / "annotation-sources.json", ann)

    backbone = load(HERE / "construct/mrna/backbone.json")
    configure_identifiers(backbone, "STRAND15")
    backbone["protein_length_aa"] = len(prot_seq)
    backbone["protein_file_sha256"] = sha(prot_path)
    backbone["protein_components_sha256"] = sha(comp_path)
    if recorded_coding_sha256 is None:
        existing_bb = HERE / "strand-strict-15/construct/mrna/backbone.json"
        if existing_bb.is_file():
            recorded_coding_sha256 = load(existing_bb)["optimization"]["recorded_coding_sequence_sha256"]
        else:
            settings = backbone["optimization"]
            fixed = backbone["fixed_elements"]
            tmp_out = mrna_dir / ".tmp-coding.txt"
            subprocess.run([
                sys.executable, str(mrna_dir / "optimizer.py"), "--protein", "-i", prot_seq,
                "--output", str(tmp_out), "--organism", settings["species_table"],
                "--cutoff", str(settings["codon_cutoff"]), "--seed", str(settings["seed"]),
                "--fold-iterations", str(settings["structure_iterations_per_pass"]),
                "--utr5-rna", fixed["utr5_dna"].replace("T", "U"),
                "--fixed-suffix-rna", (fixed["stop_block_dna"] + fixed["utr3_body_dna"]).replace("T", "U"),
                "--accessibility-trials", str(settings["accessibility_trials"]),
                "--sequences", *settings["avoid_sequences_rna"],
            ], check=True, capture_output=True, text=True)
            recorded_coding_sha256 = hashlib.sha256(tmp_out.read_text().strip().encode("ascii")).hexdigest()
            tmp_out.unlink(missing_ok=True)
    backbone["optimization"]["recorded_coding_sequence_sha256"] = recorded_coding_sha256
    save(mrna_dir / "backbone.json", backbone)

    subprocess.run(
        [sys.executable, str(mrna_dir / "run.py"), "--output", str(mrna_dir)],
        check=True, capture_output=True, text=True,
    )
    rna_seg = str(SeqIO.read(mrna_dir / "mrna-segments.fasta", "fasta").seq)
    return longest_duplex_metrics(rna_seg)


def build():
    write(OUT, panel="strand-strict-core")
    selected16 = load(HERE / "selection/selected.json")
    selected15 = load(OUT / "selection/selected.json")
    kctd_idx = next(i for i, r in enumerate(selected16) if r["gene"] == "KCTD19")

    def remap(idx):
        if idx == kctd_idx:
            return None
        if idx == len(selected16):
            return len(selected15)
        return idx if idx < kctd_idx else idx - 1

    edges16 = read_edges(HERE / "construct/protein/edge-costs.tsv")
    edges15 = []
    for edge in edges16:
        li, ri = remap(edge["left_index"]), remap(edge["right_index"])
        if li is not None and ri is not None:
            edges15.append(dict(edge, left_index=li, right_index=ri))
    kept_ids = {edge["edge_id"] for edge in edges15}
    triplets16 = read_triplets(HERE / "construct/protein/processing-triplets.tsv")
    triplets15 = []
    for trip in triplets16:
        ri = remap(trip["region_index"])
        if ri is not None and trip["incoming_edge_id"] in kept_ids and trip["outgoing_edge_id"] in kept_ids:
            triplets15.append(dict(trip, region_index=ri))

    path, solver_info = solve(len(selected15), edges15, triplets15, seconds=180)
    expected = {ph["objective"]: ph["value"] for ph in solver_info["phases"]}
    if objective_values(edges15, path, triplets15) != expected:
        raise ValueError("15-region contingency path does not match optimal objective values")

    recorded16 = load(HERE / "construct/protein/construct.json")
    alleles = load(HERE / "exports/source-case.json")["hla_class_I"]
    values = {
        (row["peptide"], row["allele"]): (float(row["affinity_nm"]), float(row["affinity_percentile"]))
        for row in iter_rows(HERE / "construct/all-junction-affinity.tsv.gz")
    }
    emit_assembly(OUT, selected15, path, solver_info, recorded16["SEC"], recorded16["MITD"], alleles, values)
    seq15, _, joins15 = assemble(path, selected15, recorded16["SEC"], recorded16["MITD"])
    imm15 = OUT / "immunogenicity"
    imm15.mkdir(parents=True, exist_ok=True)
    c2_ctx15 = {(r["Peptide"], r["Context"]): r for r in rows(imm15 / "mixmhc2-context.tsv")}
    n2_15 = {r["Peptide"]: r for r in rows(imm15 / "mixmhc2-no-context.tsv")}
    c2_junctions15, c2_summary15 = assess_assembled_class2(seq15, joins15, c2_ctx15, n2_15)
    table(imm15 / "class2-junction-scores.tsv", c2_junctions15)
    bounds15 = rows(OUT / "construct/protein/boundary-predictions.tsv")
    mix = {r["Peptide"]: r for r in rows(HERE / "immunogenicity/mixmhcpred3.tsv")}
    prime = {r["Peptide"]: r for r in rows(HERE / "immunogenicity/prime21.tsv")}
    junc_summary15 = {
        "count_unit": "Placed peptide-HLA pairs; repeated peptide occurrences count separately.",
        "MHCflurry_affinity_rank_le_0_5pct": sum(float(r["affinity_percentile"]) <= 0.5 for r in bounds15),
        "MixMHCpred_rank_le_0_5pct": sum(
            float(mix[r["peptide"]]["%Rank_" + r["allele"][4:].replace("*", "").replace(":", "")]) <= 0.5
            for r in bounds15
        ),
        "PRIME_rank_le_0_5pct": sum(
            float(prime[r["peptide"]]["%Rank_" + r["allele"][4:].replace("*", "").replace(":", "")]) <= 0.5
            for r in bounds15
        ),
        **c2_summary15,
        "full_chain_assembled_15mers_verified": len(seq15) - 14,
        "substituted_flank_contexts_used": 0,
        "scope": "Exact full-chain Class-I and Class-II evaluation of the 481-aa 15-region fallback; '-' only at actual termini, zero substituted flanking contexts.",
    }
    save(imm15 / "junction-multi-model-summary.json", junc_summary15)
    duplex15 = encode_strict_mrna(OUT, selected15, path, triplets15, alleles, values)
    mrna15 = OUT / "construct/mrna"
    ver15 = load(mrna15 / "verification.json")
    imm15_files = [
        "mixmhc2-context.tsv", "mixmhc2-no-context.tsv",
        "class2-junction-scores.tsv", "junction-multi-model-summary.json",
    ]
    for extra_imm in ("cassette-class1-presentation.tsv", "placed-class1-leads.tsv"):
        if resolve_path(imm15 / extra_imm).is_file():
            imm15_files.append(extra_imm)
    save(OUT / "contingency-receipt.json", {
        "construct_variant": ver15["construct_variant"],
        "export_identifiers": ver15["export_identifiers"],
        "purpose": "Pre-solved KCTD19-excised 15-region contingency construct (protein + known-segment mRNA) requiring >=1 coding-strand mutant RNA fragment",
        "excised_gene": "KCTD19",
        "excised_reason": "0 coding-strand ('-') and 2 opposite-strand ('+') mutant RNA fragments in RNA_TN25-112103.bam (somatic DNA-supported, clonality unestablished, coding-RNA status unconfirmed)",
        "selected_regions": len(selected15),
        "protein_length_aa": len(seq15),
        "sense_coding_nt": ver15["lengths"]["sense_coding_nt"],
        "known_rna_segment_nt": ver15["lengths"]["known_rna_segments_nt"][0],
        "known_template_segment_bp": ver15["lengths"]["known_template_segments_bp"][0],
        "final_paired_bracket_run_nt": ver15["final_paired_bracket_run_nt"],
        "known_rna_mfe_duplex_metrics": duplex15,
        "chosen_edge_ids": [edge["edge_id"] for edge in path],
        "solver_objectives": expected,
        "assembled_junction_multi_model_summary": junc_summary15,
        "upstream_junction_check_status": ver15["upstream_junction_check_status"],
        "immunogenicity_sha256": {
            name: sha(imm15 / name) for name in imm15_files
        },
        "mrna_verification_sha256": sha(mrna15 / "verification.json"),
        "clinical_release_ready": False,
    })
    return expected


if __name__ == "__main__":
    print(build())
