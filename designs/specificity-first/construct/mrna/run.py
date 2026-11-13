"""Encode the fixed routed protein with PopRNA and annotate known RNA/DNA segments.

The supplied heuristic uses a pinned human codon table, synonymous motif/repeat
repair, alternative-frame stops and ViennaRNA MFE structure heuristics. PolyA and
any downstream extension remain unspecified. Reproduction starts no cloud job.
"""
import argparse
import csv
import hashlib
import json
import math
import re
import subprocess
import sys
import tempfile
from importlib.metadata import version
from itertools import groupby
from pathlib import Path

import RNA as vienna
import python_codon_tables as pct
from Bio import SeqIO
from Bio.Seq import Seq
from genbank import build_record, load_inputs, verify_export

HERE = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sequence_hash(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def read_rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def hits(sequence, pattern):
    return [i for i in range(len(sequence) - len(pattern) + 1)
            if sequence.startswith(pattern, i)]


def fasta(path, name, sequence):
    path.write_text(">" + name + "\n" + "".join(
        sequence[i:i + 80] + "\n" for i in range(0, len(sequence), 80)))


def metrics(coding_dna, protein, table, settings):
    weights = [table[aa][coding_dna[3*i:3*i+3]] / max(table[aa].values())
               for i, aa in enumerate(protein)]
    gc = lambda text: sum(base in "GC" for base in text) / len(text)
    windows = [gc(coding_dna[i:i+100]) for i in range(len(coding_dna) - 99)]
    longest = max((len(list(group)), base) for base, group in groupby(coding_dna))
    rna = coding_dna.replace("T", "U")
    structure, energy = vienna.fold(rna)
    words = {}
    for i in range(len(rna)-11):
        word = rna[i:i+12]
        words[word] = words.get(word, 0)+1
    frames = {}
    for offset in (1, 2):
        triplets = [rna[i:i+3] for i in range(offset, len(rna)-2, 3)]
        runs, run = [], 0
        for codon in triplets:
            if codon in ("UAA", "UAG", "UGA"):
                runs.append(run)
                run = 0
            else:
                run += 1
        runs.append(run)
        frames[str(offset)] = dict(stop_codons=sum(c in ("UAA", "UAG", "UGA") for c in triplets),
                                   longest_stop_free_run_codons=max(runs))
    return {
        "cai_relative_to_frozen_human_table": round(math.exp(sum(map(math.log, weights)) / len(weights)), 8),
        "codons_with_maximum_table_frequency": sum(abs(w - 1.0) < 1e-12 for w in weights),
        "total_sense_codons": len(weights),
        "coding_gc_percent": round(100 * gc(coding_dna), 4),
        "coding_100nt_window_gc_percent_min": round(100 * min(windows), 4),
        "coding_100nt_window_gc_percent_max": round(100 * max(windows), 4),
        "coding_longest_homopolymer": {"base": longest[1], "length": longest[0]},
        "rare_codons_below_cutoff": sum(table[aa][coding_dna[3*i:3*i+3]] < settings["codon_cutoff"]
                                        for i, aa in enumerate(protein)),
        "coding_forbidden_motif_occurrences": {motif: len(hits(rna, motif)) for motif in settings["avoid_sequences_rna"]},
        "repeated_12nt_pairs": sum(n*(n-1)//2 for n in words.values()),
        "coding_mfe_kcal_per_mol": round(float(energy), 2),
        "coding_longest_paired_bracket_run": max((len(m.group()) for m in re.finditer("[()]+", structure)), default=0),
        "alternative_reading_frames": frames,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE,
                        help="Output directory; inputs always come from this design directory.")
    args = parser.parse_args()
    output = args.output.resolve()
    cfg = json.loads((HERE / "backbone.json").read_text())
    annotation_sources, selected, architecture, annotation_paths = load_inputs(HERE)
    table_path = HERE / "human-codon-usage.json"
    table_record = json.loads(table_path.read_text())
    table = table_record["codon_usage"]
    require(table == pct.get_codons_table(cfg["optimization"]["species_table"]),
            "Installed human codon table differs from the frozen table.")
    vienna.cvar.temperature = cfg["optimization"]["folding_temperature_celsius"]
    require(file_hash(table_path) == cfg["optimization"]["codon_table_file_sha256"],
            "Frozen codon table checksum mismatch.")
    packages = {}
    for line in (HERE / "requirements.txt").read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        name, expected = line.split("==")
        packages[name] = version(name)
        require(packages[name] == expected, f"Install pinned {name}=={expected}.")

    parent = (HERE / cfg["protein_design"]).resolve()
    protein_path, parts_path = parent / "cassette.fasta", parent / "components.tsv"
    require(file_hash(protein_path) == cfg["protein_file_sha256"], "Protein input checksum mismatch.")
    require(file_hash(parts_path) == cfg["protein_components_sha256"], "Protein component checksum mismatch.")
    protein = str(SeqIO.read(protein_path, "fasta").seq)
    parent_parts = read_rows(parts_path)
    cursor = 0
    for part in parent_parts:
        start, end = int(part["start_aa"]), int(part["end_aa"])
        require(start == cursor and protein[start:end] == part["aa"], "Protein component mismatch.")
        cursor = end
    antigen_parts = [p for p in parent_parts if p["kind"] == "antigen"]
    require(cursor == len(protein) == cfg["protein_length_aa"] and protein.startswith("M"), "Expected the pinned routed protein.")
    require(len(antigen_parts) == len({p["region_id"] for p in antigen_parts}) == len(selected),
            "Expected every distinct selected antigen region.")

    check_path = (HERE / cfg["junction_verification"]).resolve()
    check = json.loads(check_path.read_text())
    require(check["status"] == "FRESH_FULL_JUNCTION_CHECK_PASSED" and
            check["protein_sha256"] == file_hash(protein_path) and
            check["components_sha256"] == file_hash(parts_path),
            "Fresh junction verification does not match the protein input.")

    fixed = cfg["fixed_elements"]
    context = cfg["template_context"]
    utr5, utr3, stops = fixed["utr5_dna"], fixed["utr3_body_dna"], fixed["stop_block_dna"]
    extension = fixed["post_poly_a_extension"]
    poly_a = fixed["poly_a"]
    require(poly_a["status"] == "unspecified" and poly_a["length_nt"] is None
            and poly_a["sequence_dna"] is None, "This design leaves PolyA unspecified.")
    require(extension["status"] == "unspecified" and extension["length_nt"] is None
            and extension["sequence_dna"] is None, "The post-PolyA extension is unspecified.")
    before = context["prefix_before_transcribed_sequence_dna"]
    promoter = context["promoter_dna"]
    require(stops == "TAATAGTGA", "Unexpected stop block.")
    require((len(utr5), len(utr3)) == (63, 71), "Unexpected UTR dimensions.")
    require(before.endswith(promoter), "T7 promoter context mismatch.")
    require(all(set(x) <= set("ACGT") for x in [utr5, utr3, stops, before]),
            "Non-DNA character in known elements.")

    # Optimise the known upstream segment only. No artificial junction is made
    # across the unspecified PolyA, and its boundaries cannot be assessed.
    coding_start = len(before) + len(utr5)
    coding_end = coding_start + 3 * len(protein)
    settings = cfg["optimization"]
    with tempfile.TemporaryDirectory(prefix="yash-poprna-") as temporary:
        candidate = Path(temporary) / "coding.txt"
        command = [sys.executable, str(HERE / "optimizer.py"), "--protein", "-i", protein,
                   "--output", str(candidate), "--organism", settings["species_table"],
                   "--cutoff", str(settings["codon_cutoff"]), "--seed", str(settings["seed"]),
                   "--fold-iterations", str(settings["structure_iterations_per_pass"]),
                   "--sequences", *settings["avoid_sequences_rna"]]
        result = subprocess.run(command, text=True, capture_output=True)
        require(result.returncode == 0, "PopRNA failed: " + result.stderr + result.stdout)
        coding = candidate.read_text().strip()
        optimizer_log = result.stdout
    require(set(coding) <= set("ACGT") and len(coding) == 3*len(protein), "Invalid PopRNA coding sequence.")
    # Design-specific pin: null until the first deterministic run of THIS
    # protein, then frozen by pin_coding_sequence.py and re-checked on reruns.
    if settings.get("recorded_coding_sequence_sha256"):
        require(sequence_hash(coding) == settings["recorded_coding_sequence_sha256"],
                "PopRNA output differs from the recorded design.")
    fixed_suffix = stops + utr3
    template_before = before + utr5 + coding + fixed_suffix
    rna_before = (utr5 + coding + fixed_suffix).replace("T", "U")
    cds = coding + stops[:3]
    require(str(Seq(coding).translate()) == protein and str(Seq(cds).translate(to_stop=True)) == protein,
            "Protein translation changed.")
    require(template_before[:coding_start] == before + utr5
            and template_before[coding_end:] == fixed_suffix, "Known fixed context changed.")
    require(template_before[len(before):].replace("T", "U") == rna_before, "Upstream RNA/template mismatch.")
    require(not hits(template_before, "TCTAGA"), "Unexpected upstream XbaI site.")

    # Coordinates cover the known upstream segment only. No downstream
    # sequence or restriction boundary is inferred from the example backbone.
    components = []
    def add(kind, name, sequence, rna_start, aa="", region_id=""):
        components.append(dict(
            kind=kind, name=name, rna_start_0based=rna_start,
            rna_end_exclusive=rna_start+len(sequence),
            template_start_0based=len(before)+rna_start,
            template_end_exclusive=len(before)+rna_start+len(sequence),
            length_nt=len(sequence), dna=sequence, aa=aa, region_id=region_id))
    add("5UTR_Kozak", "Foundry_5UTR_Kozak", utr5, 0)
    for part in parent_parts:
        a, b = int(part["start_aa"]), int(part["end_aa"])
        nt = coding[3*a:3*b]
        require(str(Seq(nt).translate()) == part["aa"], "Translated component mismatch.")
        add(part["kind"], part["name"], nt, len(utr5)+3*a, part["aa"], part["region_id"])
    offset = len(utr5) + len(coding)
    for kind, name, sequence in [
        ("stop_codon", "primary_stop_TAA", stops[:3]),
        ("additional_stops", "TAG_TGA", stops[3:]),
        ("3UTR", "Foundry_3UTR", utr3),
        ("polyA_unspecified", "PolyA", ""),
    ]:
        add(kind, name, sequence, offset)
        offset += len(sequence)
    require(offset == len(rna_before), "Known component partition mismatch.")
    require("".join(p["dna"] for p in components) == rna_before.replace("U", "T"),
            "Known component sequences mismatch.")

    # Public coordinates are relative to the NAMED segment. PolyA has no
    # assigned coordinates, sequence or length; it is not a zero-length tail.
    public_components = []
    for part in components:
        row = dict(part)
        unknown = part["kind"] == "polyA_unspecified"
        row["segment"] = "unspecified" if unknown else "before_PolyA"
        if unknown:
            for key in ("rna_start_0based", "rna_end_exclusive", "template_start_0based",
                        "template_end_exclusive", "length_nt"):
                row[key] = None
        public_components.append(row)

    output.mkdir(parents=True, exist_ok=True)
    fasta(output / "cds.fasta", "Yash_SEC_antigens_MITD_CDS_DNA_including_first_stop", cds)
    def segments(filename, alphabet, sequences):
        text = ""
        for name, sequence in sequences:
            text += f">Yash_{alphabet}_{name}_partial_PolyA_unspecified\n"
            text += "".join(sequence[i:i+80]+"\n" for i in range(0, len(sequence), 80))
        (output / filename).write_text(text)
    segments("mrna-segments.fasta", "RNA", [("before_PolyA", rna_before)])
    segments("template-segments.fasta", "DNA", [("before_PolyA", template_before)])
    (output / "rna-layout.txt").write_text(
        f"5′ UTR → SEC → {len(selected)} antigens + linkers → MITD → stops → 3′ UTR → PolyA\n\n"
        "PolyA: AAAAA…AAAAA\n"
        "The A-count is illustrative. PolyA length, sequence and format are unspecified.\n"
        "Any extension after PolyA depends on its architecture and is unspecified.\n"
        "mrna-segments.fasta and template-segments.fasta contain the known upstream segment only.\n"
        "These are partial sequences, not complete molecules.\n")
    with (output / "components.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(public_components[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(public_components)
    records = build_record(template_before, protein, components, cfg,
                           annotation_sources, selected, architecture)
    SeqIO.write(records, output / "template-segments.gb", "genbank")
    annotation_check = verify_export(output / "template-segments.gb", records, components,
                                     protein, (rna_before,), selected)

    after_metrics = metrics(coding, protein, table, settings)
    require(after_metrics["rare_codons_below_cutoff"] == 0, "Excluded rare codon remains.")
    require(not any(after_metrics["coding_forbidden_motif_occurrences"].values()), "Forbidden coding motif remains.")
    (output / "optimizer.log").write_text(optimizer_log)
    input_hashes = {name: file_hash(HERE / name) for name in
                    ["backbone.json", "human-codon-usage.json", "requirements.txt", "run.py", "optimizer.py",
                     "genbank.py", "annotation-sources.json", *annotation_paths]}
    input_hashes.update({cfg["protein_design"] + "/cassette.fasta": file_hash(protein_path),
                         cfg["protein_design"] + "/components.tsv": file_hash(parts_path),
                         cfg["junction_verification"]: file_hash(check_path)})
    lengths = {"protein_aa": len(protein), "sense_coding_nt": len(coding),
               "cds_including_first_stop_nt": len(cds), "coding_plus_three_stops_nt": len(coding)+len(stops),
               "transcribed_rna_nt": None, "supplied_template_context_bp": None,
               "known_rna_segments_nt": [len(rna_before)],
               "known_template_segments_bp": [len(template_before)]}
    manifest = {
        "name": cfg["name"], "scope": "Known upstream RNA/DNA segments; PolyA and any downstream extension are unspecified.",
        "complete_molecular_sequence": False, "polyA": poly_a, "lengths": lengths,
        "post_poly_a_extension": extension,
        "segment_order": ["before_PolyA", "PolyA (unspecified)"],
        "coordinate_convention": "zero-based, end-exclusive within the named RNA/template segment; PolyA coordinates and length are blank",
        "optimization": cfg["optimization"], "packages": packages,
        "optimization_metrics": after_metrics,
        "heuristic_targets": {
            "paired_bracket_run_target_nt": settings["paired_bracket_run_target_nt"],
            "final_paired_bracket_run_nt": after_metrics["coding_longest_paired_bracket_run"],
            "folding_target_retained": after_metrics["coding_longest_paired_bracket_run"] <= settings["paired_bracket_run_target_nt"],
            "repeat_target_met": after_metrics["repeated_12nt_pairs"] == 0,
            "note": "Later synonymous-codon and repeat stages can change the folding heuristic. These are not guarantees of expression or stability."
        },
        "fixed_elements": fixed, "source": cfg["source"], "input_sha256": input_hashes,
        "sequence_sha256": {"rna_before_PolyA": sequence_hash(rna_before),
                            "template_before_PolyA": sequence_hash(template_before),
                            "protein": sequence_hash(protein)},
        "antigen_regions": len(antigen_parts), "protein_sequence_changed": False,
        "cap_chemistry": None, "nucleotide_modifications": None,
        "raw_examples_required_for_rerun": False, "genbank_annotation": annotation_check,
        "assessment_scope": {"new_HLA_predictions": False, "RNA_structure_prediction": True,
                             "expression_prediction": False, "PolyA_and_its_boundaries_assessed": False,
                             "experimental_validation": "not supplied", "template_topology": context["topology"]},
    }
    (output / "construct.json").write_text(json.dumps(manifest, indent=2) + "\n")
    files = ["cds.fasta", "mrna-segments.fasta", "template-segments.fasta", "components.tsv",
             "template-segments.gb", "rna-layout.txt", "construct.json", "optimizer.log"]
    verification = {
        "status": "KNOWN_NUCLEOTIDE_SEGMENTS_VERIFIED_POLYA_UNSPECIFIED",
        "scope": "Known segments and coding constraints only; not a complete RNA/template sequence.",
        "complete_construct_verified": False, "polyA_status": "unspecified",
        "exact_parent_protein_translation": True, "all_antigen_regions_preserved": True,
        "all_internal_linker_amino_acids_preserved": True, "SEC_and_MITD_preserved": True,
        "known_UTRs_and_upstream_context_unchanged": True, "three_stops_TAA_TAG_TGA_preserved": True,
        "post_poly_a_extension_status": "unspecified", "PolyA_and_its_boundaries_assessed": False,
        "internal_XbaI_sites": 0, "whole_template_XbaI_sites": None,
        "all_constraints_pass": True,
        "coding_motif_constraints_pass": True, "rare_codon_cutoff_pass": True,
        "folding_heuristic_target_retained": after_metrics["coding_longest_paired_bracket_run"] <= settings["paired_bracket_run_target_nt"],
        "final_paired_bracket_run_nt": after_metrics["coding_longest_paired_bracket_run"],
        "repeated_12nt_pairs": after_metrics["repeated_12nt_pairs"],
        "lengths": lengths, "genbank_annotation": annotation_check,
        "input_sha256": input_hashes, "output_sha256": {name: file_hash(output / name) for name in files},
    }
    (output / "verification.json").write_text(json.dumps(verification, indent=2) + "\n")
    print(json.dumps({"status": verification["status"], "lengths": lengths,
                      "optimized_coding_metrics": after_metrics, "output": str(output)}, indent=2))


if __name__ == "__main__":
    main()
