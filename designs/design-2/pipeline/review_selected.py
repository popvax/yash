"""Render the selected-region evidence and explicit reviewer observations.

This reads completed results. It neither changes selection rules nor executes
new variant calling, RNA reconstruction or HLA inference.
"""
from yash_vaccine.paths import resolve_input
from pathlib import Path
import csv
import hashlib
import json

CODE = Path(__file__).resolve().parent
HERE = CODE.parent

# These are reviewer observations about these exact sequences, not new filters.
NOTES = {
    "KLF5": "Substantial DNA and RNA support, with 50 reference-only normal fragments. Its large supplied expression value raises the exported comparison score, but all 12 eligible mutation regions are included regardless of that rank. This is not a measured peptide-abundance ranking.",
    "SF1": "Substantial RNA allele and longer-context support. DNA and RNA alternate fractions are similar. The normal sample now supports the somatic interpretation that was unresolved without it.",
    "AQP3": "Newly selected C11Y. DNA evidence is substantial, but the RNA estimate is only 2 alternate of 4 fragments; both support the longer context. The apparent 50% RNA fraction is poorly sampled. No extra three-fragment rule was introduced.",
    "RXRA": "Good local DNA/RNA support. RNA alternate fraction is higher than DNA model fraction, which can reflect allele-specific expression or sampling; the discrepancy is not itself evidence of an invalid sequence.",
    "HSPG2": "Strong DNA evidence and adequate normal reference support. RNA is limited to 5 alternate of 9 fragments, with 4 supporting the longer context. Rare population observation is retained, rather than labelled absent.",
    "TGFBR2": "Newly selected frameshift region. Eighty-five alternate RNA fragments, including 67 for the longer context, give substantially stronger sequence support than the other newly selected region AQP3. Transcript-dependent protein numbering is not a separate genomic mutation. This is independent paired rediscovery, not a claim that the clinical frameshift was previously unknown.",
    "PPIC": "The RNA-derived window also retains the previously reviewed N190S background difference. The focal P187S is supported somatic here. The best peptide ATDGHDRSL includes the focal change and ends before that background difference. Do not describe the complete window as having exactly one difference from reference or claim that every background residue has been separately classified.",
    "KCTD19": "Only two alternate RNA fragment names out of seven; both support the longer context. The matched normal now resolves focal somatic support, but does not improve that small RNA sample.",
    "VPS11": "Seven alternate RNA fragments and four for the longer context provide modest sequence support. No exact match in the searched normal/reference sequence models was found for the scored mutant peptides; that is a scoped computational negative.",
    "ARHGAP40": "Nine alternate RNA fragments and eight for the longer context. Rare population observations remain below the policy cutoff, while direct normal reference evidence supports the focal somatic interpretation.",
    "PLEC": "Eight alternate RNA fragments support the longer context. The best target prediction is approximately 413 nM, relatively near the 500-nM gate. It qualifies computationally without establishing tumour-cell display or a T-cell response.",
    "SYT12": "Only two alternate RNA fragments out of eight; both support the longer context. Adequate normal evidence supports the somatic call but does not eliminate uncertainty from sparse RNA sampling.",
    "F5": "RNA reads support the join and extra S-to-F change that makes the scored peptides differ from the main reference. Recorded support is 5 fragments at the changed base and 2 across the full 25-aa window. The earlier targeted check reported 5 / 1,441 RNA fragments at the base and 0 alternate / 1,300 reference tumour-DNA fragments.",
    "MUC5B": "Two full-window RNA fragments support the join and extra S-to-I change that makes the scored peptides differ from the main reference. Low-complexity mapping and molecule independence remain unresolved. The recorded event/region tie-break is retained.",
    "MUC6": "RNA reads support the join and extra T-to-M change: 9 fragments at the changed base and 3 across the full window. Two full-window fragments agree with the nominated strand; one opposes it. Earlier targeted counts were 9 / 269 RNA fragments and 8 alternate / 5,583 reference tumour-DNA fragments. The selected sequence also matches alternate-reference translation ENSP00000487059.1. At the tested position (chr11:1017576), Yash's normal DNA showed {normal_ref} fragments with the usual base and {normal_alt} with the methionine-producing change.",
    "MUC3A": "Three full-window RNA fragments support the join and extra I-to-F change that makes the scored peptides differ from the main reference. All three agree with the nominated strand. The repetitive mucin context still complicates alignment attribution.",
    "HLA-B": "The first 15 residues ITQRKWEAAREAEQR occur in compatible HLA-C reference proteins. They are not the scored tumour-specific target: scoring uses peptides overlapping the nominated junction at window residues 18 (one based). Best target EQRSLPGGHVR is approximately 66 nM on A*33:03 and has no match in the searched personal-normal models. Two full-window fragments support a hypothesis in a highly homologous HLA region, not a validated novel HLA protein.",
    "PRAME": "Intentional reference-protein target. The recorded regional audit has at least two matching fragments per codon and two covering YLHARLREL, but no single merged pair spans all 25 residues. No PRAME IHC or patient tumour-cell peptide measurement was found in that review. Do not use the transcript-level 85-fragment count as peptide support.",
    "SPATA31A6": "Weakest local CTA sequence evidence: 24 of 25 codons have no qualifying reference-matching RNA fragment in the recorded regional audit; neither the full region nor the best peptide is covered. The best peptide is shared with other SPATA31A-family reference proteins. Positive supplied gene expression does not resolve this regional/isoform attribution. Kept under the existing CTA policy, not newly confirmed by normal DNA.",
    "LRRC36": "Recorded regional evidence is sparse: one matching fragment at 24 codons, none at one; one supports SYAGKIHSI and none establishes the whole 25-aa region. The source testis label is not evidence of testis-exclusive expression.",
    "YBX2": "One fragment pair matches the full selected window, two match VAPPPMVAEI, and each codon has at least two matching fragment names in the recorded regional check. This is local RNA support, not protein presentation. Normal-tissue restriction remains a separate question.",
    "TEKT5": "Two fragment pairs match the full selected window and TIDDTLQTL in the recorded regional check. The literal source normal-tissue label Actin is not an organ and cannot support a tissue-restriction claim; it is preserved as a provenance error. Normal protein/tissue validation remains unresolved.",
    "CDH1--NRG1": "Alternative predicted frameshift-tail context at the recorded breakpoint, distinct from the explicitly included in-frame driver junction. Exact-breakpoint review reports 62 junction-plus-spanning units and 473 unfused left-partner units; weighting uses 62/(62+473) times supplied CDH1 expression. Other generic baseline columns are not this denominator. Translation of this alternative context and detailed saved artifact metrics remain unresolved.",
    "CDH1::NRG1": "Explicit in-frame driver-junction region. The recorded exact-breakpoint call reports 200 junction reads and one spanning fragment. Its best non-reference prediction is approximately 508.95 nM, an explicit inclusion exception above 500 nM. The NAT motif is a potential N-linked glycosylation sequon; occupancy and its effect on processing/presentation were not modelled. This region and the alternative tail are two hypotheses from the same gene pair, not two independent drivers.",
}


def load(name):
    return json.loads((resolve_input(HERE / name)).read_text())


def main():
    selected = load("selection/selected.json")
    splice_check = load("review/splice-sequence-check.json")
    assert {r["gene"] for r in splice_check["comparisons"]} == {"F5", "MUC5B", "MUC3A", "MUC6"}
    assert all(r["reverted_primary_nonself_peptides"] == 0 for r in splice_check["comparisons"])
    normal = splice_check["MUC6_normal_base"]
    assert normal["position_1based"] == 1017576
    notes = dict(NOTES)
    notes["MUC6"] = notes["MUC6"].format(
        normal_ref=normal["base_counts"][normal["reference_base"]],
        normal_alt=normal["base_counts"][normal["methionine_associated_alternate_base"]])
    dna = {r["variant_id"]: r for r in load("exports/variant-evidence.json")}
    rna = {r["variant_id"]: r for r in load("exports/rna-rna-evidence.json")}
    proteins = {r["protein_id"]: r for r in load("exports/rna-primary-proteins.json")}
    per_hla = list(csv.DictReader((resolve_input(HERE / "exports/rank-best-peptide-per-HLA.tsv")).open(), delimiter="\t"))
    assert len(selected) == len(NOTES) == 24 and {r["gene"] for r in selected} == set(NOTES)
    rows = []
    lines = ["# Review of all 24 selected regions", "",
             "Reviewed against the completed matched-normal outputs. This is a computational evidence review, not an experimental or clinical validation. No new selection rule is introduced here.", "",
             "The ten retained mutation regions now have direct normal-reference evidence and paired native PASS calls. Seven old targets are present in white-blood-cell DNA and are removed. AQP3 and TGFBR2 are added. All twelve RNA-source sequences remain unchanged; their earlier regional reviews remain relevant.", "",
             "Counts use fragment names, not UMI-confirmed independent molecules. Mutation RNA ALT/total describes local allele support; the longer protein-context count is different. Supplied expression retains the vendor's transformation and is not asserted to be calibrated TPM.", "",
             "All 12 mutation regions meeting the evidence, expression-availability, binding and construct-compatibility criteria are included. Their expression-weighted scores order comparison tables, not a top-N inclusion list or the assembled construct. Expression magnitude and RNA ALT fraction do not remove lower-ranked qualifying mutation regions. Splice and CTA selections still use expression-weighted ranking.", "",
             "| ID | Source | Gene / change | Best affinity (nM) | Normal REF / ALT | RNA ALT / total | Longer-context fragments | Supplied gene value |", "|---|---|---|---:|---:|---:|---:|---:|"]
    for region in selected:
        gene = region["gene"]
        facts = dict(antigen_id=region["antigen_id"], region_id=region["region_id"], source=region["source"],
                     gene=gene, protein_change=region.get("protein_change", ""), sequence=region["sequence"],
                     best_affinity_nm=float(region["best_target_affinity_nm"]),
                     max_KD19=region["hydrophobicity"]["max_mean_kd"], reviewer_observation=notes[gene])
        if region["source"] == "MUTATION":
            v = dna[region["variant_id"]]; r = rna[region["variant_id"]]; p = proteins[region["protein_id"]]
            assert v["primary_DNA_eligible"] and not v["manual_review_flags"]
            n = v["normal_evidence"]
            assert n["classification"] == "NORMAL_REFERENCE_SUPPORTED"
            best = min((x for x in per_hla if x["variant_id"] == v["variant_id"]), key=lambda x: float(x["affinity_nm"]))
            facts.update(normal_REF=n["read_evidence"]["fragment_counts"].get("REF", 0), normal_ALT=n["read_evidence"]["fragment_counts"].get("ALT", 0),
                         tumour_DNA_DP=v["DNA_DP"], tumour_DNA_ALT=v["DNA_AD"][1], tumour_model_AF=v["DNA_model_AF"], GERMQ=v["GERMQ"],
                         RNA_ALT=r["rna_alt_fragments"], RNA_total=r["rna_total_fragments"], protein_context_fragments=p["supporting_fragments"],
                         supplied_gene_value=region["supplied_gene_value"], best_peptide=best["peptide"], best_allele=best["allele"],
                         normal_status=n["classification"], population_status=v["population_status"], population_global_AF=v["QC_PASS_combined_global_AF"])
        else:
            facts.update(supplied_gene_value=(region.get("supplied_transformed_gene_tpm", region.get("expression_proxy_left", "")) if region["source"] != "FUSION" else ""),
                         expression_proxy=region.get("expression_proxy_left", ""),
                         RNA_full_window_fragments=region.get("fragment_count", ""), normal_detection=region.get("normal_detection", ""),
                         strand_agreeing=region.get("strand_agreeing_fragments", ""), strand_opposing=region.get("strand_opposing_fragments", ""),
                         motif_agrees_with_gene=region.get("motif_agrees_with_gene", ""))
        rows.append(facts)
        counts = (f"{facts['normal_REF']} / {facts['normal_ALT']} | {facts['RNA_ALT']} / {facts['RNA_total']} | {facts['protein_context_fragments']}"
                  if facts['source'] == 'MUTATION' else "— | — | —")
        expression = f"{float(facts['supplied_gene_value']):.3f}" if facts['supplied_gene_value'] not in [None, ""] else "—"
        lines.append(f"| {facts['antigen_id']} | {facts['source']} | {gene} {facts['protein_change']} | {facts['best_affinity_nm']:.2f} | {counts} | {expression} |")
    lines += ["", "## Individual assessments", ""]
    for facts in rows:
        lines += [f"### {facts['antigen_id']} · {facts['gene']} {facts['protein_change']}".rstrip(), "",
                  f"`{facts['sequence']}` · `{facts['region_id']}` · max KD19 {facts['max_KD19']:.3f}.", "", facts["reviewer_observation"], ""]
        if facts['source'] == 'SPLICE':
            lines += [f"Recorded complete-window fragments: {facts['RNA_full_window_fragments']}; strand agreement/opposition {facts['strand_agreeing']}/{facts['strand_opposing']}. The splice motif disagrees with the nominated gene direction. The selected window passed the recorded read-orientation rule and annotated upstream-frame reconstruction. ND means not detected at these coordinates in the two searched normal catalogs, not absent from every healthy tissue. New normal DNA does not settle RNA direction, splicing or translation.", ""]
    lines += ["## Construct-level findings", "",
              "The 724-aa protein contains each selected region exactly once, 43 linker residues, a 26-aa SEC element and a 55-aa MITD. Fresh inference covered all 6,600 order/linker options (578,112 peptide/HLA queries), and a separate final-junction run covered 6,132 queries. All four recorded optimizer objectives were proven optimal for this finite search and its unwanted-junction binding objective. This does not optimize cleavage or intended antigen presentation.", "",
              "The final junction scan still has one row at affinity percentile ≤0.5%, 26 at ≤2%, and 24 at affinity ≤500 nM; none is ≤50 nM. Passing verification means agreement with the objective and independent scan, not absence of predicted junction binders. There are no high-hydropathy windows confined to antigens/linkers; SEC/MITD overlaps remain separately reported and are not topology-validated.", "",
              "The fixed protein was synonymously encoded with PopRNA (human usage, 0.10 cutoff, seed 0). Translation, forbidden motifs, codon cutoff and all 64 GenBank feature annotations pass their hard checks. Coding GC is 56.08%. The final folding heuristic is 21 consecutive paired brackets versus a target of 14, so that heuristic target is **not retained** after the sequential repairs. Repeated 12-nt pairs are zero. These structure/codon heuristics are not expression or stability measurements. PolyA and the downstream extension remain unspecified.", "",
              "## Evidence", "",
              "- [All 17 original mutation reassessments](final/current-17-final-audit.tsv), including seven not retained.",
              "- [Selected regions](../selection/selected.tsv), [DNA evidence](../inputs/mutations/variant-evidence.tsv), [RNA evidence](../inputs/rna/rna-evidence.tsv), and [per-allele best mutation peptides](../inputs/predictions/mutations/best-peptide-per-HLA.tsv).",
              "- [Splice sequence comparisons and focal MUC6 normal-DNA evidence](splice-sequence-check.json), reproduced with [the check script](../pipeline/check_splice_sequences.py).",
              "- [Earlier regional observations for unchanged RNA sources](../../initial/antigen-review.md): retained as historical measurements, without importing its old tumour-only somatic assignments.",
              "- [Fresh final junction verification](../construct/final-check.json) and [nucleotide/GenBank verification](../construct/mrna/verification.json).", ""]
    out = HERE / "review"; out.mkdir(exist_ok=True)
    (out / "antigen-review.md").write_text("\n".join(lines))
    (out / "antigen-review.json").write_text(json.dumps(rows, indent=2) + "\n")
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with (out / "antigen-review.tsv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    hashes = {str(p.relative_to(HERE)): hashlib.sha256(resolve_input(p).read_bytes()).hexdigest() for p in
              [HERE / "selection/selected.json", HERE / "exports/variant-evidence.json", HERE / "exports/rna-rna-evidence.json",
               HERE / "exports/rna-primary-proteins.json", HERE / "exports/rank-best-peptide-per-HLA.tsv",
               HERE / "review/splice-sequence-check.json", Path(__file__)]}
    (out / "antigen-review-inputs.json").write_text(json.dumps(hashes, indent=2) + "\n")
    print("Reviewed", len(rows), "selected regions; selection unchanged")


if __name__ == "__main__":
    main()
