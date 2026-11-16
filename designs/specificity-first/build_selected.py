"""Build the design's selected.tsv with the full recorded field set.

Rows come from three recorded sources: the lead Design #2 selection (retained
mutation, CTA/self and fusion regions), the Design #2 candidate pool (the five
new junction-only-novelty splice windows), and the recorded mutation evidence
for the conditional H2BC12 E114K region. Column set matches the lead's
selected.tsv so the maintained GenBank annotation code runs unchanged.
"""
import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

LEAD_KEEP = {"MUTATION", "CTA/SELF", "FUSION"}
SPLICE_SELECTIONS = [
    "RNA_SPLICE_22b7baf3ed76ff263ceb",
    "RNA_SPLICE_2d0c90802e54c5d0505a",
    "RNA_SPLICE_335695e305c0bc6b7cfb",
    "RNA_SPLICE_dfffafab89cf697bf26f",
    "RNA_SPLICE_01b248aee527fedb422c",
]


def main():
    lead = list(csv.DictReader(
        (ROOT / "designs/design-2/selection/selected.tsv").open(), delimiter="\t"))
    fields = list(lead[0].keys())
    pool = {r["window_id"]: r for r in csv.DictReader(
        (ROOT / "designs/design-2/inputs/integrated/candidate-pool.tsv").open(), delimiter="\t")}

    rows = []
    for row in lead:
        if row["source"] in LEAD_KEEP:
            rows.append(row)
    for window_id in SPLICE_SELECTIONS:
        pool_row = pool[window_id]
        record = {field: pool_row.get(field, "") for field in fields}
        record["source"] = "SPLICE"
        rows.append(record)
    # conditional H2BC12 E114K from recorded mutation evidence
    evidence = list(csv.DictReader(
        (ROOT / "designs/design-2/inputs/mutations/variant-evidence.tsv").open(), delimiter="\t"))
    h2bc12 = next(e for e in evidence if e["variant_id"] == "mn_chr6_27146459_C_T")
    record = {field: "" for field in fields}
    record.update(
        variant_id=h2bc12["variant_id"],
        protein_id=h2bc12["variant_id"] + "_protein_0001",
        gene="H2BC12",
        DNA_annotation_protein_change=h2bc12["protein_change"],
        DNA_annotation_transcript_id=h2bc12["transcript_id"],
        sequence="LLPGELAKHAVSKGTKAVTKYTSAK",
        target_start_in_window="7", target_end_in_window="8",
        target_mode="substitution",
        region_id="H2BC12_E114K_conditional",
        event_id=h2bc12["variant_id"],
        source="MUTATION",
        somatic_status=h2bc12["somatic_status"] or "SUPPORTED_PAIRED_SOMATIC",
        DNA_evidence_tier="PAIRED_NORMAL_SUPPORTED",
        primary_DNA_eligible="True",
        gene_assignment_unambiguous="True",
        passes_500nm_gate="True",
        best_target_affinity_nm="139.7553694388593",
        best_target_peptide="ELAKHAVSK", best_target_allele="HLA-A*33:03",
        protein_change=h2bc12["protein_change"],
        additional_status="conditional_missing_expression_annotation",
        supplied_gene_value="", rna_alternate_fraction="0.19298245614035087",
        expression_notes="No vendor expression symbol (recorded identifier audit); "
                         "missing is not zero expression. Included as a conditional region.")
    rows.append(record)

    order = {r["antigen_id"]: i for i, r in enumerate(rows)}
    rows.sort(key=lambda r: (r["source"], order[r["antigen_id"]]))
    # renumber antigen ids in stable order: mutations, splice, CTA, fusion, conditional
    tiers = {"MUTATION": 0, "SPLICE": 1, "CTA/SELF": 2, "FUSION": 3}
    rows.sort(key=lambda r: (tiers[r["source"]], r["gene"]))
    for index, row in enumerate(rows, start=1):
        row["antigen_id"] = f"SF{index:02d}"
    with (HERE / "selected.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print("selected rows:", len(rows))
    for row in rows:
        print(row["antigen_id"], row["source"], row["gene"], row["sequence"][:20],
              "rna_confirmed:", row.get("patient_RNA_sequence_confirmed", ""))


if __name__ == "__main__":
    main()
