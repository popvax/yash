# Extracted native lenstools_prepare_cta_annotation_ref; only input path uses argv.
import csv
from collections import defaultdict

import sys
annot_path = sys.argv[1]
collapsed_path = "cta_annotation_ref.collapsed.tsv"
report_path = "cta_annotation_ref.duplicate_transcripts.tsv"

missing_values = {"", "NA", "N/A", "None", "none", "nan"}
union_columns = {
    "gene_detectable_normal_tissues": ",",
    "gene_main_subcellular_location": ";",
    "gene_addtl_subcellular_location": ";",
}

with open(annot_path, newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    fieldnames = reader.fieldnames or []
    if "transcript" not in fieldnames:
        raise SystemExit("CTA annotation reference is missing required column: transcript")
    rows = list(reader)

groups = defaultdict(list)
for row in rows:
    groups[row.get("transcript", "")].append(row)

def split_values(value, sep):
    if value is None:
        return []
    value = value.strip()
    if value in missing_values:
        return []
    return [item.strip() for item in value.split(sep) if item.strip() and item.strip() not in missing_values]

def union_values(rows, column, sep):
    values = set()
    for row in rows:
        values.update(split_values(row.get(column, ""), sep))
    if values:
        return sep.join(sorted(values))
    for row in rows:
        value = row.get(column, "")
        if value and value not in missing_values:
            return value
    return "N/A"

collapsed_rows = []
duplicate_reports = []
for transcript, tx_rows in groups.items():
    base = dict(tx_rows[0])
    for column, sep in union_columns.items():
        if column in fieldnames:
            base[column] = union_values(tx_rows, column, sep)
    collapsed_rows.append(base)

    if len(tx_rows) > 1:
        conflicting_non_union = []
        for column in fieldnames:
            if column in union_columns:
                continue
            values = {row.get(column, "") for row in tx_rows}
            if len(values) > 1:
                conflicting_non_union.append(column)
        gene_names = sorted({row.get("gene_name", "") for row in tx_rows if row.get("gene_name", "")})
        duplicate_reports.append({
            "transcript": transcript,
            "row_count": str(len(tx_rows)),
            "unique_row_count": str(len({tuple(row.get(column, "") for column in fieldnames) for row in tx_rows})),
            "gene_names": ",".join(gene_names),
            "unioned_columns": ",".join(column for column in union_columns if column in fieldnames),
            "conflicting_non_union_columns": ",".join(conflicting_non_union),
        })

    # lenstools strips each line before splitting on tabs, which drops trailing empty
    # fields and can break annotation lookups by column index. Normalize missing cells.
    for column in fieldnames:
        value = base.get(column, "")
        if value is None or value in missing_values:
            base[column] = "N/A"

with open(collapsed_path, "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", extrasaction="ignore")
    writer.writeheader()
    writer.writerows(collapsed_rows)

with open(report_path, "w", newline="") as handle:
    report_fields = ["transcript", "row_count", "unique_row_count", "gene_names", "unioned_columns", "conflicting_non_union_columns"]
    writer = csv.DictWriter(handle, fieldnames=report_fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(duplicate_reports)
