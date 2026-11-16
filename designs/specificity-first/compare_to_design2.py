"""Compare this design's selected regions with lead Design #2.

The repository's `yash-vaccine compare` helper only supports recipe-stage
outputs, so this script compares the two independent designs' selected.tsv
tables directly: shared regions (matched on gene + source + sequence),
regions removed, regions added, and per-source counts. Output:
comparison-vs-design-2/comparison.json (+ .tsv).
"""
import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
LEAD = HERE.parent / "design-2" / "selection" / "selected.tsv"


def rows(path):
    return list(csv.DictReader(path.open(newline=""), delimiter="\t"))


def key(row):
    return (row["gene"], row["source"], row["sequence"])


def main():
    lead, new = rows(LEAD), rows(HERE / "selected.tsv")
    lead_keys = {key(r): r for r in lead}
    new_keys = {key(r): r for r in new}
    shared = sorted(set(lead_keys) & set(new_keys))
    removed = sorted(set(lead_keys) - set(new_keys))
    added = sorted(set(new_keys) - set(lead_keys))

    def summary(keys, table):
        out = []
        for gene, source, sequence in keys:
            row = table[(gene, source, sequence)]
            out.append(dict(
                gene=gene, source=source, region_id=row["region_id"],
                antigen_id=row.get("antigen_id", ""),
                sequence=sequence,
                best_target_affinity_nm=row.get("best_target_affinity_nm", ""),
                additional_status=row.get("additional_status", ""),
                native_event_review_flags=row.get("native_event_review_flags", "")))
        return out

    counts = lambda table: {
        source: sum(1 for r in table if r["source"] == source)
        for source in sorted({r["source"] for r in table})}
    result = dict(
        lead_design="designs/design-2 (24 regions, 724 aa protein)",
        this_design="designs/specificity-first (25 regions, 748 aa protein)",
        shared_regions=summary(shared, new_keys),
        removed_from_lead=summary(removed, lead_keys),
        added_here=summary(added, new_keys),
        region_counts=dict(lead=counts(lead), specificity_first=counts(new)),
        protein=dict(lead_aa=724, specificity_first_aa=748))
    out = HERE / "comparison-vs-design-2"
    out.mkdir(exist_ok=True)
    (out / "comparison.json").write_text(json.dumps(result, indent=2) + "\n")
    with (out / "comparison.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=[
            "change", "gene", "source", "region_id", "antigen_id", "sequence",
            "best_target_affinity_nm", "additional_status", "native_event_review_flags"],
            delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for label, entries in [("shared", result["shared_regions"]),
                               ("removed_from_lead", result["removed_from_lead"]),
                               ("added_here", result["added_here"])]:
            for entry in entries:
                writer.writerow({"change": label, **entry})
    print(json.dumps({k: result[k] for k in ("region_counts", "protein")}, indent=2))
    print("shared:", len(shared), "removed:", len(removed), "added:", len(added))


if __name__ == "__main__":
    main()
