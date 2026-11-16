"""Build the construct region table from the design's selected.tsv (single
source of truth for antigen identity, order-independent)."""
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main():
    selected = list(csv.DictReader((HERE / "selected.tsv").open(), delimiter="\t"))
    fields = ["antigen_id", "gene", "source", "evidence_tier", "sequence",
              "target_start_in_window", "target_end_in_window", "region_id",
              "window_id", "event_id", "note"]
    with (HERE / "regions.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in selected:
            writer.writerow(dict(
                antigen_id=row["antigen_id"], gene=row["gene"], source=row["source"],
                evidence_tier=row.get("additional_status") or (
                    "junction_only_novelty" if row["source"] == "SPLICE"
                    else "paired_DNA_supported" if row["source"] == "MUTATION"
                    else "recorded_RNA_source"),
                sequence=row["sequence"],
                target_start_in_window=row["target_start_in_window"] or "0",
                target_end_in_window=row["target_end_in_window"] or "25",
                region_id=row["region_id"], window_id=row.get("window_id", ""),
                event_id=row.get("event_id", ""),
                note=f"protein_change={row.get('protein_change', '')}"))
    print("regions:", len(selected))
    assert len({r["sequence"] for r in selected}) == len(selected)
    assert all(len(r["sequence"]) == 25 for r in selected)


if __name__ == "__main__":
    main()
