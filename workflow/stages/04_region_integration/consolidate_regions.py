"""Select the highest-affinity eligible 25-aa window per event/gene."""
import collections
import csv
import json
from pathlib import Path
import sys


def read(root, name):
    with (root / name).open() as source:
        return list(csv.DictReader(source, delimiter="\t"))


def main(root):
    root = Path(root)
    windows = read(root, "all-scored-windows.tsv")
    best = {}
    for row in windows:
        if row["primary_window_eligible"] != "True":
            continue
        key = row["source"], row["event_id"]
        if key not in best or (-float(row["native_affinity_sum"]), row["window_id"]) < (-float(best[key]["native_affinity_sum"]), best[key]["window_id"]):
            best[key] = row
    mutations = read(root, "mutation-regions.tsv")
    assert len(mutations) == 33
    rows = [dict(row, current_construct_region=False) for _, row in sorted(best.items())] + mutations
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with (root / "compact-pool.tsv").open("w") as destination:
        writer = csv.DictWriter(destination, columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    summary = dict(status="COMPLETE", scored_window_placements=len(windows),
        compact_rows=len(rows), compact_sources=dict(collections.Counter(row["source"] for row in rows)),
        mutation_regions_unchanged=33, window_selection="Highest native affinity sum; window ID resolves ties")
    for name in ["summary.json", "status.json"]:
        (root / name).write_text(json.dumps(summary, indent=2) + "\n")
    return summary


if __name__ == "__main__":
    print(json.dumps(main(sys.argv[1]), indent=2))
