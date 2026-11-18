"""Strict, unpadded 19-aa hydropathy screen for a routed protein construct."""
from collections import Counter

WINDOW = 19
KD_TENTHS = dict(zip(
    "ARNDCQEGHILKMFPSTWYV",
    [18, -45, -35, -35, 25, -35, -35, -4, -32, 45, 38, -39, 19, 28, -16, -8, -7, -9, -13, 42],
))
THRESHOLD_TENTHS = 16


def scan(sequence):
    if not isinstance(sequence, str) or not sequence or set(sequence) - KD_TENTHS.keys():
        raise ValueError("Expected a nonempty uppercase standard-amino-acid sequence")
    values = [KD_TENTHS[aa] for aa in sequence]
    return [
        {
            "start_aa": start,
            "end_aa": start + WINDOW,
            "sequence": sequence[start:start + WINDOW],
            "mean_kd": sum(values[start:start + WINDOW]) / (10 * WINDOW),
            "above_threshold": sum(values[start:start + WINDOW]) > THRESHOLD_TENTHS * WINDOW,
        }
        for start in range(max(0, len(sequence) - WINDOW + 1))
    ]


def summarize(sequence):
    windows = scan(sequence)
    if not windows:
        return {
            "status": "NOT_SCORABLE_SHORTER_THAN_19_AA",
            "max_mean_kd": None,
            "windows_scanned": 0,
            "above_threshold_windows": 0,
        }
    best = max(windows, key=lambda row: row["mean_kd"])
    return {
        "status": "EXCLUDE_HYDROPHOBIC" if any(
            row["above_threshold"] for row in windows
        ) else "PASS",
        "max_mean_kd": best["mean_kd"],
        "max_start_aa": best["start_aa"],
        "max_end_aa": best["end_aa"],
        "max_window": best["sequence"],
        "windows_scanned": len(windows),
        "above_threshold_windows": sum(row["above_threshold"] for row in windows),
    }


def protein_screen(sequence, parts):
    """Scan complete construct windows and report short antigens as not scorable."""
    cursor = 0
    routing = []
    antigen_status = []
    for part in parts:
        start, end = int(part["start_aa"]), int(part["end_aa"])
        direct = part["kind"] == "linker" and part["aa"] == "" and end == start
        if start != cursor or (end <= start and not direct) or sequence[start:end] != part["aa"]:
            raise ValueError("Protein components must tile the sequence exactly")
        cursor = end
        if part["kind"] == "antigen":
            individual = scan(part["aa"])
            antigen_status.append({
                "name": part["name"],
                "length_aa": len(part["aa"]),
                "status": (
                    "NOT_SCORABLE_SHORTER_THAN_19_AA"
                    if len(part["aa"]) < WINDOW
                    else "EXCLUDE_HYDROPHOBIC"
                    if any(row["above_threshold"] for row in individual)
                    else "PASS"
                ),
                "windows_scanned": len(individual),
            })
        if part["kind"] in ("signal_peptide", "mitd"):
            if (part["kind"] == "signal_peptide" and start != 0) or (
                part["kind"] == "mitd" and end != len(sequence)
            ):
                raise ValueError("Only terminal routing elements may be exempted")
            routing.append((start, end, part["kind"]))
    if cursor != len(sequence):
        raise ValueError("Protein components do not cover the full sequence")
    if len({kind for _, _, kind in routing}) != len(routing):
        raise ValueError("Duplicate routing element")

    windows = scan(sequence)
    for row in windows:
        start, end = row["start_aa"], row["end_aa"]
        touched = [kind for left, right, kind in routing if start < right and end > left]
        contained = any(left <= start and end <= right for left, right, _ in routing)
        row.update(
            classification=(
                "INTENTIONAL_ROUTING_DOMAIN" if contained
                else "ROUTING_BOUNDARY_REVIEW" if touched
                else "ANTIGEN_LINKER_CASSETTE"
            ),
            routing_overlap=";".join(touched),
        )

    internal = [row for row in windows if row["classification"] == "ANTIGEN_LINKER_CASSETTE"]
    unexpected = [row for row in internal if row["above_threshold"]]
    excluded_antigens = [
        row for row in antigen_status if row["status"] == "EXCLUDE_HYDROPHOBIC"
    ]
    report = {
        "status": "PASS" if not unexpected and not excluded_antigens else "FAIL",
        "policy": {
            "scale": "Kyte-Doolittle (1982)",
            "window_aa": WINDOW,
            "exclude_mean_strictly_greater_than": 1.6,
            "edge_padding": False,
        },
        "protein_length_aa": len(sequence),
        "windows_scanned": len(windows),
        "antigen_linker_windows": len(internal),
        "unexpected_high_windows": len(unexpected),
        "max_antigen_linker_mean_kd": max(
            (row["mean_kd"] for row in internal), default=None
        ),
        "all_antigen_regions_individually_pass": all(
            row["status"] == "PASS" for row in antigen_status
        ),
        "antigen_regions": antigen_status,
        "short_antigen_policy": (
            "Segments shorter than 19 aa are not individually scorable; "
            "no sequence is padded. Complete windows in the assembled construct are screened."
        ),
        "routing_boundary_policy": (
            "All routing-overlap scores are exported for review; exemptions do not "
            "establish cleavage, membrane topology, solubility or biological safety."
        ),
        "high_window_classifications": dict(Counter(
            row["classification"] for row in windows if row["above_threshold"]
        )),
    }
    return windows, report
