"""Strict parsing and complete-record comparison for the numerical appendix."""
import csv
import gzip
import hashlib
import json
import math
import os
import re

INPUTS = ("selection/selected.tsv", "immunogenicity/class1-counterfactual-contrasts.tsv",
          "immunogenicity/mixmhcpred3.tsv", "safety/offtarget-homologs.tsv.gz",
          "construct/all-junction-affinity.tsv.gz")
SHORT = {"HLA-A*01:01": "A0101", "HLA-A*33:03": "A3303", "HLA-B*08:01": "B0801",
         "HLA-B*15:08": "B1508", "HLA-C*01:02": "C0102", "HLA-C*07:02": "C0702"}
GUARDANT_COVERED = {"RXRA", "TGFBR2", "CDH1::NRG1"}
PINNED_MODEL_AGGREGATE = "fae75536be69324a9ea891b2e91c8d346b0f5d48e0a7fe073db7ee699807e7c9"
PINNED_INPUT_DIGESTS = {
    "construct/all-junction-affinity.tsv.gz": "c843a7adf31698f7c10b161cc94207a691a5aafbb79d4273fe45dbbe77da96c6",
    "immunogenicity/class1-counterfactual-contrasts.tsv": "671bf3aaf8fb5745bb0472f9b452ac3111cb9d73ae72c728c6c50659e38cf650",
    "immunogenicity/mixmhcpred3.tsv": "3d96a99f5a99cc5e7bb9b7d3bb418329db468dbda01c1c6721ef72bb75621124",
    "safety/offtarget-homologs.tsv.gz": "79ac1fca416697d7beb8cd2ac5fa0be120dc14536a3318b8557be98a64bbe3d8",
    "selection/selected.tsv": "5a6a6b5aa4059b0abf8428ad6d81fcab9ea117d6e995cfa5786c345e21a90233"}
SOURCE_LABEL = "designs/design-4 (unmodified)"
LOOKUP_METHOD = ("exact (peptide, allele) key; source 1 = construct/all-junction-affinity.tsv.gz, "
                 "source 2 = shared_allele_homolog_nm from safety/offtarget-homologs.tsv.gz")
ENCODED = ("gene", "mutant_peptide", "wildtype_peptide", "allele", "mutant_affinity_nm", "wildtype_affinity_nm",
           "dai_wt_over_mut", "is_design4_lead", "binds_le_500nM", "encoded_in_design4_window")
REGION = ("gene", "window_aa", "design4_lead", "design4_lead_allele", "design4_lead_nm", "design4_lead_dai",
          "best_dai_in_window_le_500nM", "peptide_achieving_best_dai", "allele_achieving_best_dai",
          "higher_dai_choice_exists", "choice_kind", "n_mutation_spanning_pep_allele_dai_le_1_and_le500nM",
          "guardant360_covered", "dai_defined")
RESIDUAL = ("gene", "peptide", "allele", "mutant_affinity_nm", "dai", "encoded_window_aa", "note")
TIER = ("query_peptide", "homolog_peptide", "shared_allele", "mhcflurry_query_nm", "mhcflurry_homolog_nm",
        "mhcflurry_relative", "mhcflurry_query_percentile", "mhcflurry_homolog_percentile", "mixmhcpred_query_rank",
        "mixmhcpred_homolog_rank", "both_exceed_500nM_mhcflurry", "both_le_2pct_mixmhcpred",
        "both_le_2pct_mhcflurry_percentile", "homolog_genes")
ASSET_NAMES = {"allele_sequences.csv", "info.txt", "manifest.csv", "percent_ranks.json",
               "pseudosequences.mhcflurry.39aa.csv"} | {
    "weights_PAN-CLASS1-1-" + s + ".npz" for s in (
        "023e475ab69a6b2e", "21872705eecbad28", "5c21c209aa71dcbe", "5e8cfd8467919f6a",
        "a1852d5c4dd657ba", "b8cfcfd3e5868e3d", "c734d88715df1670", "e7ade34181bb9774",
        "edb7484c5406eaca", "f38e10af5757cedc", "f54dfdc7a7f34534")}


class Invalid(ValueError):
    """A missing, malformed or inconsistent declared input/output."""


def require(condition, message):
    if not condition:
        raise Invalid(message)


def _resolve(path):
    s = str(path)
    if not os.path.exists(s) and not s.endswith(".gz") and os.path.exists(s + ".gz"):
        return s + ".gz", True
    return s, False


def sha256(path):
    resolved, is_gz_fallback = _resolve(path)
    if is_gz_fallback:
        with gzip.open(resolved, "rb") as f:
            return hashlib.file_digest(f, "sha256").hexdigest()
    with open(resolved, "rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def open_text_auto(path):
    resolved, _ = _resolve(path)
    if resolved.endswith(".gz"):
        return gzip.open(resolved, "rt", encoding="utf-8", newline="")
    return open(resolved, "r", encoding="utf-8", newline="")


def keys(value, expected, label):
    require(type(value) is dict and set(value) == set(expected), f"{label}: expected exactly fields {sorted(expected)}")


def equal(actual, expected, label):
    """Compare JSON recursively without treating a Boolean as a count."""
    require(type(actual) is type(expected), f"{label}: wrong type")
    if type(expected) is dict:
        keys(actual, expected, label)
        for k, v in expected.items():
            equal(actual[k], v, f"{label}.{k}")
    elif type(expected) is list:
        require(len(actual) == len(expected), f"{label}: wrong list length")
        for i, (a, e) in enumerate(zip(actual, expected)):
            equal(a, e, f"{label}[{i}]")
    else:
        require(actual == expected, f"{label}: expected {expected!r}, got {actual!r}")


def number(value, label, *, blank=False, positive=False, percentile=False, integer=False):
    if blank and value == "":
        return None
    require(type(value) is str and value.strip() == value and value != "", f"{label}: numeric cell required")
    if integer:
        require(re.fullmatch(r"0|[1-9][0-9]*", value) is not None, f"{label}: nonnegative integer required")
    try:
        result = float(value)
    except (ValueError, OverflowError):
        raise Invalid(f"{label}: invalid number {value!r}") from None
    require(math.isfinite(result), f"{label}: non-finite value")
    require(not positive or result > 0, f"{label}: positive value required")
    require(not percentile or 0 <= result <= 100, f"{label}: percentile/rank outside [0, 100]")
    return result


def boolean(value, label):
    require(value in ("True", "False"), f"{label}: expected True or False")
    return value == "True"


def residual_low_dai(dai, mutant_nm):
    """Apply the residual cutoff to unrounded source values, before display rounding."""
    return dai <= 1.0 and mutant_nm <= 500.0


def read_json(path):
    def pairs(values):
        result = {}
        for k, v in values:
            require(k not in result, f"{path}: duplicate JSON key {k}")
            result[k] = v
        return result
    def constant(value):
        raise Invalid(f"{path}: non-finite JSON number {value}")
    resolved, _ = _resolve(path)
    opener = gzip.open if resolved.endswith(".gz") else open
    with opener(resolved, "rt") as f:
        return json.load(f, object_pairs_hook=pairs, parse_constant=constant)


def read_tsv(path, schema=None, required=(), comments=False):
    resolved, _ = _resolve(path)
    opener = gzip.open if resolved.endswith(".gz") else open
    with opener(resolved, "rt", newline="") as f:
        reader = csv.DictReader((x for x in f if not comments or not x.startswith("#")), delimiter="\t", strict=True)
        header = reader.fieldnames or []
        require(len(set(header)) == len(header), f"{path}: duplicate header")
        if schema is not None:
            require(header == list(schema), f"{path}: wrong table schema")
        require(set(required) <= set(header), f"{path}: required source columns missing")
        rows = list(reader)
    for i, row in enumerate(rows, 2):
        require(None not in row and all(v is not None for v in row.values()), f"{path}:{i}: wrong field count")
    return rows


def unique(rows, fields, label):
    result = {}
    for row in rows:
        key = tuple(row[f] for f in fields)
        require(key not in result, f"{label}: duplicate key {key}")
        result[key] = row
    return result


def compare_records(path, schema, expected, key_fields, numeric=(), integers=(), booleans=()):
    rows = read_tsv(path, schema)
    actual = unique(rows, key_fields, path.name)
    want = unique(expected, key_fields, "derived " + path.name)
    require(set(actual) == set(want), f"{path.name}: missing or unexpected row keys")
    for key, row in actual.items():
        exp = want[key]
        for field in schema:
            label = f"{path.name} {key} {field}"
            if field in numeric or field in integers:
                got = number(row[field], label, blank=exp[field] is None,
                             positive=field in numeric, integer=field in integers)
                require((got is None and exp[field] is None) or
                        (got is not None and exp[field] is not None and abs(got - exp[field]) <= 1e-8),
                        f"{label}: differs from source")
            elif field in booleans:
                equal(boolean(row[field], label), exp[field], label)
            else:
                equal(row[field], exp[field], label)
    return rows


def aggregate(files):
    return hashlib.sha256("\n".join(f"{k}:{v}" for k, v in sorted(files.items())).encode()).hexdigest()


def hash_assets(directory):
    return {os.path.relpath(os.path.join(root, name), directory): sha256(os.path.join(root, name))
            for root, _, names in os.walk(directory) for name in names if not name.endswith(".csv.bz2")}


def validate_asset_receipt(assets):
    keys(assets, ("status", "resolved_affinity_dir", "loaded_explicitly_from_this_dir", "asset_count",
                  "aggregate_sha256", "percent_ranks_json_sha256", "manifest_csv_sha256",
                  "excluded_from_hash", "per_file_sha256"), "model assets")
    equal(assets["status"], "hashed", "model assets.status")
    equal(assets["loaded_explicitly_from_this_dir"], True, "model assets.loaded_explicitly_from_this_dir")
    require(type(assets["resolved_affinity_dir"]) is str and os.path.isabs(assets["resolved_affinity_dir"]),
            "model assets.resolved_affinity_dir: absolute generation path required")
    equal(assets["asset_count"], len(ASSET_NAMES), "model assets.asset_count")
    equal(assets["excluded_from_hash"], "*.csv.bz2 (training data)", "model assets.excluded_from_hash")
    files = assets["per_file_sha256"]
    keys(files, ASSET_NAMES, "model assets.per_file_sha256")
    require(all(type(v) is str and re.fullmatch(r"[0-9a-f]{64}", v) for v in files.values()), "model assets: invalid per-file SHA-256")
    equal(assets["aggregate_sha256"], aggregate(files), "model assets.aggregate_sha256")
    equal(assets["aggregate_sha256"], PINNED_MODEL_AGGREGATE, "pinned model asset identity")
    equal(assets["percent_ranks_json_sha256"], files["percent_ranks.json"], "model assets.percent_ranks_json_sha256")
    equal(assets["manifest_csv_sha256"], files["manifest.csv"], "model assets.manifest_csv_sha256")


# Canonical descriptive metadata; these constants are not read back from receipts.
DIGEST_NOTE = 'All five Design 4 inputs read by this appendix. These hashes do not verify other Design 4 files or authenticate specimen custody.'
FRESH_NOTE = 'mhcflurry_* columns are fresh model output. All 46 affinity cells are corroborated against stored Design 4 sources on exact (peptide, allele) keys: 26 from construct/all-junction-affinity.tsv.gz and 20 from shared_allele_homolog_nm in safety/offtarget-homologs.tsv.gz. Ordinary verification repeats this lookup, checks every table field against its applicable source or derivation, and verifies table integrity. Percentile values are checked for range and consistency; --recompute independently reruns the model to check every affinity and percentile value.'
INTERPRETATION = "The >500 nM group is defined by an MHCflurry absolute-nM threshold. Design 4's predicate is not nM-ONLY: every pair in that group has BOTH MixMHCpred ranks <=2% on the shared allele, so all of them satisfy the original criterion. Different cutoffs on MHCflurry can classify a pair differently; the comparison of MHCflurry and MixMHCpred ordering is reported separately. Neither observation establishes that a flag is wrong, allele occupancy, cross-reactivity or harm. The groups also overlap: the >500 nM pairs are a subset of the stronger/weaker split, not a third disjoint category."
CROSS_MODEL_NOTE = 'Of the pairs where both MixMHCpred ranks are available, this many order the query and homolog in the opposite direction from MHCflurry. That is a genuine cross-model ordering difference. It quantifies no clinical consequence.'
SHARED_PREDICATE = 'MHCflurry nM <= 500 OR MHCflurry percentile <= 2 OR MixMHCpred rank <= 2, per peptide on the same allele. A disjunction, not an nM-only rule.'
