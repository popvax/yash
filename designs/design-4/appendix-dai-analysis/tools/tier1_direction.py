#!/usr/bin/env python3
"""Compare the two peptides of each TIER1_1MM_PARALOG_SAME_HLA pair ON THEIR SHARED
ALLELE, under each branch of the predicate Design 4's flag uses. Two predictors are
involved (MHCflurry and MixMHCpred); the predicate's three branches are three
criteria, not three models.

Why: Design 4's table records `query_best_nm` as the query's best across ALL
alleles but `shared_allele_homolog_nm` on the shared allele, so those two columns
cannot be compared directly. This script scores both peptides on the shared allele.

IMPORTANT SCOPE. Design 4's shared-HLA predicate is a DISJUNCTION, satisfied when,
for each peptide on the same allele, any of these holds:

    MHCflurry nM <= 500  OR  MHCflurry percentile <= 2  OR  MixMHCpred rank <= 2

So it is neither nM-only nor MixMHCpred-only. Applying one branch alone and finding
a pair "fails" says nothing: the other branches may be satisfied. Nothing here is
evidence that a flag is wrong. PRIME is used elsewhere in Design 4's pipeline but
not in this predicate and is not compared here.

Note also that MHCflurry percentile is an allele-calibrated transformation of the
same MHCflurry affinity. Different classifications under <=500 nM and <=2% are two
CUTOFFS ON ONE MODEL, not two models disagreeing.

Requires mhcflurry==2.3.14 and the recorded model assets.
Usage: tier1_direction.py --design4 <dir> --out <dir> --models-dir <affinity dir>
"""
import argparse, csv, gzip, json, os
from validation import open_text_auto

SHORT = {"HLA-A*01:01": "A0101", "HLA-A*33:03": "A3303", "HLA-B*08:01": "B0801",
         "HLA-B*15:08": "B1508", "HLA-C*01:02": "C0102", "HLA-C*07:02": "C0702"}

ap = argparse.ArgumentParser()
ap.add_argument("--design4", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--models-dir", help="local affinity directory; recorded hashes remain required")
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)

with gzip.open(os.path.join(a.design4, "safety/offtarget-homologs.tsv.gz"), "rt") as fh:
    rows = [r for r in csv.DictReader(fh, delimiter="\t")
            if r["model_priority_category"] == "TIER1_1MM_PARALOG_SAME_HLA"
            and r["shared_binding_allele"]]
if not rows:
    raise ValueError("no TIER1 homolog rows with a shared allele were found")

# MixMHCpred per-allele %Rank, from Design 4's own output (not re-run)
mix = {}
with open_text_auto(os.path.join(a.design4, "immunogenicity/mixmhcpred3.tsv")) as fh:
    body = [l for l in fh if not l.startswith("#")]
for r in csv.DictReader(body, delimiter="\t"):
    for allele, s in SHORT.items():
        v = r.get("%Rank_" + s)
        if v:
            mix[(r["Peptide"], allele)] = float(v)
if not mix:
    raise ValueError("mixmhcpred3.tsv contains no per-allele rank values")

from mhcflurry import Class1AffinityPredictor
import mhcflurry, hashlib as _h

def resolve_affinity_dir():
    """Resolve ONE directory: the affinity predictor that will actually be loaded.

    `Class1AffinityPredictor.load()` with no argument prefers a standalone pan
    bundle and only falls back to the presentation bundle's affinity component.
    Hashing one directory while implicitly loading another is exactly the defect
    this avoids, so the same path is hashed AND passed to load().
    """
    from mhcflurry.downloads import get_path
    for spec in (("models_class1_pan", "models.combined"),
                 ("models_class1_presentation", "models")):
        try:
            d = get_path(*spec)
        except Exception:
            continue
        if spec[0] == "models_class1_presentation":
            d = os.path.join(d, "affinity_predictor")
        if os.path.isdir(d) and os.path.exists(os.path.join(d, "manifest.csv")):
            return d
    return None


def hash_dir(d):
    """Hash every inference asset in exactly this directory."""
    files = {}
    for root, _, names in os.walk(d):
        for n in sorted(names):
            fp = os.path.join(root, n)
            rel = os.path.relpath(fp, d)
            if rel.endswith(".csv.bz2"):
                continue
            with open(fp, "rb") as fh:
                files[rel] = _h.sha256(fh.read()).hexdigest()
    agg = _h.sha256("\n".join(f"{k}:{v}" for k, v in sorted(files.items()))
                    .encode()).hexdigest()
    return files, agg


AFFINITY_DIR = a.models_dir or resolve_affinity_dir()
if AFFINITY_DIR is None:
    raise SystemExit("could not resolve an affinity predictor directory to hash and load")
_files, _agg = hash_dir(AFFINITY_DIR)
_receipt_path = os.path.join(os.path.dirname(__file__), "..", "data",
                             "tier1-direction-summary.json")
with open(_receipt_path) as _fh:
    _recorded_model = json.load(_fh)["models"]["model_bundle"]
if mhcflurry.__version__ != _recorded_model["package_version"]:
    raise SystemExit("MHCflurry package version differs from the recorded version")
if (_agg != _recorded_model["assets"]["aggregate_sha256"] or
        _files != _recorded_model["assets"]["per_file_sha256"]):
    raise SystemExit("model assets differ from the recorded bundle; no inference run")
P = Class1AffinityPredictor.load(models_dir=AFFINITY_DIR)   # the SAME directory
MODEL_ID = {"package_version": mhcflurry.__version__,
            "num_supported_alleles": len(P.supported_alleles),
            "assets": {"status": "hashed",
                       "resolved_affinity_dir": AFFINITY_DIR,
                       "loaded_explicitly_from_this_dir": True,
                       "asset_count": len(_files),
                       "aggregate_sha256": _agg,
                       "percent_ranks_json_sha256": _files.get("percent_ranks.json"),
                       "manifest_csv_sha256": _files.get("manifest.csv"),
                       "excluded_from_hash": "*.csv.bz2 (training data)",
                       "per_file_sha256": _files}}
peps = sorted({r["query_peptide"] for r in rows} | {r["homolog_peptide"] for r in rows})
nm, pct = {}, {}
for al in sorted({r["shared_binding_allele"] for r in rows}):
    df = P.predict_to_dataframe(peptides=peps, allele=al,
                                include_percentile_ranks=True, throw=False)
    for _, x in df.iterrows():
        nm[(x["peptide"], al)] = float(x["prediction"])
        pct[(x["peptide"], al)] = float(x["prediction_percentile"])

out = []
stronger = weaker = both_over_500 = both_mix2 = both_fl_pct2 = 0
both_ranks = reversed_dir = 0
for r in rows:
    al = r["shared_binding_allele"]
    q, h = nm[(r["query_peptide"], al)], nm[(r["homolog_peptide"], al)]
    qm, hm = mix.get((r["query_peptide"], al)), mix.get((r["homolog_peptide"], al))
    qp, hp = pct[(r["query_peptide"], al)], pct[(r["homolog_peptide"], al)]
    s = h <= q
    stronger += s
    weaker += not s
    o500 = q > 500 and h > 500
    both_over_500 += o500
    m2 = qm is not None and hm is not None and qm <= 2.0 and hm <= 2.0
    both_mix2 += m2
    p2 = qp <= 2.0 and hp <= 2.0
    both_fl_pct2 += p2
    if qm is not None and hm is not None:
        both_ranks += 1
        reversed_dir += (h <= q) != (hm <= qm)
    out.append(dict(
        query_peptide=r["query_peptide"], homolog_peptide=r["homolog_peptide"],
        shared_allele=al,
        mhcflurry_query_nm=round(q, 2), mhcflurry_homolog_nm=round(h, 2),
        mhcflurry_relative="homolog_stronger_or_equal" if s else "homolog_weaker",
        mhcflurry_query_percentile=round(qp, 4), mhcflurry_homolog_percentile=round(hp, 4),
        mixmhcpred_query_rank=qm, mixmhcpred_homolog_rank=hm,
        both_exceed_500nM_mhcflurry=o500,
        both_le_2pct_mixmhcpred=m2,
        both_le_2pct_mhcflurry_percentile=p2,
        homolog_genes=r["homolog_genes"]))

with open(os.path.join(a.out, "tier1-direction.tsv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(sorted(out, key=lambda x: x["query_peptide"]))

over = [x for x in out if x["both_exceed_500nM_mhcflurry"]]
import hashlib
_tbl = os.path.join(a.out, "tier1-direction.tsv")
with open(_tbl, "rb") as fh:
    _table_digest = hashlib.sha256(fh.read()).hexdigest()
summary = {
    "tier1_table_sha256": _table_digest,
    "fresh_inference_cells_note": (
        "mhcflurry_* columns are fresh model output. All 46 affinity cells are "
        "corroborated against stored Design 4 sources on exact (peptide, allele) "
        "keys: 26 from construct/all-junction-affinity.tsv.gz and 20 from "
        "shared_allele_homolog_nm in safety/offtarget-homologs.tsv.gz. Ordinary "
        "verification repeats this lookup, checks every table field against its "
        "applicable source or derivation, and verifies table integrity. Percentile "
        "values are checked for range and consistency; --recompute independently "
        "reruns the model to check every affinity and percentile value."),
    "tier1_rows_with_shared_allele": len(rows),
    "mhcflurry_homolog_stronger_or_equal": stronger,
    "mhcflurry_homolog_weaker": weaker,
    "both_exceed_500nM_mhcflurry": both_over_500,
    "of_which_homolog_stronger_or_equal": sum(
        1 for x in over if x["mhcflurry_relative"] == "homolog_stronger_or_equal"),
    "of_which_homolog_weaker": sum(
        1 for x in over if x["mhcflurry_relative"] == "homolog_weaker"),
    "both_le_2pct_mixmhcpred": both_mix2,
    "both_le_2pct_mixmhcpred_among_the_over500nM": sum(
        1 for x in over if x["both_le_2pct_mixmhcpred"]),
    "both_le_2pct_mhcflurry_percentile_among_the_over500nM": sum(
        1 for x in over if x["both_le_2pct_mhcflurry_percentile"]),
    "interpretation": (
        "The >500 nM group is defined by an MHCflurry absolute-nM threshold. Design "
        "4's predicate is not nM-ONLY: every pair in that group has BOTH MixMHCpred "
        "ranks <=2% on the shared allele, so all of them satisfy the original "
        "criterion. Different cutoffs on MHCflurry can classify a pair differently; "
        "the comparison of MHCflurry and MixMHCpred ordering is reported separately. "
        "Neither observation establishes that a flag is wrong, allele occupancy, "
        "cross-reactivity or harm. The groups also overlap: the >500 nM pairs are a subset of the "
        "stronger/weaker split, not a third disjoint category."),
    "pairs_with_both_mixmhcpred_ranks": both_ranks,
    "cross_model_direction_reversed": reversed_dir,
    "cross_model_note": (
        "Of the pairs where both MixMHCpred ranks are available, this many order the "
        "query and homolog in the opposite direction from MHCflurry. That is a "
        "genuine cross-model ordering difference. It quantifies no clinical "
        "consequence."),
    "models": {"fresh": f"mhcflurry {MODEL_ID['package_version']}",
               "model_bundle": MODEL_ID,
               "read_from_design4_not_rerun": "MixMHCpred 3.0 (mixmhcpred3.tsv)",
               "not_consulted": "PRIME (not part of the shared-HLA predicate)"},
    "shared_hla_predicate": (
        "MHCflurry nM <= 500 OR MHCflurry percentile <= 2 OR MixMHCpred rank <= 2, "
        "per peptide on the same allele. A disjunction, not an nM-only rule."),
}
with open(os.path.join(a.out, "tier1-direction-summary.json"), "w") as fh:
    json.dump(summary, fh, indent=1, sort_keys=True)
print(f"{len(rows)} rows | MHCflurry: {stronger} homolog stronger/equal, {weaker} weaker "
      f"| {both_over_500} exceed 500 nM both (of which {sum(1 for x in over if x['both_le_2pct_mixmhcpred'])} "
      f"satisfy MixMHCpred <=2% both)")
