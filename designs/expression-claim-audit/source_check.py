"""First local source check against the bundled vendor expression table."""

import argparse
import hashlib
import json
from decimal import Decimal, InvalidOperation, ROUND_HALF_UP, localcontext
from pathlib import Path


SOURCE_RELATIVE = Path("data/expression/vendor-transformed-expression.tsv")
ROOT = next(
    (parent for parent in Path(__file__).resolve().parents if (parent / SOURCE_RELATIVE).is_file()),
    Path(__file__).resolve().parent,
)
SOURCE = ROOT / SOURCE_RELATIVE
EXPECTED_HASH = "cc45875929d679620c4fee759e06f6664a7e478234d8403b33d72b34043b98b9"
PRACTICE_SOURCE = Path(__file__).with_name("blank-zero-practice.tsv")
PRACTICE_HASH = "afdb1d21ee706f39e0057d20c86d7971276a912f2566d19af8f52dc84870fa15"
# Names from the repository's H2BC12-expression-identifier-audit.json.
H2BC12_NAMES = ("H2BC12", "H2BFT", "HIST1H2BK", "H2BFAiii")
SOURCE_SCALE = "vendor_transformed"
MAX_ROUNDING_PRECISION = 1000


def numeric_value(text):
    try:
        value = Decimal(text)
    except InvalidOperation as error:
        raise ValueError("{!r} is not a valid number.".format(text)) from error
    if not value.is_finite() or value < 0:
        raise ValueError("Values must be finite and nonnegative.")
    return value


def rounded_value(value, decimal_places):
    precision = max(28, len(value.as_tuple().digits) + max(0, value.adjusted()) + decimal_places + 2)
    if precision > MAX_ROUNDING_PRECISION:
        raise ValueError("Rounding request exceeds the 1000-digit working-precision limit.")
    # Increase precision so rounding large or very precise source values remains exact.
    with localcontext() as context:
        context.prec = precision
        unit = Decimal((0, (1,), -decimal_places))
        return value.quantize(unit, rounding=ROUND_HALF_UP)


def check(gene, sample, claim_text, scale=SOURCE_SCALE, practice=False, decimal_places=None):
    claimed_gene = gene
    claimed_sample = sample
    claimed_scale = scale
    gene = gene.strip()
    sample = sample.strip()
    scale = scale.strip()
    if not gene or not sample or not scale:
        raise ValueError("Gene, sample and scale must contain a non-whitespace character.")
    claim = numeric_value(claim_text)
    if decimal_places is not None:
        if decimal_places < 0:
            raise ValueError("Decimal places must be nonnegative.")
        if claim != rounded_value(claim, decimal_places):
            raise ValueError("Claimed number has more precision than the declared decimal places.")
    source = SOURCE
    expected_hash = EXPECTED_HASH
    if practice:
        source = PRACTICE_SOURCE
        expected_hash = PRACTICE_HASH
    source_bytes = source.read_bytes()
    fingerprint = hashlib.sha256(source_bytes).hexdigest()
    if fingerprint != expected_hash:
        raise ValueError("Source bytes changed, possibly including line endings. Review or restore the pinned source before proceeding.")

    if source_bytes.startswith(b"\xef\xbb\xbf"):
        raise ValueError("Source has a UTF-8 BOM; this table format requires UTF-8 without a BOM.")
    lines = [line.removesuffix("\r") for line in source_bytes.decode("utf-8").split("\n")]
    sample_headers = [line for line in lines if line.startswith("##SampleID=")]
    if len(sample_headers) != 1 or lines.count("##Gene\tTPM") != 1:
        raise ValueError("Expected sample identity and column header were not found.")
    source_sample = sample_headers[0].split("=", 1)[1]
    if not source_sample or source_sample != source_sample.strip():
        raise ValueError("Source sample ID is blank or has surrounding whitespace.")
    column_header_index = lines.index("##Gene\tTPM")
    if any(line and not line.startswith("#") for line in lines[:column_header_index]):
        raise ValueError("Source gene column header must precede all data rows.")
    if lines.index(sample_headers[0]) > column_header_index:
        raise ValueError("Source sample header must precede the gene column header.")
    rows = {}
    for line_number, line in enumerate(lines, start=1):
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 2:
            raise ValueError("Source line {}: expected two tab-separated fields.".format(line_number))
        name, text = fields
        if not name or name != name.strip() or name in rows:
            raise ValueError("Source line {}: blank, padded or duplicate gene name.".format(line_number))
        if text.strip():
            try:
                numeric_value(text)
            except ValueError as error:
                raise ValueError("Source line {} ({}): {}".format(line_number, name, error)) from error
        rows[name] = {"value": text, "line": line_number}
    names = H2BC12_NAMES if gene == "H2BC12" else (gene,)
    report = {
        "source": str(source.relative_to(ROOT)),
        "synthetic_practice": practice,
        "source_sha256": fingerprint,
        "source_sample": source_sample,
        "claimed_sample": claimed_sample,
        "sample": sample,
        "gene": gene,
        "claimed_gene": claimed_gene,
        "claimed_value": claim_text,
        "parsed_claim_value": str(claim),
        "measurement_scale": "supplied vendor-transformed gene value",
        "claimed_scale": claimed_scale,
        "scale": scale,
        "decimal_places": decimal_places,
        "comparison_rule": "exact decimal equality",
        "explanation_review": "needs human review",
    }
    if decimal_places is not None:
        report["comparison_rule"] = "source rounded with ROUND_HALF_UP to declared decimal places"
    if sample != source_sample:
        report.update(status="gene_or_sample_differs", reason="Sample ID differs from the file header.")
        return report
    if scale != SOURCE_SCALE:
        report.update(
            status="insufficient_information",
            reason="Claimed scale has no documented conversion to the supplied vendor-transformed scale.",
        )
        return report

    matches = [dict(gene=name, **rows[name]) for name in names if name in rows]
    report.update(checked_names=list(names), matching_rows=matches)
    if not matches:
        status = "insufficient_information"
        reason = "Exact name absent. Other gene identifiers have not been resolved."
        if gene == "H2BC12":
            status = "measurement_missing_from_source"
            reason = "No value found under the four names documented in the existing H2BC12 audit."
        report.update(status=status, reason=reason)
        return report

    present = [row for row in matches if row["value"].strip()]
    if not present:
        report.update(
            status="measurement_missing_from_source",
            reason="Matching gene rows exist, but their value fields are blank.",
        )
        return report
    if len(present) != len(matches):
        report.update(
            status="insufficient_information",
            reason="Checked names include both blank and numeric values; review their relationship.",
        )
        return report
    values = {numeric_value(row["value"]) for row in present}
    if len(values) != 1:
        report.update(status="insufficient_information", reason="Checked names have conflicting values.")
        return report
    comparison_value = next(iter(values))
    if decimal_places is not None:
        comparison_value = rounded_value(comparison_value, decimal_places)
    report["source_value_for_comparison"] = str(comparison_value)
    status = "value_differs"
    reason = "Compared exact numbers on the source's supplied scale."
    if claim == comparison_value:
        status = "value_matches"
        if decimal_places is not None:
            status = "value_matches_after_rounding"
    if decimal_places is not None:
        reason = "Compared with the source rounded to {} decimal places using ROUND_HALF_UP.".format(decimal_places)
    report.update(status=status, reason=reason)
    return report


def unique_claim_fields(pairs):
    fields = {}
    for key, value in pairs:
        if key in fields:
            raise ValueError("Duplicate claim field: {}".format(key))
        fields[key] = value
    return fields


class JsonNumber(str):
    """Keep the literal JSON number for reporting and exact Decimal parsing."""


def read_claim(path):
    claim = json.loads(
        Path(path).read_text(encoding="utf-8"),
        object_pairs_hook=unique_claim_fields,
        parse_float=JsonNumber,
        parse_int=JsonNumber,
        parse_constant=numeric_value,
    )
    required = {"gene", "sample", "value", "scale"}
    if not isinstance(claim, dict) or not required.issubset(claim):
        raise ValueError("Claim file requires gene, sample, value and scale fields.")
    if set(claim) - required - {"decimal_places"}:
        raise ValueError("Claim file contains an unknown field.")
    if any(not isinstance(claim[key], str) or not claim[key].strip() for key in ("gene", "sample", "scale")):
        raise ValueError("Gene, sample and scale must be nonempty strings.")
    if not isinstance(claim["value"], str) or not claim["value"].strip():
        raise ValueError("Value must be a nonempty numeric string or a JSON number.")
    places = claim.get("decimal_places")
    if isinstance(places, JsonNumber) and places.lstrip("-").isdigit():
        places = int(places)
        claim["decimal_places"] = places
    if places is not None and (type(places) is not int or places < 0):
        raise ValueError("Decimal places must be a nonnegative integer or null.")
    return claim


class ClaimArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        raise ValueError(message)


class SingleValue(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if getattr(namespace, self.dest, None) is not None:
            parser.error("Option repeated: {}".format(option_string))
        setattr(namespace, self.dest, values)


def main():
    parser = ClaimArgumentParser(description=__doc__, allow_abbrev=False)
    parser.add_argument("--claim", action=SingleValue, help="Read gene, sample, value and scale from a JSON claim file.")
    parser.add_argument("--gene", action=SingleValue, help="Source name with surrounding whitespace trimmed, or H2BC12 with documented aliases.")
    parser.add_argument("--sample", action=SingleValue)
    parser.add_argument("--value", action=SingleValue, help="Claimed number on the vendor-transformed scale.")
    parser.add_argument(
        "--decimal-places",
        action=SingleValue,
        type=int,
        help="Explicit rounding precision. Default is exact comparison. Halfway values round up.",
    )
    parser.add_argument(
        "--practice",
        action="store_true",
        help="Use the bundled fictional blank/zero practice table instead of patient data.",
    )
    parser.add_argument(
        "--scale",
        action=SingleValue,
        help="Claimed measurement scale. Default: vendor_transformed. Other scales remain unresolved.",
    )
    try:
        args = parser.parse_args()
        if args.claim:
            overrides = (args.gene, args.sample, args.value, args.scale, args.decimal_places)
            if any(value is not None for value in overrides):
                raise ValueError("Use a claim file or individual claim options, without mixing them.")
            claim = read_claim(args.claim)
            report = check(
                claim["gene"], claim["sample"], claim["value"], claim["scale"],
                args.practice, claim.get("decimal_places"),
            )
        else:
            if args.gene is None or args.sample is None or args.value is None:
                raise ValueError("Supply --claim or all of --gene, --sample and --value.")
            scale = args.scale if args.scale is not None else SOURCE_SCALE
            report = check(args.gene, args.sample, args.value, scale, args.practice, args.decimal_places)
    except (OSError, UnicodeError, ValueError, InvalidOperation, RecursionError) as error:
        print(json.dumps({"status": "input_error", "reason": str(error)}, indent=2))
        return 2
    print(json.dumps(report, indent=2))
    return 0 if report["status"] in {"value_matches", "value_matches_after_rounding"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
