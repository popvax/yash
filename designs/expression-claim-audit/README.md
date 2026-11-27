# Expression-claim audit

Check whether a structured numerical expression claim agrees with the source table for
the stated gene and sample. This independent computational analysis accompanies issue #10.
It uses the existing Design 3 context and expression records. A match establishes source
agreement under the declared comparison policy.

A gene is an instruction in a cell, like a recipe. RNA is a working copy. Expression
describes the amount of those copies. The source contains vendor-transformed values,
which this tool compares on their supplied scale.

## Inputs and dependencies

Python 3.9 or newer and its standard library are sufficient for these independent scripts.
Local execution used Python 3.9.13. The maintained workflow has separate Python 3.11/3.12
requirements. Those interpreters were unavailable for local execution of this contribution.

- [Expression source](../../data/expression/vendor-transformed-expression.tsv), identified
  by exact bytes and sample header in [source-identity.json](source-identity.json).
- [Existing H2BC12 identifier audit](../design-2/review/H2BC12-expression-identifier-audit.json).
- [Design 3 methods](../design-3/methods.md) document the missing-expression exception.
- `fixtures/`, `claims/` and `blank-zero-practice.tsv` are newly authored fictional examples.
  The alias fixture tests the code path using invented values and sample identifiers.

The expression file is already bundled. The commands below require no downloads or cloud
account. The checker performs no prediction inference or training.

## Reproduce

From the repository root, compare a fictional rounded claim:

```text
python designs/expression-claim-audit/fixture_check.py --claim designs/expression-claim-audit/claims/rounded-number.json
```

Replay all 36 fictional examples, including relevant reported-field checks:

```text
python designs/expression-claim-audit/run.py --output work/expression-claim-audit-replay.json
```

Also exercise six read-only checks against the existing expression table:

```text
python designs/expression-claim-audit/run.py --with-source --output work/expression-claim-audit-source-replay.json
```

Choose a fresh output filename on each run. Source integration constructs claims from
the pinned table and records statuses only. It checks file access and comparison behavior.
Its expected answers are developer-authored. The H2BC12 claim of 6 is deliberately invented.
The integration replay verifies the expression file's bytes, while the earlier identifier
audit is referenced and its full byte record is not verified by this command.

To check a new source claim, supply `--gene`, `--sample`, `--value` and optionally `--scale`
to `source_check.py`, or use `--claim` with JSON fields `gene`, `sample`, `value`, `scale`
and optional `decimal_places`. Saved claims must declare the scale. Individual options
default to `vendor_transformed`. The fixture launcher selects only declared fictional profiles.

## Rules and results

| Outcome | Meaning | Single-claim exit code |
|---|---|---:|
| `value_matches` | Exact numerical agreement | 0 |
| `value_matches_after_rounding` | Agreement after explicitly requested HALF_UP rounding | 0 |
| `value_differs` | Comparable source value disagrees | 1 |
| `gene_or_sample_differs` | The sample differs from the source header | 1 |
| `measurement_missing_from_source` | A matching field is blank or the bounded H2BC12 lookup finds no value | 1 |
| `insufficient_information` | Unresolved name, scale or conflicting evidence | 1 |
| `input_error` | Malformed claim or source | 2 |

`run.py` reports whether demonstration outcomes match expectations. Its exit 0 means the
examples behaved as specified, including examples expected to reject a claim.

The prepared local run observed 36/36 expected synthetic outcomes and 6/6 source-integration
outcomes. See [validation-summary.json](validation-summary.json) for the recorded code and
input identities. Full receipts are written to the chosen `work/` path on replay.

Rules include blank versus numeric zero, exact decimals, optional HALF_UP decimal-place
rounding and a 1000-digit intermediate working-precision limit. Gene/sample/scale padding
is trimmed with originals retained. Source identifiers must be nonempty and unpadded,
headers must precede data, and malformed or negative values invalidate the whole table.
The format is UTF-8 without a BOM, using physical LF or CRLF lines.

JSON number tokens retain their written form. `parsed_claim_value` exposes their Decimal
meaning. Decimal text syntax includes equivalent Unicode digits, underscores and signed
zero. Original source cells and parsed/rounded comparison values remain distinguishable.
Identifier case and internal characters are significant. Unsupported names stay unresolved.

## Review status and limitations

AI assisted with implementation, case authoring, review and documentation. Two external
AI reviews prompted fixes that were reproduced locally. The contribution is submitted
with developer checks and explicit limitations. The contributor reports completing their
own self-review. Independent human/domain validation and held-out evaluation have not been performed.

The examples are development material and reuse cases seen during implementation. Their
counts do not estimate real-world accuracy or maintainer time savings. Written explanations,
general gene-identifier resolution and unit conversions remain outside automatic approval.
HALF_UP is this tool's declared policy and needs review when another reporting convention applies.
The H2BC12 name mapping is reused from the recorded audit and does not establish biological
activity or measurement equivalence independently. Predictions are estimates when cited in
the source project; this checker makes no conclusions from prediction models.

Fingerprints identify exact bytes and detect accidental changes. They do not authenticate
mutable inputs against a malicious editor. A changed line ending can change a fingerprint.
Scoped attributes preserve this contribution's checker, fixture and recorded JSON bytes in Git.

## Contribution and sharing scope

The repository invites independent analyses under `designs/`. Existing expression data
and audit evidence are referenced in place. Published source tables, original measurements,
models, raw reads, review transcripts and private development files are outside this addition.
The supplied fixtures and claims are fictional. Source integration emits status summaries.

No general repository license was found. The contribution guide invites code/data use for
contributions, while the third-party notices explicitly leave licensing of repository code,
patient data and designs unassigned. This analysis does not assign licenses to those inputs.
Broader redistribution, commercial reuse or training rights remain NOT VERIFIED.
