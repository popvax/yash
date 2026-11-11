# Recorded execution source

recorded-source-2026-11-09.tar.gz preserves 187 project-source/documentation files
needed for the recorded handoff. The companion JSON binds every member to its
original size/SHA-256 and records the archive's current hash.

Three full upstream RNA modules are **not bundled**. Their original target
hashes remain under `external_files`; provider pointers, original hashes and
minimal reconstruction patches are in
[upstream-sources.json](../../workflow/stages/02_lens/upstream-sources.json).
[External-resource instructions](../../docs/external-resources.md) explain how
to obtain and prepare them directly from the providers under their terms.
The pre-exclusion archive hash is retained as provenance, not a download link.

This snapshot is separate from maintained src/ and current docs. Historical
scientific measurements and sequences retain their original checksums.
`RecordedSource` verifies every bundled member, validates the external identity
records against native ledgers, and reports unbundled source separately. It does
not claim to have byte-verified provider code that is not present.

The lead replay uses the recorded project optimizer and edge costs, then current
assembly/verification code, to reconstruct the design. It does not run the three
upstream RNA modules or repeat model inference. A raw-source rerun must obtain
those upstream dependencies; an explicit read of missing source fails with
instructions. Applicable notices accompany the minimal patches and historical
reference identities in [THIRD_PARTY_NOTICES.md](../../THIRD_PARTY_NOTICES.md).
