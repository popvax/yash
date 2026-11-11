# Minimal changes to provider source

These zero-context patches retain only the changes required for the recorded
RNA analysis. Full LENSTools/MSNAF Python modules are not distributed here or in
the source snapshot. Original and target SHA-256 values, provider image names
and extraction paths are in ../upstream-sources.json.

The annotation patch changes header compatibility, the splice patch changes
patient-RNA review handling, and the MSNAF patch selects 8–11-aa peptide lengths.
The patches reconstruct the original executed bytes, without adding new source
headers to those historical outputs. See [third-party notices](../../../../THIRD_PARTY_NOTICES.md)
for the retained MSNAF MIT notice.

Obtain originals directly from the named providers under their terms, then run
../prepare_upstream.py as described in [external resources](../../../../docs/external-resources.md).
Nothing here downloads a model, accepts a license or starts compute.
