# Download metadata

Start with [the data guide](../README.md) and [file catalog](../files.tsv).

- `access.json`: current download endpoint.
- `public-object-manifest.json`: generated download metadata for the public files.
- `public-access.json`: recorded endpoint verification.
- `external-resources.json`: identities and provider links for upstream resources.

`python tools/catalog.py --write` regenerates the local manifest and marked
README tables from the file catalog, preserving published file identities.
[Original inventories and transfer records](../../provenance/data-catalog/README.md)
are separate from the current dataset menu.
