# Original data inventories

The current dataset menu is [data/README.md](../../data/README.md), backed by
[data/files.tsv](../../data/files.tsv). These files retain earlier transfer and
inventory records with their original scope and bytes:

- [Original object inventory](original-object-manifest.json): 460 private archive
  objects, including references and software that are obtained upstream today.
- [Original input/source inventory](original-sources.json).
- [Matched-normal transfer](matched-normal-transfer.json),
  [transfer verification](matched-normal-transfer-verification.json) and
  [download verification](matched-normal-download-verification.json).
- [Published release identity](published-release.json): the metadata and key-set
  identity of the unchanged 145-object public release. Its native object records
  are retained in the original inventory.

`python tools/catalog.py` checks the current catalog and generated download
manifest against those published identities. `--write` updates local generated
tables and JSON; it does not publish anything or change access.
