# Download data

Start with [the data guide](../../data/README.md) for direct links, sizes and
descriptions. Downloads need no account. Reports and HLA summaries are already
in the repository.

## Choose a dataset

From an installed checkout:

```sh
yash-vaccine data list
yash-vaccine data list --dataset tumour-rna
yash-vaccine data restore --dataset tumour-rna --destination ../yash-data
```

Restore prints the selected files' total size and starts no download. Add
`--execute` to download. Dataset selectors can be repeated. Bundled datasets such as `hla` and
`clinical-reports` are copied from the checkout into the destination without a
network request; their repository-relative paths and checksums are retained.

For a small example requiring only Python and this checkout:

```sh
python tools/r2/restore.py --dataset expression --destination ../yash-data --execute
```

Select individual file IDs from [data/files.tsv](../../data/files.tsv), for example:

```sh
python tools/r2/restore.py --file paired-variants-pass --file paired-variants-pass-index \
  --destination ../yash-data --execute
```

`--key` still accepts original storage keys. Legacy `--category` selections are
supported; `raw` includes the 143-GB normal BAM. Use dataset or file selectors
when you want a smaller download.

## Verification and resume

Each downloaded file is checked against its complete SHA-256 before restoration.
Rerun the command to resume a saved partial download. The existing prefix is
rehashed, the remaining byte range is requested and the final checksum must match.
Completed cache files are verified and reused. A conflicting existing file stops
restoration; a corrupt partial file is not installed.

Files retain their original `/data`-relative paths so recorded upstream commands
can use a restored directory mounted at `/data`. Cache hard links avoid a second
copy of large files. No Modal account is needed to download or replay a design.

`--anonymous` forces public HTTPS when private credentials are also present.
`--public-base-url` overrides the endpoint; HTTP is accepted only on loopback
for local checks. Neither option publishes data.

## Catalog maintenance

[data/files.tsv](../../data/files.tsv) supplies file and dataset selection.
`python tools/catalog.py` checks local hashes, published identities and generated
views; `--write` updates local README tables and download JSON. The published
object keys and bytes are locked by their original release records.

[Maintainer tools](public/README.md) cover private transfer, credentials, endpoint
verification and deployment. These are separate from downloading patient data.
