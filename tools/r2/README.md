# Download and restore data

The selected release contains 145 patient-data and case-analysis objects in a
private R2 bucket. A public read-only endpoint makes those selected files
available without credentials. [access.json](../../data/catalog/access.json)
records the URL and current availability. When the public route is disabled,
maintainers can use authorized private S3 access.

Local design replay uses the included evidence and needs no downloads or Modal.
Third-party software, models and references come from their
[providers](../../docs/external-resources.md).

## Public-release downloads

From the repository root, inspect an offline plan:

```sh
python tools/r2/restore.py --destination ../yash-data --category raw
```

Select an individual file and add `--execute` to download:

```sh
python tools/r2/restore.py --destination ../yash-data \
  --key raw/TN25-112103/vendor-transformed-expression.tsv --execute
```

A fresh user's download needs no account or credentials. `--anonymous` forces
public HTTPS even if private S3 credentials are present. `--public-base-url` is
an explicit host override; plain HTTP is accepted only on loopback for private
verification. Neither option publishes anything.

The normal BAM is approximately 143 GB. Original normal FASTQs are excluded;
RNA FASTQs are included. Inspect the byte total before selecting an entire
category. See the [normal-data guide](../../data/matched-normal/README.md).

Downloads are checked against the full SHA-256 manifest before restoration.
Rerun the same command to resume a saved partial download: the saved prefix is
rehashed, the missing suffix is requested by byte range, and the complete hash
must match. A server that fails to confirm the requested range is rejected
before appending. Completed cache files are reused after verification. A corrupt
partial file is not installed; remove the named invalid cache file before retrying.

Files restore under their original `/data`-relative paths. A differing existing
file stops restoration. Hard links avoid duplicating large BAMs in the cache.

## Authorized private access

For private administration:

```sh
python -m pip install '.[storage]'
```

Obtain read-only bucket credentials from the owner. Keep them outside Git in
`~/.config/yash-r2.env` with mode `0600`, or use environment variables:

```text
R2_ACCESS_KEY_ID=...
R2_SECRET_ACCESS_KEY=...
```

The bucket and S3 endpoint default to access.json. Use
`--credentials-file /private/path/to/settings.env` for another credential file.
No credentials are required for the public-release route.

## Verify the release endpoint

This bounded anonymous check validates every allowed HEAD response,
rejects unlisted keys, downloads three small files with full SHA-256 checks and
checks three normal-BAM byte ranges. It does not re-download the large BAMs:

```sh
python tools/r2/public/verify_access.py --output work/public-access-check.json
```

For a full dataset stream/hash check, which reads approximately 160 GB:

```sh
python tools/r2/verify_remote.py --anonymous --output work/full-public-get.json
```

`--resume` can reuse completed full-GET receipts only when keys, sizes and hashes
match. Keep fresh checks distinct from reused receipts.

## Maintain transfers and publishing

`launch_outputs_upload.py --plan PLAN` plans a private Modal-to-R2 transfer;
`--execute` launches it. Credentials enter through an ephemeral Modal secret.
`finish_outputs_upload.py` collects receipts; `publish_outputs_catalog.py`
updates private metadata without exposing objects publicly.

The [release Worker](public/README.md) serves only the selected allowlist.
Its ordinary configuration is a disabled fallback; the public profile serves
the selected files. The legacy endpoint remains disabled.
