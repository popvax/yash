# Selected-data endpoint

**https://yash-data.samarth-3f0.workers.dev/INDEX.json** lists the public data.
[access.json](../../../data/catalog/access.json) records the endpoint;
[public-access.json](../../../data/catalog/public-access.json) records the
anonymous verification results.

The Worker serves exactly the 145 patient-data and case-analysis objects in the
[manifest](../../../data/catalog/public-object-manifest.json), plus an index and
manifest. It supports streaming, GET/HEAD and single byte ranges. The bucket
remains private. Unlisted keys, writes, original normal FASTQs and third-party
software/model/reference mirrors are rejected. The legacy endpoint is disabled.

## Download and verify

Use the [download guide](../README.md). No credentials are required, and interrupted
downloads resume with full SHA-256 verification before restoration.

```sh
python tools/r2/public/verify_access.py --output work/public-access-check.json
```

This checks every allowlisted HEAD response, denied keys, three small complete
files and three normal-BAM ranges. It does not re-download all large files.

## Maintain the endpoint

Regenerate the local allowlist with `python tools/r2/public/build_catalog.py --execute`.
Only the selected patient/results manifest is used; generation publishes nothing.

`python tools/r2/public/launch.py` prints a plan. With `--publish`, it deploys
`wrangler.launch.jsonc`, verifies anonymous HTTPS access and records availability.
Failed verification triggers deployment of the disabled fallback. Preview URLs
remain disabled. Deployment changes require an explicit user instruction.

To withdraw public downloads without deleting private data:

```sh
npx wrangler@4.157.0 deploy --config tools/r2/public/wrangler.jsonc
```

## Private administration

Ordinary downloads need no credentials. Authorized maintainers can install
`python -m pip install '.[storage]'` and keep private bucket credentials outside
Git, in `~/.config/yash-r2.env` with mode `0600`, or environment variables.
The expected names are `R2_ACCESS_KEY_ID` and `R2_SECRET_ACCESS_KEY`;
`--credentials-file` selects another private settings file. Never commit them.

`launch_outputs_upload.py` plans a private Modal-to-R2 transfer and requires
`--execute` to run. `finish_outputs_upload.py` collects receipts;
`publish_outputs_catalog.py` updates private metadata after full-object checks.
Local transfer records go under `provenance/data-catalog/transfers/`; the original
inventories remain intact. Public dataset selection is maintained separately in
`data/files.tsv`. No transfer or deployment is performed by catalog generation.
