# Data catalog and access

This repository includes the candidate tables, analysis code and portable
checkpoints. The retained catalog lists 145 patient-data and case-analysis objects.
Third-party tools, trained models and general references are obtained from their
providers, not mirrored by this project.

We are sharing these data to support the lead matched-normal Design #2, the
earlier tumour-only Design #1, and independent approaches to vaccine design for Yash. See the
[directions to explore](../../README.md#other-directions-to-explore).

The selected files are available without credentials through a public read-only
endpoint. The bucket itself remains private. Check [access.json](access.json) for
the URL and [public-access.json](public-access.json) for the verification receipt.
Object keys and hashes remain unchanged.

[File manifest](public-object-manifest.json) · [Download and restore instructions](../../tools/r2/README.md) · [Access configuration](access.json)

[Retained-source notices](../../THIRD_PARTY_NOTICES.md) ·
[Required credits](../../ACKNOWLEDGEMENTS.md)

The retained case catalog does not contain third-party code, pretrained weights,
training datasets, reference libraries or compiled containers. The
[provider guide](../../docs/external-resources.md) and
[external identities](external-resources.json) identify required upstream files.
The original [native inventory](object-manifest.json) is retained as provenance,
not as a promise that its withdrawn reference/software keys are public downloads.

## Download

[Commands and access configuration](../../tools/r2/README.md) cover individual files,
archive restoration and private credentials. A plan works offline:

```sh
python tools/r2/restore.py --destination ../yash-data --category raw
```

Add --execute to download, or use --key for an individual object below.

## Clinical reports

The [report directory](../reports/) contains the Stanford pathology/IHC
and STAMP, Caris and Guardant PDFs.
All four are included when you clone this repository and can also be downloaded
individually from GitHub.

## Raw sequencing and supplied values

| File | Size |
|---|---:|
| DNA_TN25-112103.bam | 13.26 GB |
| DNA_TN25-112103.bam.bai | 4.4 MB |
| DNA_TN25-112103.vcf | 487.5 kB |
| RNA_TN25-112103.bam | 1.32 GB |
| RNA_TN25-112103.bam.bai | 3.0 MB |
| RNA_TN25-112103_S45.R1.fastq.gz | 678.5 MB |
| RNA_TN25-112103_S45.R2.fastq.gz | 703.0 MB |
| vendor-transformed-expression.tsv | 461.1 kB |

## Matched-normal DNA and Design #2 outputs

The [normal-data guide](../matched-normal/README.md) lists the final **143.29-GB
normal BAM**, its index, normal coding GVCF, complete filtered/PASS paired calls,
QC and the supporting Design #2 outputs. [Its manifest](../matched-normal/manifest.json)
records all new objects and hashes. The original normal FASTQs are excluded.
The existing RNA FASTQs listed above are included in the selected download catalog.

Use `--category raw` for sequencing inputs (including the normal BAM) and
`--category analysis` for the new analysis outputs. Individual `--key` selections
avoid downloading a whole category.

## References and software

Use the [upstream acquisition guide](../../docs/external-resources.md). Original
versions, hashes and intended paths are in [external-resources.json](external-resources.json).
Only explicitly needed inputs are prerequisites for raw-source reruns; the
recorded local design replay needs no third-party reference or model download.

Object keys, sizes and hashes of the retained patient/results files are unchanged.
Private archival objects are not deleted or exposed. The default Worker profile
keeps downloads disabled; the separate launch profile enables only the selected
allowlist.
