# Matched-normal DNA and Design #2 data

Yash's normal-cell reference is genomic DNA from white blood cells, sample
`TN25-112103-N`. The final aligned reads and reusable analysis outputs are in the
same R2 bucket as the original tumour/RNA inputs. Small candidate and evidence
tables remain in [Design #2](../../designs/design-2/README.md), the lead design.

The original normal R1/R2 FASTQs are excluded from this distribution. The final
BAM retains all 1,212,214,833 alignment records from the normal alignment; duplicate
reads are marked rather than removed. The earlier pre-calibration and sort BAMs
are processing intermediates and are not duplicated here.

## Download individual files

Sizes are decimal GB/MB. The [manifest](manifest.json) records every object's exact
size, source path, SHA-256 and S3 location. The
[downloader](../../tools/r2/README.md) supports credential-free public reads
and authorized private S3 administration; see the
[current access configuration](../catalog/access.json).
The bucket stays private and original normal FASTQs remain excluded. General
reference assemblies, model weights and software are obtained directly from
[their providers](../../docs/external-resources.md), not this data endpoint.

| File | Size | Contents |
|---|---:|---|
| Final normal BAM | 143.29 GB | Coordinate-sorted, duplicate-marked, BQSR-applied normal reads. |
| Normal BAM index | 9.06 MB | Random access to the BAM. |
| Normal coding GVCF | 30.30 MB | Native normal genotypes and reference-confidence blocks in the assessed coding/flanking intervals. |
| GVCF index | 556 kB | Random access to the GVCF. |
| All paired alleles | 9.95 MB | Normalized paired tumour–normal calls, including rejected alleles and effective per-allele filter flags. |
| All-allele index | 715 kB | Random access to all paired calls. |
| Paired PASS alleles | 357 kB | The 5,210 native-PASS alleles; downstream candidate eligibility applies additional rules. |
| PASS index | 48 kB | Random access to PASS calls. |
| Native filtered VCF | 11.13 MB | Native GATK filtering output before normalized unique-allele consolidation. |
| Read-quality report | 496 kB | QC-only original-read assessment; trimming and filtering were disabled. |
| Normal callability summary | 456 B | Exact assessed scope and callable-base counts. |

The catalog also includes called-normal alleles, bounded personal-normal sequence
models and unresolved cases, RNA reconstruction, rankings, source manifests,
all-option/final-junction predictions, final sequences and native completion
receipts. Original caller filters and unsuccessful candidate outcomes are retained.
[Transfer verification](verification.json) records full-object R2 GET checksums.
The [download check](download-verification.json) records a successful restore
through the existing downloader, indexed normal-GVCF access, paired sample labels
and the 5,210-record paired PASS count.

## Restore with checksum verification

Run from the repository root. The default is a plan; add `--execute` to download:

```sh
python tools/r2/restore.py --destination ../yash-data \
  --key raw/TN25-112103/normal/normal.recal.bam \
  --key raw/TN25-112103/normal/normal.recal.bam.bai
```

For just normal genotypes:

```sh
python tools/r2/restore.py --destination ../yash-data \
  --key analysis/TN25-112103/design-2/normal-coding/gather/normal.coding.g.vcf.gz \
  --key analysis/TN25-112103/design-2/normal-coding/gather/normal.coding.g.vcf.gz.tbi
```

`--category raw` now includes the large normal BAM. `--category analysis` includes
all Design #2 output objects. Inspect the printed byte total before executing.
Restoration preserves original `/data`-relative paths, so mount the destination
as `/data` when reusing recorded processing commands.

## Scope and reference

The BAM is the full normal alignment. **Normal GVCF genotyping covers primary
nuclear annotated CDS/start/stop intervals plus 100-bp flanks**; it is not a
whole-genome genotype callset. Within that 76,021,670-base scope, 97.23% meets
native GQ≥20 and depth≥20. Target-specific normal read/genotype evidence remains
separate from this overall coverage summary.

The analysis reference FASTA, fingerprint map and BQSR known-sites resources
are external prerequisites. Obtain references from [their providers](../../docs/external-resources.md)
and match the recorded identities; the project does not publicly mirror them.

The analysis reference SHA-256 is
`c3e61ad6d770db22a6f27fdbcd3dba93d0d80f649018c8887ec803f12ea90450`.
The original Caris custom masks are unavailable; this is the compatible analysis
reference described in [Design #2 methods](../../designs/design-2/methods.md).
The normal sequencing instrument is unconfirmed and recorded as unknown.
The paired callset uses the original tumour BAM with this matched-normal BAM;
its native filter PASS and the subsequent research-selection policy are distinct.
