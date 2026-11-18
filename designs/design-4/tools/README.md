# Optional raw-read and predictor verification tools

## 1. Guardant360 ctDNA locus pileups and MHCflurry Class-I predictions (`audit_ctdna_loci.py`)

`audit_ctdna_loci.py` computes the 3-BAM fragment/strand pileups, paired Mutect2 VCF lookup and MHCflurry 2.3.14 (`2.3.0` models) 6-allele Class-I affinity predictions for all 13 nucleotide-resolved Guardant360 alterations (`exports/ctdna-locus-pileups.tsv`, `exports/ctdna-mhcflurry-predictions.tsv`, `exports/ctdna-audit-receipt.json`). Offline verification runs without BAMs or predictor weights; `--execute` recomputes both tables from local BAMs/VCF and MHCflurry models.

## 2. Optional RNA/normal concordance check (`verify_rna_normal_identity.py`)

`verify_rna_normal_identity.py` compares a local RNA BAM with the supplied normal coding gVCF at common-marker positions. It reports aggregate counts; detailed genotypes require the separate `--marker-output` option. It makes no network requests or uploads.

From the `tools/` directory, install `identity-requirements.txt` in a separate environment. Supply the RNA BAM and index, normal coding gVCF and GATK GetPileupSummaries table:

```sh
python verify_rna_normal_identity.py \
  --rna-bam /local/path/RNA.bam \
  --rna-index /local/path/RNA.bam.bai \
  --normal-gvcf /local/path/normal-coding.g.vcf \
  --normal-pileups /local/path/normal-pileups.table \
  --output /local/path/new-aggregate.json
```

The input identities and recorded result are in [identity-concordance.json](../review/identity-concordance.json). The table is a GATK population-marker input, not a Picard haplotype map. The checker selects normal heterozygous or homozygous-alternate SNPs, so its opposite-homozygote denominator covers normal 1/1 versus RNA 0/0 only. Threshold matches are not an identity probability, clinical authentication or exhaustive contamination assessment. The main offline replay does not rerun this optional raw-read check.
