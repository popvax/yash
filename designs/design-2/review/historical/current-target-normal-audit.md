> **Historical checkpoint, superseded.** “Pending” below describes the early normal-only audit. See the [completed paired reassessment](../final/report.md) for current findings.

# Normal-DNA audit of the 17 recorded mutation targets

All 17 positions had usable normal evidence: complete native genotypes with GQ 99 and at least 20 qualifying fragment names. Seven target alleles are present in white-blood-cell DNA; ten have reference-only support in this audit. This is an independently nominated normal check, not the new paired-discovery result.

White-blood-cell presence is not automatically an inherited-disease finding. The full paired call set, broader normal genotype assessment, fresh RNA reconstruction and patient-specific self-sequence review remain required. The recorded initial design is unchanged.

| Gene | Recorded protein change | Normal ALT / REF fragments | Normal GT | GQ | Normal finding |
|---|---|---:|---|---:|---|
| KLF5 | D418N | 0 / 50 | 0/0 | 99 | Reference supported; paired result pending |
| SF1 | R218Q | 0 / 57 | 0/0 | 99 | Reference supported; paired result pending |
| KMT2D | G3189R | 26 / 35 | 0/1 | 99 | Present in normal |
| RXRA | A241S | 0 / 42 | 0/0 | 99 | Reference supported; paired result pending |
| HSPG2 | R3952L | 0 / 46 | 0/0 | 99 | Reference supported; paired result pending |
| PARP3 | K37Q | 23 / 29 | 0/1 | 99 | Present in normal |
| CAST | I309V | 22 / 15 | 0/1 | 99 | Present in normal |
| CCNI | P285S | 17 / 23 | 0/1 | 99 | Present in normal |
| PPIC | P187S | 0 / 47 | 0/0 | 99 | Reference supported; paired result pending |
| PPP2R1A | K280N | 26 / 22 | 0/1 | 99 | Present in normal |
| KCTD19 | G920R | 0 / 49 | 0/0 | 99 | Reference supported; paired result pending |
| VPS11 | D151A | 0 / 32 | 0/0 | 99 | Reference supported; paired result pending |
| ARHGAP40 | R438W | 0 / 49 | 0/0 | 99 | Reference supported; paired result pending |
| NRXN2 | D203H | 23 / 33 | 0/1 | 99 | Present in normal |
| PLEC | R2171W | 0 / 26 | 0/0 | 99 | Reference supported; paired result pending |
| DOCK8 | G1251R | 17 / 18 | 0/1 | 99 | Present in normal |
| SYT12 | P391L | 0 / 56 | 0/0 | 99 | Reference supported; paired result pending |

The read count uses primary, nonduplicate, vendor-QC-passing alignments with mapping quality ≥20 and focal base quality ≥20, after BQSR. Paired reads sharing a name count once. Disagreeing overlapping mates remain separate. Counts are supporting fragment names, not UMI-established independent molecules.

Native HaplotypeCaller genotyped the explicitly nominated alleles with ±100-bp intervals and `EMIT_ALL_CONFIDENT_SITES`. The final report also checks every original target against the broader coding-region GVCF and local read evidence, including targets absent from new Mutect2 output. Normal-reference support alone is not proof of somatic origin.

- [Measured normal evidence](current-mutation-audit.tsv)
- [Complete record and genomic alleles](current-mutation-audit.json)
- [Native completion receipt](current-target-normal-audit-completion.json)
