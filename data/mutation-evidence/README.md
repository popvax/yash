# Segmented tumour-only mutation evidence

This evidence supports independent tumour DNA discovery, segmented native
FilterMutectCalls filtering, Isovar RNA reconstruction and class-I affinity scoring.

- [Policy](policy.json): native filtering and additional research gates.
- [DNA measurements](annotate/variant-evidence.tsv): 1,060 protein-changing native-PASS variants.
- [Eligibility decisions](eligibility/decisions.tsv) and [verification](eligibility/verification.json):
  the research policy leaves 105 primary DNA candidates.
- [RNA outcomes](rna/rna-evidence.tsv): one outcome for each primary DNA candidate.
- [Protein contexts](rna/primary-proteins.tsv): 23 usable RNA-supported contexts.
- [Ranked regions](rank/segmented-primary-ranked-regions.tsv): 18 qualifying mutation regions.
- [All scored windows](rank/all-scored-windows.tsv), [affinities](rank/peptide-HLA-affinity.tsv)
  and [reference matches](rank/reference-matches.json): independent score-recount inputs.
- [Native commands and statistics](filter/segmented/) and [source locations](source.json):
  provenance for the native calls and exported measurements.

Native GATK mapping/base-quality filtering remains in effect. The additional
research policy uses depth, alternate support, variant fraction, GERMQ, population
frequency and context-review rules. Ancestry-frequency exclusion has no minimum
sample-size requirement. No matched-normal somatic confirmation is available.

From the repository root, reproduce the DNA decisions with:

```sh
python -m yash_vaccine.mutation_policy --output work/mutation-eligibility
```

The complete `yash-vaccine reproduce` command also verifies RNA-supported window
scores and the final selected design. These commands use recorded evidence;
they do not repeat native calling, RNA reconstruction or model inference.

To export a completed segmented source run into this portable evidence layout:

```sh
python tools/export_mutation_evidence.py --source /path/to/completed-run --output work/mutation-evidence
```

The export selects native-PASS coding calls and the primary DNA pool's measured
RNA outcomes and predictions. It does not generate new measurements. Raw-read
execution contracts are in [upstream documentation](../../docs/upstream.md).
