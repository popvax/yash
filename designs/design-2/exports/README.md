# Lead-design evidence and replay inputs

These tables are recorded inputs/evidence, rather than another selected design.
The canonical outputs are [selection/](../selection/) and [construct/](../construct/).

- `variant-evidence.*`: paired DNA, population and focal normal evidence.
- `rna-*`: RNA reconstruction/support, contexts and unsuccessful outcomes.
- `rank-*`: recorded mutation affinity, window scores and per-allele peptides. Weighted ranks are descriptive; all eligible mutation regions were retained.
- `integrated-candidate-pool.tsv`: source-specific eligible regions before final selection. Offline replay applies the recorded gates and source counts here.
- `integrated-*`: RNA-source rescoring, personal-self matches, selection decisions and comparison with Design #1.

Original stage completion receipts live once under [execution/](../execution/README.md).
Repeated `integrated-selected` copies were removed; use
[the canonical selected table](../selection/selected.tsv).
