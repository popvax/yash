# Executed independent DNA/RNA pipeline sources

These workers implement independent tumour-only discovery, population
assessment, native segmented filtering, Isovar reconstruction and native Vaxrank
class-I window scoring. They are case-specific Modal workers; several execute
at import. Do not bulk-import them. Annotation workers apply the research
policy without extra median MQ/BQ cutoffs or a minimum ancestry-group allele count.
Recorded upstream receipts identify the code and policy executed for their data;
the original annotation policy is preserved with that receipt. Offline
eligibility is independently calculated from the recorded measurements by
`src/yash_vaccine/mutation_policy.py`.

| Location | Purpose |
|---|---|
| calling/preflight.py, prepare_plan.py, resplit_plan.py | Inspect inputs and prepare 143 chromosome-wide intervals |
| calling/population.py, call.py, gather.py, support.py | Population resource, independent Mutect2 discovery, sample merge and segmentation/orientation assets |
| review/setup.py, population.py | Map every emitted ALT and update global population priors from both gnomAD cohorts |
| review/filter.py | Native FilterMutectCalls with the segmented model |
| review/gather.py, annotate.py | Annotate protein-changing alleles and apply explicit DNA/population/context gates |
| review/rna.py | Native Isovar RNA/protein reconstruction |
| review/rank.py | Native affinity scoring, reference checks and all-window ranking |
| review/audit.py | Scientific consistency audit |

## Execution contract

`calling/run.py` and `review/run.py` are explicit Modal controllers. Without
`--execute` they print plans. `review/probe.py --execute` reads official source
HTTP metadata; it does not upload patient data. Controllers retain the recorded
remote result paths and image IDs. Scientific workers retain their executed source. Controllers resolve the current
checkout location; recorded remote paths and identifiers remain traceable.

```sh
python workflow/stages/01_mutation_calling/discovery/calling/run.py population --output work/discovery-calling
python workflow/stages/01_mutation_calling/discovery/review/run.py segmented --output work/discovery-review
```

Raw execution requires the same authenticated Modal volume, staged tumour DNA/RNA,
compatible GRCh38, gene annotation, population assets, completed calling outputs,
merged statistics, orientation model and matching saved images. No automatic
fresh-account/raw-read deployment is asserted. Adapt the pinned workspace roots
and rebuild matching runtimes deliberately for another environment. Some supporting
workers are invoked by their controllers within that environment, not as local CLI tools.

The portable `yash-vaccine reproduce` command verifies the initial design using
included evidence without these cloud prerequisites. Exact final-design inputs,
settings and scope are in [methods](../../../../designs/initial/methods.md) and
[data evidence](../../../../data/mutation-evidence/README.md).
