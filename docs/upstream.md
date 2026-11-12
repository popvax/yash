# Upstream discovery and native scoring

Obtain third-party tools, pretrained weights and general references directly
from their [providers](external-resources.md). Large patient inputs and case outputs
use the [documented download route](../tools/r2/README.md); current availability
is recorded in data/catalog/access.json. Full upstream RNA modules are not
bundled; minimal patches plus original/target hashes reconstruct provider source
locally. The offline design replay does not require these dependencies.

## Lead Design #2: paired tumour–normal workflow

The current lead uses paired Mutect2 discovery, measured normal genotypes/reads,
bounded personal-normal sequence models, tumour-RNA reconstruction and HLA-I
affinity screening. Start with [Design #2 methods](../designs/design-2/methods.md)
and its [controller/runtime guide](../designs/design-2/pipeline/README.md).
Those case-specific controllers retain the original Modal volume/image/run paths.
New raw-read processing and new inference require staged inputs and explicit
`--execute`; a fresh-account cloud deployment has not been verified.

For an offline replay of the completed lead, use `yash-vaccine reproduce`.
The `upstream dna/rna` CLI adapters below are explicitly **Design #1** helpers;
they do not silently launch the matched-normal branch.

## Design #1: recorded tumour-only workflow

Initial-design mutation selection uses independent tumour-only Mutect2 discovery,
segmented native filtering with current population evidence, Isovar protein
reconstruction and class-I affinity scoring. No Caris VCF nomination is used.
[Methods](../designs/initial/methods.md) state the exact settings and research gates.

## DNA and RNA

The [pipeline source](../workflow/stages/01_mutation_calling/discovery/README.md)
map each stage, its Modal runtime and its input contract. [Included evidence](../designs/initial/inputs/mutations/README.md)
contains decisions, RNA contexts, model predictions and all-window scores.
The [research policy](../designs/initial/inputs/mutations/policy.json) is applied
to those measurements by `src/yash_vaccine/mutation_policy.py`: native quality
filters remain, with no extra median MQ/BQ gate or ancestry sample-size minimum.
Execution receipts remain available alongside the measured evidence.
The offline [mutation replay](../src/yash_vaccine/mutations.py) recounts window
scores and regenerates the 18-row segmented checkpoint. All scores and original
DNA/RNA measurements remain available. The separate construct-compatibility
screen excludes PTCH1, leaving 17 mutation regions in the initial design.


Remote helpers require authenticated Modal, matching saved images and staged
BAMs/references in the recorded volume. Their case-specific controllers print
plans without `--execute`. They are distinct from the default local replay;
an anonymous fresh-account, raw-read-to-construct cloud run has not been verified.
Do not import upstream worker modules to inspect them: some run at import.

## LENS and integrated regions

The RNA-only LENS source stages supply splice, fusion and CTA/self evidence.
Their [canonical runner](../workflow/stages/02_lens/canonical_run.sh), expression
adapter, native review, exact-breakpoint fusion check and GENCODE frame
reconstruction are described in the initial-design methods. Source eligibility
and exact-self policy precede shared MHCflurry class-I scoring.

```sh
python workflow/stages/02_lens/reproduce.py --output-root /data/results/TN25-112103/lens-NEW
python workflow/stages/04_region_integration/reproduce.py --output-root /data/results/TN25-112103/replays/regions-NEW --mutation-regions /path/to/ranked-regions.tsv --lens-root /data/results/TN25-112103/lens-NEW
```

Inspect these plans before adding `--execute`. Restore the documented original
inputs first. Upstream results need their own provenance and checkpoints; cached
predictions cannot be relabelled as new inference.

## Patient-RNA splice reconstruction

Start with all 4,929 junctions passing both native MSNAF normal-control screens,
before reference translation or HLA filtering. Reconstruct actual patient RNA
query bases in verified GENCODE 37 upstream coding frames. Require a complete
75-nt haplotype in at least two qualifying fragment names, consistent overlapping
mates and a high-quality, stop-free observed path from the coding anchor to the
join. Do not fill missing bases from the reference or infer translation beyond
an observed stop codon.

Apply native LENS event/peptide review, compare join-overlapping peptides against
the full native Ensembl 103 reference inventory and run fresh six-allele
MHCflurry affinity predictions. Native Vaxrank affinity sums are multiplied by
the square root of supplied gene expression, with no splice fraction. Rank
windows, events and genes; the selection takes one region per splice gene.

[Stage instructions](../workflow/stages/04_region_integration/patient_splice/README.md)
contain explicit raw-read/GPU commands, runtime/input contracts and output names.
[Included read and prediction evidence](../data/shared-rna/splicing/README.md) supports
an offline audit of every window and score. Local replay checks the recorded
independent BAM verification; it does not reread the BAM or rerun the model.

```sh
python workflow/stages/04_region_integration/patient_splice/replay.py \
  --results data/shared-rna/splicing --output work/splice-verification.json
```

### Additional normal-reference and limited strand checks

The existing splice workflow also has an added admission stage: coordinate-aware
normal lookup across gene assignments and both strands, normal counts by source
tissue category, and read-pair orientation of complete-window witnesses. Motifs
are report-only. Current automatic selection requires ND in both searched normal
catalogs: any exact coordinate entry, including an opposite-strand or all-zero
entry, is held for review. The admission policy is recorded separately from the
unchanged count-measurement receipts. This stage adds no distinct-pattern or substitution filter and
reuses unchanged affinity/expression scores. Frozen Design #1 inputs stay intact.

See the [stage guide](../workflow/stages/04_region_integration/patient_splice/CHECKS.md)
and [completed analysis](../data/shared-rna/splicing/normal-screen/README.md). To replay it
offline from the included normal/strand measurements:

```sh
python workflow/stages/04_region_integration/patient_splice/checks.py apply \
  --source data/shared-rna/splicing --output data/shared-rna/splicing/normal-screen
python workflow/stages/04_region_integration/patient_splice/checks.py verify \
  --source data/shared-rna/splicing --output data/shared-rna/splicing/normal-screen
```

Mutations use gene expression times the local RNA alternate-fragment fraction;
splices use gene expression alone. Fusion weighting uses the documented
left-partner proxy. Missing values remain missing.

## Fresh construct inference

[Assembly](../designs/initial/construct/README.md#run-new-model-inference) describes
explicit GPU inference for junction options and independent selected-junction
verification. `yash-vaccine reproduce designs/initial` instead uses their saved verified predictions.
