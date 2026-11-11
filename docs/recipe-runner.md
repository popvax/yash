# Optional: vary Design #1's selection recipe

**This helper uses Design #1's tumour-only evidence, not the matched-normal lead
Design #2.** It varies the earlier selection and internal cassette, before
SEC/MITD and RNA encoding. It does not adapt Design #2 when you change a quota.

To develop an alternative to the lead, create a plain workspace with
`yash-vaccine init my-candidate` and use
[Design #2's methods and evidence](../designs/design-2/README.md) in your own
script. Independent methods and formats are welcome; [usage](usage.md) explains
that route. Continue below only to vary the earlier Design #1 recipe.

Run commands from the repository root, with Python 3.11 or 3.12 and the environment
installed as shown in [the README](../README.md). `python -m yash_vaccine` is an
alternative to the installed `yash-vaccine` command.

## Reproduce Design #1

```sh
yash-vaccine reproduce designs/initial
```

This explicit Design #1 command rebuilds selection, joint SEC/MITD assembly, junction checks
and RNA encoding in a temporary directory. Install with
`python -m pip install -e '.[mrna]'`. It compares the final 852-aa protein and
known RNA/DNA segments before unspecified PolyA with the recorded outputs.
All solver phases must be optimal. The temporary files are removed afterward.
To keep them in a new directory:

```sh
yash-vaccine reproduce designs/initial --output work/reproduced-initial
```

Local reproduction uses the included evidence and predictions. It does not
rerun the tumour reads, call new mutations or run new HLA model inference.

## Copy the Design #1 recipe

```sh
yash-vaccine init fewer-mutations --from designs/initial/selection
```

This copies the selection/internal-cassette stage into
designs/fewer-mutations/design.json and README.md. It produces a cassette before
SEC/MITD and RNA encoding. Changing the antigen set needs matching predictions
and checks before using the later stages; the fixed final-design scripts are not
a generic runner for changed selections. Describe your
question in the README. In design.json, change `selection.quotas.MUTATION` from
17 to 16, then run:

```sh
yash-vaccine run designs/fewer-mutations/design.json
yash-vaccine compare designs/initial/selection designs/fewer-mutations --output work/fewer-mutations-comparison
```

The example selects 28 regions and builds an internal cassette using the included
prediction cache. The comparison explains which region was removed and how the
cassette and boundary scores changed. Files remain in the named design directory;
temporary upstream calculations are discarded after successful export.

The runner refuses to replace existing results. To try another change, create
another named design, or use --output with a new work/ directory. To check an
already completed candidate:

```sh
yash-vaccine reproduce designs/fewer-mutations
```

For the recorded stage alone:

```sh
yash-vaccine reproduce designs/initial/selection
```

The compatibility shorthand `--from designs/initial` also resolves to this
selection recipe. Comparisons against `designs/initial` use its intermediate
recipe-stage cassette and explicitly label that scope.

## What to edit

| Recipe field | Meaning |
|---|---|
| name, description | Identify the candidate and the question you are asking |
| selection.quotas | Mutation, splice, CTA/self and fusion region counts; 0 omits a source |
| selection.splice_frame_policy | Requires patient_RNA_annotated_frame, including full-window RNA support and the stop-free upstream path |
| selection.splice_one_region_per_gene | Whether splice selection is restricted to distinct genes |
| selection.hydrophobicity_policy | `kyte_doolittle_19_gt_1_6`; region and internal-join compatibility gate |
| linkers.choices | Named linker sequences; an empty string is a direct join |
| optimization.seconds_per_phase | Time limit for each of the four optimization phases |

Paths in a recipe resolve relative to design.json. Evidence and cache references
include SHA-256 hashes. `init` rebases paths automatically so they still point to
the shared files. The exported recipe records any explicitly supplied new cache.

The evidence reference identifies the recorded sequence/score bundle. The current
local adapters consume this case bundle; changing its path, HLA types, expression
model or window length is not an alternative selection setting. Such changes need
new upstream evidence and a deliberately versioned adapter/bundle. See
[upstream stages](upstream.md) and
[Design #1 methods](../designs/initial/methods.md).

## Read the outputs

| File | What it tells you |
|---|---|
| README.md | Your question and notes, or a short generated summary |
| design.json | Complete selection/linker recipe and shared evidence references |
| selected.tsv, selected.fasta | Selected antigen regions and sequences |
| candidate-pool.tsv | Pool after annotated-frame filtering; no sequence-identity filter |
| selection-audit.tsv, splice-frame-summary.json | Why regions were selected or excluded |
| cassette.fasta | Assembled protein cassette |
| components.tsv, junctions.tsv | Antigen/linker coordinates and chosen joins |
| boundary-predictions.tsv | Peptide/HLA predictions at the chosen boundaries |
| construct.json | Optimization objectives, linker composition and solver results |
| verification.json | Whole-design checks and artifact hashes |
| run.json | Code hashes, revision and local environment |

FASTA, TSV and TXT outputs are plain text. JSON is used for structured recipes,
checks and provenance. Coordinate intervals are zero-based, end-exclusive.

The detailed construct-verification.json separates `required_peptides` and
`required_peptide_hla_pairs` for this search from `cache_peptides` and
`cache_peptide_hla_pairs` available in storage. The same counts appear under
`predictions` in whole-design verification.json. Required queries cover potential
joins/linkers, start-M peptides and internal antigen peptides; they are not a
count of antigens in the final cassette.

## Compare candidates

```sh
yash-vaccine compare designs/initial/selection designs/fewer-mutations
```

Comparisons match selected.tsv rows by `region_id`. The same MUC3A region moves
from AG17 to AG16 in the example because one mutation was removed. AG01-style
`antigen_id` labels are local to a selection; they connect that selection to its
construct coordinates and must not be used to join separate designs.

--output saves a readable comparison.txt plus regions.tsv and comparison.json.
The comparison records changes to evidence, case context, prediction cache,
linker choices and optimization strategy so metrics can be interpreted in context.
It does not rank experimental or clinical vaccine effectiveness.

## New linkers or missing predictions

A new design can reuse the existing cache only when every required peptide and
patient HLA allele is covered. Missing predictions stop construction before
result export. To inspect the required prediction inputs:

```sh
yash-vaccine prepare designs/fewer-mutations/design.json --output work/my-queries
yash-vaccine predict work/my-queries --output work/my-predictions
```

Preparation writes peptide/HLA inputs, the selected regions and linker choices.
Prediction defaults to a plan and starts no job. To execute, install the optional
Modal dependency and authenticate your own account first:

```sh
python -m pip install -e '.[modal]'
yash-vaccine predict work/my-queries --output work/my-predictions --execute
yash-vaccine run designs/fewer-mutations/design.json --prediction-cache work/my-predictions/construct-predictions.npz --output work/new-prediction-design
```

The recorded A10 runtime requires the model files and saved volume described in
the upstream guide. A new output directory and explicit --execute are required.
The shared frozen cache is never replaced.

## Get source data

```sh
yash-vaccine data list
yash-vaccine data restore --category raw --destination ../yash-data
```

Listing and planning work offline. Add --execute to restore and --key for an exact
object. Downloads require authorized private S3 credentials; the anonymous endpoint
is disabled. See [data](../data/README.md) and
[download instructions](../tools/r2/README.md). These commands do not change access.

## Recompute upstream evidence

Inspect the recorded-input plans without starting jobs or creating results:

```sh
yash-vaccine upstream dna --output work/dna-run
yash-vaccine upstream rna --output work/rna-run
```

For the independent tumour-only discovery and segmented-filtering inputs, use [the upstream guide](upstream.md).

## Share a design

Create a branch, give the design a descriptive name, and include its directory in
a pull request. Keep your reasoning in README.md: what changed, why, the commands
you ran, and what you learned. The recipe and result files travel together.
Provide access instructions and hashes for any extra evidence/predictions; keep
large working files out of Git. Different antigen-selection or scoring approaches
are welcome; label their scope and keep Design #1 as a reference.
