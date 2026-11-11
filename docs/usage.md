# Use the data and build another approach

Design #2 is the matched-normal lead; Design #1 was generated without matched-normal
DNA and remains available for comparison. Independent methods and formats are welcome.

## Verify the lead

With Python 3.11 or 3.12, install from the repository root:

```sh
python -m pip install -e '.[mrna]'
yash-vaccine verify
```

This read-only command checks recorded evidence, native completion receipts,
source snapshots, selected antigens, protein and GenBank annotations. It starts
no jobs and downloads nothing. `make verify` uses the same lead-design check.

## Replay the lead locally

```sh
yash-vaccine reproduce --output work/design-2-replay
```

This reselects regions from the included candidate pool, re-solves the recorded
joint order/linker cost graph, independently recounts chosen junction predictions
and reruns synonymous encoding. Selected TSV/FASTA, protein/components/order and
nucleotide outputs must match byte for byte. It uses recorded model measurements;
it does not repeat variant calling, raw-read processing or HLA model inference.

Omit `--output` to use a temporary directory, removed after verification. Existing
output directories are never overwritten. `make reproduce` also targets the lead.

```sh
yash-vaccine reproduce designs/design-2 --output work/explicit-lead
python designs/design-2/run.py --output work/scripted-lead
```

## Earlier design and comparisons

```sh
yash-vaccine reproduce designs/initial --output work/design-1-replay
yash-vaccine compare designs/initial designs/design-2 --output work/design-comparison
```

The first replays Design #1's full 29-region/852-aa protein and known nucleotide
segments. The second compares the two complete proteins by `region_id`, reporting
22 retained, two added and seven removed regions. `AG01` labels are local to each
selection. DNA/self evidence and prediction queries differ between these designs;
the comparison is not a controlled measurement of vaccine performance.

Both designs leave PolyA and downstream extensions unspecified. No complete
molecular RNA/template length is assigned.

## Start independent work

```sh
yash-vaccine init my-candidate
python designs/my-candidate/run.py
```

`init` creates a README and starter script, not an analysis or vaccine. Replace
the script with your work, or document another language/workflow. You can create
the directory manually instead. No design.json, quota, HLA class, window length
or construct format is required for an independent contribution.

Describe your question, methods, inputs, dependencies, commands, results and
unfinished work in that directory. Keep useful plain-text, TSV, FASTA or other
outputs alongside it. A design can be an antigen list or analysis without a cassette.

To change the lead Design #2, start your script from its
[candidate pool](../designs/design-2/exports/integrated-candidate-pool.tsv),
[selected regions](../designs/design-2/selection/selected.tsv) and
[methods](../designs/design-2/methods.md). The recorded lead replay checks a fixed
design; changing its inputs is not a supported way to create another design.
Document any new predictions or assembly checks your changes require.

The optional [Design #1 recipe helper](recipe-runner.md) varies the earlier
tumour-only selection and internal cassette. `init --from` does not copy the
matched-normal lead or add its SEC/MITD and RNA stages.

## Download source data

```sh
yash-vaccine data list
python tools/r2/restore.py --destination ../yash-data --category raw
```

The second command prints an offline plan. Add `--execute` to download, or use
`--key` for individual files. Public downloads need no credentials. The endpoint and availability are
recorded in data/catalog/access.json. The raw category includes the large
matched-normal BAM; inspect the byte total. Original normal FASTQs are excluded.

[Data](../data/README.md), [normal-data scope](../data/matched-normal/README.md)
and [download instructions](../tools/r2/README.md) explain files, hashes and access.
New discovery/model runs have separate [runtime contracts](upstream.md).
Their third-party dependencies are obtained directly from
[upstream providers](external-resources.md), not this project's data endpoint.
Verification reports three upstream source files as hash references when they
have not been obtained and byte-checked locally; offline replay does not run them.

## Contribute

Fork the repository, create a branch and add your work under `designs/`. Include
antigens/rationale, sequences where produced, supporting results, code, dependencies
and reproduction commands. Label unspecified sequence components explicitly.
When comparing with our approach, use the lead Design #2 unless studying Design #1.
See [CONTRIBUTING](../CONTRIBUTING.md).
