# Contribute a vaccine design for Yash

Have a vaccine design for Yash, or an analysis that could improve one? Submit a
pull request. You can use our data and code, bring your own methods, or combine
them. There is no required pipeline, antigen quota, HLA class or construct format.

## Start a design

1. Fork this repository and create a branch for your work.
2. Add a directory under `designs/`, such as `designs/your-design/`. You can
   create it yourself, or run `yash-vaccine init your-design` for a README and
   starter script. See [usage](docs/usage.md) for the optional helpers.
3. Keep your notes, code and useful results together in that directory.
4. Open a pull request from your branch to this repository's `main` branch.

Use [Design #2](designs/design-2/README.md) as the baseline unless you are studying
the earlier design. To change its antigen selection, start with its
[candidate pool](designs/design-2/exports/integrated-candidate-pool.tsv),
[selected regions](designs/design-2/selection/selected.tsv) and
[methods](designs/design-2/methods.md), then implement your changes in your own
script. `yash-vaccine init your-design` creates that workspace without copying
either design's analysis. The optional `--from` recipe helper uses **Design #1's
tumour-only evidence and internal cassette**, not the matched-normal lead.

## What to include

- **Design and rationale:** the question, your selected antigens and why you
  chose them, and what differs from the lead Design #2. Explain the selection and scoring
  methods you actually used.
- **Sequences and results:** include the antigen table, predictions or other
  results supporting your choices. Include protein and mRNA sequence files if
  produced. Label unspecified components explicitly; do not fill them with
  invented bases or present sequence sketches as complete molecules.
- **Reproduction instructions:** list inputs and their sources, tool/model
  versions, dependencies and commands. Include your script or workflow. For
  additional large inputs, provide a location and checksum instead of adding
  raw reads or model weights to Git.
- **Interpretation:** describe the results and the scope of the evidence.
  Identify remaining questions, unfinished work and any changes to assumptions
  such as HLA typing or expression weighting. Distinguish reused predictions,
  new inference and experimental findings.

Your README should let another researcher or their agent understand what you
did and run it. Simple text, TSV and FASTA are welcome; use other formats when
they suit the work. A useful analysis or antigen list can be submitted before
you have assembled a complete construct.

## Describe the pull request

Summarise your approach, link to the design directory, give the command to run,
and identify the main result files. Explain shared-code changes separately.
When reusing our antigens, retain stable `region_id` values so comparisons are
clear; `AG01`-style identifiers belong to one selection.

Keep both recorded designs and original shared evidence as references. Put
alternative results in your own directory and temporary files in ignored
`work/`. Corrections to the shared evidence or implementation are welcome;
explain the issue and what you changed.

Ideas without a design or code are welcome too: open an issue or email
yash@rnafoundry.com.
