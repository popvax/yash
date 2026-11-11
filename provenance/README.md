# Provenance records

Evidence lives with the data or result it describes:

- [Lead Design #2](../designs/design-2/README.md) records matched-normal evidence, selection and the final protein/nucleotide segments. [Native receipts](../designs/design-2/execution/README.md) preserve execution-time facts; [source snapshots](snapshots/README.md) keep historical code/docs separate from maintained code.

- [Design #1](../designs/initial/README.md) records selection, protein assembly,
  junction checks and RNA encoding. Its [manifest](../designs/initial/candidate.json)
  pins the files needed to reproduce the complete result.
- [Shared evidence manifest](../data/checkpoints/manifest.json) pins the portable
  analysis checkpoints.
- [Patient-RNA splice evidence](../data/splice-evidence/README.md) records all
  normal-filtered events, observed sequence/path/frame witnesses, native QC,
  complete reference matches, fresh affinities and independent BAM checks.
  [Offline verification](../designs/initial/selection/splice-rna-verification.json)
  recounts all scored windows. [ND-only normal screening](../data/splice-evidence/normal-screen/README.md) supplies the final splice admission and gene ranking.
- [Model and runtime records](../data/checkpoints/provenance/) identify the
  prediction environment and reference files.
- [Terminal inference](../designs/initial/construct/terminal-inference/README.md)
  and [selected-junction inference](../designs/initial/construct/selected-inference/)
  retain the executed requests, predictions and model checks.
- [Data catalog](../data/catalog/README.md) provides access instructions and
  checksums for the sequencing data and case-analysis files; general references
  and models come from their upstream providers.

Independent designs should document their inputs, methods and verification in
their own directories.
