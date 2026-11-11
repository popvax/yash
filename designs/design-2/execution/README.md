# Native execution evidence

These are historical records of the completed matched-normal analysis, not live
job status. Original commands, worker exits, bytes, hashes and remote paths are
preserved. Current results are indexed in [the design guide](../README.md).

- `stage-receipts.json` and stage directories: original completion receipts.
- `artifacts.json`: original remote/local artifact names and checksums.
- `native-execution.json`: observed successful worker exits.
- `completion-audit.json`: completion-time audit; fields such as `pushed` refer to that run, not today’s GitHub state.
- `frozen-inputs.json`: the original reused-evidence/source manifest.
- `original-scientific-files.json`: original Design #1/shared-file ledger.
- `local-paths.json`: maps old local names to canonical files after consolidation. It does not rewrite native remote receipts.
- `source-adaptations/`: records how the executed branch adapted earlier workers.

Duplicate export copies of receipts and selected results have been removed.
Selection is canonical in `../selection/`; final junction verification is
`../construct/final-check.json`; native completion receipts are canonical here.

[Recorded-source snapshots](../../../provenance/snapshots/README.md) preserve old
code/docs. Current helpers can be edited and tested without pretending that the
edited code was used for the original experiment. Scientific measurement and
sequence hashes are still checked against their canonical files.
