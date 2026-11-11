# Recorded context and regression definitions

These files support reproduction and the optional recipe runner. Independent
designs can use their own scripts and configuration; see docs/usage.md.

- case.json records the HLA-I types, window length and expression policies that
  underlie the included evidence.
- frozen-replay.json records the Design #1 annotated-frame recipe’s regression expectations.
- linkers.json is the original panel used by low-level prediction/runtime adapters.
  Named designs record their actual linker choices directly in their own recipe.

Changes to these context/regression definitions require rebuilt evidence and an
explicitly reviewed migration. Other designs can use their own methods, inputs and configuration.
