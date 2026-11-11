# Designs and independent analyses

| Design | DNA evidence | Regions | Protein | Status |
|---|---|---:|---:|---|
| [design-2/](design-2/README.md) | Paired tumour and matched-normal white-blood-cell DNA | 24 | 724 aa | Lead |
| [initial/](initial/README.md) | Tumour-only, without matched-normal DNA | 29 | 852 aa | Design #1, retained for comparison |

Both designs include every mutation region passing their recorded eligibility
and construct-compatibility gates. Each contains its own methods, selected
antigens, sequences, evidence and reproduction instructions. PolyA remains
unspecified in both; neither assigns a complete RNA/template length.

```sh
yash-vaccine verify
yash-vaccine reproduce --output work/design-2-replay
yash-vaccine compare designs/initial designs/design-2
yash-vaccine init my-candidate
```

The default verification/replay target is Design #2. Explicitly request
`designs/initial` to reproduce Design #1. Independent work can use any methods,
code and formats; no recipe, quota, HLA class or construct is prescribed.
Describe inputs, commands and findings in your design directory.

See [usage](../docs/usage.md), [contributing](../CONTRIBUTING.md) or the
[optional Design #1 recipe helper](../docs/recipe-runner.md).
