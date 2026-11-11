# Hydrophobicity compatibility screen

The initial design rejects a nominated antigen region if **any complete 19-aa
window has a mean Kyte–Doolittle hydropathy strictly greater than 1.6**. All residues
have equal weight, no terminal padding is used, and threshold decisions use exact
integer tenths before any display rounding. Unknown residues and sequences shorter
than 19 aa stop this helper rather than receiving a false PASS.

This is a conservative design choice based on the membrane-like-segment heuristic
in [Kyte and Doolittle (1982)](https://doi.org/10.1016/0022-2836(82)90515-0).
The [ExPASy scale table](https://web.expasy.org/protscale/pscale/Hydropath.Doolittle.html)
pins the residue values. The cutoff is **not** a calibrated membrane-insertion
probability or a validated vaccine-selection threshold: hydrophobic soluble
segments and less hydrophobic membrane spans exist. Passing does not establish
solubility, trafficking, processing or presentation.

## Where the rule applies

1. Screen all **223 source-nominated 25-aa regions**, before source quotas. Keep
   the input pool and DNA/RNA evidence unchanged. Exclude one PTCH1 region;
   the other 222 pass. The selected set is **29 regions**: 17 mutation, 5 splice,
   5 CTA/self and 2 fusion.
2. Screen each candidate internal antigen–linker–antigen sequence before the
   order/linker solver. Incompatible joins are infeasible, rather than penalized
   against stronger HLA binding. The optional cassette helper also screens its
   initiating methionine.
3. Independently rescan the entire assembled protein. Require zero high-scoring
   windows outside the fixed terminal SEC/MITD elements. Export every unmasked
   score. Distinguish windows wholly inside an intentional routing element from
   windows crossing a routing boundary. Boundary overlaps are retained for review;
   their topology is not validated by this screen. Each antigen must separately
   pass, so overlapping a routing element cannot hide a hydrophobic antigen.

The selected-region gate operates **after the source adapters nominate their
windows**. It does not automatically search for lower-ranked replacement windows.
For the excluded PTCH1 V1083M event, all **25** recorded 25-aa windows also exceed
the threshold (minimum maximum-KD19 1.747; chosen window 2.121), so shifting within
that recorded context does not produce a passing alternative.

This is independent of somatic filtering. PTCH1 remains a DNA/RNA-supported,
tumour-only mutation hypothesis in the shared measurements. The hydropathy rule
does not establish whether that allele is inherited.

## Code and outputs

- [Maintained implementation](../../../../src/yash_vaccine/hydrophobicity.py)
- [All nominated regions](../../selection/hydrophobicity-screen.tsv) and
  [selection summary](../../selection/hydrophobicity-summary.json)
- [Full-protein window scores](../../construct/protein/hydrophobicity-windows.tsv)
  and [full-protein summary](../../construct/protein/hydrophobicity.json)
- [Alternate PTCH1 windows](excluded-mutation-windows.tsv)
- [Independent scale/calculation verification](verification.json)

Reproduce the complete candidate with `yash-vaccine reproduce`. To independently
recount all 1,561 nominated-region windows and 834 full-protein windows using Biopython's scale table and
`ProteinAnalysis.gravy`, and inspect recorded alternatives:

```sh
python designs/initial/review/hydrophobicity/run.py --output work/hydropathy-review
```

This uses included evidence only, starts no cloud jobs and preserves the candidate.
