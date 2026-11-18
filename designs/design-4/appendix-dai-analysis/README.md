# Differential agretopicity analysis of Design 4

This research appendix compares predicted Class-I binding of mutant peptides with their matched wild-type peptides. It proposes no change to Design 4's selection, protein sequence or nucleotide sequence. It is located inside `designs/design-4/appendix-dai-analysis/`; its files are included in the design's manifest.

For the overall proposal and reading order, start with the [Design 4 reading guide](../README.md#reading-guide). The separate 15-region construct belongs to Design 4; this appendix proposes neither construct.

**DAI = predicted wild-type nM ÷ predicted mutant nM**, on the same HLA allele. Lower nM means stronger predicted binding. DAI above 1 means a predicted increase in binding; 1 means no change; below 1 means a decrease. A higher DAI does not establish a better vaccine target.

## Findings

* Among the **14 missense regions** with applicable comparisons, **10** have another eligible peptide–HLA pair with higher DAI than the existing lead pair: **8** use a different peptide and **2** use the same peptide on another allele. The cassette contains 16 regions in total.
* **18 peptide–HLA pairs, representing 17 distinct peptides**, have DAI ≤1 and predicted mutant affinity ≤500 nM. Their sequences remain encoded regardless of which peptide is labelled the lead.
* Two selected 9-mers place a changed basic residue at the C terminus and have large predicted binding ratios on HLA-A*33:03. These predictions do not resolve the absence of coding-strand mutant RNA support for KCTD19 in the analysed specimen.
* Re-examining **23 TIER1 homolog pairs** preserves their original flags. MHCflurry and MixMHCpred give opposite relative ordering for **2 of 21** pairs with complete values from both models. No clinical consequence follows from that comparison.

See [FINDINGS.md](FINDINGS.md) for the figures and source links, and [LIMITS.md](LIMITS.md) for their interpretation. These calculations establish neither tumour presentation nor T-cell recognition, safety or clinical benefit.

## Data and scope

| File | Contents |
|---|---|
| `data/encoded-peptide-dai.tsv.gz` | 3,024 single-amino-acid WT-reversion peptide–HLA comparisons |
| `data/region-dai-summary.tsv` | Existing lead pair, maximum eligible DAI pair and window length |
| `data/residual-low-dai-peptides.tsv` | The 18 pairs meeting the descriptive DAI ≤1 / mutant ≤500 nM filters |
| `data/tier1-direction.tsv` | Both members of each TIER1 pair scored on its recorded shared allele |
| `data/cell-corroboration.json` | Stored-source corroboration of the 46 paired affinity cells |
| `data/receipt.json` | DAI-table reconstruction receipt; its no-inference flag applies to `tools/build_tables.py` |

The DAI table covers mutation-spanning missense comparisons, not every peptide in the cassette. Unchanged flanking peptides have identical WT counterparts. TGFBR2 frameshift and CDH1::NRG1 fusion regions have no applicable single-amino-acid reversion in this analysis and are marked `dai_defined=False`.

## Reproduce

```sh
python tools/verify.py --design4 ..
python tools/verify.py --design4 .. --recompute --models-dir /path/to/affinity_predictor
```

Use Python 3.11 or later. `tools/build_tables.py` uses the standard library. Fresh inference requires `mhcflurry==2.3.14` and the recorded model assets. Ordinary verification checks the five named Design 4 inputs, derives the published numerical findings and reruns the stored-source lookup. `--recompute` additionally checks local model-asset digests and repeats the affinity and percentile predictions. Agreement verifies computational reproduction, not biological accuracy. Files and analyses outside that scope are not authenticated by this verifier.

The model component is `models/affinity_predictor` from the official [OpenVax 20260928 presentation-model archive](https://github.com/openvax/mhcflurry/releases/download/2.3.0/models_class1_presentation.20260928.tar.bz2). Its individual hashes and aggregate identity are recorded in `data/tier1-direction-summary.json`; model files are not included here.

To regenerate the tables, work in a separate copy of this appendix:

```sh
python tools/build_tables.py --design4 .. --out data
python tools/tier1_direction.py --design4 .. --out data --models-dir /path/to/affinity_predictor
python tools/corroborate_cells.py --design4 .. --table data/tier1-direction.tsv --out data/cell-corroboration.json
python tools/verify.py --design4 .. --recompute --models-dir /path/to/affinity_predictor
```

Fresh model outputs can differ in their last displayed digits across numerical environments. Replay allows 0.05 nM for affinity and 0.001 percentile points for percentile values; it does not require byte-identical inference output. Regeneration updates the table digest and corroboration receipt together.

The [source provenance](PROVENANCE.md) distinguishes upstream patient evidence from derived research tables. The [correction record](CORRECTIONS.md) describes the withdrawn redesign and subsequent changes.

Produced by an AI agent under non-expert human guidance. This is research analysis, not medical advice or a recommendation to manufacture or administer a product. It establishes no safety, immunogenicity or clinical benefit.
