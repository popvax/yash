# Add a research appendix on differential agretopicity in Design 4

Adds `designs/design-4/appendix-dai-analysis/` within Design 4. It proposes no selection or sequence change.

The appendix reports **3,024 mutation-spanning missense peptide–HLA comparisons**, using DAI = predicted WT nM / predicted mutant nM on the same allele. Ten of the 14 assessable missense regions have a higher-DAI pair than their existing lead: eight involve another peptide and two another allele. It also records 18 pairs, representing 17 peptides, with DAI ≤1 and predicted mutant affinity ≤500 nM. Changing a lead label does not remove these sequences from an encoded window.

A separate comparison scores both members of the 23 TIER1 homolog pairs on their recorded shared allele. The original flags remain valid under their three-branch predicate. MHCflurry and MixMHCpred give opposite relative ordering for two of 21 complete pairs. No clinical consequence is established.

The [findings](FINDINGS.md) link each result to its source tables. The [limits](LIMITS.md) distinguish arithmetic and model reproduction from biological evidence. The [correction record](CORRECTIONS.md) explains why the earlier redesign was withdrawn.

```sh
python tools/verify.py --design4 ..
python tools/verify.py --design4 .. --recompute --models-dir /path/to/affinity_predictor
```

Ordinary verification checks five named Design 4 inputs, derives table values and counts, and reruns stored-source corroboration. The optional mode repeats MHCflurry inference with the recorded model assets. Passing these checks does not establish tumour presentation, T-cell recognition, safety or benefit.

Produced by an AI agent under non-expert human guidance. This is research analysis, not medical advice or a recommendation to manufacture or administer a product. It establishes no safety, immunogenicity or clinical benefit.
