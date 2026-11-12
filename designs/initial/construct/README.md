# Joint protein assembly and RNA encoding

The fixed 29 antigen regions are assembled between patient-supported SEC and
MITD elements. Order and linkers are jointly optimized across **30 joins**.
The result is one **852-aa protein**: 26 SEC + 725 antigen + 46 linker + 55 MITD.
PolyA remains unspecified. [Full methods](../methods.md#7-add-patient-supported-routing-sequences-and-optimize-all-joins)
explain the objective and limitations.

## Search and result

There are 11 linker choices, **8,932 internal pair/linker options** and **638 terminal
options**: 9,570 candidate edges. A separate [hydropathy gate](../review/hydrophobicity/README.md)
rejects high-scoring internal joins before the solver. No remaining pair/linker
option fails in this candidate set, but the constraint remains enforced.

The verified internal-cassette cache contains 135,344 peptides against six HLA-I
alleles; the 29-region helper requires 127,115 of them. These are cached predictions,
not newly inferred pairs for the smaller selection. Fresh terminal inference
covers **14,040 peptides × six = 84,240 pairs**. The pinned MHCflurry affinity model
is used; processing/presentation is not scored.

OR-Tools CP-SAT minimizes rank ≤0.5% junction hits, then rank ≤2% hits,
then maximizes the worst junction rank, then minimizes added residues. All four
phases are proven optimal in the fixed feasible search space. The result has 0 hits
at rank ≤0.5%, 33 at rank ≤2% and 38 at affinity ≤500 nM.
A full-protein scan is checked against **7,224 fresh selected-junction predictions**.
Position repeats and allele rows are counted explicitly.

The unmasked hydropathy scan covers **834 windows**. All **753 windows wholly
outside routing elements pass**, with maximum KD19 **1.000**. Ten high-scoring windows
are wholly in MITD, and six overlap its boundary. Those are reported separately;
the heuristic does not establish their topology or prove proper trafficking.

## Files

| File | Contents |
|---|---|
| [architecture.json](architecture.json) | Patient HLA-derived routing sequences and source evidence |
| [protein/cassette.fasta](protein/cassette.fasta) | 852-aa protein |
| [protein/antigen-order.tsv](protein/antigen-order.tsv) | Order and stable region IDs |
| [protein/components.tsv](protein/components.tsv) | Component coordinates |
| [protein/junctions.tsv](protein/junctions.tsv) | All 30 joins |
| [protein/construct.json](protein/construct.json) | Objective, solver proof, hashes and model provenance |
| [protein/edge-costs.tsv](protein/edge-costs.tsv) | All 9,570 edge costs |
| [terminal-inference/](terminal-inference/README.md) | Terminal requests and fresh predictions |
| [selected-inference/](selected-inference) | Independent selected-junction predictions |
| [final-check.json](final-check.json) | Full-chain checks |
| [protein/hydrophobicity-windows.tsv](protein/hydrophobicity-windows.tsv) | Every unmasked 19-aa score and routing-overlap classification |
| [protein/hydrophobicity.json](protein/hydrophobicity.json) | Compatibility gate and summary |
| [mrna/](mrna/README.md) | Nucleotide segments and annotations |

## Reproduce

Install `python -m pip install -e '.[mrna]'` from the repository root.
`yash-vaccine reproduce` checks every stage without cloud jobs. Individual stages:

```sh
python designs/initial/construct/run.py --output work/replayed-protein
python designs/initial/construct/verify.py --protein work/replayed-protein --output work/replayed-check.json
python designs/initial/construct/mrna/run.py --output work/replayed-mrna
```

The standalone nucleotide builder reads this design's pinned protein and final check.
Full-design reproduction chains the rebuilt artifacts and compares all hashes.

## Run new model inference

The explicit launchers require authenticated Modal access to the recorded CUDA
image and reference volume, or matching replacements. Defaults use an A10,
four CPUs and 16 GiB. Missing predictions stop assembly; nothing is silently dropped.
Commands without `--execute` print a plan.

```sh
python tools/launch_predictions.py --output work/internal-predictions --execute
python designs/initial/construct/infer.py --output work/terminal-predictions --execute
python designs/initial/construct/run.py --terminal-predictions work/terminal-predictions --output work/new-protein
python designs/initial/construct/infer.py --request work/new-protein/verification-request.json --output work/selected-predictions --execute
python designs/initial/construct/verify.py --protein work/new-protein --inference work/selected-predictions --output work/new-check.json
```

The full routed runner reads the included internal cache. New antigen sets need
an explicitly updated cache, request manifests and documented design. These
commands are case-specific, rather than an unrestricted constructor.
