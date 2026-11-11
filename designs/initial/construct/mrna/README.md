# RNA encoding for Design #1

PopRNA encodes the verified **852-aa protein** using human synonymous codons,
with seed **0**, a **0.10** codon-frequency cutoff and ViennaRNA **2.7.2** folding
at **37 °C**. The 29 antigen regions, their order, linkers and patient-supported
SEC/MITD amino-acid sequences remain fixed.

The same 5′ leader/Kozak, three-stop block and 3′ UTR are used as in the RNA
Foundry backbone. **PolyA and any downstream extension are unspecified.** The
FASTA and GenBank files contain known upstream segments, without assigning
complete RNA/template lengths.

## Method and checks

The runnable [optimizer.py](optimizer.py) adapts the user-supplied
`PopRNAOptimizer_final.py`; its supplied-source SHA-256 and exact settings are
recorded in [backbone.json](backbone.json). It uses sequential motif repair,
first-codon selection, alternative-frame stops, paired-bracket-run reduction,
codon-distribution balancing and repeat removal. Its implementation corrects
input/output handling, substring and boundary checks, finite synonymous repair,
rare-codon exclusion, paired/unpaired codon classification and repeat bounds.
The pinned ViennaRNA Python API replaces the fragile RNAfold subprocess parser.

[run.py](run.py) executes that optimizer, independently checks exact protein and
component translations, preserves the known fixed context and exports the RNA,
DNA and annotations. The human table in [human-codon-usage.json](human-codon-usage.json)
is pinned to python-codon-tables **0.1.18** and checked against the installed data.
The coding sequence excludes `UUUU`, `UUUC`, `UCUAGA`, `GGUCUC`, `GAGACC`,
`AAAAAAAA`, `GGGGGGGG` and `CCCCCCCC`. These are coding-region constraints;
fixed UTRs are not modified to satisfy them.

The final coding sequence has **56.2989% GC**, standard **CAI 0.77396802**, no
excluded rare codons, no forbidden coding motifs, no repeated 12-nt pairs and a
longest homopolymer of **7 nt**. Alternative frame +1 has 41 stop codons and a
74-codon longest stop-free stretch; frame +2 has 63 stops and a 54-codon stretch.
These are computed sequence properties, not expression measurements.

Folding is a bounded heuristic. Its final paired-bracket run is **28 nt**, so the
14-nt target is **not retained after the later stages**. The bracket-run metric
is not a rigorous continuous-stem or ensemble-accessibility measure. Mandatory
translation, fixed-context, motif and codon-cutoff checks are separate from this
heuristic target. No global optimum, measured expression or stability result is
claimed.

## Outputs

- [cds.fasta](cds.fasta): **2,559-nt** coding DNA, including its first stop.
- [mrna-segments.fasta](mrna-segments.fasta): the **2,699-nt** known RNA segment before PolyA.
- [template-segments.fasta](template-segments.fasta): the **4,640-bp** known upstream DNA context.
- [template-segments.gb](template-segments.gb): checked annotations of that partial DNA record.
- [components.tsv](components.tsv): identities, known sequences and segment-relative coordinates.
- [rna-layout.txt](rna-layout.txt): component order and illustrative `AAAAA…AAAAA` PolyA sketch.
- [construct.json](construct.json) and [verification.json](verification.json): settings, metrics, checks and hashes.
- [optimizer.log](optimizer.log): deterministic stage output.
- [execution/modal-state.json](execution/modal-state.json): the recorded Modal execution and its persistent-volume receipt.

The PolyA sketch assigns neither sequence nor length. FASTA records contain only
known nucleotides. PolyA coordinates, sequence and length are blank, rather than
zero. No downstream suffix or template endpoint is inferred.

## GenBank annotations

[genbank.py](genbank.py) derives features from the selected-region records,
patient-supported routing elements and pinned [backbone sources](annotation-sources.json).
It annotates the full vaccine CDS, T7 promoter, UTRs, SEC/MITD, all 29 regions,
numbered joins, stops and sequence-checked backbone features. Stable region and
event IDs, source flags and mutation-status notes are retained. Splice notes
retain the GENCODE 37 supporting transcripts.

The record is explicitly partial. PolyA is marked as an unknown-length boundary;
no feature may span its unspecified sequence. Complete circular topology and
experimental cleavage boundaries are not asserted. The export is parsed back to
check sequence, feature locations, identifiers, translations and RNA segments.

## Reproduce

Use Python 3.11 or 3.12 and `python -m pip install -e '.[mrna]'` from the repo root:

```sh
python designs/initial/construct/mrna/run.py --output work/replayed-mrna
```

`yash-vaccine reproduce` rebuilds all stages and checks the recorded outputs. It
starts no cloud jobs and downloads no reads or prediction models.

For an explicit Modal rerun, install the Modal SDK in the launcher environment:

```sh
python designs/initial/construct/mrna/modal_run.py --output work/mrna-modal
python designs/initial/construct/mrna/modal_run.py --execute --output work/mrna-modal
```

The first command only plans. Execution uses 2 CPUs, 8 GiB and no GPU, writes to a
new directory on the existing research Volume, collects outputs and checks that
the verification receipt persists after container termination. Choose new output
directories. Cap chemistry and nucleotide modifications remain unspecified.

Maintainers recording a newly assembled protein can explicitly add
`--record-new-encoding --execute`. The launcher first runs PopRNA on Modal, checks
translation, pins the resulting coding-sequence hash in backbone.json, and then
runs the strict checked exporter. Normal reproduction always requires that pinned
hash. [Recording tool](../../../../tools/record_initial.py) refreshes the complete
initial-design manifests only after all stages pass; it does not replace checks.
