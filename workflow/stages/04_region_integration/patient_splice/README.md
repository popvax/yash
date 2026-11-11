# Patient-RNA splice reconstruction

This stage reconstructs the patient sequence **before** exact-self exclusion,
binding prediction or ranking. It starts from normal-filtered RNA junctions,
rather than an HLA-filtered peptide report. Reference annotations supply the
reading frame; actual RNA query bases supply the sequence.

## What it does

1. Start with every junction passing both recorded MSNAF `maxmin` normal-control
   decisions: **4,929 UIDs**, including those whose reference translation did not
   produce a previously admitted peptide. Recover coordinates with the original
   MSNAF functions and cross-check all 3,656 recorded translated junctions.
2. Verify each usable GENCODE 37 transcript's reference CDS translation against
   its documented protein. The reference supplies a coding-frame anchor, rather
   than replacement bases for the patient.
3. Find primary, nonduplicate, QC-pass RNA alignments with MAPQ >=20 and NH=1
   when present. Read actual query bases across the exact CIGAR junction. Retain
   insertions; do not treat skipped introns as query bases or use soft clips.
   Anchor the frame in at least 10 upstream annotated CDS bases and propagate it
   along the query sequence. Complete observed codons between the upstream
   anchor and the join must be high-quality and stop-free. A candidate after an
   observed stop codon cannot qualify by preserving its modulo-three frame.
   A consistently overlapping mate may extend the
   sequence; at least 10 aligned overlap bases must establish a single offset.
4. Enumerate every fully observed 25-aa window crossing the join. All 75 bases
   must have quality >=20. At least two distinct fragment names must support
   the same complete nucleotide haplotype. Keep alternative supported haplotypes
   separately. Never fill gaps with reference bases or combine different
   fragments into a consensus. These are fragment names, not UMI molecule counts.
   High-quality overlapping mates must agree within the window and the observed
   upstream coding-anchor path. If the same
   nucleotide sequence has multiple alignment paths, record each path with its
   own supporting fragment names instead of assigning one representative path
   to every fragment.
5. Apply native LENS event metrics, flags, penalties and tier thresholds. Apply
   its peptide coding-kmer support and cysteine rules to the **patient** coding
   sequences. The reference-only peptide-consistency test is replaced by direct
   full-window RNA support and the annotated coding-frame anchor. Preserve the
   native treatment of splice motifs; no additional motif-direction filter is
   introduced.
6. Compare the actual join-overlapping 8–11-aa peptides against the full native
   Vaxrank Ensembl 103 annotated reference inventory. Revalidate all 104,806
   protein-source associations against the current GTF and protein FASTA. Check
   the multi-pattern search against exhaustive k-mer enumeration for each of
   its 89,878 protein sequences. This accelerates the same exact-substring
   policy without a biological filter change. Exact matches cannot contribute
   to target scoring. Run fresh MHCflurry 2.3.14 affinity predictions
   for Yash's six class-I alleles, using the hash-verified recorded models.
7. Use native Vaxrank affinity contributions. Admission requires a nonself
   peptide with native RNA review PASS and affinity <=500 nM; weaker overlapping
   peptides can still contribute under the native scoring defaults. Multiply the
   region score by `sqrt(supplied_transformed_gene_tpm)`. Use no splice fraction
   or raw-read abundance multiplier. Rank all windows, then the best window per
   event and the best event per gene. Select the top five positive-expression
   genes without an identity filter.
8. Independently reopen the BAM and verify **every recorded window and every
   named supporting fragment**, including the stop-free upstream coding path,
   using pysam's query/reference-position arrays,
   rather than the reconstruction worker's CIGAR path builder.

## Run it

The default commands print plans. `--execute` starts explicit cloud work or
downloads the completed results. Use an authenticated Modal account with the
recorded `research-lens-workspace` volume, the pinned image, BAM, reference and
LENS source files. Input hashes and inference receipts are recorded with results.
This is separate from the offline replay of the recorded initial design.

Obtain the complete GENCODE 37 `gencode.v37.pc_translations.fa.gz`, then use a new
result directory under `/data/results/TN25-112103/`:

```sh
python workflow/stages/04_region_integration/patient_splice/run.py reconstruct \
  --remote-root /data/results/TN25-112103/patient-splices-NEW \
  --output work/patient-splices-NEW \
  --shards 4 \
  --proteins /path/to/gencode.v37.pc_translations.fa.gz --execute
```

Check the specific sandbox handles with `status --execute`, using the same
`--remote-root` and `--output`. After all four reconstruction shards complete
successfully, merge them. The merge proves every input event was processed
exactly once. Then run each stage, waiting for completion before a dependent one:

```sh
python workflow/stages/04_region_integration/patient_splice/run.py merge --shards 4 \
  --remote-root /data/results/TN25-112103/patient-splices-NEW \
  --output work/patient-splices-NEW --execute
python workflow/stages/04_region_integration/patient_splice/run.py review \
  --remote-root /data/results/TN25-112103/patient-splices-NEW \
  --output work/patient-splices-NEW --execute
python workflow/stages/04_region_integration/patient_splice/run.py reference \
  --remote-root /data/results/TN25-112103/patient-splices-NEW \
  --output work/patient-splices-NEW --execute
python workflow/stages/04_region_integration/patient_splice/run.py predict-score \
  --remote-root /data/results/TN25-112103/patient-splices-NEW \
  --output work/patient-splices-NEW --execute
python workflow/stages/04_region_integration/patient_splice/run.py verify \
  --remote-root /data/results/TN25-112103/patient-splices-NEW \
  --output work/patient-splices-NEW --execute
python workflow/stages/04_region_integration/patient_splice/run.py collect \
  --remote-root /data/results/TN25-112103/patient-splices-NEW \
  --output work/patient-splices-NEW --execute
```

`verify` depends only on reconstruction, so it may run alongside later scoring.
Each launch records its sandbox ID and uploaded source hashes locally. The
same stage cannot be launched again while its sandbox is active. Completed
stages are collected rather than silently repeated. Use a new result directory
for reconstruction with a different policy; --scoring-name permits an explicit
new scoring attempt on unchanged read evidence.

After collecting a completed run, independently recount its scores and ranking
offline, using the installed repository dependencies:

```sh
python workflow/stages/04_region_integration/patient_splice/replay.py \
  --results work/patient-splices-NEW \
  --output work/patient-splices-NEW/offline-verification.json
```

This checks the input partition, nucleotide translations, per-path fragment
witnesses, stop-free upstream paths, independent BAM-check receipt, exact-self decisions, six-allele query
coverage, native logistic affinity contributions, expression multiplication and
gene selection. It reuses the saved model outputs and does not reread the BAM.

## Outputs

| Directory | Contents |
|---|---|
| reconstruction/ | Every input event and outcome, verified reference transcripts, patient haplotype windows with fragment names and frame anchors, native event evidence and source hashes |
| review/ | Native peptide flags and support for actual patient coding sequences |
| reference-check/ | Complete reference-inventory validation and two independent exact-match scans |
| scoring/ | Fresh affinities, exact-reference matches, every scored window, event/gene rankings, top-five selection and model/scoring receipts |
| read-verification/ | Independent BAM verification of all reconstructed windows |

A supported local RNA sequence does not establish a complete translated isoform,
tumour specificity or presentation. Total gene expression remains a ranking proxy
for the splice product. Normal-control screening and native read review retain
their existing limitations in repeats, strand interpretation and bulk RNA.
This worker assesses linear CIGAR splice gaps and contiguous intron boundaries;
ordinary intron-skipping alignments cannot establish circular/back-splice joins.
Non-linear or unresolved coordinates remain visible in the input/outcome audit.
