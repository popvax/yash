# Yash's HLA-derived routing sequences

The SEC and MITD **amino-acid sequences** were selected from Yash's called
HLA-C*01:02 and HLA-B*15:08 alleles and checked against his tumour DNA and RNA
reads. The construct's nucleotide encoding is subsequently codon-optimised;
it does not copy his native codons.

| Element | Native protein interval | Sequence support |
|---|---|---|
| SEC | HLA-C*01:02 residues 1–26 | DNA and RNA consensus agree at all 78 bases; every codon has 87–533 DNA and 11–86 RNA supporting fragments. |
| MITD | HLA-B*15:08 residues 308–362 | DNA and RNA consensus agree at all 165 bases; every codon has 11–84 RNA supporting fragments. |

The SEC sequence is `MRVMAPRTLILLLSGALALTETWACS`. Its first 24 residues are the
annotated signal peptide; CS are the two following native HLA-C residues.
This is a 26-aa protein segment, rather than an experimentally assigned cleavage
boundary for the vaccine construct.

SEC position 25 is cysteine: 87 DNA fragments support C with no other amino acid
among qualifying full-codon observations. RNA has 67 C, two G and two F
observations. An independent exact-match search found five RNA fragments with
the complete 78-nt sequence at base quality ≥30. The proposed sequence has
support across multiple distinct fragment-coordinate patterns.

The 55-aa MITD matches the HLA-B*15:08 C-terminal segment. Positions 329 (T)
and 349 (S) distinguish it from Yash's other called HLA-B*08:01 reference:
DNA supports them with 288 and 39 fragments, and RNA with 74 and 23.
Eight RNA fragments contain an exact high-quality 63-nt stretch linking both
positions. No fragment covers the complete 165-nt MITD stretch at once;
verification uses overlapping fragments and the linked distinguishing residues.

HLA-B codon 338 crosses an intron. No qualifying DNA fragment observes all three
bases together, although each base has at least 171 supporting observations.
RNA directly supports the spliced glycine codon with 67 fragments. This is not
a missing amino acid or an alternative call. Minor discordant observations at
other positions are retained in the tables.

## Files

- [typing-comparison.json](typing-comparison.json): independent DNA/RNA T1K calls.
- [read-evidence.json](read-evidence.json): thresholds, input hashes and summaries.
- [codon-evidence.tsv](codon-evidence.tsv): every codon's observations, including
  alternative amino acids, stricter mapping counts and distinct coordinates.
- [nucleotide-evidence.tsv](nucleotide-evidence.tsv): every base's observations.
- [RNA reference details](rna-target-reference-details.json) and
  [DNA reference details](dna-target-reference-details.json): sequences and
  allele-reference coding coordinates.
- [audit_reads.py](audit_reads.py): analysis code for the retained extracted HLA reads.
- [alignment-tool-versions.txt](alignment-tool-versions.txt): BWA/samtools versions.

The allele references are from the official IPD–IMGT/HLA database; protein
comparisons use commit `d972ae287948dc2e5468a40ff3b3d27227541189`:
[HLA-C](https://github.com/ANHIG/IMGTHLA/blob/d972ae287948dc2e5468a40ff3b3d27227541189/fasta/C_prot.fasta),
[HLA-B](https://github.com/ANHIG/IMGTHLA/blob/d972ae287948dc2e5468a40ff3b3d27227541189/fasta/B_prot.fasta).
The DNA caller's C*01:02:01:23 representative and RNA caller's C*01:02:01:01
representative have the same coding sequence. The audit does not establish a
clinically validated four-field genotype.

## How the read check works

The audit remaps 23,858 extracted RNA fragments and 150,223 DNA fragments against
63 competing called HLA/related allele references. It requires base quality ≥30
and at least one mate with mapping quality ≥20. Mapped mates must agree on the
reference and face opposite directions. Overlapping mates count once; conflicting
high-quality calls are excluded. All three codon bases must be observed to count
an amino acid. Requiring the supporting read itself to have mapping quality ≥20
still gives at least 11 RNA fragments for every target codon.

These are fragment counts, not UMI-confirmed molecule counts. Exact substring
checks provide additional alignment-independent evidence; uniqueness is checked
within the called-reference panel, rather than every possible human transcript.
No matched-normal sample was used. This establishes agreement with the sequenced
tumour sample; cleavage, trafficking and presentation remain experimental questions.

For offline inspection, the tables and reference details are self-contained.
To repeat the raw-read audit, use the T1K 9f09f4dccee4e7501b98a4e822943c0fb26a4e25
HLA extraction products retained on the `research-lens-workspace` Modal volume,
under `/data/results/TN25-112103/hla-class2-20261028/`, then run:

```sh
python audit_reads.py
```

Dependencies: BWA 0.7.17, samtools 1.16.1, Biopython 1.86 and pysam 0.23.3.
The script checks CDS translation, remaps the paired reads and exports all
evidence tables under `/data/results/TN25-112103/trafficking-domain-audit-20261102/`.
It starts no cloud instance or downloads implicitly. Original sequencing files
are described in [data access](../../../../data/README.md).

The local design rebuild checks the architecture against both typing and read
evidence before assembly; it does not rerun raw-read analysis.
