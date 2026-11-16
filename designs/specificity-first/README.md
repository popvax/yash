# Design "specificity-first" — independent candidate addressing issue #1

This design responds to [issue #1](https://github.com/popvax/yash/issues/1): four of the
lead Design #2 splice regions (F5, MUC5B, MUC3A, MUC6) owe their peptide novelty to
**extra substitutions inside the window**, not to the splice junction itself. Reverting a
single residue per region zeroes every nonself junction-overlapping 8–11-mer, so those
regions would stop being tumor-specific if the substituted junction peptide is not
presented.

This design keeps the same evidence base and pipeline conventions as Design #2 but makes
**splice-junction specificity the primary selection criterion**:

- **25 regions of 25 amino acids: 13 mutation, 5 splice, 5 CTA/self and 2 CDH1::NRG1
  fusion regions**, assembled into a **748-aa protein** (Design #2: 24 regions, 724 aa).
- All four contested splice windows are **replaced** with splice windows from the same
  recorded 80-gene pool whose junction-overlapping peptides remain nonself **from the
  junction alone** — no extra substitution is required (verified by reverting every
  detected extra substitution and re-testing all 8–11-mers).
- **H2BC12 E114K** is added as a 13th mutation region, labeled
  `conditional_missing_expression_annotation`: it was excluded from Design #2 only for a
  missing symbol in the expression table, while carrying paired somatic-DNA and
  tumor-RNA support and a predicted ~21× mutant/reference binding gain on A*33:03.
- Every region carries a per-allele dependence analysis so the cost of losing any single
  HLA allele is explicit.

Everything here is **computational prediction** (MHCflurry 2.3.14 pan-allele affinity
weights), not measured presentation or immunogenicity. No wet-lab claims are made.

## What changed relative to Design #2

See `comparison-vs-design-2/comparison.tsv` (regenerate with `compare_to_design2.py`):

- **20 regions shared** with Design #2 (12 of 12 mutations kept, 5 CTA/self, 2 fusion,
  and the HLA-B splice region, which passes the same attribution test — it was not
  evaluated in issue #1).
- **4 splice regions removed**: F5, MUC5B, MUC3A, MUC6 (substitution-dependent novelty).
- **5 regions added**: F5 (different window, junction-only novelty), HLA-B retained,
  HSF1, KMT2D, TSPY2 (all junction-only novelty), plus H2BC12 E114K (conditional).

## Antigen regions

Full evidence per region: `selected.tsv` (one row per region, lead column schema).
Affinities are fresh MHCflurry predictions of the best 8–11-mer × allele in the window.

| Antigen | Gene | Source | Tier / status | Best affinity (nM) | Alleles ≤500 nM | Single-allele at 500 nM |
|---|---|---|---|---|---|---|
| SF01 | AQP3 | MUTATION | paired_DNA_supported | 29.4 | 5 | no |
| SF02 | ARHGAP40 | MUTATION | paired_DNA_supported | 138.2 | 3 | no |
| SF03 | H2BC12 | MUTATION | **conditional** (missing expression annotation) | 139.8 | 1 | **yes (A*33:03)** |
| SF04 | HSPG2 | MUTATION | paired_DNA_supported | 81.1 | 1 | **yes** |
| SF05 | KCTD19 | MUTATION | paired_DNA_supported | 36.9 | 4 | no |
| SF06 | KLF5 | MUTATION | paired_DNA_supported | 37.0 | 3 | no |
| SF07 | PLEC | MUTATION | paired_DNA_supported | 413.3 | 1 | **yes** |
| SF08 | PPIC | MUTATION | paired_DNA_supported | 36.5 | 3 | no |
| SF09 | RXRA | MUTATION | paired_DNA_supported | 58.1 | 4 | no |
| SF10 | SF1 | MUTATION | paired_DNA_supported | 41.5 | 3 | no |
| SF11 | SYT12 | MUTATION | paired_DNA_supported | 243.7 | 1 | **yes** |
| SF12 | TGFBR2 | MUTATION | paired_DNA_supported | 38.0 | 1 | **yes** |
| SF13 | VPS11 | MUTATION | paired_DNA_supported | 115.6 | 3 | no |
| SF14 | F5 | SPLICE | junction_only_novelty, PASS | 14.2 | 5 | no |
| SF15 | HLA-B | SPLICE | junction_only_novelty, PASS | 66.2 | 2 | no |
| SF16 | HSF1 | SPLICE | junction_only_novelty, PASS | 44.5 | 3 | no |
| SF17 | KMT2D | SPLICE | junction_only_novelty, PASS (READ_END_BIAS flag inherited) | 85.6 | 4 | no |
| SF18 | TSPY2 | SPLICE | junction_only_novelty, PASS | 25.7 | 4 | no |
| SF19 | LRRC36 | CTA/SELF | recorded_RNA_source | 21.0 | 5 | no |
| SF20 | PRAME | CTA/SELF | recorded_RNA_source | 15.3 | 6 | no |
| SF21 | SPATA31A6 | CTA/SELF | recorded_RNA_source | 33.6 | 5 | no |
| SF22 | TEKT5 | CTA/SELF | recorded_RNA_source | 23.7 | 6 | no |
| SF23 | YBX2 | CTA/SELF | recorded_RNA_source | 29.8 | 4 | no |
| SF24 | CDH1--NRG1 | FUSION | recorded_RNA_source (encoded N) | 31.7 | 6 | no |
| SF25 | CDH1::NRG1 | FUSION | recorded_RNA_source | 60.4 | 1 | **yes** |

All five splice regions have recorded normal-screen PASS, patient-RNA-confirmed sequence
(`patient_RNA_sequence_confirmed=True`), ≥2 supporting RNA fragments in a stop-free
annotated upstream frame, and were independently re-checked against the patient's
genotyped normal GVCF (`normal-model-check.json`: no nonself peptide found in the local
annotated normal model).

## Methods

Replay everything with `python designs/specificity-first/run.py` (offline; prerequisites
below). Stages, in order:

1. **Specificity attribution screen** (`splice_specificity_screen.py`): for all 37 splice
   windows in the Design #2 candidate pool, align both junction sides to their reference
   proteins, enumerate every 8–11-mer overlapping the junction codon, test exact matches
   against the 85,968 primary protein-coding Ensembl 103 translations, then re-test after
   reverting detected extra substitutions. Result: 22 junction-only-novelty vs 15
   substitution-dependent windows; the four contested genes reproduce issue #1's
   counterfactual exactly (26→0, 22→0, 18→0, 23→0 nonself peptides). Output:
   `splice-specificity-attribution.tsv`, `splice-specificity-summary.json`.
2. **Selection** (`build_selected.py`, `build_regions.py`): keep all 12 recorded
   mutation regions, 5 CTA/self and 2 fusion regions; select 5 junction-only-novelty
   splice windows ranked by best target affinity, RNA support, expression and allele
   coverage. H2BC12 E114K is added with an explicit conditional label.
3. **Independent normal check** (`normal_model_check.py`): rebuild CDS sequences from the
   annotated Ensembl 103 GTF, apply the patient's genotyped normal variants
   (MAX_COMBINATIONS=4096), and confirm no selected nonself peptide appears in any
   annotated normal possibility.
4. **Fresh predictions** (`build_prediction_cache.py`): 96,973 unique peptides × 6
   recorded HLA-I alleles (A*01:01, A*33:03, B*08:01, B*15:08, C*01:02, C*07:02),
   MHCflurry 2.3.14 with the 2.3.0 pan-allele affinity weights.
5. **HLA-robustness pass** (`hla_robustness.py`): eligible-allele counts at
   50/150/500/1000 nM per region (`hla-allele-dependence.tsv`), 64-subset lost-allele
   sensitivity (`lost-allele-sensitivity.tsv`), and paired mutant/reference
   discrimination for the conditional candidates (`conditional-candidates.tsv`).
6. **Assembly** (`build_terminal_cache.py`): pinned recorded routing elements (SEC =
   HLA-C*01:02 residues 1–26; MITD = HLA-B*15:08 residues 308–362), joint CP-SAT
   optimization of antigen order and G/GS linkers (lexicographic: strong-binder count,
   then mild-binder count, then worst junction percentile, then added residues; seed 31).
   All four phases finished OPTIMAL: 2 strong, 27 mild, 42 added linker residues.
   Hydrophobicity screen PASS (`construct/hydrophobicity.json`; 0 unexpected
   high-Kyte–Doolittle windows; 8 routing-boundary review windows at the fixed
   SEC/MITD edges, same policy as the lead).
7. **Fresh full-junction check** (`fresh_junction_check.py`): re-predicts all 4,476
   boundary peptide × allele pairs and compares with the optimization cache
   (max relative affinity difference 9.5e-08; threshold classifications identical).
   Output: `construct/final-check.json`.
8. **RNA encoding** (`construct/mrna/run.py`): the recorded PopRNA optimizer (pinned
   ViennaRNA 2.7.2, python-codon-tables 0.1.18, biopython 1.86, numpy 1.26.4 in a
   dedicated venv) encodes the fixed protein synonymously with the frozen human codon
   table and motif/repeat avoidance. Verified: exact protein translation, zero rare
   codons below cutoff, zero forbidden motifs, zero repeated 12-nt pairs, CAI 0.777,
   no XbaI site. The longest paired-bracket run is 24 nt against the 14-nt heuristic
   target — the lead design also misses this target (21 nt); it is recorded honestly,
   not gated. **PolyA and any downstream extension remain unspecified**, as in the lead;
   only the known upstream RNA/DNA segment is exported
   (`verification.json`: `KNOWN_NUCLEOTIDE_SEGMENTS_VERIFIED_POLYA_UNSPECIFIED`).
   The coding sequence is pinned by SHA-256 and reproduced bit-identically on rerun.
9. **Comparison** (`compare_to_design2.py`): region-level diff against Design #2.

The repository's `yash-vaccine verify` command still passes for the untouched lead
design; `yash-vaccine compare` only supports recipe-stage outputs, so
`compare_to_design2.py` compares the two independent designs' tables directly.

### Prerequisites for replay

- MHCflurry 2.3.14 with the 2.3.0 pan-allele affinity weights
  (`models_class1_pan.selected.20260928.tar.bz2`) extracted under
  `~/.local/share/mhcflurry/`.
- Ensembl 103 GRCh38 proteome (SHA-256 `f8df3d7c…0505`) and GTF
  (`55bbf597…150a`) under `/workspace/ref/`. GENCODE 37 was unreachable at design
  time; the Ensembl 103 GTF (same release as the recorded proteome) was used instead.
- The recorded normal coding GVCF (`data/` download) for the independent normal check.
- A dedicated venv for the RNA stage: `numpy==1.26.4 ViennaRNA==2.7.2
  python-codon-tables==0.1.18 biopython==1.86` (the pin conflicts with the newer numpy
  used by the prediction stages).

### Not included in Git

`prediction-cache.npz` (14 MB, SHA-256 `52e73fecc15326fdb8acc4fcb5d35ce1b22571406cb3f8d8e55411088738822b`)
is a derived intermediate —
96,973 unique peptides × 6 alleles of fresh MHCflurry predictions used by the assembly
stage. Regenerate it with `python designs/specificity-first/build_prediction_cache.py`
(same MHCflurry version and weights as above); the fresh full-junction check then
confirms the rebuilt cache matches the recorded boundary predictions.

## Conditional and lower-tier components (explicit labels)

- **SF03 H2BC12 E114K** — `conditional_missing_expression_annotation`. Strong paired
  somatic + RNA support; predicted ELAKHAVSK/A*33:03 = 139.8 nM vs ELAKHAVSE =
  2998.1 nM (~21.5× gain). Expression stays `UNAVAILABLE` (missing annotation, not zero
  expression). Single-allele dependent at 500 nM.
- **FES R804W short context** — NOT in the cassette. The construct pipeline requires
  exactly 25-aa regions; the 15-aa context `EQCWAYEPGQWPSFS` is documented in
  `conditional-candidates.tsv` as a validation candidate (predicted AYEPGQWPSF/C*07:02
  = 94.5 nM). No flanking residues were invented to force it into a 25-aa format.
- **SF17 KMT2D** — inherits `native_event_review_flags=READ_END_BIAS` from the recorded
  evidence.
- **SF24 CDH1--NRG1** — no recorded target interval; robustness analysis uses the full
  window. Encoded as N (observed), not the hypothetical D variant.

## Limitations

- Binding predictions do not establish presentation or immunogenicity; no mass-spec or
  T-cell evidence is supplied.
- The specificity screen covers the 37-window Design #2 splice pool; the full 80-gene
  pool was not exhausted.
- The independent normal check tests nonself peptides against annotated CDS possibilities
  derived from genotyped variants; it does not re-align raw normal reads (the 143 GB BAM
  was not downloaded).
- Clinical status is unknown to this design (latest public report October 2025).
- PolyA, cap chemistry, nucleotide modifications and any downstream extension are
  unspecified; the exported GenBank record is a partial construct by design.
