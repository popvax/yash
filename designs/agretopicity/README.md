# Agretopicity for the substitution regions of Design #2

This is an analysis, not a competing construct. It adds one reported quantity to
regions Design #2 already selected, and it changes nothing about the selection.

## Design and rationale

### The question

Design #2 ranks and admits mutation regions by class I binding affinity. Affinity
answers "will this peptide be presented". It does not answer "has the immune
system already been trained to ignore something that looks just like it".

Agretopicity answers the second question. LENS computes it natively and this
repository did not carry it over, so the quantity is missing for every region
even though everything needed to compute it is already here.

A worked example of why the two can disagree. A peptide binding at 20 nM whose
wild-type counterpart binds at 25 nM is a poor target, because the thymus
negatively selected T cells against a peptide that is almost the same and is
presented just as well. A peptide binding at 200 nM whose wild-type counterpart
binds at 5,000 nM is a better target, because the 25-fold gain in presentation is
new. Affinity ranks the first peptide ten times higher. Agretopicity separates
them, 0.80 against 0.04.

### Direction, which is the one thing to get right

This follows LENS exactly. In `lenstools/src/metadata.py`,
`add_agreto_metadata()` computes `mhcflurry_mt / mhcflurry_wt` on nM affinities,
so **lower is better**. That is the reciprocal of the DAI convention used in much
of the literature, where the ratio is wild-type over mutant and higher is better.
A reimplementation that follows the papers rather than the tool inverts every
ranking, so the column here is named `agretopicity_mt_over_wt` and nothing is
abbreviated to "agretopicity" alone.

### What differs from Design #2

Nothing in the design. No region is added, removed, reordered or rescored, no
sequence changes, and no recorded artifact or hash changes. This directory adds
two tables and a run record.

### Method

1. Take the 12 selected MUTATION regions from
   [`../design-2/inputs/integrated/candidate-pool.tsv`](../design-2/inputs/integrated/candidate-pool.tsv).
   11 are single-residue substitutions. TGFBR2 `p.N515fs` is a frameshift and is
   out of scope, because a frameshift has no wild-type counterpart peptide.
2. Build the wild-type window by restoring the reference residue at
   `target_start_in_window`. The reference residue is already stated in the DNA
   annotation protein change, so for KLF5 `p.D418N` the asparagine at window
   index 12 becomes aspartate again. The script asserts that the alternate
   residue is where the recorded window says it is before touching anything.
3. Enumerate every 8, 9, 10 and 11mer of the window that covers the substituted
   position. That gives 408 mutant and wild-type pairs across the 11 regions.
4. Predict both sides with the same model, MHCflurry 2.3.14 loaded from the
   `affinity_predictor` component of the `models_class1_presentation` 2.3.0
   bundle, across the patient's six class I alleles. 816 unique peptides, 4,896
   peptide and allele predictions, 9.4 seconds.
5. Divide, and report. No threshold is applied.

LENS obtains its wild-type peptide by BLASTing mutant peptides against the
antigen.garnish database and taking the best hit's subject sequence. Restoring
the annotated reference residue is more direct for a substitution and needs no
external database, and it cannot return a wrong-gene homologue. It only works for
substitutions, which is why the scope here is 11 regions rather than the 19 LENS
would cover.

### Results

Ordered by the agretopicity of each region's tightest mutant binder. MT and WT
are affinities in nM for that same peptide and allele.

| Region | Design #2 label | Gene | Change | Best mutant peptide | Allele | MT nM | WT nM | Agretopicity | Lowest in region |
|---|---|---|---|---|---|---:|---:|---:|---:|
| `REG_9dde223e9bf58444173b` | AG08 | KCTD19 | p.G920R | `RSVSYSILR` | HLA-A*33:03 | 36.9 | 11691.6 | **0.003** | 0.003 |
| `REG_9954d136ed6d4ba9f759` | AG05 | HSPG2 | p.R3952L | `LTNTHHELL` | HLA-C*01:02 | 81.1 | 14358.9 | **0.006** | 0.006 |
| `REG_f93fe5a69c414329a31d` | AG12 | SYT12 | p.P391L | `NVGHVIIGL` | HLA-C*01:02 | 243.7 | 23261.9 | **0.010** | 0.010 |
| `REG_c310fc71d7b0fb823535` | AG03 | AQP3 | p.C11Y | `SRYGEMLHI` | HLA-C*07:02 | 29.4 | 110.5 | **0.266** | 0.007 |
| `REG_1dfce9e906fa21c20196` | AG02 | SF1 | p.R218Q | `IYNSEGKQL` | HLA-C*07:02 | 41.5 | 72.8 | **0.569** | 0.083 |
| `REG_8687690f602539a7060e` | AG04 | RXRA | p.A241S | `SVEPKTETY` | HLA-B*15:08 | 58.1 | 94.0 | **0.618** | 0.251 |
| `REG_e16d2a87273645d8d7b4` | AG01 | KLF5 | p.D418N | `WRFARSNEL` | HLA-C*07:02 | 37.0 | 52.0 | **0.712** | 0.296 |
| `REG_e6aeb364b2cc7e909162` | AG07 | PPIC | p.P187S | `ATDGHDRSL` | HLA-C*01:02 | 36.5 | 49.9 | **0.732** | 0.187 |
| `REG_7176562878cc11bcbfaa` | AG10 | ARHGAP40 | p.R438W | `LHQGWPPKL` | HLA-C*07:02 | 138.2 | 89.7 | **1.541** | 0.260 |
| `REG_f0a90a31c4c6e8899224` | AG11 | PLEC | p.R2171W | `ERAEQESAWQL` | HLA-C*07:02 | 413.3 | 243.7 | **1.696** | 0.231 |
| `REG_7bb9945d495d7fb8ce3f` | AG09 | VPS11 | p.D151A | `FTAGSVTL` | HLA-C*01:02 | 115.6 | 51.3 | **2.252** | 0.295 |

The AG labels are carried only so the table can be read against the recorded
cassette. This analysis keys on `region_id`, since AG numbering belongs to one
selection.

Three observations, in order of how much they would change a decision.

**Three regions create presentation where there was none.** KCTD19 p.G920R binds
at 36.9 nM against 11,691.6 nM for its wild-type counterpart, a 317-fold gain at
a position no T cell repertoire was tolerized against. HSPG2 p.R3952L and SYT12
p.P391L are the same pattern, with wild-type counterparts at 14,358.9 nM and
23,261.9 nM. On a quality view these are the strongest substitution antigens in
the construct.

**Three regions bind worse than the self peptide they came from.** VPS11 p.D151A
binds at 115.6 nM while its own wild-type binds at 51.3 nM, so the better
presented of the two peptides is the one the patient has been tolerant to since
birth. ARHGAP40 p.R438W and PLEC p.R2171W are in the same position. All three
passed the 500 nM gate and look reasonable by affinity alone.

**Affinity rank and agretopicity rank disagree in a way that matters.** SYT12 at
243.7 nM and PLEC at 413.3 nM are the two weakest admitted mutation binders in
Design #2, and they sit at opposite ends of this table, 0.010 against 1.696.
Affinity treats them as near neighbours.

The `Lowest in region` column is there to prevent an over-reading of the main
column. AQP3 p.C11Y has an unremarkable 0.266 for its tightest binder, but some
other peptide and allele pair in the same region reaches 0.007. A region can
contain a highly distinct epitope that is not its best binder.

## Sequences and results

No construct was assembled and no sequence was changed, so there is no
`cassette.fasta` or mRNA here. The mutant and wild-type peptides are columns in
the tables below rather than separate FASTA files.

- [`agretopicity.tsv`](agretopicity.tsv), 2,448 rows. One row per
  mutant and wild-type peptide pair per allele, with both affinities, both
  percentiles and the ratio.
- [`region-summary.tsv`](region-summary.tsv), 11 rows. The table
  above, plus each region's median agretopicity and the recorded best affinity
  for comparison.
- [`agretopicity-run.json`](agretopicity-run.json). Scope,
  alleles, counts, versions, the model bundle check and the parity result.

Every number here is new inference. Nothing is a relabelled cached prediction.
There is no experimental result and no raw read processing.

## Reproduction instructions

### Inputs, all already in the repository

| Input | Path |
|---|---|
| Selected mutation regions and windows | `designs/design-2/inputs/integrated/candidate-pool.tsv` |
| Region identifiers and cassette labels | `designs/design-2/selection/selected.json`, `designs/design-2/selection/selected.fasta` |
| Patient class I alleles | `config/case.json` |
| Recorded affinities, used only as a parity check | `designs/design-2/inputs/predictions/mutations/best-peptide-per-HLA.tsv` |
| Pinned model manifest | `designs/initial/inputs/checkpoints/provenance/construct-model-files.json` |

### The one external download

```sh
mhcflurry downloads fetch models_class1_presentation --release 2.3.0
mhcflurry downloads path  models_class1_presentation --release 2.3.0
```

Direct URL, from `mhcflurry/downloads.yml` at tag 2.3.14:
`https://github.com/openvax/mhcflurry/releases/download/2.3.0/models_class1_presentation.20260928.tar.bz2`,
160,356,140 bytes compressed.

The script hashes all 44 files of the extracted bundle and compares them against
`construct-model-files.json`, the manifest this repository already pins. They are
identical, in order, byte for byte, and the run record states that as
`models_match_pinned_bundle`. So this analysis uses the same weights the recorded
predictions used, and nothing new is pinned.

### Dependencies

```
mhcflurry 2.3.14
numpy     2.5.4
torch     2.15.0, CPU
Python    3.12.3
```

Those are not the versions this repository recorded. `construct-inference.json`
records torch 2.14.1+cu130 and numpy 2.4.6 on a GPU. The parity check exists
because of that difference, see below.

### Commands

```sh
pip install 'mhcflurry==2.3.14'
MODELS="$(mhcflurry downloads path models_class1_presentation --release 2.3.0)/models/affinity_predictor"
python designs/agretopicity/run.py --models "$MODELS"
```

The script writes its three output files into `designs/agretopicity/`. It reads nothing outside
the repository and the model bundle, starts no remote job and downloads nothing
itself.

### Parity against the recorded predictions

The run compares its own predictions against every row of
`best-peptide-per-HLA.tsv` for the 11 variants in scope, 66 peptide and allele
pairs. Maximum relative affinity difference **4.1e-07**, with no row falling
outside the scan. Per region, the tightest mutant binder found here agrees with
the recorded `best_target_affinity_nm` to within 3.5e-07.

In other words, the mutant side of every ratio reproduces this repository's own
numbers on different hardware and different library versions. Only the wild-type
side is new.

## Interpretation

### What the evidence covers

11 of the 24 regions in Design #2. Class I only. Predicted affinity for both
sides of each ratio, from the same model and the same bundle. No claim about
presentation, cleavage, T-cell recognition or clinical benefit.

### What it does not cover

- Splice, fusion and CTA/self regions have no wild-type counterpart that this
  method can construct. LENS reports agretopicity for SNVs, indels, fusions and
  CTA/self, and not for splice, so a fusion extension would be possible with more
  care about which side of the junction the comparison is drawn from. It is not
  attempted here.
- TGFBR2 p.N515fs is excluded for the same structural reason. A frameshift tail
  is not a substitution of one residue.
- The wild-type window is reconstructed from the annotated reference residue. It
  is not read off the patient's own normal sequencing, so a germline variant in
  the same window would not be reflected. For these 11 regions the reference
  residue is the one the annotation states, and nothing here checks it against
  matched-normal reads.
- Agretopicity is one half of the Łuksza and Balachandran neoantigen quality
  model. The other half, foreignness relative to a pathogen epitope set, is not
  computed. Quality in that model is a product of the two, so this is a component
  and not a quality score.
- No threshold. A grep across the LENS repositories for an agretopicity threshold
  returns nothing, and `prioritize_peptides` in `lenstools/src/pep_utils.py`
  builds its priority score from cancer cell fraction, read support and binding
  affinity, with agretopicity absent. It is reported for the analyst. Applying a
  cutoff here would invent a decision rule that its own authors declined to make.

### Open questions I cannot answer from this evidence

1. Does a high agretopicity, meaning a mutant that binds worse than its own
   wild-type, actually predict a dead antigen in a vaccine context, or only in
   the immunoediting setting where the quality model was fitted? The Łuksza and
   Balachandran work is observational in resected pancreatic cancer. It is the
   only outcome data in this disease, and it is not a vaccine trial.
2. For the three regions with agretopicity above 1, is the right response to drop
   them, keep them for antigen spreading, or ignore the metric? This analysis
   does not answer that, and it deliberately does not act on it.
3. VPS11 p.D151A at 2.252 is also the region whose recorded best affinity is
   115.6 nM. Would a quality-aware selection have reached for a different 25mer
   window over the same variant? The window choice is upstream of everything
   here.

### Unfinished work

A fusion extension, a foreignness term to pair with this one, and matched-normal
confirmation of the reference residue. None is started.
