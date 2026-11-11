# Supplementary mutation sequence review

This review covers the **18 upstream mutation candidates** before hydropathy
screening. The construct includes 17 after PTCH1 is excluded for compatibility. The
review adds observations and individual interpretation; it does not introduce
new eligibility gates or replace native GATK filtering or Isovar reconstruction.

[Summary](summary.tsv) distinguishes DNA alternate reads/depth, alternate RNA
fragment names, the union denominator and fragments supporting the longer Isovar
context. These are different measurements. [Results](result.json) record the
fresh independent DNA/RNA base check, keyed by stable variant and region IDs.

Every reviewed SNV has alternate bases in both datasets and both read
orientations under the supplementary quality policy. The eight regions in the supplementary reference-context comparison differ
from a compatible annotated reference by exactly their intended substitution. KCTD19, NRXN2 and SYT12 each have only two
alternate-supporting RNA fragment names and two longer-context fragments. These
small samples are retained explicitly. No molecular independence is inferred
without UMIs, and no somatic status or peptide presentation is established.

The existing PPIC window includes its documented N190S background difference.
The [reference-context comparison](reference-context-check.json) records the
local reference sequences, compatible transcript IDs and every difference.
GENCODE protein references are annotations, not normal-tissue measurements.

From the repository root, with staged BAMs in the recorded Modal volume:

```sh
python designs/initial/review/mutation-read-check/launch.py --output work/selected-mutation-check --execute
```

Omit `--execute` for a plan. The launcher uses the request's explicit BAM paths,
does read-only locus checks, records the executed source/request hashes,
terminates its sandbox and verifies the persistent-volume result. It downloads
only the small result. [Request](request.json), [worker](check.py) and
[execution receipt](modal-state.json) document this run.

To repeat the local reference comparison, download the reference identified and
hashed in its result, then run:

```sh
python designs/initial/review/mutation-read-check/reference_context.py --reference-proteins /path/to/gencode.v37.pc_translations.fa.gz --output work/reference-context-check.json
```

This is an ungapped comparison within supporting transcripts allowing at most
two substitutions. It stops when no compatible context is found; it is not a
general protein-alignment or de novo transcript-reconstruction method.

This raw-read review covers the 18 upstream mutation candidates. The separate
hydropathy compatibility gate excludes PTCH1 from assembly; its original measurements
and this execution receipt remain unchanged. The construct contains 17 mutation regions.
