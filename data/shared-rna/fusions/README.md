# CDH1::NRG1 fusion evidence

The initial design includes two different RNA-derived hypotheses from this gene
pair. They have distinct CDH1 breakpoints and must not share support fractions.

| Product | GRCh38 breakpoints | Selected 25-aa region | Best target affinity |
|---|---|---|---:|
| Predicted frameshift tail | chr16:68815759:+ → chr8:32727949:+ | EPDTFMEQKITYIYIHHWDKPSCKM | 70.20 nM |
| In-frame junction | chr16:68819425:+ → chr8:32727949:+ | STYTALIIATDNATSTSTTGTSHLV | 508.95 nM |

The frameshift-context window comes from the native LENS/STAR-Fusion pool.
FusionInspector confirmation uses the exact breakpoint pair. Its recorded local
support is 62 junction-plus-spanning units and 473 unfused left-partner units;
the ranking proxy uses the supplied left-partner gene value times 62/(62+473).
The source window scores and reference-filter policy are retained in
[the LENS checkpoint](../checkpoints/lens-scored-windows.tsv.gz).

The in-frame source call reports 200 junction reads and one spanning fragment.
The [request](inframe-request.json) records its source protein, exact breakpoints,
join codon and six HLA-I alleles. The [binding result](inframe-binding-review.json)
contains fresh, model-file-verified MHCflurry 2.3.14 predictions for 38 short
peptides across six alleles. Eight reference-matching peptides are excluded.
ATDNATSTST is the strongest remaining target at 508.95 nM, rank 0.7111%, for
HLA-A*01:01. The [selection record](inframe-selection.json) pins the inputs.

The in-frame driver junction is explicitly included despite being just above the
general 500 nM gate. It receives no invented expression-weighted score or RNA
fraction borrowed from the frameshift event. The replay adapter is
[candidates.py](../../../src/yash_vaccine/candidates.py). The executed binding worker
is [included](inframe_binding.py); original raw-read evidence and references are
listed in the [data catalog](../../catalog/README.md).

Read/support counts are not UMI molecule counts. RNA support and predicted HLA
binding do not establish full-length protein translation or presentation.
