#!/usr/bin/env bash
set -euo pipefail
WORK=/work/lens-expression-aware
export APPTAINER_BINDPATH=/data,/work
mkdir -p "$WORK/annotation-test"
cat > "$WORK/annotation-test/in.tsv" <<'TSV'
identity	peptide	allele	pos
cta1	AAAAAAAAA	HLA-A*01:01	0
fusion1	CCCCCCCCC	HLA-A*01:01	0
splice1	DDDDDDDDD	HLA-A*01:01	0
TSV
cat > "$WORK/annotation-test/in.fa" <<'FA'
>cta1
;antigen_source:CTA/SELF transcript_id:ENST1.1 pep_context:AAAAAAAAA
AAAAAAAAA
>fusion1
;antigen_source:FUSION fusion_left_breakpoint:chr1:100:+ pep_context:CCCCCCCCC
CCCCCCCCC
>splice1
;antigen_source:SPLICE splice_coords:chr1:100-200(+) pep_context:DDDDDDDDD
DDDDDDDDD
FA
singularity exec --cleanenv "$WORK/containers/lenstools.sif" python3 /opt/lenstools/lenstools.py annotate-pmhcs \
 -s "$WORK/annotation-test/in.tsv" -f "$WORK/annotation-test/in.fa" -o "$WORK/annotation-test/original.tsv"
singularity exec --cleanenv --bind "$WORK/metadata.annotation-header-fix.py:/opt/lenstools/src/metadata.py:ro" \
 "$WORK/containers/lenstools.sif" python3 /opt/lenstools/lenstools.py annotate-pmhcs \
 -s "$WORK/annotation-test/in.tsv" -f "$WORK/annotation-test/in.fa" -o "$WORK/annotation-test/patched.tsv"
cmp "$WORK/annotation-test/original.tsv" "$WORK/annotation-test/patched.tsv"
echo 'Native annotation header patch: byte-identical output for all three sources.'
