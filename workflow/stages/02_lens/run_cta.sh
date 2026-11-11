#!/usr/bin/env bash
set -euo pipefail
ROOT=${LENS_RUN_ROOT:-/data/results/TN25-112103/lens-expression-aware-20261030}
REF=/data/raft/references/homo_sapiens
WORK=/work/lens-expression-aware
OLD=/data/results/TN25-112103/original-fastqs/outputs/samples/TN25_research/TN25-112103/ar-TN25-112103
export APPTAINER_BINDPATH=/data,/work
mkdir -p "$ROOT/cta" "$ROOT/source"
cp "$WORK/prepare_expression.py" "$ROOT/source/"
python3 "$WORK/prepare_expression.py" \
 --expression /data/inputs/TN25-112103/vendor-transformed-expression.tsv \
 --quants "$OLD/salmon_aln_quant/TN25_research-TN25-112103-ar-TN25-112103.quant.sf" \
 --gtf "$REF/annot/gencode.v37.annotation.with.hervs.gtf" \
 --genes "$REF/cta_self/cta.pirlygenes.06APR2026.gene_list" --out "$ROOT/cta"
lens() { singularity exec --cleanenv "$WORK/containers/lenstools.sif" python3 /opt/lenstools/lenstools.py "$@"; }
lens filter-expressed-self-genes -q "$ROOT/cta/expression-adapted.quant.sf" \
 -m SuppliedGeneExpression -t 0 -g "$REF/cta_self/cta.pirlygenes.06APR2026.gene_list" \
 -f "$ROOT/cta/cta.selected.gtf" -o "$ROOT/cta/expressed-selfs.txt"
lens get-expressed-selfs-bed -e "$ROOT/cta/expressed-selfs.txt" \
 -g "$ROOT/cta/cta.selected.gtf" -o "$ROOT/cta/cds.regions.tsv"
awk '{print $1 ":" $2 "-" $3}' "$ROOT/cta/cds.regions.tsv" | sort -u > "$ROOT/cta/cds.regions.txt"
mapfile -t regions < "$ROOT/cta/cds.regions.txt"
singularity exec --cleanenv "$WORK/containers/starfusion.sif" samtools faidx \
 "$REF/fasta/Homo_sapiens.assembly38.fa" "${regions[@]}" > "$ROOT/cta/cds.fa"
lens make-self-antigen-peptides -s "$ROOT/cta/cds.fa" -g "$ROOT/cta/cta.selected.gtf" \
 -e "$ROOT/cta/expressed-selfs.txt" -o "$ROOT/cta/cta.pep.fa" -n "$ROOT/cta/cta.nt.fa"
