#!/usr/bin/env bash
# Run native pre-HLA source tools from preserved original-FASTQ checkpoints.
set -euo pipefail
ROOT=${LENS_RUN_ROOT:-/data/results/TN25-112103/lens-expression-aware-20261030}
OLD=/data/results/TN25-112103/original-fastqs/outputs
REF=/data/raft/references/homo_sapiens
WORK=/work/lens-expression-aware
export APPTAINER_BINDPATH=/data,/work
mkdir -p "$ROOT"/{splice,fusion,logs}
cd "$WORK"
case "${1:?splice or fusion}" in
splice)
  singularity exec --cleanenv --pwd "$WORK" \
    --bind "$WORK/msnaf-core.8-11.py:/usr/local/lib/python3.10/site-packages/msnaf/core.py:ro" \
    "$WORK/containers/msnaf1.2.1.sif" msnaf \
    --counts "$OLD/samples/multisample/msnaf/counts.msnaf.input.txt" \
    --refs "$REF/splice_variants/snaf/snaf-data" \
    --genome-fasta "$REF/fasta/Homo_sapiens.assembly38.fa" \
    --output "$ROOT/splice/msnaf_peptides.csv" --species human \
    --filter-mode maxmin --t-min 20 --n-max 3
  ;;
fusion)
  singularity exec --cleanenv --pwd "$WORK" "$WORK/containers/starfusion.sif" \
    /usr/local/lib/STAR-Fusion/FusionInspector/FusionInspector \
    --fusions "$OLD/samples/TN25_research/TN25-112103/ar-TN25-112103/starfusion/TN25_research-TN25-112103-ar-TN25-112103.star-fusion.fusion_predictions.abridged.coding_effect.tsv" \
    --left_fq /data/inputs/TN25-112103/original/RNA_TN25-112103_S45.R1.fastq.gz \
    --right_fq /data/inputs/TN25-112103/original/RNA_TN25-112103_S45.R2.fastq.gz \
    --genome_lib_dir "$REF/fusions/starfusion/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play" \
    --output_dir "$ROOT/fusion/fusioninspector" --CPU 6 \
    --annotate --examine_coding_effect --write_intermediate_results \
    --STAR_xtra_params "--outTmpDir $WORK/fusion-star-tmp"
  ;;
*) exit 2;;
esac
