#!/usr/bin/env bash
# ==============================================================================
# 02_alignment_and_quant.sh: STAR Alignment and RSEM Quantification
# Reference: Macaca mulatta (rheMac10)
# ==============================================================================
set -euo pipefail

CLEAN_DIR="data/cleaned_fastq"
STAR_INDEX="reference/star_rheMac10_index"
RSEM_REF="reference/rsem_rheMac10_ref/rheMac10"
ALIGN_DIR="aligned"
QUANT_DIR="rsem_genes_results"

mkdir -p "$ALIGN_DIR" "$QUANT_DIR"

for fq in "$CLEAN_DIR"/*_clean.fastq.gz; do
    [ -f "$fq" ] || continue
    sample=$(basename "$fq" _clean.fastq.gz)

    STAR \
        --runThreadN 8 \
        --genomeDir "$STAR_INDEX" \
        --readFilesIn "$fq" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${ALIGN_DIR}/${sample}_" \
        --outSAMtype BAM Unsorted \
        --quantMode TranscriptomeSAM

    rsem-calculate-expression \
        --bam --no-bam-output -p 8 \
        "${ALIGN_DIR}/${sample}_Aligned.toTranscriptome.out.bam" \
        "$RSEM_REF" \
        "${QUANT_DIR}/${sample}"
done
