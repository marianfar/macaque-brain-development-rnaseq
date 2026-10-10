#!/usr/bin/env bash
# ==============================================================================
# 01_qc_and_trim.sh: Quality Control and Read Preprocessing
# ==============================================================================
set -euo pipefail

RAW_DIR="data/raw_fastq"
QC_RAW_DIR="qc_raw"
CLEAN_DIR="data/cleaned_fastq"
QC_CLEAN_DIR="qc_cleaned"

mkdir -p "$QC_RAW_DIR" "$CLEAN_DIR" "$QC_CLEAN_DIR"

# 1. Initial FastQC across raw FASTQs
fastqc "$RAW_DIR"/*.fastq.gz -o "$QC_RAW_DIR"
multiqc "$QC_RAW_DIR" -o "$QC_RAW_DIR"

# 2. Standardize Phred64 (Illumina 1.5) to Phred33
for fq in "$RAW_DIR"/*_phred64.fastq.gz; do
    [ -f "$fq" ] || continue
    base=$(basename "$fq" _phred64.fastq.gz)
    seqtk seq -Q64 -V "$fq" | gzip > "${RAW_DIR}/${base}_phred33.fastq.gz"
done

# 3. Trimming and filtering with fastp
for fq in "$RAW_DIR"/*_phred33.fastq.gz; do
    [ -f "$fq" ] || continue
    base=$(basename "$fq" _phred33.fastq.gz)
    fastp \
        -i "$fq" \
        -o "${CLEAN_DIR}/${base}_clean.fastq.gz" \
        --html "${QC_CLEAN_DIR}/${base}_fastp.html" \
        --json "${QC_CLEAN_DIR}/${base}_fastp.json" \
        -w 8 -q 30 -u 30 \
        --cut_front --cut_tail \
        --cut_window_size 4 --cut_mean_quality 30 \
        --length_required 25
done

# 4. Final post-trim MultiQC
fastqc "$CLEAN_DIR"/*.fastq.gz -o "$QC_CLEAN_DIR"
multiqc "$QC_CLEAN_DIR" -o "$QC_CLEAN_DIR"
