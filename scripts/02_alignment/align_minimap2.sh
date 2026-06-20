#!/bin/bash

set -euo pipefail

# Align Oxford Nanopore reads to the P. falciparum 3D7 reference
# genome using minimap2 and generate BAM quality statistics.

PROJECT_DIR="."

FASTQ_DIR="${PROJECT_DIR}/trimmed_fastq"
REF="${PROJECT_DIR}/PlasmoDB-68_Pfalciparum3D7_Genome.fasta"
BED="${PROJECT_DIR}/amplicons.bed"

BAM_DIR="${PROJECT_DIR}/bam_files"
DEPTH_DIR="${PROJECT_DIR}/depth_files"
SUMMARY_DIR="${PROJECT_DIR}/mapping_summaries"

mkdir -p "${BAM_DIR}"
mkdir -p "${DEPTH_DIR}"
mkdir -p "${SUMMARY_DIR}"

COVERAGE_MATRIX="${PROJECT_DIR}/coverage_matrix.tsv"

echo -e "Sample\tMSP2\tCSP\tMDR1\tCRT\tK13" > "${COVERAGE_MATRIX}"

for fq in "${FASTQ_DIR}"/*.fastq.gz
do

    sample=$(basename "${fq}" _trimmed.fastq.gz)

    echo "Processing ${sample}"

    minimap2 \
        -t 8 \
        -ax map-ont \
        --secondary=no \
        -p 0.9 \
        "${REF}" \
        "${fq}" \
    | samtools view -b -F 4 \
    | samtools sort -o "${BAM_DIR}/${sample}.bam"

    samtools index "${BAM_DIR}/${sample}.bam"

    samtools depth -a "${BAM_DIR}/${sample}.bam" \
        > "${DEPTH_DIR}/${sample}.depth.txt"

    samtools flagstat "${BAM_DIR}/${sample}.bam" \
        > "${SUMMARY_DIR}/${sample}.flagstat.txt"

    samtools idxstats "${BAM_DIR}/${sample}.bam" \
        > "${SUMMARY_DIR}/${sample}.idxstats.txt"

    MSP2=$(bedtools coverage -a "${BED}" -b "${BAM_DIR}/${sample}.bam" -mean | awk '$4=="MSP2"{print $NF}')
    CSP=$(bedtools coverage -a "${BED}" -b "${BAM_DIR}/${sample}.bam" -mean | awk '$4=="CSP"{print $NF}')
    MDR1=$(bedtools coverage -a "${BED}" -b "${BAM_DIR}/${sample}.bam" -mean | awk '$4=="MDR1"{print $NF}')
    CRT=$(bedtools coverage -a "${BED}" -b "${BAM_DIR}/${sample}.bam" -mean | awk '$4=="CRT"{print $NF}')
    K13=$(bedtools coverage -a "${BED}" -b "${BAM_DIR}/${sample}.bam" -mean | awk '$4=="K13"{print $NF}')

    echo -e "${sample}\t${MSP2}\t${CSP}\t${MDR1}\t${CRT}\t${K13}" \
        >> "${COVERAGE_MATRIX}"

done

echo "Alignment completed."
