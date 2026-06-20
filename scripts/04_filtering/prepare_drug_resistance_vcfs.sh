#!/bin/bash

set -euo pipefail

# Prepare annotated VCF files for drug resistance marker analysis.

PROJECT_DIR="."

CLAIR3_DIR="${PROJECT_DIR}/clair3_output"
VCF_DIR="${PROJECT_DIR}/all_vcfs"
FILTERED_DIR="${PROJECT_DIR}/filtered_vcfs"
DRUG_DIR="${PROJECT_DIR}/drug_marker_vcfs"
ANNOTATED_DIR="${PROJECT_DIR}/annotated_vcfs"

REGIONS="${PROJECT_DIR}/regions.txt"
SNPEFF="${PROJECT_DIR}/snpEff/snpEff.jar"

mkdir -p "${VCF_DIR}"
mkdir -p "${FILTERED_DIR}"
mkdir -p "${DRUG_DIR}"
mkdir -p "${ANNOTATED_DIR}"

##################################################
# Collect Clair3 VCF files
##################################################

for sample in "${CLAIR3_DIR}"/*_clair3
do

    mv "${sample}/merge_output.vcf.gz" \
       "${VCF_DIR}/$(basename "${sample}").vcf.gz"

done

##################################################
# Filter variants and retain SNPs
##################################################

for file in "${VCF_DIR}"/*.vcf.gz
do

    sample=$(basename "${file}" .vcf.gz)

    bcftools view \
        -e 'QUAL <= 15 || FMT/DP < 10' \
        "${file}" \
    | bcftools view \
        --types snps \
        -Oz \
        -o "${FILTERED_DIR}/${sample}_snps.vcf.gz"

done

##################################################
# Extract drug resistance genes
##################################################

for file in "${FILTERED_DIR}"/*.vcf.gz
do

    sample=$(basename "${file}" _snps.vcf.gz)

    bcftools view \
        -R "${REGIONS}" \
        "${file}" \
        -o "${DRUG_DIR}/${sample}_DRG.vcf"

done

##################################################
# Annotate variants
##################################################

for file in "${DRUG_DIR}"/*.vcf
do

    sample=$(basename "${file}" .vcf)

    java -Xmx4g \
        -jar "${SNPEFF}" \
        -v Pf3D7 \
        "${file}" \
        -stats "${ANNOTATED_DIR}/${sample}.html" \
        > "${ANNOTATED_DIR}/${sample}.ann.vcf"

done

echo "Drug resistance VCF preparation completed."

