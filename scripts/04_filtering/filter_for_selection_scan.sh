#!/usr/bin/env bash
# filter_for_selection_scan.sh — second, stricter filtering pass for the allele-frequency and FST analyses
#
# Inputs:  ${PROCESSED_DIR}/population_filtered.vcf.gz   from filter_population_vcf.sh
# Outputs: ${DATA_ROOT}/selection_track/step1..step4.vcf.gz   intermediates
#          ${PROCESSED_DIR}/selection_scan_filtered.vcf.gz     final call set
#
# Usage:   bash scripts/04_filtering/filter_for_selection_scan.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config.sh"

# --- filtering parameters, as used for the published analysis -------------------
MIN_DP=10              # FORMAT/DP, minimum genotype depth
MIN_QUAL=30            # QUAL, minimum site-level variant quality
SNP_MISSING_FIRST=0.2  # drop SNPs missing in more than this fraction of samples
SAMPLE_MISSING=0.2     # drop samples missing more than this fraction of SNPs
SNP_MISSING_FINAL=0.1  # final, stricter SNP missingness limit

# A minor-allele-frequency filter was considered but not applied, so that rare
# variants showing large temporal shifts remain visible to the scan:
# APPLY_MAF=0 ; MIN_MAF=0.01

for TOOL in bcftools vcftools; do
    command -v "${TOOL}" >/dev/null 2>&1 || {
        echo "ERROR: ${TOOL} not found on PATH." >&2; exit 1; }
done

INPUT_VCF="${PROCESSED_DIR}/population_filtered.vcf.gz"
[ -f "${INPUT_VCF}" ] || {
    echo "ERROR: input VCF not found: ${INPUT_VCF}" >&2
    echo "       Run filter_population_vcf.sh first." >&2; exit 1; }

WORKDIR="${DATA_ROOT}/selection_track"
mkdir -p "${WORKDIR}" "${PROCESSED_DIR}"
cd "${WORKDIR}"

count_vcf() {
    local VCF="$1"
    printf "%s variants, %s samples" \
        "$(bcftools view -H "${VCF}" | wc -l | tr -d ' ')" \
        "$(bcftools query -l "${VCF}" | wc -l | tr -d ' ')"
}

echo "============================================================"
echo "Input    : ${INPUT_VCF}"
echo "Starting : $(count_vcf "${INPUT_VCF}")"
echo "============================================================"

bcftools index --force "${INPUT_VCF}"

# --- 1. depth and site quality ------------------------------------------------------
echo ""
echo "--- 1. FORMAT/DP >= ${MIN_DP} and QUAL >= ${MIN_QUAL} ---"
bcftools view -i "FORMAT/DP>=${MIN_DP} && QUAL>=${MIN_QUAL}" "${INPUT_VCF}" \
    -Oz -o step1.vcf.gz
bcftools index --force step1.vcf.gz
echo "    $(count_vcf step1.vcf.gz)"

# --- 2. first SNP missingness pass ------------------------------------------------------
echo ""
echo "--- 2. SNP missingness < ${SNP_MISSING_FIRST} ---"
bcftools view -i "F_MISSING<${SNP_MISSING_FIRST}" step1.vcf.gz -Oz -o step2.vcf.gz
bcftools index --force step2.vcf.gz
echo "    $(count_vcf step2.vcf.gz)"

# --- 3. sample missingness ------------------------------------------------------------------
# Samples are assessed only after the worst SNPs have gone, so that a sample is
# not condemned by sites that are about to be dropped anyway.
echo ""
echo "--- 3. dropping samples missing more than ${SAMPLE_MISSING} ---"
vcftools --gzvcf step2.vcf.gz --missing-indv --out miss 2>/dev/null
awk -v t="${SAMPLE_MISSING}" 'NR > 1 && $5 > t {print $1}' miss.imiss > bad_samples.txt
echo "    dropping $(wc -l < bad_samples.txt | tr -d ' ') samples"

bcftools view -S ^bad_samples.txt step2.vcf.gz --force-samples -Oz -o step3.vcf.gz
bcftools index --force step3.vcf.gz
echo "    $(count_vcf step3.vcf.gz)"

# --- 4. second, stricter SNP missingness pass --------------------------------------------------
# Repeated after the sample removal, because dropping poorly covered samples
# improves the missingness of every remaining site.
echo ""
echo "--- 4. SNP missingness < ${SNP_MISSING_FINAL} ---"
bcftools view -i "F_MISSING<${SNP_MISSING_FINAL}" step3.vcf.gz -Oz -o step4.vcf.gz
bcftools index --force step4.vcf.gz
echo "    $(count_vcf step4.vcf.gz)"

# --- 5. biallelic SNPs only ------------------------------------------------------------------------
echo ""
echo "--- 5. biallelic SNPs ---"
FINAL_VCF="${PROCESSED_DIR}/selection_scan_filtered.vcf.gz"
bcftools view -m2 -M2 -v snps step4.vcf.gz -Oz -o "${FINAL_VCF}"
bcftools index --force --tbi "${FINAL_VCF}"
echo "    $(count_vcf "${FINAL_VCF}")"

echo ""
echo "============================================================"
echo "Final output : ${FINAL_VCF}"
echo "Samples dropped at step 3 are listed in ${WORKDIR}/bad_samples.txt"
echo "Next         : Rscript scripts/05_downstream_analysis/selection_scan.R"
echo "============================================================"
