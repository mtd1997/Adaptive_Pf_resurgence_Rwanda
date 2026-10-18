#!/usr/bin/env bash
# =============================================================================
# prepare_drug_resistance_vcfs.sh — filter and annotate VCFs for the
#                                   drug-resistance marker analysis
#
# WHAT IT DOES
#   For each sample, takes the per-sample VCF from Clair3 and prepares it for
#   drug-resistance genotyping, in four steps:
#
#     1. collect   copy the Clair3 VCF out to a working folder
#     2. filter    keep high-quality SNPs (QUAL > 15, depth > 10)
#     3. subset    restrict to the drug-resistance regions
#     4. annotate  add functional consequences with SnpEff
#
#   This is the drug-resistance track. The population-genetics track starts
#   from the gVCFs instead; see combine_and_genotype_gvcfs.sh.
#
# INPUTS
#   ${VCF_DIR}/<SampleID>_clair3/merge_output.vcf.gz   from stage 03
#   ${REGIONS_FILE}                                     drug-resistance regions
#   SnpEff database Pf3D7 (see scripts/04_filtering/README.md)
#
# OUTPUTS
#   ${DATA_ROOT}/drug_track/01_collected/<SampleID>.vcf.gz
#   ${DATA_ROOT}/drug_track/02_filtered/<SampleID>_snps.vcf.gz
#   ${DATA_ROOT}/drug_track/03_regions/<SampleID>_DRG.vcf.gz
#   ${ANNOTATED_DIR}/<SampleID>.ann.vcf      input to stage 05
#   ${LOG_DIR}/drug_vcf_failures.txt
#
# REQUIRES
#   bcftools 1.23, SnpEff 5.2f, Java
#
# USAGE
#   bash scripts/04_filtering/prepare_drug_resistance_vcfs.sh
#
#   Samples already annotated are skipped, so the script can be re-run safely.
#
# -----------------------------------------------------------------------------
# A NOTE ON REFERENCE CALLS                           <<< CHECK THIS SETTING
#
#   Estimating wild-type frequency requires distinguishing two situations that
#   look identical in an ordinary VCF: a position confidently called as
#   reference, and a position with no data at all. Clair3's output records
#   confident reference calls, but `bcftools view --types snps` discards them,
#   because a reference call is not a SNP.
#
#   KEEP_REF_CALLS=1  keeps reference calls alongside SNPs. This matches the
#                     Methods, which state that reference genotype calls were
#                     preserved to allow estimation of wild-type frequencies.
#   KEEP_REF_CALLS=0  keeps SNPs only, as in the original script.
#
#   The two give different denominators for every wild-type frequency in the
#   paper, so this must be set to whichever was actually used.
# -----------------------------------------------------------------------------
#
# Part of: Diallo et al., Nature Communications (2026)
#          https://github.com/mtd1997/Adaptive_Pf_resurgence_Rwanda
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config.sh"

KEEP_REF_CALLS="${KEEP_REF_CALLS:-1}"

# --- checks -------------------------------------------------------------------
for TOOL in bcftools java; do
    command -v "${TOOL}" >/dev/null 2>&1 || {
        echo "ERROR: ${TOOL} not found on PATH." >&2; exit 1; }
done

[ -f "${SNPEFF_JAR}" ] || {
    echo "ERROR: SnpEff jar not found: ${SNPEFF_JAR}" >&2
    echo "       Set SNPEFF_JAR in config.sh." >&2; exit 1; }

[ -f "${REGIONS_FILE}" ] || {
    echo "ERROR: regions file not found: ${REGIONS_FILE}" >&2; exit 1; }

shopt -s nullglob
CLAIR3_DIRS=("${VCF_DIR}"/*_clair3)
shopt -u nullglob

[ "${#CLAIR3_DIRS[@]}" -gt 0 ] || {
    echo "ERROR: no *_clair3 folders in ${VCF_DIR}" >&2
    echo "       Run stage 03 first." >&2; exit 1; }

COLLECTED_DIR="${DATA_ROOT}/drug_track/01_collected"
FILTERED_DIR="${DATA_ROOT}/drug_track/02_filtered"
DRUG_DIR="${DATA_ROOT}/drug_track/03_regions"
mkdir -p "${COLLECTED_DIR}" "${FILTERED_DIR}" "${DRUG_DIR}" "${ANNOTATED_DIR}" "${LOG_DIR}"

FAILURE_LOG="${LOG_DIR}/drug_vcf_failures.txt"
: > "${FAILURE_LOG}"

if [ "${KEEP_REF_CALLS}" -eq 1 ]; then
    TYPE_FILTER="--types snps,ref"
    TYPE_NOTE="SNPs and confident reference calls"
else
    TYPE_FILTER="--types snps"
    TYPE_NOTE="SNPs only"
fi

echo "============================================================"
echo "Samples       : ${#CLAIR3_DIRS[@]}"
echo "Retaining     : ${TYPE_NOTE}"
echo "Quality filter: QUAL > 15 and FMT/DP > 10"
echo "Regions       : ${REGIONS_FILE}"
echo "Annotated out : ${ANNOTATED_DIR}"
echo "============================================================"

INDEX=0
DONE=0
SKIPPED=0
FAILED=0

for CLAIR3_SUBDIR in "${CLAIR3_DIRS[@]}"; do

    SAMPLE="$(basename "${CLAIR3_SUBDIR}" _clair3)"
    INDEX=$((INDEX + 1))

    SOURCE_VCF="${CLAIR3_SUBDIR}/merge_output.vcf.gz"
    COLLECTED="${COLLECTED_DIR}/${SAMPLE}.vcf.gz"
    FILTERED="${FILTERED_DIR}/${SAMPLE}_snps.vcf.gz"
    REGIONAL="${DRUG_DIR}/${SAMPLE}_DRG.vcf.gz"
    ANNOTATED="${ANNOTATED_DIR}/${SAMPLE}.ann.vcf"

    if [ -s "${ANNOTATED}" ]; then
        echo "  [${INDEX}/${#CLAIR3_DIRS[@]}] ${SAMPLE}: already annotated, skipping"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    if [ ! -f "${SOURCE_VCF}" ]; then
        echo "    MISSING: ${SOURCE_VCF}" >&2
        echo "${SAMPLE}	missing merge_output.vcf.gz" >> "${FAILURE_LOG}"
        FAILED=$((FAILED + 1))
        continue
    fi

    echo "  [${INDEX}/${#CLAIR3_DIRS[@]}] ${SAMPLE}"

    if ! {
        # 1. Collect. Copied, not moved: moving would destroy the Clair3
        #    output, break this script on a second run, and defeat the
        #    completion check that lets stage 03 resume.
        cp "${SOURCE_VCF}" "${COLLECTED}"
        bcftools index --force --tbi "${COLLECTED}"

        # 2. Filter on quality and depth.
        #    QUAL > 15 and FMT/DP > 10, as stated in the Methods. Note the
        #    original used 'FMT/DP < 10' as the exclusion, which keeps DP = 10;
        #    'FMT/DP <= 10' is what "depth > 10" means.
        bcftools view -e 'QUAL <= 15 || FMT/DP <= 10' "${COLLECTED}" -Ou \
        | bcftools view ${TYPE_FILTER} -Oz -o "${FILTERED}"
        bcftools index --force --tbi "${FILTERED}"

        # 3. Restrict to the drug-resistance regions.
        #    -R needs an index on the input, which is why each step above
        #    indexes what it writes.
        bcftools view -R "${REGIONS_FILE}" "${FILTERED}" -Oz -o "${REGIONAL}"
        bcftools index --force --tbi "${REGIONAL}"

        # 4. Annotate.
        java -Xmx4g -jar "${SNPEFF_JAR}" \
            -v "${SNPEFF_DB}" \
            -stats "${ANNOTATED_DIR}/${SAMPLE}_snpeff.html" \
            "${REGIONAL}" \
            > "${ANNOTATED}"
    } 2>> "${LOG_DIR}/drug_vcf_${SAMPLE}.log"; then
        echo "    FAILED: ${SAMPLE} (see ${LOG_DIR}/drug_vcf_${SAMPLE}.log)" >&2
        echo "${SAMPLE}	processing error" >> "${FAILURE_LOG}"
        rm -f "${ANNOTATED}"
        FAILED=$((FAILED + 1))
        continue
    fi

    DONE=$((DONE + 1))

done

echo "============================================================"
echo "Annotated : ${DONE}"
echo "Skipped   : ${SKIPPED}"
echo "Failed    : ${FAILED}"
[ "${FAILED}" -gt 0 ] && echo "            see ${FAILURE_LOG}"
echo ""
echo "Next      : python3 scripts/05_downstream_analysis/extract_genotypes.py"
echo "============================================================"
