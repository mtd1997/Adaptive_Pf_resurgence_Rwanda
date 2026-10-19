#!/usr/bin/env bash
# combine_and_genotype_gvcfs.sh — joint genotyping across all samples
#
# Inputs:  ${VCF_DIR}/<SampleID>_clair3/merge_output.gvcf.gz   from stage 03
#          ${REF_GENOME} and its .fai and .dict
#          ${AMPLICON_BED}  used to build the interval list, so genotyping is
#          restricted to the seven targeted loci
# Outputs: ${DATA_ROOT}/population_track/gvcf/<SampleID>.gvcf.gz   collected, indexed
#          ${DATA_ROOT}/population_track/sample_map.txt
#          ${DATA_ROOT}/population_track/intervals.list
#          ${DATA_ROOT}/population_track/genomicsdb/                GenomicsDB workspace
#          ${PROCESSED_DIR}/joint_genotyped.vcf.gz                  input to filtering
#
# Usage:   bash scripts/04_filtering/combine_and_genotype_gvcfs.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config.sh"

# --- checks ---------------------------------------------------------------------
for TOOL in gatk bcftools samtools; do
    command -v "${TOOL}" >/dev/null 2>&1 || {
        echo "ERROR: ${TOOL} not found on PATH." >&2; exit 1; }
done

[ -f "${REF_GENOME}" ] || { echo "ERROR: reference not found: ${REF_GENOME}" >&2; exit 1; }
[ -f "${AMPLICON_BED}" ] || { echo "ERROR: BED not found: ${AMPLICON_BED}" >&2; exit 1; }

# GATK needs a sequence dictionary as well as the .fai
if [ ! -f "${REF_GENOME}.fai" ]; then
    echo "Creating reference index..."
    samtools faidx "${REF_GENOME}"
fi
REF_DICT="${REF_GENOME%.fasta}.dict"
if [ ! -f "${REF_DICT}" ]; then
    echo "Creating sequence dictionary..."
    gatk CreateSequenceDictionary -R "${REF_GENOME}"
fi

POP_DIR="${DATA_ROOT}/population_track"
GVCF_COLLECTED="${POP_DIR}/gvcf"
SAMPLE_MAP="${POP_DIR}/sample_map.txt"
INTERVALS="${POP_DIR}/intervals.list"
DATABASE="${POP_DIR}/genomicsdb"
OUTPUT="${PROCESSED_DIR}/joint_genotyped.vcf.gz"

mkdir -p "${GVCF_COLLECTED}" "${PROCESSED_DIR}" "${LOG_DIR}"

# --- 1. collect and index the gVCFs ------------------------------------------------
# Clair3 writes every sample's gVCF under the same name, merge_output.gvcf.gz,
# inside a per-sample folder. They are copied out under their sample names so
# that GenomicsDBImport can refer to them individually.
echo "--- collecting gVCFs ---"

shopt -s nullglob
CLAIR3_DIRS=("${VCF_DIR}"/*_clair3)
shopt -u nullglob

[ "${#CLAIR3_DIRS[@]}" -gt 0 ] || {
    echo "ERROR: no *_clair3 folders in ${VCF_DIR}" >&2; exit 1; }

COLLECTED=0
MISSING=0
for CLAIR3_SUBDIR in "${CLAIR3_DIRS[@]}"; do
    SAMPLE="$(basename "${CLAIR3_SUBDIR}" _clair3)"
    SOURCE="${CLAIR3_SUBDIR}/merge_output.gvcf.gz"
    TARGET="${GVCF_COLLECTED}/${SAMPLE}.gvcf.gz"

    if [ ! -f "${SOURCE}" ]; then
        echo "  missing gVCF for ${SAMPLE}" >&2
        MISSING=$((MISSING + 1))
        continue
    fi

    [ -f "${TARGET}" ] || cp "${SOURCE}" "${TARGET}"
    # GenomicsDBImport requires an index for every input gVCF
    [ -f "${TARGET}.tbi" ] || bcftools index --tbi "${TARGET}"
    COLLECTED=$((COLLECTED + 1))
done

echo "  collected ${COLLECTED} gVCFs (${MISSING} missing)"
[ "${COLLECTED}" -gt 0 ] || { echo "ERROR: nothing to genotype." >&2; exit 1; }

# --- 2. build the sample map -----------------------------------------------------------
# Two tab-separated columns: sample name, then path to that sample's gVCF.
echo "--- building sample map ---"
: > "${SAMPLE_MAP}"
for GVCF in "${GVCF_COLLECTED}"/*.gvcf.gz; do
    SAMPLE="$(basename "${GVCF}" .gvcf.gz)"
    printf "%s\t%s\n" "${SAMPLE}" "$(cd "$(dirname "${GVCF}")" && pwd)/$(basename "${GVCF}")" \
        >> "${SAMPLE_MAP}"
done
echo "  $(wc -l < "${SAMPLE_MAP}" | tr -d ' ') samples in ${SAMPLE_MAP}"

# Verify the names in the map match the names inside the gVCF headers, before
# spending hours on an import that would fail at the end.
echo "--- verifying sample names against gVCF headers ---"
MISMATCHES=0
while IFS=$'\t' read -r NAME PATH_TO_GVCF; do
    HEADER_NAME="$(bcftools query -l "${PATH_TO_GVCF}" | head -1)"
    if [ "${NAME}" != "${HEADER_NAME}" ]; then
        echo "  MISMATCH: map says '${NAME}', gVCF header says '${HEADER_NAME}'" >&2
        MISMATCHES=$((MISMATCHES + 1))
    fi
done < "${SAMPLE_MAP}"

if [ "${MISMATCHES}" -gt 0 ]; then
    echo "ERROR: ${MISMATCHES} sample-name mismatches; GenomicsDBImport would fail." >&2
    echo "       Check --sample_name in stage 03." >&2
    exit 1
fi
echo "  all names match"

# --- 3. build the interval list from the amplicon BED ------------------------------------
# Genotyping is restricted to the seven targeted loci. Deriving the intervals
# from the same BED used for alignment keeps the two in step; a separate
# hand-written intervals.list would be one more file to drift out of date.
awk '{printf "%s:%d-%d\n", $1, $2 + 1, $3}' "${AMPLICON_BED}" > "${INTERVALS}"
echo "--- intervals ---"
sed 's/^/  /' "${INTERVALS}"

# --- 4. import into GenomicsDB -------------------------------------------------------------
# GenomicsDBImport refuses to write into an existing workspace, so a previous
# attempt must be cleared before retrying.
if [ -d "${DATABASE}" ]; then
    echo "--- removing previous GenomicsDB workspace ---"
    rm -rf "${DATABASE}"
fi

echo "--- GenomicsDBImport ---"
# --batch-size keeps the number of simultaneously open files bounded. Importing
# over a thousand samples at once can exhaust the open-file limit.
gatk GenomicsDBImport \
    --genomicsdb-workspace-path "${DATABASE}" \
    --sample-name-map "${SAMPLE_MAP}" \
    --reader-threads 5 \
    --batch-size 50 \
    -L "${INTERVALS}" \
    2>&1 | tee "${LOG_DIR}/genomicsdbimport.log"

# --- 5. joint genotyping ----------------------------------------------------------------------
echo "--- GenotypeGVCFs ---"
gatk GenotypeGVCFs \
    -R "${REF_GENOME}" \
    -V "gendb://${DATABASE}" \
    -O "${OUTPUT}" \
    -L "${INTERVALS}" \
    2>&1 | tee "${LOG_DIR}/genotypegvcfs.log"

bcftools index --force --tbi "${OUTPUT}"

echo "============================================================"
echo "Samples genotyped : $(bcftools query -l "${OUTPUT}" | wc -l | tr -d ' ')"
echo "Sites             : $(bcftools index --nrecords "${OUTPUT}")"
echo "Output            : ${OUTPUT}"
echo "Next              : bash scripts/04_filtering/filter_population_vcf.sh"
echo "============================================================"
