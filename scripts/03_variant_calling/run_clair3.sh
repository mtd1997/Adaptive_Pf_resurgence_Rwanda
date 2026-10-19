#!/usr/bin/env bash
# run_clair3.sh — haploid variant calling with Clair3
#
# Inputs:  ${BAM_DIR}      sorted, indexed BAMs from stage 02
#          ${REF_GENOME}   P. falciparum 3D7 v3 reference, PlasmoDB release 68
# Outputs: ${VCF_DIR}/<SampleID>_clair3/merge_output.vcf.gz
#          ${VCF_DIR}/<SampleID>_clair3/merge_output.gvcf.gz
#          ${LOG_DIR}/clair3_failures.txt   samples that failed, if any
#
# Usage:   bash scripts/03_variant_calling/run_clair3.sh

set -euo pipefail

# --- load shared paths ----------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config.sh"

# --- Docker needs absolute paths for -v mounts -------------------------------------
# A relative path such as ./bam_files is silently misread by Docker, so every
# mounted location is resolved to an absolute path first.
abs_dir() { cd "$1" 2>/dev/null && pwd; }
abs_file() { printf "%s/%s\n" "$(cd "$(dirname "$1")" && pwd)" "$(basename "$1")"; }

# --- checks -------------------------------------------------------------------------
if ! command -v docker >/dev/null 2>&1; then
    echo "ERROR: docker not found on PATH." >&2
    echo "       Clair3 is run from its official container image." >&2
    exit 1
fi

if ! docker info >/dev/null 2>&1; then
    echo "ERROR: the Docker daemon is not running." >&2
    echo "       Start Docker Desktop and try again." >&2
    exit 1
fi

if [ ! -f "${REF_GENOME}" ]; then
    echo "ERROR: reference genome not found: ${REF_GENOME}" >&2
    exit 1
fi

if [ ! -f "${REF_GENOME}.fai" ]; then
    echo "ERROR: reference index not found: ${REF_GENOME}.fai" >&2
    echo "       Clair3 requires it. Create it with:" >&2
    echo "         samtools faidx ${REF_GENOME}" >&2
    exit 1
fi

if [ ! -d "${BAM_DIR}" ]; then
    echo "ERROR: BAM folder not found: ${BAM_DIR}" >&2
    echo "       Run stage 02 first." >&2
    exit 1
fi

shopt -s nullglob
BAM_FILES=("${BAM_DIR}"/*.bam)
shopt -u nullglob

if [ "${#BAM_FILES[@]}" -eq 0 ]; then
    echo "ERROR: no BAM files in ${BAM_DIR}" >&2
    exit 1
fi

# --- warn if the model does not match the basecalling model --------------------------
MODEL_NAME="$(basename "${CLAIR3_MODEL}")"
case "${MODEL_NAME}" in
    *_sup_*) : ;;   # matches this study's SUP basecalling
    *_hac_*)
        echo "WARNING: CLAIR3_MODEL is a HAC model (${MODEL_NAME}), but the reads" >&2
        echo "         were re-basecalled with SUP. See the note in this script's" >&2
        echo "         header. Continuing in 5 seconds; Ctrl+C to stop." >&2
        sleep 5
        ;;
    *)
        echo "NOTE: could not tell from its name whether ${MODEL_NAME} matches" >&2
        echo "      the SUP basecalling used in this study." >&2
        ;;
esac

mkdir -p "${VCF_DIR}" "${LOG_DIR}"

REF_ABS="$(abs_file "${REF_GENOME}")"
REF_DIR_ABS="$(dirname "${REF_ABS}")"
BAM_DIR_ABS="$(abs_dir "${BAM_DIR}")"
VCF_DIR_ABS="$(abs_dir "${VCF_DIR}")"

# This study used a rerio model, which lives outside the container and so must
# be mounted into it. A model bundled in the image needs no mount.
MODEL_MOUNT=()
if [ -n "${CLAIR3_MODEL_DIR:-}" ]; then
    MODEL_DIR_ABS="$(abs_dir "${CLAIR3_MODEL_DIR}")"
    if [ -z "${MODEL_DIR_ABS}" ]; then
        echo "ERROR: CLAIR3_MODEL_DIR not found: ${CLAIR3_MODEL_DIR}" >&2
        echo "       Download the rerio models (see this script's header)." >&2
        exit 1
    fi
    MODEL_PATH="${MODEL_DIR_ABS}/${MODEL_NAME}"
    if [ ! -d "${MODEL_PATH}" ]; then
        echo "ERROR: model not found: ${MODEL_PATH}" >&2
        echo "       Available in ${MODEL_DIR_ABS}:" >&2
        ls -1 "${MODEL_DIR_ABS}" 2>/dev/null | sed 's/^/         /' >&2
        exit 1
    fi
    MODEL_MOUNT=(-v "${MODEL_DIR_ABS}:${MODEL_DIR_ABS}")
else
    MODEL_PATH="${CLAIR3_MODEL}"
fi

FAILURE_LOG="${LOG_DIR}/clair3_failures.txt"
: > "${FAILURE_LOG}"

echo "============================================================"
echo "Clair3 image : ${CLAIR3_IMAGE}"
echo "Model        : ${MODEL_PATH}"
echo "Reference    : $(basename "${REF_ABS}")"
echo "Samples      : ${#BAM_FILES[@]}"
echo "Output       : ${VCF_DIR_ABS}"
echo "Threads      : ${THREADS}"
echo "============================================================"

# --- call variants per sample -----------------------------------------------------------
INDEX=0
CALLED=0
SKIPPED=0
FAILED=0

for BAM in "${BAM_FILES[@]}"; do

    SAMPLE="$(basename "${BAM}" .bam)"
    INDEX=$((INDEX + 1))
    OUTDIR="${VCF_DIR_ABS}/${SAMPLE}_clair3"

    # Completion is judged by the output file, not by the folder existing. An
    # interrupted run leaves a folder behind, and skipping on the folder alone
    # would silently pass over a sample that never finished.
    if [ -f "${OUTDIR}/merge_output.vcf.gz" ]; then
        echo "  [${INDEX}/${#BAM_FILES[@]}] ${SAMPLE}: already called, skipping"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    rm -rf "${OUTDIR}"
    mkdir -p "${OUTDIR}"

    echo "  [${INDEX}/${#BAM_FILES[@]}] ${SAMPLE}"

    if ! docker run --rm \
        -v "${BAM_DIR_ABS}:${BAM_DIR_ABS}" \
        -v "${VCF_DIR_ABS}:${VCF_DIR_ABS}" \
        -v "${REF_DIR_ABS}:${REF_DIR_ABS}" \
        "${MODEL_MOUNT[@]}" \
        "${CLAIR3_IMAGE}" \
        /opt/bin/run_clair3.sh \
            --bam_fn="${BAM_DIR_ABS}/${SAMPLE}.bam" \
            --ref_fn="${REF_ABS}" \
            --threads="${THREADS}" \
            --platform="ont" \
            --model_path="${MODEL_PATH}" \
            --output="${OUTDIR}" \
            --sample_name="${SAMPLE}" \
            --haploid_sensitive \
            --include_all_ctgs \
            --gvcf \
            --var_pct_full=1 \
            --ref_pct_full=1 \
            --no_phasing_for_fa \
            --use_whatshap_for_final_output_phasing \
            --remove_intermediate_dir \
        > "${LOG_DIR}/clair3_${SAMPLE}.log" 2>&1
    then
        echo "    FAILED: ${SAMPLE} (see ${LOG_DIR}/clair3_${SAMPLE}.log)" >&2
        echo "${SAMPLE}" >> "${FAILURE_LOG}"
        FAILED=$((FAILED + 1))
        continue
    fi

    CALLED=$((CALLED + 1))

done

# --- summary -------------------------------------------------------------------------------
echo "============================================================"
echo "Called  : ${CALLED}"
echo "Skipped : ${SKIPPED} (already complete)"
echo "Failed  : ${FAILED}"
[ "${FAILED}" -gt 0 ] && echo "          see ${FAILURE_LOG}"
echo ""
echo "Next    : bash scripts/04_filtering/prepare_drug_resistance_vcfs.sh"
echo "          bash scripts/04_filtering/combine_and_genotype_gvcfs.sh"
echo "============================================================"
