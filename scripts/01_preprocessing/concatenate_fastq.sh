#!/usr/bin/env bash
# =============================================================================
# concatenate_fastq.sh — merge per-barcode FASTQ files from one sequencing run
#
# WHAT IT DOES
#   MinKNOW writes many small FASTQ files per barcode (one every few minutes
#   during the run). This merges each barcode's files into a single FASTQ,
#   tagged with the run identifier.
#
#   The run tag matters: barcode numbers restart at 01 on every run, so
#   barcode01 from run 1 and barcode01 from run 2 are different samples. Merging
#   without a run tag would silently overwrite one with the other.
#
# INPUTS
#   $1  RUN_ID   short identifier for the run, e.g. run1
#   $2  RUN_DIR  folder holding the demultiplexed barcode folders from stage 00,
#                i.e. the folder that contains barcode01/, barcode02/, ...
#
# OUTPUTS
#   ${MERGED_DIR}/<barcode>_<RUN_ID>.fastq.gz   one file per barcode
#   These names must match the `barcode` column of the sample mapping table
#   used by rename_fastq.R in the next step.
#
# REQUIRES
#   bash, gzip (both standard)
#
# USAGE
#   bash scripts/01_preprocessing/concatenate_fastq.sh run1 /path/to/run1/fastq_pass
#
#   Run once per sequencing run, changing RUN_ID each time.
#
# NOTE ON `cat` WITH .gz FILES
#   Concatenating gzip files with `cat` is valid: the gzip format allows
#   streams to be joined end to end, and decompressors read the result as one
#   continuous file. No decompression or recompression is needed.
#
# Part of: Diallo et al., Nature Communications (2026)
#          https://github.com/mtd1997/Adaptive_Pf_resurgence_Rwanda
# =============================================================================

set -euo pipefail

# --- load shared paths -------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config.sh"

# --- arguments ----------------------------------------------------------------
if [ "$#" -ne 2 ]; then
    echo "ERROR: wrong number of arguments." >&2
    echo "" >&2
    echo "Usage: bash $0 <RUN_ID> <RUN_DIR>" >&2
    echo "  RUN_ID   short run identifier, e.g. run1" >&2
    echo "  RUN_DIR  folder containing barcode01/, barcode02/, ..." >&2
    echo "" >&2
    echo "Example:" >&2
    echo "  bash $0 run1 \${DATA_ROOT}/sequencing/run1/fastq_pass" >&2
    exit 1
fi

RUN_ID="$1"
RUN_DIR="$2"

# --- check the inputs exist before doing any work -------------------------------
if [ ! -d "${RUN_DIR}" ]; then
    echo "ERROR: run folder not found: ${RUN_DIR}" >&2
    exit 1
fi

shopt -s nullglob
BARCODE_DIRS=("${RUN_DIR}"/barcode*)
shopt -u nullglob

if [ "${#BARCODE_DIRS[@]}" -eq 0 ]; then
    echo "ERROR: no barcode* folders found in ${RUN_DIR}" >&2
    echo "       Expected the demultiplexed output of stage 00." >&2
    exit 1
fi

mkdir -p "${MERGED_DIR}"

echo "============================================================"
echo "Run identifier : ${RUN_ID}"
echo "Input folder   : ${RUN_DIR}"
echo "Output folder  : ${MERGED_DIR}"
echo "Barcodes found : ${#BARCODE_DIRS[@]}"
echo "============================================================"

# --- merge each barcode ----------------------------------------------------------
MERGED_COUNT=0
SKIPPED_COUNT=0

for BARCODE_DIR in "${BARCODE_DIRS[@]}"; do

    BARCODE="$(basename "${BARCODE_DIR}")"
    OUTPUT_FILE="${MERGED_DIR}/${BARCODE}_${RUN_ID}.fastq.gz"

    # count the input files for this barcode
    shopt -s nullglob
    PARTS=("${BARCODE_DIR}"/*.fastq.gz)
    shopt -u nullglob

    if [ "${#PARTS[@]}" -eq 0 ]; then
        echo "  ${BARCODE}: no FASTQ files, skipped"
        SKIPPED_COUNT=$((SKIPPED_COUNT + 1))
        continue
    fi

    # refuse to overwrite: a clash means the same RUN_ID was used twice
    if [ -f "${OUTPUT_FILE}" ]; then
        echo "ERROR: ${OUTPUT_FILE} already exists." >&2
        echo "       Has this run already been processed, or is RUN_ID reused?" >&2
        echo "       Delete the file or choose a different RUN_ID." >&2
        exit 1
    fi

    cat "${PARTS[@]}" > "${OUTPUT_FILE}"

    SIZE="$(du -h "${OUTPUT_FILE}" | cut -f1)"
    echo "  ${BARCODE}: ${#PARTS[@]} files -> ${BARCODE}_${RUN_ID}.fastq.gz (${SIZE})"
    MERGED_COUNT=$((MERGED_COUNT + 1))

done

echo "============================================================"
echo "Merged : ${MERGED_COUNT} barcodes"
echo "Skipped: ${SKIPPED_COUNT} empty barcodes"
echo "Next   : Rscript scripts/01_preprocessing/rename_fastq.R"
echo "============================================================"
