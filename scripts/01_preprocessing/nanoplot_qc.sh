#!/usr/bin/env bash
# =============================================================================
# nanoplot_qc.sh — read-level quality control with NanoPlot
#
# WHAT IT DOES
#   Runs NanoPlot on the renamed per-sample FASTQ files to summarise read
#   length, read quality and yield, so that poor-quality samples are visible
#   before any time is spent on alignment and variant calling.
#
#   Two kinds of report are produced:
#     1. one report per sample, for spotting individual failures
#     2. one combined report across all samples, for the run as a whole
#
#   pycoQC was run alongside NanoPlot on the sequencing summary files; it is
#   not scripted here because it reads MinKNOW's own summary output rather
#   than the FASTQ files.
#
# INPUTS
#   ${MERGED_DIR}  per-sample FASTQ files, named <SampleID>.fastq.gz
#
# OUTPUTS
#   ${QC_DIR}/per_sample/<SampleID>/   one NanoPlot report per sample
#   ${QC_DIR}/combined/                one report across all samples
#   ${QC_DIR}/nanoplot_summary.tsv     one row per sample: reads, bases, N50, mean Q
#
# REQUIRES
#   NanoPlot 1.46.2
#
# USAGE
#   bash scripts/01_preprocessing/nanoplot_qc.sh
#
#   Per-sample reports only (faster):
#     bash scripts/01_preprocessing/nanoplot_qc.sh --no-combined
#
# Part of: Diallo et al., Nature Communications (2026)
#          https://github.com/mtd1997/Adaptive_Pf_resurgence_Rwanda
# =============================================================================

set -euo pipefail

# --- load shared paths --------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config.sh"

RUN_COMBINED=1
if [ "${1:-}" = "--no-combined" ]; then
    RUN_COMBINED=0
fi

# --- check NanoPlot is available -------------------------------------------------
if ! command -v NanoPlot >/dev/null 2>&1; then
    echo "ERROR: NanoPlot not found on PATH." >&2
    echo "       Install it with:  conda env create -f environment.yml" >&2
    echo "       then:             conda activate pf-rwanda" >&2
    exit 1
fi

NANOPLOT_VERSION="$(NanoPlot --version 2>&1 | head -1)"

# --- check inputs exist -----------------------------------------------------------
if [ ! -d "${MERGED_DIR}" ]; then
    echo "ERROR: FASTQ folder not found: ${MERGED_DIR}" >&2
    exit 1
fi

shopt -s nullglob
FASTQ_FILES=("${MERGED_DIR}"/*.fastq.gz)
shopt -u nullglob

if [ "${#FASTQ_FILES[@]}" -eq 0 ]; then
    echo "ERROR: no .fastq.gz files in ${MERGED_DIR}" >&2
    echo "       Run concatenate_fastq.sh and rename_fastq.R first." >&2
    exit 1
fi

mkdir -p "${QC_DIR}/per_sample"

echo "============================================================"
echo "NanoPlot      : ${NANOPLOT_VERSION}"
echo "Input folder  : ${MERGED_DIR}"
echo "Samples       : ${#FASTQ_FILES[@]}"
echo "Output folder : ${QC_DIR}"
echo "Threads       : ${THREADS}"
echo "============================================================"

# --- per-sample reports ------------------------------------------------------------
SUMMARY="${QC_DIR}/nanoplot_summary.tsv"
printf "SampleID\treads\tbases\tread_N50\tmean_qual\n" > "${SUMMARY}"

INDEX=0
for FASTQ in "${FASTQ_FILES[@]}"; do

    SAMPLE="$(basename "${FASTQ}" .fastq.gz)"
    INDEX=$((INDEX + 1))
    OUTDIR="${QC_DIR}/per_sample/${SAMPLE}"

    # skip samples already done, so an interrupted run can be resumed
    if [ -f "${OUTDIR}/NanoStats.txt" ]; then
        echo "  [${INDEX}/${#FASTQ_FILES[@]}] ${SAMPLE}: already done, skipping"
    else
        echo "  [${INDEX}/${#FASTQ_FILES[@]}] ${SAMPLE}"
        NanoPlot \
            --fastq "${FASTQ}" \
            --outdir "${OUTDIR}" \
            --prefix "${SAMPLE}_" \
            --threads "${THREADS}" \
            --plots hex dot \
            --N50 \
            --tsv_stats \
            > /dev/null 2>&1 || {
                echo "    WARNING: NanoPlot failed for ${SAMPLE}, continuing" >&2
                continue
            }
    fi

    # pull the headline numbers into one table for quick comparison
    STATS="${OUTDIR}/${SAMPLE}_NanoStats.txt"
    [ -f "${STATS}" ] || STATS="${OUTDIR}/NanoStats.txt"

    if [ -f "${STATS}" ]; then
        READS=$(awk -F: '/Number of reads/   {gsub(/[ ,]/,"",$2); print $2}' "${STATS}" | head -1)
        BASES=$(awk -F: '/Total bases/       {gsub(/[ ,]/,"",$2); print $2}' "${STATS}" | head -1)
        N50=$(  awk -F: '/Read length N50/   {gsub(/[ ,]/,"",$2); print $2}' "${STATS}" | head -1)
        MEANQ=$(awk -F: '/Mean read quality/ {gsub(/[ ,]/,"",$2); print $2}' "${STATS}" | head -1)
        printf "%s\t%s\t%s\t%s\t%s\n" \
            "${SAMPLE}" "${READS:-NA}" "${BASES:-NA}" "${N50:-NA}" "${MEANQ:-NA}" \
            >> "${SUMMARY}"
    fi

done

# --- combined report ----------------------------------------------------------------
if [ "${RUN_COMBINED}" -eq 1 ]; then
    echo ""
    echo "Building combined report across all samples..."
    NanoPlot \
        --fastq "${FASTQ_FILES[@]}" \
        --outdir "${QC_DIR}/combined" \
        --threads "${THREADS}" \
        --plots hex dot \
        --N50 \
        --tsv_stats \
        > /dev/null 2>&1 || \
        echo "  WARNING: combined report failed (too many input files?); per-sample reports are unaffected" >&2
fi

echo "============================================================"
echo "Per-sample reports : ${QC_DIR}/per_sample/"
echo "Summary table      : ${SUMMARY}"
[ "${RUN_COMBINED}" -eq 1 ] && echo "Combined report    : ${QC_DIR}/combined/"
echo "Next               : bash scripts/02_alignment/align_minimap2.sh"
echo "============================================================"
