#!/usr/bin/env bash
# align_minimap2.sh — align reads to the 3D7 reference and summarise coverage
#
# Inputs:  ${MERGED_DIR}   per-sample FASTQ files named <SampleID>.fastq.gz (stage 01)
#          ${REF_GENOME}   P. falciparum 3D7 v3 reference, PlasmoDB release 68
#          ${AMPLICON_BED} BED file of the targeted amplicons. Column 4 must hold the
#          gene name; the coverage table takes its columns from this
# Outputs: ${BAM_DIR}/<SampleID>.bam(.bai)        sorted, indexed alignments
#          ${QC_DIR}/depth/<SampleID>.depth.txt   per-base depth within the amplicons
#          ${QC_DIR}/mapping/<SampleID>.flagstat.txt
#          ${QC_DIR}/mapping/<SampleID>.idxstats.txt
#          ${TABLE_DIR}/coverage_matrix.tsv       mean depth per sample per amplicon
#          ${TABLE_DIR}/mapping_summary.tsv       reads, mapped reads, mapped % per sample
#
# Usage:   bash scripts/02_alignment/align_minimap2.sh

set -euo pipefail

# --- load shared paths ---------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config.sh"

# --- check the tools are available ------------------------------------------------
for TOOL in minimap2 samtools bedtools; do
    if ! command -v "${TOOL}" >/dev/null 2>&1; then
        echo "ERROR: ${TOOL} not found on PATH." >&2
        echo "       conda env create -f environment.yml && conda activate pf-rwanda" >&2
        exit 1
    fi
done

# --- check the inputs exist --------------------------------------------------------
for FILE in "${REF_GENOME}" "${AMPLICON_BED}"; do
    if [ ! -f "${FILE}" ]; then
        echo "ERROR: required file not found: ${FILE}" >&2
        echo "       Check the paths set in config.sh." >&2
        exit 1
    fi
done

if [ ! -d "${MERGED_DIR}" ]; then
    echo "ERROR: FASTQ folder not found: ${MERGED_DIR}" >&2
    exit 1
fi

shopt -s nullglob
FASTQ_FILES=("${MERGED_DIR}"/*.fastq.gz)
shopt -u nullglob

if [ "${#FASTQ_FILES[@]}" -eq 0 ]; then
    echo "ERROR: no .fastq.gz files in ${MERGED_DIR}" >&2
    echo "       Run stage 01 first." >&2
    exit 1
fi

DEPTH_DIR="${QC_DIR}/depth"
MAPPING_DIR="${QC_DIR}/mapping"
mkdir -p "${BAM_DIR}" "${DEPTH_DIR}" "${MAPPING_DIR}" "${TABLE_DIR}" "${LOG_DIR}"

# --- build the minimap2 index once -------------------------------------------------
# Without this, minimap2 re-indexes the reference for every sample, which for a
# study of this size costs hours of pointless work.
REF_INDEX="${REF_GENOME%.fasta}.mmi"
if [ ! -f "${REF_INDEX}" ]; then
    echo "Building minimap2 index (first run only)..."
    minimap2 -x map-ont -d "${REF_INDEX}" "${REF_GENOME}"
fi

# --- read the amplicon names from the BED file --------------------------------------
# bedtools coverage reports intervals in the order they appear in the BED file,
# so the header built here and the values collected below always line up.
GENE_NAMES="$(awk 'NF >= 4 {printf "\t%s", $4}' "${AMPLICON_BED}")"
N_GENES="$(awk 'NF >= 4' "${AMPLICON_BED}" | wc -l | tr -d ' ')"

if [ "${N_GENES}" -eq 0 ]; then
    echo "ERROR: no usable intervals in ${AMPLICON_BED}" >&2
    echo "       Expected at least 4 columns, with the gene name in column 4." >&2
    exit 1
fi

COVERAGE_MATRIX="${TABLE_DIR}/coverage_matrix.tsv"
MAPPING_SUMMARY="${TABLE_DIR}/mapping_summary.tsv"
FAILURE_LOG="${LOG_DIR}/align_failures.txt"

printf "Sample%s\n" "${GENE_NAMES}"                           > "${COVERAGE_MATRIX}"
printf "Sample\ttotal_reads\tmapped_reads\tmapped_percent\n"  > "${MAPPING_SUMMARY}"
: > "${FAILURE_LOG}"

echo "============================================================"
echo "Reference  : $(basename "${REF_GENOME}")"
echo "Amplicons  : ${N_GENES} ($(echo "${GENE_NAMES}" | tr '\t' ' ' | sed 's/^ //'))"
echo "Samples    : ${#FASTQ_FILES[@]}"
echo "BAM output : ${BAM_DIR}"
echo "Threads    : ${THREADS}"
echo "============================================================"

# --- align each sample ---------------------------------------------------------------
INDEX=0
ALIGNED=0
SKIPPED=0
FAILED=0

for FASTQ in "${FASTQ_FILES[@]}"; do

    # Stage 01 names files <SampleID>.fastq.gz, so stripping .fastq.gz gives the
    # sample ID directly.
    SAMPLE="$(basename "${FASTQ}" .fastq.gz)"
    INDEX=$((INDEX + 1))
    BAM="${BAM_DIR}/${SAMPLE}.bam"

    if [ -f "${BAM}" ] && [ -f "${BAM}.bai" ]; then
        echo "  [${INDEX}/${#FASTQ_FILES[@]}] ${SAMPLE}: already aligned, skipping"
        SKIPPED=$((SKIPPED + 1))
    else
        echo "  [${INDEX}/${#FASTQ_FILES[@]}] ${SAMPLE}"

        # A single failed sample should not abort a run of over a thousand, so
        # errors are caught here and recorded rather than stopping the loop.
        if ! {
            minimap2 -t "${THREADS}" -ax map-ont --secondary=no -p 0.9 \
                     "${REF_INDEX}" "${FASTQ}" 2>/dev/null \
            | samtools view -b -F 4 - \
            | samtools sort -@ "${THREADS}" -o "${BAM}" - 2>/dev/null \
            && samtools index "${BAM}"
        }; then
            echo "    FAILED: ${SAMPLE}" >&2
            echo "${SAMPLE}" >> "${FAILURE_LOG}"
            rm -f "${BAM}" "${BAM}.bai"
            FAILED=$((FAILED + 1))
            continue
        fi
        ALIGNED=$((ALIGNED + 1))
    fi

    # --- per-sample statistics ---------------------------------------------------
    samtools flagstat "${BAM}" > "${MAPPING_DIR}/${SAMPLE}.flagstat.txt"
    samtools idxstats "${BAM}" > "${MAPPING_DIR}/${SAMPLE}.idxstats.txt"

    # Depth is restricted to the amplicon intervals. Using -a across the whole
    # 23 Mb genome would write mostly zeros, at a few hundred MB per sample.
    samtools depth -a -b "${AMPLICON_BED}" "${BAM}" \
        > "${DEPTH_DIR}/${SAMPLE}.depth.txt"

    # --- mean depth per amplicon -------------------------------------------------
    # One bedtools call covers every amplicon. The original script called
    # bedtools once per gene, re-reading the whole BAM each time.
    COVERAGE_ROW="$(bedtools coverage -a "${AMPLICON_BED}" -b "${BAM}" -mean \
                    | awk '{printf "\t%s", $NF}')"
    printf "%s%s\n" "${SAMPLE}" "${COVERAGE_ROW}" >> "${COVERAGE_MATRIX}"

    # --- mapped read counts --------------------------------------------------------
    awk -v s="${SAMPLE}" '
        /in total/        { total  = $1 }
        / mapped \(/      { if (mapped == "") mapped = $1 }
        END {
            pct = (total > 0) ? (100 * mapped / total) : 0
            printf "%s\t%s\t%s\t%.2f\n", s, total, mapped, pct
        }' "${MAPPING_DIR}/${SAMPLE}.flagstat.txt" >> "${MAPPING_SUMMARY}"

done

# --- summary ---------------------------------------------------------------------------
echo "============================================================"
echo "Aligned : ${ALIGNED}"
echo "Skipped : ${SKIPPED} (already had a BAM)"
echo "Failed  : ${FAILED}"
[ "${FAILED}" -gt 0 ] && echo "          see ${FAILURE_LOG}"
echo ""
echo "Coverage matrix : ${COVERAGE_MATRIX}"
echo "Mapping summary : ${MAPPING_SUMMARY}"
echo "Next            : bash scripts/03_variant_calling/run_clair3.sh"
echo "============================================================"
