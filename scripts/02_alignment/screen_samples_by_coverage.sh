#!/usr/bin/env bash
# screen_samples_by_coverage.sh — sample eligibility screen on amplicon coverage
#
# Inputs:  ${BAM_DIR}       sorted, indexed BAMs from align_minimap2.sh
#          ${AMPLICON_BED}  the seven targeted amplicons
# Outputs: ${TABLE_DIR}/amplicon_read_counts.tsv   reads per sample per amplicon
#          ${TABLE_DIR}/amplicon_mean_depth.tsv    mean depth per sample per amplicon
#          ${TABLE_DIR}/samples_passing.txt        one sample ID per line
#          ${TABLE_DIR}/samples_failing.tsv        excluded samples and which amplicons failed
#
# Usage:   bash scripts/02_alignment/screen_samples_by_coverage.sh
#          bash scripts/02_alignment/screen_samples_by_coverage.sh --min-reads 50

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config.sh"

MIN_READS=50
if [ "${1:-}" = "--min-reads" ] && [ -n "${2:-}" ]; then
    MIN_READS="$2"
fi

for TOOL in samtools bedtools; do
    command -v "${TOOL}" >/dev/null 2>&1 || {
        echo "ERROR: ${TOOL} not found on PATH." >&2; exit 1; }
done

[ -f "${AMPLICON_BED}" ] || { echo "ERROR: BED not found: ${AMPLICON_BED}" >&2; exit 1; }
[ -d "${BAM_DIR}" ]      || { echo "ERROR: BAM folder not found: ${BAM_DIR}" >&2; exit 1; }

shopt -s nullglob
BAM_FILES=("${BAM_DIR}"/*.bam)
shopt -u nullglob
[ "${#BAM_FILES[@]}" -gt 0 ] || { echo "ERROR: no BAMs in ${BAM_DIR}" >&2; exit 1; }

mkdir -p "${TABLE_DIR}"

GENE_NAMES="$(awk 'NF >= 4 {printf "\t%s", $4}' "${AMPLICON_BED}")"
N_GENES="$(awk 'NF >= 4' "${AMPLICON_BED}" | wc -l | tr -d ' ')"

COUNTS="${TABLE_DIR}/amplicon_read_counts.tsv"
DEPTHS="${TABLE_DIR}/amplicon_mean_depth.tsv"
PASSING="${TABLE_DIR}/samples_passing.txt"
FAILING="${TABLE_DIR}/samples_failing.tsv"

printf "Sample%s\n" "${GENE_NAMES}" > "${COUNTS}"
printf "Sample%s\n" "${GENE_NAMES}" > "${DEPTHS}"
: > "${PASSING}"
printf "Sample\tfailing_amplicons\tread_counts\n" > "${FAILING}"

echo "============================================================"
echo "Samples      : ${#BAM_FILES[@]}"
echo "Amplicons    : ${N_GENES}"
echo "Criterion    : >= ${MIN_READS} reads in EVERY amplicon"
echo "============================================================"

INDEX=0
PASSED=0
FAILED=0

for BAM in "${BAM_FILES[@]}"; do

    SAMPLE="$(basename "${BAM}" .bam)"
    INDEX=$((INDEX + 1))

    # bedtools coverage reports intervals in BED order. Column 4 of the output
    # beyond the BED fields is the count of overlapping reads; -mean instead
    # gives mean depth. Both are collected, so the screen can be stated either
    # as a read count or as a depth.
    COUNT_ROW="$(bedtools coverage -a "${AMPLICON_BED}" -b "${BAM}" -counts \
                 | awk '{printf "\t%s", $NF}')"
    DEPTH_ROW="$(bedtools coverage -a "${AMPLICON_BED}" -b "${BAM}" -mean \
                 | awk '{printf "\t%.1f", $NF}')"

    printf "%s%s\n" "${SAMPLE}" "${COUNT_ROW}" >> "${COUNTS}"
    printf "%s%s\n" "${SAMPLE}" "${DEPTH_ROW}" >> "${DEPTHS}"

    # Which amplicons, if any, fall below the threshold?
    FAILED_GENES="$(paste \
        <(awk 'NF >= 4 {print $4}' "${AMPLICON_BED}") \
        <(printf "%s" "${COUNT_ROW}" | tr '\t' '\n' | sed '/^$/d') \
        | awk -v m="${MIN_READS}" '$2 < m {printf "%s,", $1}' | sed 's/,$//')"

    if [ -z "${FAILED_GENES}" ]; then
        echo "${SAMPLE}" >> "${PASSING}"
        PASSED=$((PASSED + 1))
    else
        printf "%s\t%s\t%s\n" "${SAMPLE}" "${FAILED_GENES}" \
            "$(printf "%s" "${COUNT_ROW}" | sed 's/^\t//' | tr '\t' ',')" >> "${FAILING}"
        FAILED=$((FAILED + 1))
    fi

    [ $((INDEX % 100)) -eq 0 ] && echo "  ...${INDEX}/${#BAM_FILES[@]}"

done

echo "============================================================"
echo "Passing : ${PASSED}  -> ${PASSING}"
echo "Failing : ${FAILED}  -> ${FAILING}"
echo ""
echo "Which amplicon most often caused exclusion:"
awk -F'\t' 'NR > 1 {n = split($2, a, ","); for (i = 1; i <= n; i++) c[a[i]]++}
            END {for (g in c) printf "  %-6s %d\n", g, c[g]}' "${FAILING}" | sort -k2 -rn
echo ""
echo "VALIDATION: compare ${PASSED} with the sample count carried into joint"
echo "            genotyping for the published analysis. See the provenance"
echo "            note at the top of this script."
echo "============================================================"
