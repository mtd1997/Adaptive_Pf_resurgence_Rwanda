#!/usr/bin/env bash
################################################################################
# PER-AMPLICON COVERAGE FROM THE ALIGNMENTS
################################################################################
#
# WHAT IT DOES
#   Measures mean sequencing depth over each of the seven amplicons, for every
#   BAM file, and writes one long-format table:
#
#       sample    gene    coverage
#       RM197     MSP2    1250.18
#       RM197     CSP     515.07
#       ...
#
#   This is the input to scripts/06_figures/figure02a_coverage.R.
#
# WHY REGENERATE IT
#   The coverage table used for the original Figure 2a was produced early in
#   the project and its sample roster no longer matches the analysed set: it
#   held a laboratory control, some non-study identifiers and a batch of
#   placeholder-named samples, and was at the same time missing samples that
#   are in the final analysis. The difference could not be reconstructed.
#
#   Measuring coverage again from the current alignments gives a table whose
#   every row can be accounted for. It does not re-align anything, so it is a
#   single pass over the BAMs rather than a reprocessing.
#
# DECIMAL SEPARATOR
#   This script writes a period, e.g. 1250.18. The original table used a
#   comma, because it was generated under a locale that formats numbers that
#   way. The figure script accepts either, so both old and new tables work.
#
# RESUMING
#   Samples already present in the output are skipped, so the script can be
#   interrupted and restarted. Pass --force to measure everything again.
#
# HOW TO RUN
#   From the top of the repository:
#
#     bash scripts/02_alignment/compute_amplicon_coverage.sh
#
#   Options:
#     --bam-dir <dir>   where the BAM files are; defaults to BAM_DIR in
#                       config.sh
#     --out <file>      output table; defaults to
#                       data/processed/coverage_all_samples.tsv
#     --force           recompute every sample, ignoring existing output
#
################################################################################

set -euo pipefail


# ------------------------------------------------------------------------------
# Load the shared paths.
# ------------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

if [[ ! -f "${PROJECT_ROOT}/config.sh" ]]; then
    echo "ERROR: config.sh not found at ${PROJECT_ROOT}/config.sh" >&2
    exit 1
fi

# shellcheck source=/dev/null
source "${PROJECT_ROOT}/config.sh"


# ------------------------------------------------------------------------------
# Arguments
# ------------------------------------------------------------------------------

BAM_DIRECTORY="${BAM_DIR:-}"
OUTPUT_FILE="${PROJECT_ROOT}/data/processed/coverage_all_samples.tsv"
FORCE=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bam-dir)
            BAM_DIRECTORY="$2"
            shift 2
            ;;
        --out)
            OUTPUT_FILE="$2"
            shift 2
            ;;
        --force)
            FORCE=1
            shift
            ;;
        *)
            echo "ERROR: unknown option: $1" >&2
            echo "       Valid options: --bam-dir, --out, --force" >&2
            exit 1
            ;;
    esac
done


# ------------------------------------------------------------------------------
# Check everything before starting, and say exactly what is wrong.
# ------------------------------------------------------------------------------

if [[ -z "${BAM_DIRECTORY}" ]]; then
    echo "ERROR: no BAM directory. Either set BAM_DIR in config.sh or pass" >&2
    echo "       --bam-dir /path/to/bams" >&2
    exit 1
fi

if [[ ! -d "${BAM_DIRECTORY}" ]]; then
    echo "ERROR: BAM directory not found: ${BAM_DIRECTORY}" >&2
    exit 1
fi

AMPLICON_BED="${AMPLICON_BED:-${PROJECT_ROOT}/data/metadata/amplicons.bed}"

if [[ ! -f "${AMPLICON_BED}" ]]; then
    echo "ERROR: amplicon BED file not found: ${AMPLICON_BED}" >&2
    exit 1
fi

if ! command -v bedtools >/dev/null 2>&1; then
    echo "ERROR: bedtools was not found on PATH." >&2
    exit 1
fi


# ------------------------------------------------------------------------------
# The BED file must be tab-separated with the gene name in column 4. Spaces
# instead of tabs would silently put the coordinates in the wrong columns, so
# this is checked rather than assumed.
# ------------------------------------------------------------------------------

BED_COLUMNS="$(awk -F'\t' 'NF > 0 {print NF; exit}' "${AMPLICON_BED}")"

if [[ "${BED_COLUMNS}" -lt 4 ]]; then
    echo "ERROR: ${AMPLICON_BED} has ${BED_COLUMNS} tab-separated column(s)." >&2
    echo "       It needs at least four, with the gene name in column 4." >&2
    echo "       Check for spaces instead of tabs:" >&2
    echo "         cat -A ${AMPLICON_BED} | head -2" >&2
    echo "       Tabs appear as ^I." >&2
    exit 1
fi

AMPLICON_COUNT="$(awk -F'\t' 'NF >= 4' "${AMPLICON_BED}" | wc -l | tr -d ' ')"


# ------------------------------------------------------------------------------
# Find the BAM files.
# ------------------------------------------------------------------------------

BAM_LIST="$(mktemp)"
trap 'rm -f "${BAM_LIST}"' EXIT

find "${BAM_DIRECTORY}" -maxdepth 1 -name '*.bam' -type f \
    | sort > "${BAM_LIST}"

BAM_COUNT="$(wc -l < "${BAM_LIST}" | tr -d ' ')"

if [[ "${BAM_COUNT}" -eq 0 ]]; then
    echo "ERROR: no .bam files found in ${BAM_DIRECTORY}" >&2
    exit 1
fi

mkdir -p "$(dirname "${OUTPUT_FILE}")"

echo "BAM directory:  ${BAM_DIRECTORY}"
echo "BAM files:      ${BAM_COUNT}"
echo "Amplicon BED:   ${AMPLICON_BED}"
echo "Amplicons:      ${AMPLICON_COUNT}"
echo "Output:         ${OUTPUT_FILE}"
echo


# ------------------------------------------------------------------------------
# Work out which samples are already done, so an interrupted run resumes.
# ------------------------------------------------------------------------------

DONE_LIST="$(mktemp)"
trap 'rm -f "${BAM_LIST}" "${DONE_LIST}"' EXIT

if [[ "${FORCE}" -eq 1 ]] || [[ ! -s "${OUTPUT_FILE}" ]]; then

    printf 'sample\tgene\tcoverage\n' > "${OUTPUT_FILE}"
    : > "${DONE_LIST}"

    if [[ "${FORCE}" -eq 1 ]]; then
        echo "--force given: measuring every sample again."
        echo
    fi

else

    # A sample counts as done only if it has a row for every amplicon.
    # A sample interrupted part-way through would otherwise be left with an
    # incomplete set of rows and skipped as if it were finished.
    awk -F'\t' -v expected="${AMPLICON_COUNT}" '
        NR > 1 { count[$1]++ }
        END { for (sample in count) if (count[sample] == expected) print sample }
    ' "${OUTPUT_FILE}" | sort > "${DONE_LIST}"

    DONE_COUNT="$(wc -l < "${DONE_LIST}" | tr -d ' ')"

    # Remove the rows belonging to incomplete samples before they are
    # measured again. Without this, a retried sample's new rows are appended
    # on top of its stale ones and it ends up duplicated in the table --
    # which would plot it twice and distort the medians.

    ROWS_BEFORE="$(($(wc -l < "${OUTPUT_FILE}" | tr -d ' ') - 1))"

    CLEANED="$(mktemp)"

    awk -F'\t' '
        NR == FNR { complete[$0] = 1; next }
        FNR == 1  { print; next }
        ($1 in complete)
    ' "${DONE_LIST}" "${OUTPUT_FILE}" > "${CLEANED}"

    mv "${CLEANED}" "${OUTPUT_FILE}"

    ROWS_AFTER="$(($(wc -l < "${OUTPUT_FILE}" | tr -d ' ') - 1))"
    STALE_ROWS="$((ROWS_BEFORE - ROWS_AFTER))"

    if [[ "${DONE_COUNT}" -gt 0 ]]; then
        echo "Resuming: ${DONE_COUNT} sample(s) already measured will be skipped."
        echo "Pass --force to measure everything again."
    fi

    if [[ "${STALE_ROWS}" -gt 0 ]]; then
        echo "Discarded ${STALE_ROWS} row(s) from sample(s) that were measured"
        echo "only part-way through; those samples will be measured again."
    fi

    if [[ "${DONE_COUNT}" -gt 0 ]] || [[ "${STALE_ROWS}" -gt 0 ]]; then
        echo
    fi
fi


# ------------------------------------------------------------------------------
# Measure.
#
# bedtools coverage -mean reports the mean depth over each BED interval, one
# line per interval, in the BED's order. Column 4 is the gene name and the
# last column is the mean.
# ------------------------------------------------------------------------------

PROCESSED=0
SKIPPED=0
FAILED=0
FAILED_SAMPLES=""

while IFS= read -r BAM_PATH; do

    # Strip the .bam extension, and the .targets suffix the alignment step
    # adds, so the sample name matches the metadata and the VCF headers.
    SAMPLE="$(basename "${BAM_PATH}" .bam)"
    SAMPLE="${SAMPLE%.targets}"

    if grep -qxF "${SAMPLE}" "${DONE_LIST}" 2>/dev/null; then
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    PROCESSED=$((PROCESSED + 1))

    printf '\r  [%d/%d] %-40s' \
        "$((PROCESSED + SKIPPED))" "${BAM_COUNT}" "${SAMPLE}"

    # Write to a temporary file first. If bedtools fails part-way, the output
    # table is not left with a half-written sample in it.
    SAMPLE_OUTPUT="$(mktemp)"

    if bedtools coverage \
            -a "${AMPLICON_BED}" \
            -b "${BAM_PATH}" \
            -mean \
            2>/dev/null \
        | awk -F'\t' -v sample="${SAMPLE}" \
            'NF >= 5 { printf "%s\t%s\t%s\n", sample, $4, $NF }' \
            > "${SAMPLE_OUTPUT}"
    then

        LINES_WRITTEN="$(wc -l < "${SAMPLE_OUTPUT}" | tr -d ' ')"

        if [[ "${LINES_WRITTEN}" -eq "${AMPLICON_COUNT}" ]]; then
            cat "${SAMPLE_OUTPUT}" >> "${OUTPUT_FILE}"
        else
            FAILED=$((FAILED + 1))
            FAILED_SAMPLES="${FAILED_SAMPLES}
  ${SAMPLE}  (got ${LINES_WRITTEN} amplicon(s), expected ${AMPLICON_COUNT})"
        fi

    else
        FAILED=$((FAILED + 1))
        FAILED_SAMPLES="${FAILED_SAMPLES}
  ${SAMPLE}  (bedtools failed)"
    fi

    rm -f "${SAMPLE_OUTPUT}"

done < "${BAM_LIST}"

printf '\r%-60s\r' ''


# ------------------------------------------------------------------------------
# Report.
# ------------------------------------------------------------------------------

TOTAL_ROWS="$(($(wc -l < "${OUTPUT_FILE}" | tr -d ' ') - 1))"

TOTAL_SAMPLES="$(awk -F'\t' 'NR > 1 {print $1}' "${OUTPUT_FILE}" \
    | sort -u | wc -l | tr -d ' ')"

echo "============================================================"
echo "COVERAGE MEASURED"
echo "============================================================"
echo "  measured now:     ${PROCESSED}"
echo "  skipped:          ${SKIPPED}"
echo "  failed:           ${FAILED}"
echo
echo "  samples in table: ${TOTAL_SAMPLES}"
echo "  rows in table:    ${TOTAL_ROWS}"
echo
echo "  ${OUTPUT_FILE}"

if [[ "${FAILED}" -gt 0 ]]; then
    echo
    echo "These samples did not produce a complete set of amplicons:"
    echo "${FAILED_SAMPLES}"
    echo
    echo "Re-running will retry them, since only complete samples are skipped."
fi

echo
echo "Per-gene medians:"

# Arrays are indexed as values[gene SUBSEP position] rather than
# values[gene][position]. The two-bracket form is a GNU awk extension and
# fails on the BSD awk that macOS ships.

awk -F'\t' '
    NR > 1 && $3 != "" {
        value = $3
        gsub(/,/, ".", value)
        count[$2]++
        values[$2 SUBSEP count[$2]] = value + 0
    }
    END {
        for (gene in count) {
            n = count[gene]

            # Insertion sort. This runs once per gene over the samples, so
            # it is fast enough and needs no external sort.
            for (i = 2; i <= n; i++) {
                key = values[gene SUBSEP i]
                j = i - 1
                while (j > 0 && values[gene SUBSEP j] > key) {
                    values[gene SUBSEP (j + 1)] = values[gene SUBSEP j]
                    j--
                }
                values[gene SUBSEP (j + 1)] = key
            }

            if (n % 2) {
                median = values[gene SUBSEP ((n + 1) / 2)]
            } else {
                median = (values[gene SUBSEP (n / 2)] \
                        + values[gene SUBSEP (n / 2 + 1)]) / 2
            }

            printf "  %-8s n=%-6d median=%.1f\n", gene, n, median
        }
    }
' "${OUTPUT_FILE}" | sort

echo
echo "Next step — draw the figure:"
echo "  Rscript scripts/06_figures/figure02a_coverage.R \\"
echo "      --coverage ${OUTPUT_FILE}"
