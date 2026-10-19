#!/usr/bin/env bash
# annotate_population_vcf.sh — add functional annotations to the population VCF
#
# Inputs:  ${PROCESSED_DIR}/population_filtered.vcf.gz   from filter_population_vcf.sh
#          A SnpEff database for P. falciparum 3D7, named by SNPEFF_DB in config.sh
# Outputs: ${PROCESSED_DIR}/final_population.ann.vcf       annotated call set
#          ${RESULTS_DIR}/qc/final_population_snpeff.html  SnpEff summary report
#          ${RESULTS_DIR}/qc/final_population_snpeff.genes.txt
#
# Usage:   bash scripts/04_filtering/annotate_population_vcf.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config.sh"

# --- checks ---------------------------------------------------------------------
command -v java >/dev/null 2>&1 || {
    echo "ERROR: java not found on PATH; SnpEff needs it." >&2; exit 1; }

[ -f "${SNPEFF_JAR}" ] || {
    echo "ERROR: SnpEff jar not found: ${SNPEFF_JAR}" >&2
    echo "       Set SNPEFF_JAR in config.sh." >&2; exit 1; }

INPUT="${PROCESSED_DIR}/population_filtered.vcf.gz"
[ -f "${INPUT}" ] || {
    echo "ERROR: input VCF not found: ${INPUT}" >&2
    echo "       Run filter_population_vcf.sh first." >&2; exit 1; }

# Fail early and clearly if the database is missing, rather than letting SnpEff
# annotate everything as intergenic — which it does silently when the genome
# name is unknown, producing a file that looks fine but carries no gene calls.
SNPEFF_DIR="$(dirname "${SNPEFF_JAR}")"
if [ -d "${SNPEFF_DIR}/data" ] && [ ! -d "${SNPEFF_DIR}/data/${SNPEFF_DB}" ]; then
    echo "ERROR: SnpEff database '${SNPEFF_DB}' not found in ${SNPEFF_DIR}/data/" >&2
    echo "       Databases present:" >&2
    ls -1 "${SNPEFF_DIR}/data" 2>/dev/null | sed 's/^/         /' >&2
    echo "       See the note on building the database in this script's header." >&2
    exit 1
fi

OUTPUT="${PROCESSED_DIR}/final_population.ann.vcf"
STATS="${QC_DIR}/final_population_snpeff.html"
mkdir -p "${PROCESSED_DIR}" "${QC_DIR}"

echo "============================================================"
echo "Input    : ${INPUT}"
echo "Database : ${SNPEFF_DB}"
echo "Output   : ${OUTPUT}"
echo "============================================================"

java -Xmx4g -jar "${SNPEFF_JAR}" \
    -v "${SNPEFF_DB}" \
    -stats "${STATS}" \
    "${INPUT}" \
    > "${OUTPUT}"

# --- sanity check on the result ------------------------------------------------------
# An annotated file in which nothing carries an ANN field, or in which every
# variant is intergenic, means the database did not match the reference.
TOTAL=$(grep -vc "^#" "${OUTPUT}" || true)
ANNOTATED=$(grep -v "^#" "${OUTPUT}" | grep -c "ANN=" || true)

echo ""
echo "Variants           : ${TOTAL}"
echo "Carrying ANN       : ${ANNOTATED}"

if [ "${ANNOTATED}" -eq 0 ]; then
    echo "" >&2
    echo "WARNING: no variant received an annotation." >&2
    echo "         The database name probably does not match the reference." >&2
fi

echo ""
echo "Genes represented:"
grep -v "^#" "${OUTPUT}" \
    | grep -o 'ANN=[^;]*' \
    | cut -d'|' -f4 \
    | sort | uniq -c | sort -rn | head -10 | sed 's/^/  /'

echo ""
echo "SnpEff report : ${STATS}"
echo "Next          : Rscript scripts/05_downstream_analysis/population_genomics.R"
echo "============================================================"
