#!/usr/bin/env bash
# =============================================================================
# filter_population_vcf.sh — progressive missingness filtering of the
#                            jointly genotyped VCF
#
# WHAT IT DOES
#   Takes the jointly genotyped VCF across all seven targeted loci and filters
#   it down to the call set used for the population-structure analyses:
#   principal-coordinate analysis and the neighbour-joining tree (Figure 4).
#
#   Three steps:
#     1. genotype-level quality and depth
#     2. biallelic SNPs only
#     3. progressive missingness filtering, in seven rounds
#
# WHY PROGRESSIVE RATHER THAN A SINGLE THRESHOLD
#   Site missingness and sample missingness are entangled: a handful of badly
#   covered samples makes every site look poor, and a handful of poor sites
#   makes every sample look poor. Applying one strict threshold to either in
#   isolation therefore discards data that is in fact usable.
#
#   This script instead tightens the threshold in stages, from 0.7 down to 0.1,
#   and within each stage filters sites first, then samples. Early rounds strip
#   out the worst offenders on both axes, which improves the statistics seen by
#   the next round. By the final round at 0.1, the remaining sites and samples
#   are genuinely well covered rather than merely surviving a single cut.
#
#   Each round writes its own output and a summary, so the trade-off between
#   stringency and data retained is visible rather than assumed.
#
# INPUTS
#   ${PROCESSED_DIR}/joint_genotyped.vcf.gz    from combine_and_genotype_gvcfs.sh
#
# OUTPUTS
#   ${DATA_ROOT}/population_track/filtering/step1_quality.vcf.gz
#   ${DATA_ROOT}/population_track/filtering/step2_biallelic.vcf.gz
#   ${DATA_ROOT}/population_track/filtering/filtered_<threshold>.recode.vcf.gz
#   ${DATA_ROOT}/population_track/filtering/summary_<threshold>.txt
#   ${DATA_ROOT}/population_track/filtering/filtering_summary.tsv
#   ${PROCESSED_DIR}/population_filtered.vcf.gz   the final round, copied out
#
# REQUIRES
#   bcftools 1.23, vcftools 0.1.16, tabix/bgzip (htslib)
#
# USAGE
#   bash scripts/04_filtering/filter_population_vcf.sh
#
# DOWNSTREAM
#   This output feeds two different paths:
#     - annotate_population_vcf.sh, then population_genomics.R   (Figure 4)
#     - filter_for_selection_scan.sh, a second and stricter pass (Figure 5)
#
# NOTE ON SAMPLE COUNTS
#   This script does not remove the three isolates later found to be of
#   non-Rwandan origin. They are excluded at the R stage, when genotypes are
#   joined to the sample metadata, taking the analysed set to 1,176.
#
# Part of: Diallo et al., Nature Communications (2026)
#          https://github.com/mtd1997/Adaptive_Pf_resurgence_Rwanda
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config.sh"

# --- filtering parameters, as used for the published analysis ------------------
MIN_DP=3          # FORMAT/DP, minimum genotype depth
MIN_GQ=5          # FORMAT/GQ, minimum genotype quality
THRESHOLDS=(0.7 0.6 0.5 0.4 0.3 0.2 0.1)   # missingness, loosest to strictest

# --- checks ---------------------------------------------------------------------
for TOOL in bcftools vcftools bgzip tabix; do
    command -v "${TOOL}" >/dev/null 2>&1 || {
        echo "ERROR: ${TOOL} not found on PATH." >&2
        echo "       conda env create -f environment.yml && conda activate pf-rwanda" >&2
        exit 1; }
done

INPUT_VCF="${PROCESSED_DIR}/joint_genotyped.vcf.gz"
[ -f "${INPUT_VCF}" ] || {
    echo "ERROR: input VCF not found: ${INPUT_VCF}" >&2
    echo "       Run combine_and_genotype_gvcfs.sh first." >&2; exit 1; }

WORKDIR="${DATA_ROOT}/population_track/filtering"
mkdir -p "${WORKDIR}" "${PROCESSED_DIR}"
cd "${WORKDIR}"

SUMMARY_TABLE="${WORKDIR}/filtering_summary.tsv"
printf "stage\tvariants\tsamples\n" > "${SUMMARY_TABLE}"

count_vcf() {
    # prints "<variants>\t<samples>" for a VCF
    local VCF="$1"
    printf "%s\t%s" \
        "$(bcftools view -H "${VCF}" | wc -l | tr -d ' ')" \
        "$(bcftools query -l "${VCF}" | wc -l | tr -d ' ')"
}

echo "============================================================"
echo "Input      : ${INPUT_VCF}"
printf "Starting   : %s variants, %s samples\n" $(count_vcf "${INPUT_VCF}")
echo "Depth      : FORMAT/DP >= ${MIN_DP}"
echo "Quality    : FORMAT/GQ >= ${MIN_GQ}"
echo "Thresholds : ${THRESHOLDS[*]}"
echo "============================================================"

printf "input\t%s\n" "$(count_vcf "${INPUT_VCF}")" >> "${SUMMARY_TABLE}"

# --- step 1: genotype quality and depth -------------------------------------------
echo ""
echo "--- step 1: genotype quality and depth ---"
bcftools filter -i "FORMAT/DP>=${MIN_DP} && FORMAT/GQ>=${MIN_GQ}" "${INPUT_VCF}" \
    -Oz -o step1_quality.vcf.gz
tabix -f -p vcf step1_quality.vcf.gz
printf "step1_quality\t%s\n" "$(count_vcf step1_quality.vcf.gz)" >> "${SUMMARY_TABLE}"
printf "    %s variants, %s samples\n" $(count_vcf step1_quality.vcf.gz)

# --- step 2: biallelic SNPs only ---------------------------------------------------
# Multiallelic sites and indels are excluded: the distance measures used for
# the PCoA and the tree assume biallelic markers.
echo ""
echo "--- step 2: biallelic SNPs ---"
bcftools view -v snps -m2 -M2 step1_quality.vcf.gz -Oz -o step2_biallelic.vcf.gz
tabix -f -p vcf step2_biallelic.vcf.gz
printf "step2_biallelic\t%s\n" "$(count_vcf step2_biallelic.vcf.gz)" >> "${SUMMARY_TABLE}"
printf "    %s variants, %s samples\n" $(count_vcf step2_biallelic.vcf.gz)

# --- step 3: progressive missingness -------------------------------------------------
CURRENT_VCF="step2_biallelic.vcf.gz"

for THRESHOLD in "${THRESHOLDS[@]}"; do

    echo ""
    echo "--- missingness threshold ${THRESHOLD} ---"

    # ---- sites ----
    # vcftools writes .lmiss with F_MISS, the fraction of missing genotypes, in
    # column 6. `NR > 1` skips the header row; without it the header is parsed
    # as a data line and a bogus entry reaches the keep list.
    vcftools --gzvcf "${CURRENT_VCF}" --missing-site --out tmp_site 2>/dev/null
    awk -v t="${THRESHOLD}" 'NR > 1 && $6 < t {print $1"\t"$2}' tmp_site.lmiss > keep_sites.txt
    echo "    sites kept  : $(wc -l < keep_sites.txt | tr -d ' ')"

    vcftools --gzvcf "${CURRENT_VCF}" \
        --positions keep_sites.txt \
        --recode --recode-INFO-all \
        --out tmp_sites_filtered 2>/dev/null
    bgzip -f tmp_sites_filtered.recode.vcf
    tabix -f -p vcf tmp_sites_filtered.recode.vcf.gz

    # ---- samples ----
    # .imiss holds F_MISS in column 5.
    vcftools --gzvcf tmp_sites_filtered.recode.vcf.gz --missing-indv --out tmp_ind 2>/dev/null
    awk -v t="${THRESHOLD}" 'NR > 1 && $5 < t {print $1}' tmp_ind.imiss > keep_samples.txt
    echo "    samples kept: $(wc -l < keep_samples.txt | tr -d ' ')"

    vcftools --gzvcf tmp_sites_filtered.recode.vcf.gz \
        --keep keep_samples.txt \
        --recode --recode-INFO-all \
        --out "filtered_${THRESHOLD}" 2>/dev/null
    bgzip -f "filtered_${THRESHOLD}.recode.vcf"
    tabix -f -p vcf "filtered_${THRESHOLD}.recode.vcf.gz"

    # ---- per-round summary ----
    {
        echo "Threshold ${THRESHOLD}"
        echo "Variants:"
        bcftools view -H "filtered_${THRESHOLD}.recode.vcf.gz" | wc -l
        echo "Samples:"
        bcftools query -l "filtered_${THRESHOLD}.recode.vcf.gz" | wc -l
    } > "summary_${THRESHOLD}.txt"

    printf "missingness_%s\t%s\n" "${THRESHOLD}" \
        "$(count_vcf "filtered_${THRESHOLD}.recode.vcf.gz")" >> "${SUMMARY_TABLE}"
    printf "    -> %s variants, %s samples\n" $(count_vcf "filtered_${THRESHOLD}.recode.vcf.gz")

    CURRENT_VCF="filtered_${THRESHOLD}.recode.vcf.gz"

done

# --- publish the final round --------------------------------------------------------------
FINAL_VCF="${PROCESSED_DIR}/population_filtered.vcf.gz"
cp "${CURRENT_VCF}" "${FINAL_VCF}"
bcftools index --force --tbi "${FINAL_VCF}"

rm -f tmp_site.* tmp_ind.* tmp_sites_filtered.* keep_sites.txt keep_samples.txt

echo ""
echo "============================================================"
echo "Filtering summary:"
column -t -s$'\t' "${SUMMARY_TABLE}" | sed 's/^/  /'
echo ""
echo "Final output : ${FINAL_VCF}"
echo "Next         : bash scripts/04_filtering/annotate_population_vcf.sh   (Figure 4)"
echo "               bash scripts/04_filtering/filter_for_selection_scan.sh (Figure 5)"
echo "============================================================"
