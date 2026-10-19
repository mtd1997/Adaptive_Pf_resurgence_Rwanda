#!/usr/bin/env bash
# run_all_figures.sh — regenerate every figure and report what reproduced
#
# Inputs:  the paths set in config.R
#          $1  optional: the saved Kigali tree, needed by Figure 4d
#
# Outputs: results/figures/..., results/tables/...
#          results/logs/<script>.log   full output of each run
#
# Usage:   bash run_all_figures.sh
#          bash run_all_figures.sh /path/to/..._Kigali_tree.nwk

set -uo pipefail

cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

TREE_REFERENCE="${1:-${FIGURE04D_REFERENCE_TREE:-}}"
LOG_DIR="results/logs"
mkdir -p "${LOG_DIR}"

# Order matters where one script hands a file to another:
#   figure01 and figureS1b write panels that figureS1_assemble reads
#   figure04_ibs_distance writes the matrix figure04_ibs_pcoa reads
#   figure05_statistics writes the tables figure05_panels reads
SCRIPTS=(
  figure01_sampling_maps
  figureS1b_malaria_cases
  figureS1_assemble
  figure02a_coverage
  figure02b_gene_mutation_trends
  figure02c_resistance_categories
  figure03a_k13_lollipop
  figure03b_k13_category_pies
  figure03d_k13_district_map
  figure04_ibs_distance
  figure04_ibs_pcoa
  figure04d_distance_tree
  figure05_selection_scan_statistics
  figure05_selection_scan_panels
  figureS3_k13_comparisons
)

declare -a NAMES STATUS VERDICT
failures=0

for name in "${SCRIPTS[@]}"; do

  script="scripts/06_figures/${name}.R"
  log="${LOG_DIR}/${name}.log"

  args=()
  if [[ "${name}" == "figure04d_distance_tree" ]]; then
    if [[ -n "${TREE_REFERENCE}" && -f "${TREE_REFERENCE}" ]]; then
      args=(--reference "${TREE_REFERENCE}")
    else
      NAMES+=("${name}")
      STATUS+=("skipped")
      VERDICT+=("no reference tree given")
      printf '  %-38s skipped\n' "${name}"
      continue
    fi
  fi

  printf '  %-38s ' "${name}"

  if Rscript "${script}" "${args[@]}" > "${log}" 2>&1; then
    state="ok"
  else
    state="FAILED"
    failures=$((failures + 1))
  fi

  # Pull the strongest self-check line from the log. Patterns are tried in
  # order of strength, so a weak "PASS" never masks a definitive verdict.
  verdict=""
  for pattern in \
    'DOES NOT MATCH THE PUBLISHED PANEL' \
    'do NOT match published' \
    'does NOT match' \
    'IS REPRODUCED' \
    'Every count matches published[^.]*' \
    'Every per-year total matches[^.]*' \
    'All four group sizes match[^.]*' \
    'Figure 3b is reproduced' \
    'Using the [a-z]+ tree' \
    'the two implementations agree[^.]*' \
    'Every category matches[^.]*' \
    'PASS —[^.]*' \
    'FAIL —[^.]*'
  do
    verdict="$(grep -aoE "${pattern}" "${log}" | head -1)"
    [[ -n "${verdict}" ]] && break
  done
  [[ -z "${verdict}" ]] && verdict="(no self-check)"

  NAMES+=("${name}")
  STATUS+=("${state}")
  VERDICT+=("${verdict}")

  printf '%s\n' "${state}"
done

echo
printf '%s\n' "========================================================================"
printf '%s\n' "SUMMARY"
printf '%s\n' "========================================================================"
printf '%-38s %-9s %s\n' "script" "status" "self-check"
for i in "${!NAMES[@]}"; do
  printf '%-38s %-9s %s\n' "${NAMES[$i]}" "${STATUS[$i]}" "${VERDICT[$i]}"
done

echo
if (( failures == 0 )); then
  echo "All scripts completed. Logs in ${LOG_DIR}/"
else
  echo "${failures} script(s) failed. See ${LOG_DIR}/ for the full output."
fi

exit "${failures}"
