#!/usr/bin/env bash
################################################################################
# CAPTURE SOFTWARE VERSIONS
################################################################################
#
# WHAT IT DOES
#   Asks every tool on this machine for its version and writes the answers to
#   environment/software_versions.tsv.
#
# WHY NOT JUST WRITE THE TABLE BY HAND
#   A hand-written version list drifts away from reality, and the error is
#   invisible: nothing breaks, the Methods section is simply wrong. Reading
#   the versions off the machine that ran the analysis removes the chance to
#   mistype one.
#
#   The file it writes is also what the Methods software list should be
#   checked against before submission.
#
# WHAT IT DOES WITH TOOLS IT CANNOT FIND
#   It records them as "not found" rather than leaving a gap. A tool that is
#   missing from this machine but was used for the analysis needs its version
#   filled in by hand, and "not found" makes that obvious. Add a comment line
#   beside it saying where the version came from.
#
# HOW TO RUN
#   From the top of the repository:
#
#     bash environment/capture_versions.sh
#
#   Then read the file it writes, fill in anything marked "not found" or
#   "fill in by hand", and commit it.
#
################################################################################

set -uo pipefail
# Note: -e is deliberately omitted. Several of the version commands below
# exit non-zero even when they print a version, and one missing tool should
# not abandon the whole capture.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

OUTPUT_FILE="${SCRIPT_DIR}/software_versions.tsv"

if [[ -f "${PROJECT_ROOT}/config.sh" ]]; then
    # shellcheck source=/dev/null
    source "${PROJECT_ROOT}/config.sh"
fi


# ------------------------------------------------------------------------------
# Pull the version out of a command's output, reading it from stdin.
#
# Tools disagree about where they print the version, how they label it, and
# whether they use stdout or stderr, so both streams are searched.
#
# Simply taking the first version-shaped string does not work. Java prints a
# JAVA_TOOL_OPTIONS preamble containing a proxy address, and taking the first
# match returns "127.0.0.1" as the Java version. Other tools print an install
# path or a URL first.
#
# So the version is looked for in context before falling back:
#
#   1. version "21.0.12.1"   quoted, as Java and several JVM tools print it
#   2. version 4.6.1.0       labelled but unquoted
#   3. 2.30-r1287            bare, for tools that print the number alone
#
# Environment-variable dumps are dropped first, since they are never the
# version and frequently contain something that looks like one.
# ------------------------------------------------------------------------------

VERSION_PATTERN='[0-9]+\.[0-9]+(\.[0-9]+)*[0-9a-zA-Z._+-]*'

extract_version() {

    local text
    text="$(cat)"

    text="$(printf '%s\n' "${text}" | grep -v 'Picked up .*_OPTIONS' || true)"

    local candidate

    # 1. A quoted version following the word "version".
    candidate="$(
        printf '%s\n' "${text}" \
            | grep -oiE 'version[[:space:]]*"[^"]+"' \
            | head -1 \
            | grep -oE "${VERSION_PATTERN}" \
            | head -1
    )"

    if [[ -n "${candidate}" ]]; then
        printf '%s' "${candidate}"
        return
    fi

    # 2. An unquoted version following the word "version", with an optional
    #    colon and an optional leading "v".
    candidate="$(
        printf '%s\n' "${text}" \
            | grep -oiE "version:?[[:space:]]*v?${VERSION_PATTERN}" \
            | head -1 \
            | grep -oE "${VERSION_PATTERN}" \
            | head -1
    )"

    if [[ -n "${candidate}" ]]; then
        printf '%s' "${candidate}"
        return
    fi

    # 3. The first version-shaped string anywhere.
    printf '%s\n' "${text}" \
        | grep -oE "${VERSION_PATTERN}" \
        | head -1
}

report() {
    local name="$1"
    local version="$2"
    local source_note="$3"

    printf '%s\t%s\t%s\n' "${name}" "${version}" "${source_note}" \
        >> "${OUTPUT_FILE}"

    printf '  %-22s %-24s %s\n' "${name}" "${version}" "${source_note}"
}

# Ask a command for its version, or report that it is absent.
#   $1 tool label for the table
#   $2 executable name
#   $3.. the arguments that make it print a version
capture() {
    local label="$1"
    local executable="$2"
    shift 2

    if ! command -v "${executable}" >/dev/null 2>&1; then
        report "${label}" "not found" "not installed on this machine"
        return
    fi

    local output
    output="$("${executable}" "$@" 2>&1)"

    local version
    version="$(printf '%s' "${output}" | extract_version)"

    if [[ -z "${version}" ]]; then
        version="present, version not parsed"
    fi

    report "${label}" "${version}" "$(command -v "${executable}")"
}


# ==============================================================================
# Start the file
# ==============================================================================

{
    printf '# Software versions for the Rwanda P. falciparum amplicon analysis\n'
    printf '# Generated by environment/capture_versions.sh\n'
    printf '# Machine: %s\n' "$(uname -srm)"
    printf '# Date:    %s\n' "$(date -u '+%Y-%m-%d %H:%M:%S UTC')"
    printf '#\n'
    printf '# Anything marked "not found" or "fill in by hand" needs a value\n'
    printf '# entered manually. Add a note saying where you got it.\n'
    printf '#\n'
    printf 'TOOL\tVERSION\tSOURCE\n'
} > "${OUTPUT_FILE}"

echo "Capturing software versions"
echo

echo "Alignment and BAM handling:"
capture "minimap2"  minimap2  --version
capture "samtools"  samtools  --version
capture "bedtools"  bedtools  --version

echo
echo "Variant calling and manipulation:"
capture "bcftools"  bcftools  --version
capture "vcftools"  vcftools  --version
capture "tabix"     tabix     --version
capture "bgzip"     bgzip     --version

echo
echo "Quality control:"
capture "NanoPlot"  NanoPlot  --version

echo
echo "Runtimes:"
capture "python3"   python3   --version
capture "R"         R         --version
capture "java"      java      -version


# ------------------------------------------------------------------------------
# GATK. Launched from a jar or a wrapper script rather than being on PATH, so
# the path from config.sh is tried first.
# ------------------------------------------------------------------------------

echo
echo "GATK:"

GATK_VERSION=""
GATK_SOURCE=""

if [[ -n "${GATK_DIR:-}" ]] && [[ -x "${GATK_DIR}/gatk" ]]; then

    GATK_VERSION="$("${GATK_DIR}/gatk" --version 2>&1 | extract_version)"
    GATK_SOURCE="${GATK_DIR}/gatk"

elif command -v gatk >/dev/null 2>&1; then

    GATK_VERSION="$(gatk --version 2>&1 | extract_version)"
    GATK_SOURCE="$(command -v gatk)"

fi

if [[ -n "${GATK_VERSION}" ]]; then
    report "GATK" "${GATK_VERSION}" "${GATK_SOURCE}"
else
    report "GATK" "not found" "set GATK_DIR in config.sh"
fi


# ------------------------------------------------------------------------------
# SnpEff. Reports its version when run with no arguments, on stderr.
# ------------------------------------------------------------------------------

echo
echo "SnpEff:"

SNPEFF_VERSION=""
SNPEFF_SOURCE=""

if [[ -n "${SNPEFF_DIR:-}" ]] && [[ -f "${SNPEFF_DIR}/snpEff.jar" ]]; then

    SNPEFF_VERSION="$(
        java -jar "${SNPEFF_DIR}/snpEff.jar" -version 2>&1 | extract_version
    )"
    SNPEFF_SOURCE="${SNPEFF_DIR}/snpEff.jar"

fi

if [[ -n "${SNPEFF_VERSION}" ]]; then
    report "SnpEff" "${SNPEFF_VERSION}" "${SNPEFF_SOURCE}"
else
    report "SnpEff" "not found" "set SNPEFF_DIR in config.sh"
fi


# ------------------------------------------------------------------------------
# Clair3, which ran from a Docker image. The image tag is what pins it, so
# that is what gets recorded — not the version string inside the container,
# which lagged behind the tag. See tools/README.md.
# ------------------------------------------------------------------------------

echo
echo "Clair3:"

report "Clair3 (Docker image)" \
    "${CLAIR3_IMAGE:-not set in config.sh}" \
    "image tag; see tools/README.md"

report "Clair3 basecalling model" \
    "${CLAIR3_MODEL:-not set in config.sh}" \
    "from ONT Rerio"

if command -v docker >/dev/null 2>&1; then
    capture "Docker" docker --version
else
    report "Docker" "not found" "not installed on this machine"
fi


# ------------------------------------------------------------------------------
# Things that cannot be read off this machine: the sequencing instrument
# software, and the two tools used only through their interfaces. These are
# written as placeholders so they are visibly unfilled rather than absent.
# ------------------------------------------------------------------------------

echo
echo "Not detectable from this machine — fill these in by hand:"

report "MinKNOW"          "25.09.16"         "sequencing run metadata"
report "Dorado"           "fill in by hand"  "see MinKNOW run report"
report "PlasmoDB release" "fill in by hand"  "release the reference came from"
report "Qualimap"         "fill in by hand"  "used through its interface"
report "IGV"              "fill in by hand"  "used through its interface"


# ==============================================================================
# R packages
# ==============================================================================
#
# Written to a separate file, because there are many of them and they belong
# with the R analyses rather than in the command-line tool table.
# ==============================================================================

R_PACKAGE_FILE="${SCRIPT_DIR}/r_package_versions.tsv"

if command -v Rscript >/dev/null 2>&1; then

    echo
    echo "Capturing R package versions"

    Rscript -e '
      packages <- c(
        "vcfR", "dplyr", "tidyr", "readr", "stringr", "tibble", "purrr",
        "ggplot2", "ggrepel", "patchwork", "vegan", "ape", "poppr",
        "adegenet", "sf", "scales", "RColorBrewer"
      )

      rows <- lapply(packages, function(package_name) {
        version <- tryCatch(
          as.character(utils::packageVersion(package_name)),
          error = function(e) "not installed"
        )
        data.frame(PACKAGE = package_name, VERSION = version)
      })

      result <- do.call(rbind, rows)

      output_file <- commandArgs(trailingOnly = TRUE)[[1]]

      writeLines(
        c(
          paste("# R:", R.version.string),
          paste("# Platform:", R.version$platform),
          paste("# Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
          "#",
          "PACKAGE\tVERSION"
        ),
        output_file
      )

      write.table(
        result, output_file,
        sep = "\t", quote = FALSE,
        row.names = FALSE, col.names = FALSE,
        append = TRUE
      )
    ' "${R_PACKAGE_FILE}"

    echo "  written: ${R_PACKAGE_FILE}"

else
    echo
    echo "  Rscript not found; skipping R package versions."
fi


# ==============================================================================
# Done
# ==============================================================================

echo
echo "============================================================"
echo "WRITTEN"
echo "============================================================"
echo "  ${OUTPUT_FILE}"

if [[ -f "${R_PACKAGE_FILE}" ]]; then
    echo "  ${R_PACKAGE_FILE}"
fi

echo
MISSING_COUNT="$(
    grep -cE $'\t(not found|fill in by hand)\t' "${OUTPUT_FILE}" || true
)"

if [[ "${MISSING_COUNT}" -gt 0 ]]; then
    echo "${MISSING_COUNT} entr(ies) still need a value:"
    echo
    grep -E $'\t(not found|fill in by hand)\t' "${OUTPUT_FILE}" \
        | awk -F'\t' '{printf "  %-26s %s\n", $1, $3}'
    echo
    echo "Fill those in, then check this table against the Methods section"
    echo "before submission."
else
    echo "Every tool was detected. Check this table against the Methods"
    echo "section before submission."
fi
