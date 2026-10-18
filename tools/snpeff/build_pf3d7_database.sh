#!/usr/bin/env bash
################################################################################
# BUILD THE SNPEFF DATABASE FOR PLASMODIUM FALCIPARUM 3D7
################################################################################
#
# WHAT IT DOES
#   Builds the SnpEff genome database used to annotate variants in this
#   project. SnpEff does not ship a P. falciparum database, so one has to be
#   built from the reference genome and its annotation.
#
# A NOTE ON PROVENANCE
#   The original build was done interactively and the exact commands were not
#   recorded at the time. This script reconstructs them from the SnpEff
#   configuration entry and the database directory that the build left
#   behind, both of which survived. It is a documented reconstruction, not a
#   transcript, and that is stated plainly here.
#
#   Nothing about the annotation depends on how the database was built, only
#   on which reference and annotation went into it. Those are pinned below,
#   so a database built with this script annotates variants identically.
#
# WHAT YOU NEED FIRST
#   1. SnpEff installed, and SNPEFF_DIR in config.sh pointing at it
#   2. The P. falciparum 3D7 reference genome FASTA
#   3. The matching GFF3 annotation
#
#   Both reference files come from PlasmoDB. Use the same release for each —
#   mixing a genome from one release with an annotation from another produces
#   annotations that are subtly wrong rather than obviously broken.
#
#     https://plasmodb.org/plasmo/app/downloads/
#
#   The release used in this project is recorded in
#   ../environment/software_versions.tsv.
#
# HOW TO RUN
#   From the top of the repository:
#
#     bash tools/snpeff/build_pf3d7_database.sh
#
#   The build takes a few minutes and needs a couple of gigabytes of memory.
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
# The genome name. This is the string SnpEff uses to identify the database,
# and it must match what the annotation step passes to SnpEff. Changing it
# here means changing it in scripts/04_filtering/annotate_population_vcf.sh
# as well.
# ------------------------------------------------------------------------------

GENOME_NAME="Pf3D7"
GENOME_DESCRIPTION="Plasmodium_falciparum_3D7"


# ------------------------------------------------------------------------------
# Check that everything needed is present, and say exactly what is missing
# rather than failing part-way through the build.
# ------------------------------------------------------------------------------

if [[ -z "${SNPEFF_DIR:-}" ]]; then
    echo "ERROR: SNPEFF_DIR is not set in config.sh." >&2
    echo "       Set it to the directory containing snpEff.jar." >&2
    exit 1
fi

SNPEFF_JAR="${SNPEFF_DIR}/snpEff.jar"
SNPEFF_CONFIG="${SNPEFF_DIR}/snpEff.config"

for REQUIRED in "${SNPEFF_JAR}" "${SNPEFF_CONFIG}"; do
    if [[ ! -f "${REQUIRED}" ]]; then
        echo "ERROR: not found: ${REQUIRED}" >&2
        exit 1
    fi
done

if [[ -z "${REFERENCE_FASTA:-}" ]] || [[ ! -f "${REFERENCE_FASTA}" ]]; then
    echo "ERROR: REFERENCE_FASTA is not set in config.sh, or the file is" >&2
    echo "       missing. It should be the P. falciparum 3D7 genome FASTA." >&2
    exit 1
fi

if [[ -z "${REFERENCE_GFF:-}" ]] || [[ ! -f "${REFERENCE_GFF}" ]]; then
    echo "ERROR: REFERENCE_GFF is not set in config.sh, or the file is" >&2
    echo "       missing. It should be the GFF3 annotation matching the" >&2
    echo "       genome FASTA above, from the same PlasmoDB release." >&2
    exit 1
fi

if ! command -v java >/dev/null 2>&1; then
    echo "ERROR: java was not found on PATH. SnpEff needs a Java runtime." >&2
    exit 1
fi


# ------------------------------------------------------------------------------
# SnpEff expects its input files in data/<GENOME_NAME>/ under its own
# directory, with fixed names.
# ------------------------------------------------------------------------------

GENOME_DATA_DIR="${SNPEFF_DIR}/data/${GENOME_NAME}"

echo "SnpEff directory:  ${SNPEFF_DIR}"
echo "Genome name:       ${GENOME_NAME}"
echo "Database will be:  ${GENOME_DATA_DIR}"
echo


# ------------------------------------------------------------------------------
# A database that already exists is left alone. Rebuilding it would be a
# waste of time at best, and at worst would mask the fact that the published
# annotation came from a different build.
# ------------------------------------------------------------------------------

if [[ -f "${GENOME_DATA_DIR}/snpEffectPredictor.bin" ]]; then
    echo "The database already exists:"
    echo "  ${GENOME_DATA_DIR}/snpEffectPredictor.bin"
    echo
    echo "Nothing to do. To force a rebuild, remove that file first."
    exit 0
fi

mkdir -p "${GENOME_DATA_DIR}"


# ------------------------------------------------------------------------------
# Register the genome in snpEff.config, if it is not already there.
#
# The entry SnpEff needs is a single line:
#
#   Pf3D7.genome : Plasmodium_falciparum_3D7
#
# Appending it is idempotent: the grep guard means running this script twice
# does not add a duplicate entry.
# ------------------------------------------------------------------------------

CONFIG_ENTRY="${GENOME_NAME}.genome : ${GENOME_DESCRIPTION}"

if grep -qE "^[[:space:]]*${GENOME_NAME}\.genome[[:space:]]*:" "${SNPEFF_CONFIG}"; then

    echo "The genome is already registered in snpEff.config:"
    grep -E "^[[:space:]]*${GENOME_NAME}\.genome[[:space:]]*:" "${SNPEFF_CONFIG}"

else

    echo "Registering the genome in snpEff.config"

    # Keep a copy before editing someone else's installation.
    cp "${SNPEFF_CONFIG}" "${SNPEFF_CONFIG}.backup.$(date +%Y%m%d%H%M%S)"

    {
        printf '\n'
        printf '# Added by tools/snpeff/build_pf3d7_database.sh\n'
        printf '# Plasmodium falciparum 3D7, for the Rwanda amplicon analysis\n'
        printf '%s\n' "${CONFIG_ENTRY}"
    } >> "${SNPEFF_CONFIG}"
fi

echo


# ------------------------------------------------------------------------------
# Copy the reference files into place under the names SnpEff requires.
#
# These are copies rather than symlinks: SnpEff reads them during the build,
# and a symlink into an external drive that is later unplugged turns a
# working installation into a confusing failure.
# ------------------------------------------------------------------------------

echo "Copying the reference genome"
cp -f "${REFERENCE_FASTA}" "${GENOME_DATA_DIR}/sequences.fa"

echo "Copying the annotation"
cp -f "${REFERENCE_GFF}" "${GENOME_DATA_DIR}/genes.gff"

echo


# ------------------------------------------------------------------------------
# Build.
#
#   -gff3   the annotation is GFF3
#   -v      verbose, so a failure says which gene it failed on
#
# -noCheckCds and -noCheckProtein are deliberately NOT used. Those flags
# silence the consistency checks between the annotation and the genome, and
# those checks are exactly what catches a mismatched PlasmoDB release. If the
# build fails on them, the reference files do not agree with each other and
# the right fix is to get a matching pair.
# ------------------------------------------------------------------------------

echo "Building the database. This takes a few minutes."
echo

java -Xmx4g -jar "${SNPEFF_JAR}" \
    build \
    -gff3 \
    -v \
    -dataDir "${SNPEFF_DIR}/data" \
    -c "${SNPEFF_CONFIG}" \
    "${GENOME_NAME}"

echo


# ------------------------------------------------------------------------------
# Verify.
# ------------------------------------------------------------------------------

if [[ ! -f "${GENOME_DATA_DIR}/snpEffectPredictor.bin" ]]; then
    echo "ERROR: the build finished but snpEffectPredictor.bin was not" >&2
    echo "       created. Look back through the output above for the" >&2
    echo "       failure." >&2
    exit 1
fi

echo "============================================================"
echo "DATABASE BUILT"
echo "============================================================"
echo "  ${GENOME_DATA_DIR}/snpEffectPredictor.bin"
echo

echo "Checking that SnpEff can load it and that the target genes are present:"
echo

java -Xmx4g -jar "${SNPEFF_JAR}" \
    genes2bed \
    -c "${SNPEFF_CONFIG}" \
    -dataDir "${SNPEFF_DIR}/data" \
    "${GENOME_NAME}" \
    2>/dev/null \
    | head -5

echo
echo "The database is ready. Annotation step:"
echo "  bash scripts/04_filtering/annotate_population_vcf.sh"
