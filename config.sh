# All file paths for the bash pipeline (stages 01-04), in one place.
# Scripts source this file instead of hard-coding paths.
#   source "$(dirname "$0")/../../config.sh"

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- Reference genome: P. falciparum 3D7 v3, PlasmoDB release 68 ---
# Download from https://plasmodb.org/
REF_GENOME="${PROJECT_DIR}/reference/PlasmoDB-68_Pfalciparum3D7_Genome.fasta"

# --- Small inputs tracked in this repository ---
AMPLICON_BED="${PROJECT_DIR}/data/metadata/amplicons.bed"
REGIONS_FILE="${PROJECT_DIR}/data/metadata/regions.txt"
BARCODE_MAP="${PROJECT_DIR}/data/metadata/sample_barcode_mapping.csv"
SAMPLE_METADATA="${PROJECT_DIR}/data/metadata/samples_metadata.tsv"

# --- Sequencing data: NOT in this repository ---
# Obtain from the repository given in the paper's Data Availability statement,
# then set DATA_ROOT to wherever you downloaded it.
DATA_ROOT="${HOME}/pf_rwanda_data"
FASTQ_DIR="${DATA_ROOT}/fastq"
BAM_DIR="${DATA_ROOT}/bam"
VCF_DIR="${DATA_ROOT}/vcf"
GVCF_DIR="${DATA_ROOT}/gvcf"
ANNOTATED_DIR="${DATA_ROOT}/annotated"
FILTERED_DIR="${DATA_ROOT}/filtered"
DRUG_DIR="${DATA_ROOT}/drug_resistance"

# --- Outputs ---
PROCESSED_DIR="${PROJECT_DIR}/data/processed"
RESULTS_DIR="${PROJECT_DIR}/results"
FIGURE_DIR="${RESULTS_DIR}/figures"
TABLE_DIR="${RESULTS_DIR}/tables"

# --- Tools not on PATH ---
SNPEFF_JAR="${HOME}/tools/snpEff/snpEff.jar"
SNPEFF_DB="Pf3D7"
CLAIR3_IMAGE="hkubal/clair3:v2.0.0"

THREADS=8

mkdir -p "${PROCESSED_DIR}" "${FIGURE_DIR}" "${TABLE_DIR}"

# --- stage 01 ---
MERGED_DIR="${DATA_ROOT}/fastq_merged"
QC_DIR="${RESULTS_DIR}/qc"
LOG_DIR="${RESULTS_DIR}/logs"
