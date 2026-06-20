#!/bin/bash

set -euo pipefail

# Variant calling using Clair3 in haploid mode for
# Oxford Nanopore sequencing data.

PROJECT_DIR="."

INPUT_DIR="${PROJECT_DIR}/bam_files"
OUTPUT_DIR="${PROJECT_DIR}/clair3_output"

REF="${PROJECT_DIR}/PlasmoDB-68_Pfalciparum3D7_Genome.fasta"

MODEL_DIR="/path/to/clair3_models"
MODEL_NAME="r1041_e82_400bps_hac_v410"

THREADS=8

mkdir -p "${OUTPUT_DIR}"

for bam in "${INPUT_DIR}"/*.bam
do

    sample=$(basename "${bam}" .bam)

    outdir="${OUTPUT_DIR}/${sample}_clair3"

    echo "Processing ${sample}"

    if [[ -d "${outdir}" ]]; then
        echo "Skipping ${sample}"
        continue
    fi

    docker run --rm \
        -v "${INPUT_DIR}:${INPUT_DIR}" \
        -v "${OUTPUT_DIR}:${OUTPUT_DIR}" \
        -v "${MODEL_DIR}:${MODEL_DIR}" \
        -v "$(dirname "${REF}")":"$(dirname "${REF}")" \
        hkubal/clair3:latest \
        /opt/bin/run_clair3.sh \
            --bam_fn="${bam}" \
            --ref_fn="${REF}" \
            --threads="${THREADS}" \
            --platform="ont" \
            --model_path="${MODEL_DIR}/${MODEL_NAME}" \
            --haploid_sensitive \
            --include_all_ctgs \
            --var_pct_full=1 \
            --ref_pct_full=1 \
            --no_phasing_for_fa \
            --output="${outdir}" \
            --sample_name="${sample}" \
            --gvcf \
            --remove_intermediate_dir

done

echo "Variant calling completed."
