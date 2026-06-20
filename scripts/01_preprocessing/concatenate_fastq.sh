#!/bin/bash

set -euo pipefail

# Concatenate FASTQ files generated for each barcode
# after MinKNOW basecalling and demultiplexing.

INPUT_DIR="sequencing_run"
OUTPUT_DIR="concatenated_fastq_files"

mkdir -p "${OUTPUT_DIR}"

for BARCODE_DIR in "${INPUT_DIR}"/barcode*
do

    BARCODE=$(basename "${BARCODE_DIR}")

    echo "Processing ${BARCODE}"

    cat "${BARCODE_DIR}"/*.fastq.gz \
        > "${OUTPUT_DIR}/${BARCODE}.fastq.gz"

done

