#!/bin/bash

set -euo pipefail

FASTQ_DIR="trimmed_fastq"
OUTDIR="nanoplot_qc"

mkdir -p "${OUTDIR}"

NanoPlot \
    --fastq "${FASTQ_DIR}"/*.fastq.gz \
    --outdir "${OUTDIR}" \
    --threads 8 \
    --plots hex dot \
    --N50


