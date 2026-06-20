#!/bin/bash

set -euo pipefail

# Joint genotyping of individual gVCFs using GATK.

REF="PlasmoDB-68_Pfalciparum3D7_Genome.fasta"

GVCF_DIR="gvcf_files"

DATABASE="genomicsdb"

SAMPLE_MAP="sample_map.txt"

INTERVALS="intervals.list"

OUTPUT="joint_genotyped.vcf.gz"

##################################################
# Create sample map
##################################################

ls "${GVCF_DIR}"/*.gvcf.gz \
| sed 's/.gvcf.gz//' \
| awk '{print $1"\t"$1".gvcf.gz"}' \
> "${SAMPLE_MAP}"

##################################################
# Import gVCFs into GenomicsDB
##################################################

gatk GenomicsDBImport \
    --genomicsdb-workspace-path "${DATABASE}" \
    --sample-name-map "${SAMPLE_MAP}" \
    --reader-threads 5 \
    -L "${INTERVALS}"

##################################################
# Joint genotyping
##################################################

gatk GenotypeGVCFs \
    -R "${REF}" \
    -V gendb://"${DATABASE}" \
    -O "${OUTPUT}" \
    -L "${INTERVALS}"

echo "Joint genotyping completed."

