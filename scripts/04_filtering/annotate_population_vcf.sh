```bash
#!/bin/bash

set -euo pipefail

# Annotate the filtered population VCF using SnpEff.

INPUT="final_population.vcf.gz"

OUTPUT="final_population.ann.vcf"

STATS="final_population_snpeff.html"

SNPEFF="snpEff/snpEff.jar"

java -Xmx4g \
    -jar "${SNPEFF}" \
    -v Pf3D7 \
    "${INPUT}" \
    -stats "${STATS}" \
    > "${OUTPUT}"

echo "Population VCF annotation completed."

