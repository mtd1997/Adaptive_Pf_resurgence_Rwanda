#!/usr/bin/env python3
"""
extract_population_variant_table.py — site annotations from the population VCF

WHAT IT DOES
    Unpacks the SnpEff ANN field from the jointly genotyped, filtered and
    annotated population VCF into a flat table of variant sites.

    Unlike the drug-resistance track, this table carries no per-sample
    genotypes. It describes the sites themselves: where each variant is, which
    gene it falls in, and what it does to the protein. selection_scan.R joins
    it to per-site FST and allele-frequency values so that each point in
    Figure 5 can be labelled by gene and variant type.

INPUT
    --vcf    final_population.ann.vcf from annotate_population_vcf.sh

OUTPUT
    --out    TSV: CHROM POS REF ALT QUAL ALLELE TYPE IMPACT
                  GENE_NAME GENE_ID FEATURE_ID HGVS_CODING HGVS_PROTEIN

USAGE
    python3 scripts/05_downstream_analysis/extract_population_variant_table.py \
        --vcf data/processed/final_population.ann.vcf \
        --out data/processed/population_variant_table.tsv

Part of: Diallo et al., Nature Communications (2026)
         https://github.com/mtd1997/Adaptive_Pf_resurgence_Rwanda
"""

import argparse
import csv
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vcf_common import open_vcf, parse_ann, ANN_FIELDS   # noqa: E402

HEADER = ["CHROM", "POS", "REF", "ALT", "QUAL"] + list(ANN_FIELDS.keys())


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--vcf", default="data/processed/final_population.ann.vcf")
    parser.add_argument("--out", default="data/processed/population_variant_table.tsv")
    args = parser.parse_args()

    if not os.path.isfile(args.vcf):
        sys.exit("ERROR: VCF not found: %s\n"
                 "       Run annotate_population_vcf.sh first." % args.vcf)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    print("Input  : %s" % args.vcf)
    print("Output : %s" % args.out)

    sites = 0
    rows = 0
    sites_without_ann = 0
    genes = Counter()
    types = Counter()

    with open_vcf(args.vcf) as handle, \
            open(args.out, "w", newline="") as out_handle:

        writer = csv.DictWriter(out_handle, fieldnames=HEADER, delimiter="\t")
        writer.writeheader()

        for line in handle:
            if line.startswith("#"):
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue

            chrom, pos, ref, alt, qual, info = (
                cols[0], cols[1], cols[3], cols[4], cols[5], cols[7])
            sites += 1

            wrote_any = False
            for annotation in parse_ann(info):
                row = {"CHROM": chrom, "POS": pos, "REF": ref,
                       "ALT": alt, "QUAL": qual}
                row.update(annotation)
                writer.writerow(row)
                genes[annotation["GENE_NAME"]] += 1
                types[annotation["TYPE"]] += 1
                rows += 1
                wrote_any = True

            if not wrote_any:
                sites_without_ann += 1

    print()
    print("=" * 60)
    print("Variant sites       : %d" % sites)
    print("Annotation rows     : %d" % rows)
    print("Sites with no ANN   : %d" % sites_without_ann)
    if sites and sites_without_ann == sites:
        print("    WARNING: nothing was annotated. The SnpEff database")
        print("    probably does not match the reference.")
    print()
    print("Sites per gene:")
    for gene, count in genes.most_common(10):
        print("  %-12s %d" % (gene or "(none)", count))
    print()
    print("Variant types:")
    for variant_type, count in types.most_common(8):
        print("  %-28s %d" % (variant_type or "(none)", count))
    print("=" * 60)


if __name__ == "__main__":
    main()
