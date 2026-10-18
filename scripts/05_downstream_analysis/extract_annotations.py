#!/usr/bin/env python3
"""
extract_annotations.py — pull SnpEff annotations out of annotated VCFs

WHAT IT DOES
    Reads every annotated per-sample VCF from the drug-resistance track and
    unpacks the SnpEff ANN field into columns: the affected allele, the
    variant type, its predicted impact, the gene, the transcript, and the
    coding and protein changes in HGVS notation.

    A single variant can carry several ANN entries, one per affected
    transcript or allele, and each becomes its own row. The protein change
    (HGVS_PROTEIN) is what identifies a drug-resistance marker: K13 R561H
    appears here as p.Arg561His.

INPUT
    --input-dir   folder of <SampleID>.ann.vcf files from stage 04

OUTPUT
    --out         TSV: SAMPLE CHROM POS REF ALT ALLELE TYPE IMPACT
                       GENE_NAME GENE_ID FEATURE_ID HGVS_CODING HGVS_PROTEIN

USAGE
    python3 scripts/05_downstream_analysis/extract_annotations.py \
        --input-dir data/processed/annotated_vcfs \
        --out       data/processed/annotation_information.tsv

Part of: Diallo et al., Nature Communications (2026)
         https://github.com/mtd1997/Adaptive_Pf_resurgence_Rwanda
"""

import argparse
import csv
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vcf_common import open_vcf, find_vcfs, parse_ann, ANN_FIELDS   # noqa: E402

HEADER = ["SAMPLE", "CHROM", "POS", "REF", "ALT"] + list(ANN_FIELDS.keys())


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--input-dir", default="data/processed/annotated_vcfs")
    parser.add_argument("--out", default="data/processed/annotation_information.tsv")
    args = parser.parse_args()

    vcfs = find_vcfs(args.input_dir)
    print("Samples : %d" % len(vcfs))
    print("Output  : %s" % args.out)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    rows = 0
    records_without_ann = 0
    genes = Counter()

    with open(args.out, "w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, fieldnames=HEADER, delimiter="\t")
        writer.writeheader()

        for index, (sample, path) in enumerate(sorted(vcfs.items()), start=1):
            with open_vcf(path) as handle:
                for line in handle:
                    if line.startswith("#"):
                        continue

                    cols = line.rstrip("\n").split("\t")
                    if len(cols) < 8:
                        continue

                    chrom, pos, ref, alt, info = (
                        cols[0], cols[1], cols[3], cols[4], cols[7])

                    wrote_any = False
                    for annotation in parse_ann(info):
                        row = {"SAMPLE": sample, "CHROM": chrom, "POS": pos,
                               "REF": ref, "ALT": alt}
                        row.update(annotation)
                        writer.writerow(row)
                        genes[annotation["GENE_NAME"]] += 1
                        rows += 1
                        wrote_any = True

                    if not wrote_any:
                        records_without_ann += 1

            if index % 200 == 0:
                print("  ...%d/%d" % (index, len(vcfs)))

    print()
    print("=" * 60)
    print("Rows written            : %d" % rows)
    print("Records with no ANN     : %d" % records_without_ann)
    if rows and records_without_ann > rows:
        print("    WARNING: most records carry no annotation. Check that the")
        print("    SnpEff database matches the reference used for alignment.")
    print()
    print("Genes represented:")
    for gene, count in genes.most_common(15):
        print("  %-12s %d" % (gene or "(none)", count))
    print("=" * 60)


if __name__ == "__main__":
    main()
