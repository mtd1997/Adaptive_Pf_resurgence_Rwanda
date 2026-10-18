#!/usr/bin/env python3
"""
extract_genotypes.py — pull per-sample genotype fields out of annotated VCFs

WHAT IT DOES
    Reads every annotated per-sample VCF from the drug-resistance track and
    writes one long table holding the genotype call and its supporting
    evidence at each variant position: GT, GQ, DP, AD, AF and PL.

    This is the genotype half of the drug-marker table. The annotation half
    comes from extract_annotations.py, and build_drug_marker_table.py joins
    the two on (SAMPLE, CHROM, POS, REF, ALT).

INPUT
    --input-dir   folder of <SampleID>.ann.vcf files from stage 04
                  (gzipped files are read too)

OUTPUT
    --out         TSV: SAMPLE CHROM POS REF ALT QUAL FILTER INFO GT GQ DP AD AF PL

NOTE ON WILD-TYPE CALLS
    Only variant sites appear here, so every GT is 1. Positions where a sample
    is confidently wild-type are not in the VCF at all; they are recovered
    separately from the gVCF reference blocks by recover_reference_calls.py.
    Computing a wild-type frequency from this table alone would use the wrong
    denominator.

USAGE
    python3 scripts/05_downstream_analysis/extract_genotypes.py \
        --input-dir data/processed/annotated_vcfs \
        --out       data/processed/genotype_information.tsv

Part of: Diallo et al., Nature Communications (2026)
         https://github.com/mtd1997/Adaptive_Pf_resurgence_Rwanda
"""

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vcf_common import open_vcf, find_vcfs, parse_format   # noqa: E402

FORMAT_KEYS = ["GT", "GQ", "DP", "AD", "AF", "PL"]
HEADER = ["SAMPLE", "CHROM", "POS", "REF", "ALT",
          "QUAL", "FILTER", "INFO"] + FORMAT_KEYS


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--input-dir", default="data/processed/annotated_vcfs")
    parser.add_argument("--out", default="data/processed/genotype_information.tsv")
    parser.add_argument("--pass-only", action="store_true",
                        help="keep only records with FILTER=PASS")
    args = parser.parse_args()

    vcfs = find_vcfs(args.input_dir)
    print("Samples : %d" % len(vcfs))
    print("Output  : %s" % args.out)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    rows = 0
    skipped_malformed = 0
    skipped_filter = 0
    samples_without_records = []

    with open(args.out, "w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, fieldnames=HEADER, delimiter="\t")
        writer.writeheader()

        for index, (sample, path) in enumerate(sorted(vcfs.items()), start=1):
            sample_rows = 0

            with open_vcf(path) as handle:
                for line in handle:
                    if line.startswith("#"):
                        continue

                    cols = line.rstrip("\n").split("\t")
                    if len(cols) < 10:
                        skipped_malformed += 1
                        continue

                    (chrom, pos, _, ref, alt, qual,
                     filt, info, fmt, sample_col) = cols[:10]

                    if args.pass_only and filt != "PASS":
                        skipped_filter += 1
                        continue

                    call = parse_format(fmt, sample_col)

                    row = {
                        "SAMPLE": sample, "CHROM": chrom, "POS": pos,
                        "REF": ref, "ALT": alt, "QUAL": qual,
                        "FILTER": filt, "INFO": info,
                    }
                    for key in FORMAT_KEYS:
                        row[key] = call.get(key, ".")

                    writer.writerow(row)
                    rows += 1
                    sample_rows += 1

            if sample_rows == 0:
                samples_without_records.append(sample)

            if index % 200 == 0:
                print("  ...%d/%d" % (index, len(vcfs)))

    print()
    print("=" * 60)
    print("Rows written          : %d" % rows)
    print("Malformed lines       : %d" % skipped_malformed)
    if args.pass_only:
        print("Dropped by FILTER     : %d" % skipped_filter)
    print("Samples with no calls : %d" % len(samples_without_records))
    if samples_without_records:
        preview = ", ".join(samples_without_records[:10])
        if len(samples_without_records) > 10:
            preview += ", ..."
        print("    %s" % preview)
        print("    (expected for samples carrying no variants at these loci)")
    print("=" * 60)


if __name__ == "__main__":
    main()
