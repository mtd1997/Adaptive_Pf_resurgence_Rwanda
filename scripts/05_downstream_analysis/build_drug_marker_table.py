#!/usr/bin/env python3
"""
build_drug_marker_table.py — join genotypes and annotations into one table

WHAT IT DOES
    Combines the two tables produced by extract_genotypes.py and
    extract_annotations.py into the drug-resistance marker table used by the
    downstream R analyses.

    The join key is (SAMPLE, CHROM, POS, REF, ALT): the same sample, at the
    same position, with the same alleles. Both tables derive their sample
    names the same way, so the key is consistent.

    Rows without a protein-level change are dropped. A drug-resistance marker
    is defined by its amino-acid substitution, so intronic, intergenic and
    other non-coding annotations carry no marker information. This filter is
    deliberate and is reported in the summary below.

INPUTS
    --genotypes     genotype_information.tsv
    --annotations   annotation_information.tsv

OUTPUT
    --out           drug_marker_table.tsv

OPTIONS
    --keep-noncoding   keep rows with an empty HGVS_PROTEIN

USAGE
    python3 scripts/05_downstream_analysis/build_drug_marker_table.py \
        --genotypes   data/processed/genotype_information.tsv \
        --annotations data/processed/annotation_information.tsv \
        --out         data/processed/drug_marker_table.tsv

Part of: Diallo et al., Nature Communications (2026)
         https://github.com/mtd1997/Adaptive_Pf_resurgence_Rwanda
"""

import argparse
import csv
import os
import sys
from collections import Counter

# Column order of the output. Declared explicitly rather than derived from the
# inputs, so the layout is stable from run to run.
DESIRED_ORDER = [
    "SAMPLE", "CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", "INFO",
    "GT", "GQ", "DP", "AF", "AD", "PL",
    "ALLELE", "GENE_NAME", "GENE_ID", "HGVS_CODING", "HGVS_PROTEIN",
    "TYPE", "IMPACT", "FEATURE_ID",
]

KEY_COLUMNS = ["SAMPLE", "CHROM", "POS", "REF", "ALT"]


def read_genotypes(path):
    """Index the genotype table by join key."""
    if not os.path.isfile(path):
        sys.exit("ERROR: genotype table not found: %s\n"
                 "       Run extract_genotypes.py first." % path)

    index = {}
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [c for c in KEY_COLUMNS if c not in (reader.fieldnames or [])]
        if missing:
            sys.exit("ERROR: %s lacks column(s): %s"
                     % (path, ", ".join(missing)))
        for row in reader:
            index[tuple(row[c] for c in KEY_COLUMNS)] = row
    return index


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--genotypes", default="data/processed/genotype_information.tsv")
    parser.add_argument("--annotations", default="data/processed/annotation_information.tsv")
    parser.add_argument("--out", default="data/processed/drug_marker_table.tsv")
    parser.add_argument("--keep-noncoding", action="store_true")
    args = parser.parse_args()

    genotypes = read_genotypes(args.genotypes)
    print("Genotype records  : %d" % len(genotypes))

    if not os.path.isfile(args.annotations):
        sys.exit("ERROR: annotation table not found: %s\n"
                 "       Run extract_annotations.py first." % args.annotations)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    written = 0
    unmatched = 0
    dropped_noncoding = 0
    unmatched_samples = Counter()
    genes = Counter()

    with open(args.annotations) as ann_handle, \
            open(args.out, "w", newline="") as out_handle:

        reader = csv.DictReader(ann_handle, delimiter="\t")
        writer = csv.DictWriter(out_handle, fieldnames=DESIRED_ORDER,
                                delimiter="\t", extrasaction="ignore")
        writer.writeheader()

        for row in reader:
            key = tuple(row[c] for c in KEY_COLUMNS)

            genotype = genotypes.get(key)
            if genotype is None:
                # An annotation with no matching genotype. Counted rather than
                # silently discarded: a large number here means the two tables
                # were built from different inputs.
                unmatched += 1
                unmatched_samples[row["SAMPLE"]] += 1
                continue

            merged = genotype.copy()
            merged.update(row)

            if not merged.get("HGVS_PROTEIN", "").strip():
                dropped_noncoding += 1
                if not args.keep_noncoding:
                    continue

            writer.writerow({c: merged.get(c, "") for c in DESIRED_ORDER})
            genes[merged.get("GENE_NAME", "")] += 1
            written += 1

    print("Rows written      : %d" % written)
    print("Unmatched         : %d" % unmatched)
    if unmatched:
        print("    affecting %d sample(s), worst:" % len(unmatched_samples))
        for sample, count in unmatched_samples.most_common(5):
            print("      %-18s %d" % (sample, count))
        if unmatched > written:
            print("    WARNING: more rows failed to match than matched.")
            print("    Check that both tables came from the same VCFs.")
    print("Non-coding rows   : %d %s"
          % (dropped_noncoding, "(kept)" if args.keep_noncoding else "(dropped)"))
    print()
    print("Markers by gene:")
    for gene, count in genes.most_common(15):
        print("  %-12s %d" % (gene or "(none)", count))
    print()
    print("Output            : %s" % args.out)


if __name__ == "__main__":
    main()
