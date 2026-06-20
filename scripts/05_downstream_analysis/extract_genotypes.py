import os
import gzip
import csv

# Extract genotype and FORMAT information from annotated VCF files.

INPUT_DIR = "annotated_vcfs"
OUTPUT_FILE = "genotype_information.tsv"

FORMAT_KEYS = ["GT", "GQ", "DP", "AD", "AF", "PL"]

HEADER = [
    "SAMPLE",
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "QUAL",
    "FILTER",
    "INFO"
] + FORMAT_KEYS

with open(OUTPUT_FILE, "w", newline="") as out:

    writer = csv.DictWriter(
        out,
        fieldnames=HEADER,
        delimiter="\t"
    )

    writer.writeheader()

    for filename in os.listdir(INPUT_DIR):

        if not filename.endswith(".ann.vcf"):
            continue

        sample = filename.replace(".ann.vcf", "")

        filepath = os.path.join(INPUT_DIR, filename)

        with open(filepath) as f:

            for line in f:

                if line.startswith("#"):
                    continue

                cols = line.strip().split("\t")

                if len(cols) < 10:
                    continue

                chrom, pos, _, ref, alt, qual, fltr, info, fmt, sample_data = cols[:10]

                fmt_keys = fmt.split(":")
                fmt_values = sample_data.split(":")

                fmt_dict = dict(zip(fmt_keys, fmt_values))

                row = {

                    "SAMPLE": sample,
                    "CHROM": chrom,
                    "POS": pos,
                    "REF": ref,
                    "ALT": alt,
                    "QUAL": qual,
                    "FILTER": fltr,
                    "INFO": info

                }

                for key in FORMAT_KEYS:

                    row[key] = fmt_dict.get(key, ".")

                writer.writerow(row)

print("Genotype extraction completed.")
