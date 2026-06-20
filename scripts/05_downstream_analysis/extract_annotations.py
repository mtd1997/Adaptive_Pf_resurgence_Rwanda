import os
import csv

# Extract SnpEff annotations from annotated VCF files.

INPUT_DIR = "annotated_vcfs"

OUTPUT_FILE = "annotation_information.tsv"

FIELDS = [

    "SAMPLE",
    "CHROM",
    "POS",
    "REF",
    "ALT",

    "ALLELE",
    "TYPE",
    "IMPACT",
    "GENE_NAME",
    "GENE_ID",
    "FEATURE_ID",
    "HGVS_CODING",
    "HGVS_PROTEIN"

]

with open(OUTPUT_FILE, "w", newline="") as out:

    writer = csv.DictWriter(out, fieldnames=FIELDS, delimiter="\t")

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

                fields = line.strip().split("\t")

                chrom = fields[0]
                pos = fields[1]
                ref = fields[3]
                alt = fields[4]
                info = fields[7]

                ann = [x for x in info.split(";") if x.startswith("ANN=")]

                if not ann:
                    continue

                ann = ann[0].replace("ANN=", "").split(",")

                for item in ann:

                    parts = item.split("|")

                    if len(parts) < 11:
                        continue

                    writer.writerow({

                        "SAMPLE": sample,
                        "CHROM": chrom,
                        "POS": pos,
                        "REF": ref,
                        "ALT": alt,

                        "ALLELE": parts[0],
                        "TYPE": parts[1],
                        "IMPACT": parts[2],
                        "GENE_NAME": parts[3],
                        "GENE_ID": parts[4],
                        "FEATURE_ID": parts[6],
                        "HGVS_CODING": parts[9],
                        "HGVS_PROTEIN": parts[10]

                    })

print("Annotation extraction completed.")
