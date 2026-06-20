import csv

# Merge genotype and annotation information into a final mutation table.

GENOTYPES = "genotype_information.tsv"

ANNOTATIONS = "annotation_information.tsv"

OUTPUT = "drug_marker_table.tsv"

genotypes = {}

with open(GENOTYPES) as f:

    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        key = (
            row["SAMPLE"],
            row["CHROM"],
            row["POS"],
            row["REF"],
            row["ALT"]
        )

        genotypes[key] = row

DESIRED_ORDER = [

    "SAMPLE",
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "QUAL",
    "FILTER",
    "INFO",

    "GT",
    "GQ",
    "DP",
    "AF",
    "AD",
    "PL",

    "ALLELE",
    "GENE_NAME",
    "GENE_ID",
    "HGVS_CODING",
    "HGVS_PROTEIN",
    "TYPE",
    "IMPACT",
    "FEATURE_ID"

]

with open(ANNOTATIONS) as ann, open(OUTPUT, "w", newline="") as out:

    ann_reader = csv.DictReader(ann, delimiter="\t")

    writer = csv.DictWriter(

        out,

        fieldnames=DESIRED_ORDER,

        delimiter="\t"

    )

    writer.writeheader()

    for row in ann_reader:

        key = (

            row["SAMPLE"],
            row["CHROM"],
            row["POS"],
            row["REF"],
            row["ALT"]

        )

        if key not in genotypes:
            continue

        merged = genotypes[key].copy()

        merged.update(row)

        if merged["HGVS_PROTEIN"] == "":
            continue

        writer.writerow(

            {k: merged.get(k, "") for k in DESIRED_ORDER}

        )

print("Drug marker table generated.")
