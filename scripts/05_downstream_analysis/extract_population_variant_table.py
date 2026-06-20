import csv

INPUT = "final_population.ann.vcf"
OUTPUT = "population_variant_table.tsv"

HEADER = [
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "QUAL",
    "ALLELE",
    "TYPE",
    "IMPACT",
    "GENE_NAME",
    "GENE_ID",
    "FEATURE_ID",
    "HGVS_CODING",
    "HGVS_PROTEIN"
]

with open(INPUT) as infile, open(OUTPUT, "w", newline="") as outfile:

    writer = csv.writer(outfile, delimiter="\t")
    writer.writerow(HEADER)

    for line in infile:

        if line.startswith("#"):
            continue

        fields = line.strip().split("\t")

        chrom = fields[0]
        pos = fields[1]
        ref = fields[3]
        alt = fields[4]
        qual = fields[5]
        info = fields[7]

        ann = [x for x in info.split(";") if x.startswith("ANN=")]

        if not ann:
            continue

        annotations = ann[0].replace("ANN=", "").split(",")

        for entry in annotations:

            parts = entry.split("|")

            if len(parts) < 11:
                continue

            writer.writerow([

                chrom,
                pos,
                ref,
                alt,
                qual,

                parts[0],
                parts[1],
                parts[2],
                parts[3],
                parts[4],
                parts[6],
                parts[9],
                parts[10]

            ])

print("Population annotation table generated.")
