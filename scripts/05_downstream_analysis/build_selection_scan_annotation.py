#!/usr/bin/env python3
"""
================================================================================
BUILD THE ANNOTATION TABLE USED BY THE SELECTION SCAN
================================================================================

WHAT IT DOES
    Reads a SnpEff-annotated VCF and writes one row per variant describing
    which gene it falls in, what it does to the protein, and whether it is a
    known drug-resistance marker.

    This single script replaces a chain of three scripts (numbered 01, 02 and
    03) that produced `03_annotated_variants_snpeff.tsv` during the original
    analysis. Those scripts were not kept. The table they produced was still
    available, and its column set is fully documented by the script that read
    it, so this is a direct reimplementation rather than a guess -- but it is
    a reimplementation, and that is stated here so no reader assumes otherwise.

WHY GENES COME FROM THE BED FILE, NOT FROM SNPEFF
    SnpEff reports whichever gene name the Plasmodium falciparum 3D7
    annotation uses, and those names are not the short symbols used in the
    manuscript (MSP2, CSP, DHFR, MDR1, CRT, DHPS, K13). Assigning genes by
    position against `data/metadata/amplicons.bed` gives exactly the seven
    manuscript symbols and cannot drift. SnpEff's own gene ID and name are
    still written out, in GENE_ID and SNPEFF_GENE_NAME, so nothing is lost.

WHY THE RESISTANCE LIST IS A DATA FILE
    Which mutations count as resistance markers decides which points are
    drawn as diamonds in the selection-scan figure. That decision belongs
    somewhere a reader can inspect it, so it lives in
    `data/metadata/resistance_markers.tsv` rather than inside this code.

HOW TO CHECK IT AGAINST THE ORIGINAL
    Pass --compare-to with the annotation table you already have:

        python3 build_selection_scan_annotation.py \
            --vcf        data/processed/selection_scan_annotated.vcf.gz \
            --bed        data/metadata/amplicons.bed \
            --markers    data/metadata/resistance_markers.tsv \
            --out        results/tables/annotated_variants_snpeff_final.tsv \
            --compare-to /path/to/annotated_variants_snpeff_final.tsv

    The script then reports every variant whose gene, effect, protein change,
    resistance flag or figure label differs. If anything differs, fix it
    before regenerating the figure -- otherwise the diamonds and the labels
    will not match the published panel.

OUTPUT COLUMNS
    CHROM POS REF ALT              position and alleles, as in the VCF
    variant_id                     CHROM_POS_REF_ALT
    MUTATION_CLASS                 SNP / INSERTION / DELETION / MNP
    GENE_NAME                      manuscript symbol, from the BED file
    GENE_ID                        SnpEff gene identifier
    SNPEFF_GENE_NAME               SnpEff's own gene name
    EFFECT                         simplified: nonsynonymous / synonymous /
                                   stop_gained / frameshift / other
    DNA_CHANGE                     SnpEff HGVS.c (transcript orientation)
    REF_AA ALT_AA AA_POS           parsed from HGVS.p
    PROTEIN_CHANGE                 one-letter form, e.g. N51I
    IS_RESISTANCE                  True / False, from the marker file
    LABEL                          figure label (see labelling rule below)
    GENOMIC_LABEL                  GENE:POS_REF>ALT
    BIOLOGICAL_LABEL               GENE:PROTEIN_CHANGE where one exists
    SNPEFF_EFFECT SNPEFF_IMPACT    SnpEff's own terms, unmodified
    SNPEFF_HGVS_C SNPEFF_HGVS_P
    SNPEFF_FEATURE_ID
    SNPEFF_WARNINGS

THE LABELLING RULE
    Reproduced exactly from the original analysis:

      nonsynonymous with a protein change  ->  GENE:PROTEIN_CHANGE
      everything else                      ->  GENE:POS_REF>ALT

    This is why the K13 stop-gained variant is labelled K13:1726187_G>A in
    the figure rather than K13:Q271*. The protein form is still written to
    BIOLOGICAL_LABEL if you prefer it.

A NOTE ON CODONS
    REF_CODON and ALT_CODON are written as NA. SnpEff's ANN field does not
    carry codon sequences, so they cannot be recovered from the VCF alone.
    Nothing downstream uses them; they are kept only so the column set still
    matches what the original selection-scan script expected.

REQUIREMENTS
    Python 3.8 or newer. Standard library only -- nothing to install.
================================================================================
"""

import argparse
import csv
import gzip
import os
import re
import sys
from collections import OrderedDict


# ------------------------------------------------------------------------------
# Three-letter to one-letter amino-acid codes, plus the stop codon.
# ------------------------------------------------------------------------------

AMINO_ACIDS = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*", "Sec": "U", "Pyl": "O", "Xaa": "X",
}

# p.Asn51Ile  ->  ("Asn", "51", "Ile")
HGVS_PROTEIN = re.compile(
    r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*|Ter)"
)

# ANN subfields, in the order SnpEff writes them.
ANN_FIELDS = [
    "allele", "effect", "impact", "gene_name", "gene_id",
    "feature_type", "feature_id", "transcript_biotype",
    "rank", "hgvs_c", "hgvs_p", "cdna_pos", "cds_pos",
    "aa_pos", "distance", "warnings",
]

# SnpEff effect terms grouped into the simplified categories the original
# analysis used. Order matters: the first group that matches wins.
EFFECT_GROUPS = [
    ("stop_gained",    ("stop_gained", "stop_lost", "start_lost")),
    ("frameshift",     ("frameshift_variant",)),
    ("nonsynonymous",  ("missense_variant", "initiator_codon_variant",
                        "conservative_missense_variant",
                        "disruptive_inframe_deletion",
                        "disruptive_inframe_insertion",
                        "conservative_inframe_deletion",
                        "conservative_inframe_insertion")),
    ("synonymous",     ("synonymous_variant", "stop_retained_variant",
                        "start_retained_variant")),
]


# ------------------------------------------------------------------------------
# Small helpers
# ------------------------------------------------------------------------------

def open_maybe_gzip(path):
    """Open a plain or gzip-compressed text file."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def fail(message):
    sys.stderr.write("ERROR: %s\n" % message)
    sys.exit(1)


def classify_mutation(ref, alt):
    """SNP, MNP, INSERTION or DELETION, from the allele lengths."""
    if len(ref) == 1 and len(alt) == 1:
        return "SNP"
    if len(ref) == len(alt):
        return "MNP"
    if len(ref) < len(alt):
        return "INSERTION"
    return "DELETION"


def simplify_effect(snpeff_effect):
    """
    Collapse SnpEff's effect vocabulary into the five categories the original
    analysis used. The labelling rule tests for 'nonsynonymous' exactly, so
    this mapping has to stay stable.
    """
    terms = snpeff_effect.split("&")
    for simple_name, members in EFFECT_GROUPS:
        for term in terms:
            if term in members:
                return simple_name
    return "other"


def parse_protein_change(hgvs_p):
    """
    Turn 'p.Asn51Ile' into ('N', 'I', '51', 'N51I').
    Returns four empty strings when there is no protein-level change.
    """
    if not hgvs_p:
        return "", "", "", ""

    match = HGVS_PROTEIN.match(hgvs_p)
    if match is None:
        return "", "", "", ""

    ref_three, position, alt_three = match.groups()

    if alt_three == "*":
        alt_three = "Ter"

    ref_one = AMINO_ACIDS.get(ref_three, "")
    alt_one = AMINO_ACIDS.get(alt_three, "")

    if not ref_one or not alt_one:
        return "", "", "", ""

    return ref_one, alt_one, position, "%s%s%s" % (ref_one, position, alt_one)


# ------------------------------------------------------------------------------
# Reading the amplicon BED file
# ------------------------------------------------------------------------------

def read_amplicons(path):
    """
    Read the BED file into a list of (chromosome, start, end, gene).

    BED coordinates are half-open and zero-based: a feature covering the
    first base of a chromosome is written as 0 1. VCF positions are
    one-based. The comparison below accounts for that.
    """
    regions = []

    with open(path, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")

            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue

            fields = line.split("\t")

            if len(fields) < 4:
                fail(
                    "%s line %d has %d column(s); the BED file needs at least "
                    "four, tab-separated, with the gene name in column 4."
                    % (path, line_number, len(fields))
                )

            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                fail(
                    "%s line %d: columns 2 and 3 must be whole numbers."
                    % (path, line_number)
                )

            regions.append((fields[0], start, end, fields[3].strip()))

    if not regions:
        fail("%s contained no usable regions." % path)

    return regions


def gene_at_position(regions, chrom, pos):
    """Return the gene whose amplicon covers this one-based position."""
    for region_chrom, start, end, gene in regions:
        if region_chrom == chrom and start < pos <= end:
            return gene
    return ""


# ------------------------------------------------------------------------------
# Reading the resistance-marker whitelist
# ------------------------------------------------------------------------------

def read_markers(path):
    """
    Read the marker file into a set of (GENE, PROTEIN_CHANGE) pairs.

    Comment lines are skipped, and the header is found rather than assumed
    to be the first line, so the explanatory block at the top of the file
    does not have to be removed.
    """
    markers = set()
    header = None

    with open(path, "rt") as handle:
        for line in handle:
            line = line.rstrip("\n")

            if not line.strip() or line.startswith("#"):
                continue

            fields = [field.strip() for field in line.split("\t")]

            if header is None:
                header = fields
                for required in ("GENE", "PROTEIN_CHANGE"):
                    if required not in header:
                        fail(
                            "%s is missing the %s column."
                            % (path, required)
                        )
                continue

            row = dict(zip(header, fields))
            gene = row.get("GENE", "")
            change = row.get("PROTEIN_CHANGE", "")

            if gene and change:
                markers.add((gene, change))

    if not markers:
        fail("%s listed no markers." % path)

    return markers


# ------------------------------------------------------------------------------
# Picking one SnpEff annotation per variant
# ------------------------------------------------------------------------------

def parse_ann(info_field):
    """Pull the ANN= entries out of a VCF INFO field."""
    for entry in info_field.split(";"):
        if entry.startswith("ANN="):
            annotations = []
            for record in entry[4:].split(","):
                parts = record.split("|")
                parts += [""] * (len(ANN_FIELDS) - len(parts))
                annotations.append(dict(zip(ANN_FIELDS, parts)))
            return annotations
    return []


def choose_annotation(annotations, gene):
    """
    SnpEff reports one annotation per transcript, most consequential first.

    Prefer an annotation that has a protein-level change, because that is
    what the figure labels use. Within that preference, keep SnpEff's own
    ordering. Falling back to the first annotation reproduces SnpEff's
    default choice.
    """
    if not annotations:
        return {}

    for annotation in annotations:
        if annotation.get("hgvs_p"):
            return annotation

    return annotations[0]


# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

OUTPUT_COLUMNS = [
    "CHROM", "POS", "REF", "ALT",
    "variant_id",
    "MUTATION_CLASS",
    "GENE_NAME", "GENE_ID", "SNPEFF_GENE_NAME",
    "EFFECT",
    "DNA_CHANGE",
    "REF_CODON", "ALT_CODON",
    "REF_AA", "ALT_AA", "AA_POS",
    "PROTEIN_CHANGE",
    "IS_RESISTANCE",
    "LABEL", "GENOMIC_LABEL", "BIOLOGICAL_LABEL",
    "SNPEFF_EFFECT", "SNPEFF_IMPACT",
    "SNPEFF_HGVS_C", "SNPEFF_HGVS_P",
    "SNPEFF_FEATURE_ID", "SNPEFF_WARNINGS",
]

COMPARE_COLUMNS = [
    "GENE_NAME", "EFFECT", "PROTEIN_CHANGE", "IS_RESISTANCE", "LABEL",
]


def build_rows(vcf_path, regions, markers):
    """Walk the VCF once and build one output row per variant record."""
    rows = []
    seen_ids = {}
    variants_outside_amplicons = 0
    variants_without_ann = 0

    with open_maybe_gzip(vcf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 8:
                continue

            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3].upper()
            alt = fields[4].upper()
            info = fields[7]

            # Multi-allelic records would make a one-row-per-variant table
            # ambiguous. Stopping is safer than silently keeping the first
            # alternate allele.
            if "," in alt:
                fail(
                    "%s:%d has more than one alternate allele (%s). Split "
                    "multi-allelic records first:\n"
                    "    bcftools norm -m -any -Oz -o split.vcf.gz %s"
                    % (chrom, pos, alt, vcf_path)
                )

            variant_id = "%s_%d_%s_%s" % (chrom, pos, ref, alt)

            if variant_id in seen_ids:
                fail(
                    "%s appears more than once in the VCF. Every variant must "
                    "be unique for the downstream joins to be valid."
                    % variant_id
                )
            seen_ids[variant_id] = True

            gene = gene_at_position(regions, chrom, pos)

            if not gene:
                variants_outside_amplicons += 1

            annotations = parse_ann(info)

            if not annotations:
                variants_without_ann += 1

            annotation = choose_annotation(annotations, gene)

            snpeff_effect = annotation.get("effect", "")
            hgvs_p = annotation.get("hgvs_p", "")

            effect = simplify_effect(snpeff_effect) if snpeff_effect else ""
            ref_aa, alt_aa, aa_pos, protein_change = parse_protein_change(hgvs_p)

            genomic_label = "%s:%d_%s>%s" % (gene, pos, ref, alt)

            if protein_change:
                biological_label = "%s:%s" % (gene, protein_change)
            else:
                biological_label = genomic_label

            # The original labelling rule, reproduced exactly.
            if effect == "nonsynonymous" and protein_change:
                label = "%s:%s" % (gene, protein_change)
            else:
                label = genomic_label

            is_resistance = (gene, protein_change) in markers if protein_change \
                else False

            rows.append(OrderedDict([
                ("CHROM", chrom),
                ("POS", str(pos)),
                ("REF", ref),
                ("ALT", alt),
                ("variant_id", variant_id),
                ("MUTATION_CLASS", classify_mutation(ref, alt)),
                ("GENE_NAME", gene),
                ("GENE_ID", annotation.get("gene_id", "")),
                ("SNPEFF_GENE_NAME", annotation.get("gene_name", "")),
                ("EFFECT", effect),
                ("DNA_CHANGE", annotation.get("hgvs_c", "")),
                ("REF_CODON", "NA"),
                ("ALT_CODON", "NA"),
                ("REF_AA", ref_aa),
                ("ALT_AA", alt_aa),
                ("AA_POS", aa_pos),
                ("PROTEIN_CHANGE", protein_change),
                ("IS_RESISTANCE", "True" if is_resistance else "False"),
                ("LABEL", label),
                ("GENOMIC_LABEL", genomic_label),
                ("BIOLOGICAL_LABEL", biological_label),
                ("SNPEFF_EFFECT", snpeff_effect),
                ("SNPEFF_IMPACT", annotation.get("impact", "")),
                ("SNPEFF_HGVS_C", annotation.get("hgvs_c", "")),
                ("SNPEFF_HGVS_P", hgvs_p),
                ("SNPEFF_FEATURE_ID", annotation.get("feature_id", "")),
                ("SNPEFF_WARNINGS", annotation.get("warnings", "")),
            ]))

    return rows, variants_outside_amplicons, variants_without_ann


def compare_to_original(rows, original_path, report_path):
    """
    Check the rebuilt table against the one produced by the original
    analysis, and write out every disagreement.
    """
    original = {}

    with open(original_path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if "variant_id" not in (reader.fieldnames or []):
            fail("%s has no variant_id column to compare on." % original_path)

        for row in reader:
            original[row["variant_id"]] = row

    rebuilt = {row["variant_id"]: row for row in rows}

    differences = []

    for variant_id in sorted(set(original) | set(rebuilt)):
        old = original.get(variant_id)
        new = rebuilt.get(variant_id)

        if old is None:
            differences.append((variant_id, "present only in rebuilt", "", ""))
            continue

        if new is None:
            differences.append((variant_id, "present only in original", "", ""))
            continue

        for column in COMPARE_COLUMNS:
            if column not in old:
                continue

            old_value = (old.get(column) or "").strip()
            new_value = (new.get(column) or "").strip()

            # Normalise the boolean spellings R and Python use.
            if column == "IS_RESISTANCE":
                old_value = old_value.upper().replace("TRUE", "T") \
                                            .replace("FALSE", "F")
                new_value = new_value.upper().replace("TRUE", "T") \
                                            .replace("FALSE", "F")

            if old_value != new_value:
                differences.append((variant_id, column, old_value, new_value))

    with open(report_path, "wt") as handle:
        handle.write("variant_id\tfield\toriginal\trebuilt\n")
        for row in differences:
            handle.write("%s\t%s\t%s\t%s\n" % row)

    return len(original), len(rebuilt), differences


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Build the annotation table used by the selection scan from a "
            "SnpEff-annotated VCF."
        )
    )
    parser.add_argument(
        "--vcf", required=True,
        help="SnpEff-annotated VCF (.vcf or .vcf.gz).",
    )
    parser.add_argument(
        "--bed", required=True,
        help="Amplicon BED file; column 4 supplies the gene names.",
    )
    parser.add_argument(
        "--markers", required=True,
        help="Resistance-marker whitelist (TSV).",
    )
    parser.add_argument(
        "--out", required=True,
        help="Output annotation table (TSV).",
    )
    parser.add_argument(
        "--compare-to", default=None,
        help=(
            "Optional: an annotation table from the original analysis. Every "
            "disagreement is written next to --out as "
            "annotation_comparison_report.tsv."
        ),
    )

    args = parser.parse_args()

    for label, path in (
        ("VCF", args.vcf),
        ("BED file", args.bed),
        ("marker file", args.markers),
    ):
        if not os.path.exists(path):
            fail("%s not found: %s" % (label, path))

    output_directory = os.path.dirname(os.path.abspath(args.out))
    if output_directory:
        os.makedirs(output_directory, exist_ok=True)

    regions = read_amplicons(args.bed)
    markers = read_markers(args.markers)

    sys.stderr.write(
        "Amplicons: %d    Resistance markers: %d\n"
        % (len(regions), len(markers))
    )

    rows, outside, no_ann = build_rows(args.vcf, regions, markers)

    if not rows:
        fail("The VCF contained no variant records.")

    with open(args.out, "wt") as handle:
        handle.write("\t".join(OUTPUT_COLUMNS) + "\n")
        for row in rows:
            handle.write("\t".join(row[column] for column in OUTPUT_COLUMNS))
            handle.write("\n")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------

    resistance_count = sum(
        1 for row in rows if row["IS_RESISTANCE"] == "True"
    )

    per_gene = OrderedDict()
    for row in rows:
        gene = row["GENE_NAME"] or "(outside amplicons)"
        per_gene[gene] = per_gene.get(gene, 0) + 1

    sys.stderr.write("\nVariants written: %d\n" % len(rows))
    sys.stderr.write("Resistance-associated: %d\n" % resistance_count)

    if outside:
        sys.stderr.write(
            "Outside every amplicon: %d  (these get an empty GENE_NAME and "
            "are dropped by the figure script)\n" % outside
        )

    if no_ann:
        sys.stderr.write(
            "Without a SnpEff ANN field: %d  (was the VCF actually annotated?)\n"
            % no_ann
        )

    sys.stderr.write("\nVariants per gene:\n")
    for gene, count in per_gene.items():
        sys.stderr.write("  %-22s %d\n" % (gene, count))

    sys.stderr.write("\nWritten: %s\n" % os.path.abspath(args.out))

    # ------------------------------------------------------------------
    # Optional comparison against the original table
    # ------------------------------------------------------------------

    if args.compare_to:
        if not os.path.exists(args.compare_to):
            fail("--compare-to file not found: %s" % args.compare_to)

        report_path = os.path.join(
            output_directory or ".", "annotation_comparison_report.tsv"
        )

        n_original, n_rebuilt, differences = compare_to_original(
            rows, args.compare_to, report_path
        )

        sys.stderr.write("\n" + "-" * 70 + "\n")
        sys.stderr.write("COMPARISON WITH THE ORIGINAL ANNOTATION\n")
        sys.stderr.write("-" * 70 + "\n")
        sys.stderr.write("Rows in the original: %d\n" % n_original)
        sys.stderr.write("Rows rebuilt here:    %d\n" % n_rebuilt)
        sys.stderr.write("Disagreements:        %d\n" % len(differences))
        sys.stderr.write("Report: %s\n" % os.path.abspath(report_path))

        if differences:
            sys.stderr.write(
                "\nThe two tables do not agree. Resolve this before "
                "regenerating the figure: the resistance flags decide which "
                "points are drawn as diamonds, and the labels decide what is "
                "written next to them.\n"
            )
            sys.stderr.write("\nFirst few disagreements:\n")
            for variant_id, field, old, new in differences[:10]:
                sys.stderr.write(
                    "  %-32s %-16s original=%-14s rebuilt=%s\n"
                    % (variant_id, field, old or "(empty)", new or "(empty)")
                )
        else:
            sys.stderr.write(
                "\nThe rebuilt table matches the original exactly on gene, "
                "effect, protein change, resistance flag and label.\n"
            )


if __name__ == "__main__":
    main()
