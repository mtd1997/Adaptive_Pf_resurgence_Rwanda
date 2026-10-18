"""
vcf_common.py — shared helpers for the VCF extraction scripts

Small utilities used by extract_genotypes.py, extract_annotations.py and
extract_population_variant_table.py, so that file opening, path handling and
the SnpEff ANN layout are defined once rather than repeated in each script.

Part of: Diallo et al., Nature Communications (2026)
         https://github.com/mtd1997/Adaptive_Pf_resurgence_Rwanda
"""

import gzip
import os
import sys

# -----------------------------------------------------------------------------
# SnpEff ANN field layout.
#
# SnpEff packs its annotation into a single INFO key, ANN=, with fields
# separated by '|'. The positions below follow the SnpEff 5.x specification:
#
#   0 Allele        1 Annotation (TYPE)   2 Annotation_Impact   3 Gene_Name
#   4 Gene_ID       5 Feature_Type        6 Feature_ID          7 Transcript_BioType
#   8 Rank          9 HGVS.c             10 HGVS.p
#
# A record can carry several ANN entries separated by commas, one per
# affected transcript or allele, so each is written as its own row.
# -----------------------------------------------------------------------------
ANN_FIELDS = {
    "ALLELE": 0,
    "TYPE": 1,
    "IMPACT": 2,
    "GENE_NAME": 3,
    "GENE_ID": 4,
    "FEATURE_ID": 6,
    "HGVS_CODING": 9,
    "HGVS_PROTEIN": 10,
}

MIN_ANN_PARTS = 11


def open_vcf(path):
    """Open a VCF whether or not it is gzipped."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def find_vcfs(directory, suffixes=(".ann.vcf", ".ann.vcf.gz")):
    """Return {sample: path} for annotated VCFs in a directory.

    The sample name is the filename with the suffix removed. Keeping this in
    one place means the genotype and annotation tables always derive sample
    names the same way, so the two can be joined.
    """
    if not os.path.isdir(directory):
        sys.exit("ERROR: directory not found: %s" % directory)

    found = {}
    for filename in sorted(os.listdir(directory)):
        for suffix in suffixes:
            if filename.endswith(suffix):
                found[filename[: -len(suffix)]] = os.path.join(directory, filename)
                break

    if not found:
        sys.exit(
            "ERROR: no files ending in %s found in %s"
            % (" or ".join(suffixes), directory)
        )
    return found


def parse_ann(info):
    """Yield one dict per SnpEff ANN entry in an INFO string."""
    for entry in info.split(";"):
        if not entry.startswith("ANN="):
            continue
        for annotation in entry[len("ANN="):].split(","):
            parts = annotation.split("|")
            if len(parts) < MIN_ANN_PARTS:
                continue
            yield {name: parts[index] for name, index in ANN_FIELDS.items()}
        return


def parse_format(format_field, sample_field):
    """Pair a VCF FORMAT string with its sample column."""
    return dict(zip(format_field.split(":"), sample_field.split(":")))


def resolve(path, env_name, default):
    """Allow any path to be overridden by an environment variable.

    Lets a reader point the scripts at their own data without editing them:
        GENOTYPE_OUTPUT=/my/path/out.tsv python3 extract_genotypes.py
    """
    return os.environ.get(env_name, path if path else default)
