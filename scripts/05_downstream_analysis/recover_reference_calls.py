#!/usr/bin/env python3
"""
recover_reference_calls.py — recover wild-type calls from Clair3 gVCFs

WHAT IT DOES
    A Clair3 VCF lists only variant sites. If a sample has no record at a
    drug-resistance codon, that could mean either "confidently wild-type" or
    "no usable data", and the VCF cannot tell them apart. Every genotype in
    the VCF therefore carries GT=1; there are no GT=0 calls to count.

    The gVCF does distinguish them. It records confident reference blocks:

        Pf3D7_13_v3  1724851  .  C  <NON_REF>  0  .  END=1725959
                     GT:GQ:MIN_DP:PL  0:50:933:0,2846,28406

    one record covering positions 1724851-1725959, with GT=0, a genotype
    quality, and the minimum depth observed anywhere in the block.

    This script walks each sample's gVCF and, for every site where a variant
    was seen in any sample, classifies that sample as:

        variant    a PASS variant record at the position
        reference  covered by a reference block with GT=0 meeting the quality
                   and depth thresholds
        missing    neither: no block, GT='.', or below threshold

    The result is the denominator needed for wild-type frequency.

-------------------------------------------------------------------------------
PROVENANCE                                             <<< PLEASE READ

    The original script that performed this step was not retained. This is a
    reimplementation of the described procedure: reference calls taken from
    the gVCFs and combined with the variants extracted from the VCFs.

    VALIDATE IT BEFORE RELYING ON IT. Compare the wild-type frequencies this
    produces at the kelch13, pfcrt, pfmdr1, pfdhfr and pfdhps marker codons
    against the published values in Figure 2B and Figure 3A. Agreement
    confirms the procedure; disagreement means the original differed.
-------------------------------------------------------------------------------

INPUTS
    --gvcf-dir   folder of per-sample Clair3 output folders (<Sample>_clair3/
                 containing merge_output.gvcf.gz), or of <Sample>.gvcf.gz files
    --regions    regions file, chrom/start/end, tab-separated, '#' comments ok

OUTPUT
    --out        long-format TSV, one row per sample per site:
                 SAMPLE CHROM POS REF ALT STATE GT GQ DP AF

THRESHOLDS
    --min-gq     minimum GQ for a reference call to count    (default 20)
    --min-dp     minimum MIN_DP for a reference call to count (default 10)

    Clair3 caps the GQ of reference blocks at 50, so a threshold above 50
    would reject every reference call.

USAGE
    python3 scripts/05_downstream_analysis/recover_reference_calls.py \\
        --gvcf-dir /path/to/clair3_output \\
        --regions  data/metadata/regions.txt \\
        --out      data/processed/reference_and_variant_calls.tsv

Part of: Diallo et al., Nature Communications (2026)
         https://github.com/mtd1997/Adaptive_Pf_resurgence_Rwanda
"""

import argparse
import bisect
import csv
import glob
import gzip
import os
import sys


def open_maybe_gzip(path):
    """Open a file whether or not it is gzipped."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def read_regions(path):
    """Read a regions file into [(chrom, start, end), ...]."""
    regions = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 3:
                continue
            regions.append((fields[0], int(fields[1]), int(fields[2])))
    if not regions:
        sys.exit("ERROR: no usable regions found in %s" % path)
    return regions


def in_regions(chrom, pos, regions):
    for region_chrom, start, end in regions:
        if chrom == region_chrom and start <= pos <= end:
            return True
    return False


def find_gvcfs(gvcf_dir):
    """Locate per-sample gVCFs, returning {sample: path}.

    Handles both layouts: Clair3's own <Sample>_clair3/merge_output.gvcf.gz,
    and flat <Sample>.gvcf.gz files.
    """
    found = {}

    for path in sorted(glob.glob(os.path.join(gvcf_dir, "*_clair3", "merge_output.gvcf.gz"))):
        sample = os.path.basename(os.path.dirname(path))
        if sample.endswith("_clair3"):
            sample = sample[: -len("_clair3")]
        found[sample] = path

    if not found:
        for pattern in ("*.gvcf.gz", "*.g.vcf.gz", "*.gvcf"):
            for path in sorted(glob.glob(os.path.join(gvcf_dir, pattern))):
                sample = os.path.basename(path)
                for suffix in (".gvcf.gz", ".g.vcf.gz", ".gvcf"):
                    if sample.endswith(suffix):
                        sample = sample[: -len(suffix)]
                        break
                found[sample] = path

    return found


def parse_gvcf(path, regions):
    """Read one gVCF.

    Returns (variants, blocks):
        variants  {(chrom, pos): dict}   PASS variant records in the regions
        blocks    {chrom: (starts, records)}  reference blocks, sorted by start
    """
    variants = {}
    blocks = {}

    with open_maybe_gzip(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue

            chrom, pos_text, _, ref, alt, qual, filt, info, fmt, sample_col = fields[:10]
            pos = int(pos_text)

            format_keys = fmt.split(":")
            sample_values = sample_col.split(":")
            call = dict(zip(format_keys, sample_values))
            genotype = call.get("GT", ".")

            is_reference_block = "END=" in info and alt == "<NON_REF>"

            if is_reference_block:
                end = int(info.split("END=")[1].split(";")[0])
                # Skip blocks that do not touch any region of interest
                if not any(chrom == rc and not (end < s or pos > e)
                           for rc, s, e in regions):
                    continue
                blocks.setdefault(chrom, []).append({
                    "start": pos,
                    "end": end,
                    "GT": genotype,
                    "GQ": call.get("GQ", "."),
                    "MIN_DP": call.get("MIN_DP", "."),
                })
            else:
                if not in_regions(chrom, pos, regions):
                    continue
                # Clair3 appends <NON_REF> to the ALT list; drop it.
                alleles = [a for a in alt.split(",") if a != "<NON_REF>"]
                if not alleles:
                    continue
                variants[(chrom, pos)] = {
                    "REF": ref,
                    "ALT": ",".join(alleles),
                    "QUAL": qual,
                    "FILTER": filt,
                    "GT": genotype,
                    "GQ": call.get("GQ", "."),
                    "DP": call.get("DP", "."),
                    "AF": call.get("AF", "."),
                }

    # Sort blocks so they can be searched by position
    for chrom in blocks:
        blocks[chrom].sort(key=lambda b: b["start"])
    indexed = {c: ([b["start"] for b in bl], bl) for c, bl in blocks.items()}

    return variants, indexed


def lookup_block(blocks, chrom, pos):
    """Return the reference block covering pos, or None."""
    if chrom not in blocks:
        return None
    starts, records = blocks[chrom]
    # rightmost block whose start is <= pos
    index = bisect.bisect_right(starts, pos) - 1
    if index < 0:
        return None
    block = records[index]
    if block["start"] <= pos <= block["end"]:
        return block
    return None


def as_int(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def main():
    parser = argparse.ArgumentParser(
        description="Recover wild-type calls from Clair3 gVCFs.")
    parser.add_argument("--gvcf-dir", required=True)
    parser.add_argument("--regions", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--min-gq", type=int, default=20)
    parser.add_argument("--min-dp", type=int, default=10)
    args = parser.parse_args()

    regions = read_regions(args.regions)
    gvcfs = find_gvcfs(args.gvcf_dir)

    if not gvcfs:
        sys.exit("ERROR: no gVCFs found under %s" % args.gvcf_dir)

    print("Samples        : %d" % len(gvcfs))
    print("Regions        : %d" % len(regions))
    print("Reference call : GT=0, GQ >= %d, MIN_DP >= %d"
          % (args.min_gq, args.min_dp))
    print()

    # ---- pass 1: every site variant in at least one sample -------------------
    print("Pass 1 of 2: collecting variant sites across all samples")
    per_sample = {}
    all_sites = {}

    for index, (sample, path) in enumerate(sorted(gvcfs.items()), start=1):
        variants, blocks = parse_gvcf(path, regions)
        per_sample[sample] = (variants, blocks)
        for key, record in variants.items():
            if record["FILTER"] == "PASS":
                all_sites.setdefault(key, (record["REF"], record["ALT"]))
        if index % 100 == 0:
            print("  ...%d/%d" % (index, len(gvcfs)))

    sites = sorted(all_sites)
    print("  %d variant sites found\n" % len(sites))

    if not sites:
        sys.exit("ERROR: no PASS variants in any sample within the regions.")

    # ---- pass 2: classify every sample at every site --------------------------
    print("Pass 2 of 2: classifying each sample at each site")
    counts = {"variant": 0, "reference": 0, "missing": 0}

    with open(args.out, "w", newline="") as out_handle:
        writer = csv.writer(out_handle, delimiter="\t")
        writer.writerow(["SAMPLE", "CHROM", "POS", "REF", "ALT",
                         "STATE", "GT", "GQ", "DP", "AF"])

        for sample in sorted(per_sample):
            variants, blocks = per_sample[sample]

            for (chrom, pos) in sites:
                ref, alt = all_sites[(chrom, pos)]

                record = variants.get((chrom, pos))
                if record is not None and record["FILTER"] == "PASS":
                    state = "variant"
                    row = [state, record["GT"], record["GQ"],
                           record["DP"], record["AF"]]
                else:
                    block = lookup_block(blocks, chrom, pos)
                    if block is None:
                        state, row = "missing", ["missing", ".", ".", ".", "."]
                    else:
                        gq = as_int(block["GQ"])
                        dp = as_int(block["MIN_DP"])
                        if (block["GT"] == "0"
                                and gq is not None and gq >= args.min_gq
                                and dp is not None and dp >= args.min_dp):
                            state = "reference"
                        else:
                            state = "missing"
                        row = [state, block["GT"], block["GQ"],
                               block["MIN_DP"], "."]

                counts[state] += 1
                writer.writerow([sample, chrom, pos, ref, alt] + row)

    total = sum(counts.values())
    print()
    print("=" * 60)
    print("Calls written : %d" % total)
    for state in ("variant", "reference", "missing"):
        pct = 100.0 * counts[state] / total if total else 0.0
        print("  %-10s %9d  (%5.1f%%)" % (state, counts[state], pct))
    print()
    print("Output        : %s" % args.out)
    print()
    print("VALIDATE: compare wild-type frequencies from this table against")
    print("          the published Figure 2B and Figure 3A values before use.")
    print("=" * 60)


if __name__ == "__main__":
    main()
