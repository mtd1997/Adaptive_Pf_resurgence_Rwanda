################################################################################
# FIGURE 5 — SELECTION-SCAN STATISTICS
################################################################################
#
# WHAT IT DOES
#   Compares allele frequencies between the two sampling periods and writes
#   one row per variant with:
#
#     - the alternate-allele frequency in each period
#     - the signed and absolute change between them
#     - a per-SNP F_ST
#     - empirical 95th-percentile thresholds and the resulting outlier flags
#
#   This is the only script in the repository that computes selection
#   statistics. The figure script that follows it does no arithmetic at all:
#   it reads this script's output table and draws it.
#
# WHAT IT DOES NOT DO
#   It does not change the statistical method. The genotype encoding and the
#   F_ST formula are carried over from the original analysis unchanged, so
#   the numbers reproduce the published figure. Comments explain the choices;
#   they do not alter them.
#
# INPUTS
#   1. The filtered, annotated selection-scan VCF
#   2. The sample metadata table
#   3. The annotation table from build_selection_scan_annotation.py
#
# OUTPUTS  (all under results/tables/figure05/)
#   selection_statistics_all_variants.tsv   one row per variant — the main
#                                           table, and the Source Data for
#                                           Figure 5
#   empirical_thresholds.tsv                the two 95th-percentile cutoffs
#   fst_outliers.tsv                        top 5% by F_ST
#   delta_af_outliers.tsv                   top 5% by |ΔAF|
#   double_outliers.tsv                     outliers on both statistics
#   resistance_marker_frequencies.tsv       frequency change per known marker
#   summary_by_gene.tsv                     per-gene counts
#   sample_counts.tsv                       samples per period
#   vcf_annotation_matching.tsv             variants in one input but not the
#                                           other — should be empty
#   sessionInfo.txt                         package versions
#
# PERIOD LABELS
#   The manuscript uses "Baseline-Kigali" and "Resurgence". Older metadata
#   files use "Pre-Resurgence" and "Post-Resurgence". Both are accepted; the
#   published labels are applied in memory and the original value is kept in
#   a `period_original` column so nothing is silently rewritten.
#
# HOW TO RUN
#   From the top of the repository:
#
#     Rscript scripts/06_figures/figure05_selection_scan_statistics.R
#
#   Paths come from config.R. Nothing is hard-coded to one machine.
#
################################################################################


# ==============================================================================
# 1. SETUP
# ==============================================================================

options(error = NULL)

suppressPackageStartupMessages({
  library(vcfR)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(tibble)
})


# ------------------------------------------------------------------------------
# Locate the repository root and read the shared paths.
#
# This looks for config.R in the current directory and then in the parent
# directories above it, so the script works whether it is run from the top of
# the repository or from inside scripts/06_figures/.
# ------------------------------------------------------------------------------

find_project_root <- function() {

  candidate <- normalizePath(getwd(), mustWork = TRUE)

  for (level in seq_len(5)) {

    if (file.exists(file.path(candidate, "config.R"))) {
      return(candidate)
    }

    parent <- dirname(candidate)

    if (identical(parent, candidate)) {
      break
    }

    candidate <- parent
  }

  stop(
    "Could not find config.R. Run this script from the top of the ",
    "repository:\n  Rscript scripts/06_figures/",
    "figure05_selection_scan_statistics.R",
    call. = FALSE
  )
}

PROJECT_ROOT <- find_project_root()

source(file.path(PROJECT_ROOT, "config.R"))


# ------------------------------------------------------------------------------
# Inputs and outputs
#
# Each of these falls back to a path inside the repository if config.R does
# not define it, so the script is runnable immediately after a clean clone.
# ------------------------------------------------------------------------------

path_or_default <- function(variable_name, default_relative_path) {

  if (exists(variable_name, inherits = TRUE)) {

    value <- get(variable_name, inherits = TRUE)

    if (is.character(value) && length(value) == 1 && nzchar(value)) {
      return(value)
    }
  }

  file.path(PROJECT_ROOT, default_relative_path)
}

VCF_FILE <- path_or_default(
  "SELECTION_VCF",
  "data/processed/selection_scan_annotated.vcf.gz"
)

METADATA_FILE <- path_or_default(
  "METADATA_FILE",
  "data/metadata/samples_metadata.tsv"
)

ANNOTATION_FILE <- path_or_default(
  "SELECTION_ANNOTATION",
  "results/tables/annotated_variants_snpeff_final.tsv"
)

OUTPUT_DIR <- file.path(PROJECT_ROOT, "results", "tables", "figure05")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------------
# Check the inputs before doing any work, and say exactly which is missing.
# ------------------------------------------------------------------------------

input_files <- c(
  VCF        = VCF_FILE,
  metadata   = METADATA_FILE,
  annotation = ANNOTATION_FILE
)

missing_inputs <- input_files[!file.exists(input_files)]

if (length(missing_inputs) > 0) {
  stop(
    "These input files were not found:\n",
    paste0("  ", names(missing_inputs), ": ", missing_inputs, collapse = "\n"),
    "\n\nEdit the paths in config.R, or place the files where they are ",
    "expected.",
    call. = FALSE
  )
}

message("Inputs:")
for (name in names(input_files)) {
  message(sprintf("  %-11s %s", name, input_files[[name]]))
}


# ==============================================================================
# 2. READ THE VCF
# ==============================================================================

vcf <- vcfR::read.vcfR(VCF_FILE, verbose = FALSE)

message(sprintf(
  "\nVCF: %d variants, %d samples",
  nrow(vcf@fix),
  ncol(vcf@gt) - 1L
))


# ==============================================================================
# 3. READ AND TIDY THE METADATA
# ==============================================================================

metadata <- readr::read_tsv(METADATA_FILE, show_col_types = FALSE)

# Metadata files in this project have used both upper- and lower-case names
# for these two columns. Accept either.
if ("PERIOD" %in% names(metadata) && !"period" %in% names(metadata)) {
  metadata <- dplyr::rename(metadata, period = PERIOD)
}

if ("REGION" %in% names(metadata) && !"region" %in% names(metadata)) {
  metadata <- dplyr::rename(metadata, region = REGION)
}

required_metadata_columns <- c("SAMPLE", "period", "region")

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  names(metadata)
)

if (length(missing_metadata_columns) > 0) {
  stop(
    "The metadata is missing: ",
    paste(missing_metadata_columns, collapse = ", "),
    call. = FALSE
  )
}

metadata <- metadata %>%
  dplyr::mutate(
    dplyr::across(
      c(SAMPLE, period, region),
      ~ stringr::str_trim(as.character(.x))
    )
  ) %>%
  dplyr::filter(!is.na(SAMPLE), SAMPLE != "") %>%
  dplyr::distinct(SAMPLE, .keep_all = TRUE)


# ------------------------------------------------------------------------------
# Apply the published period labels, keeping the original values.
# ------------------------------------------------------------------------------

BASELINE_GROUP   <- "Baseline-Kigali"
RESURGENCE_GROUP <- "Resurgence"

metadata <- metadata %>%
  dplyr::mutate(
    period_original = period,
    # recode() leaves any value it does not match unchanged, so metadata
    # files that already carry the published labels pass through untouched.
    period = dplyr::recode(
      period,
      "Pre-Resurgence"  = BASELINE_GROUP,
      "Post-Resurgence" = RESURGENCE_GROUP
    ),
    group = period
  )

unexpected_periods <- setdiff(
  unique(stats::na.omit(metadata$period)),
  c(BASELINE_GROUP, RESURGENCE_GROUP)
)

if (length(unexpected_periods) > 0) {
  warning(
    "Unrecognised period labels, which will be excluded from the ",
    "comparison: ",
    paste(unexpected_periods, collapse = ", "),
    call. = FALSE
  )
}

message("\nSamples per period in the metadata:")
print(table(metadata$period, useNA = "ifany"))


# ==============================================================================
# 4. READ THE ANNOTATION
# ==============================================================================

annotation <- readr::read_tsv(ANNOTATION_FILE, show_col_types = FALSE)

required_annotation_columns <- c(
  "CHROM", "POS", "REF", "ALT",
  "GENE_NAME", "EFFECT", "PROTEIN_CHANGE",
  "IS_RESISTANCE", "LABEL", "variant_id"
)

missing_annotation_columns <- setdiff(
  required_annotation_columns,
  names(annotation)
)

if (length(missing_annotation_columns) > 0) {
  stop(
    "The annotation table is missing: ",
    paste(missing_annotation_columns, collapse = ", "),
    "\n\nRegenerate it with:\n",
    "  python3 scripts/05_downstream_analysis/",
    "build_selection_scan_annotation.py --help",
    call. = FALSE
  )
}

annotation <- annotation %>%
  dplyr::mutate(
    CHROM          = as.character(CHROM),
    POS            = as.integer(POS),
    REF            = stringr::str_to_upper(as.character(REF)),
    ALT            = stringr::str_to_upper(as.character(ALT)),
    IS_RESISTANCE  = tidyr::replace_na(as.logical(IS_RESISTANCE), FALSE),
    snp_id         = paste0(CHROM, "_", POS),
    variant_id_check = paste(CHROM, POS, REF, ALT, sep = "_")
  )

# The variant_id column is what every downstream join uses. Verify it really
# is CHROM_POS_REF_ALT rather than trusting the label.
id_mismatches <- dplyr::filter(annotation, variant_id != variant_id_check)

if (nrow(id_mismatches) > 0) {

  readr::write_tsv(
    id_mismatches,
    file.path(OUTPUT_DIR, "annotation_variant_id_mismatches.tsv")
  )

  stop(
    nrow(id_mismatches),
    " annotation rows have a variant_id that is not CHROM_POS_REF_ALT. ",
    "See annotation_variant_id_mismatches.tsv.",
    call. = FALSE
  )
}

annotation <- dplyr::select(annotation, -variant_id_check)

duplicated_annotation_ids <- annotation %>%
  dplyr::count(variant_id, name = "n") %>%
  dplyr::filter(n > 1)

if (nrow(duplicated_annotation_ids) > 0) {
  stop(
    nrow(duplicated_annotation_ids),
    " variant IDs appear more than once in the annotation. Each variant ",
    "must have exactly one annotation row.",
    call. = FALSE
  )
}


# ==============================================================================
# 5. EXTRACT GENOTYPES AND ALIGN SAMPLES
# ==============================================================================

genotype_characters <- vcfR::extract.gt(vcf, element = "GT", as.numeric = FALSE)

if (is.null(genotype_characters)) {
  stop("Could not read the GT field from the VCF.", call. = FALSE)
}

depth_matrix <- vcfR::extract.gt(vcf, element = "DP", as.numeric = TRUE)

if (is.null(depth_matrix)) {
  warning(
    "No DP field in the VCF; mean depth will be reported as NA.",
    call. = FALSE
  )
}


# ------------------------------------------------------------------------------
# Keep only samples present in both the VCF and the metadata, in the VCF's
# order. This is the step that drops the three samples found not to be of
# Rwandan origin during the metadata check.
# ------------------------------------------------------------------------------

vcf_samples      <- colnames(genotype_characters)
metadata_samples <- metadata$SAMPLE
shared_samples   <- intersect(vcf_samples, metadata_samples)

message("\nSample reconciliation:")
message(sprintf("  in the VCF:            %d", length(vcf_samples)))
message(sprintf("  in the metadata:       %d", length(metadata_samples)))
message(sprintf("  retained (in both):    %d", length(shared_samples)))
message(sprintf(
  "  VCF only, dropped:     %d",
  length(setdiff(vcf_samples, metadata_samples))
))
message(sprintf(
  "  metadata only, dropped:%d",
  length(setdiff(metadata_samples, vcf_samples))
))

if (length(shared_samples) == 0) {
  stop(
    "No sample names are shared between the VCF and the metadata. Check ",
    "that both use the same naming convention.",
    call. = FALSE
  )
}

genotype_characters <- genotype_characters[, shared_samples, drop = FALSE]

if (!is.null(depth_matrix)) {
  depth_matrix <- depth_matrix[, shared_samples, drop = FALSE]
}

metadata <- metadata[match(shared_samples, metadata$SAMPLE), , drop = FALSE]

stopifnot(identical(colnames(genotype_characters), metadata$SAMPLE))

if (!is.null(depth_matrix)) {
  stopifnot(identical(colnames(depth_matrix), metadata$SAMPLE))
}

message("\nSamples per period after reconciliation:")
print(table(metadata$group, useNA = "ifany"))


# ==============================================================================
# 6. ENCODE GENOTYPES
# ==============================================================================
#
# P. falciparum is haploid at these loci, so each genotype is a single
# allele index rather than a diploid pair.
#
#   "."  or missing  ->  NA
#   "0"              ->  0   reference allele
#   any other number ->  1   alternate allele
#
# Carried over from the original analysis unchanged.
# ==============================================================================

encode_genotypes <- function(genotype_matrix) {

  missing_codes <- c(".", "", "NA", "./.", ".|.")

  encoded <- apply(genotype_matrix, 2, function(column) {

    column[column %in% missing_codes] <- NA_character_

    numeric_column <- suppressWarnings(as.numeric(column))

    ifelse(
      is.na(numeric_column), NA_real_,
      ifelse(numeric_column == 0, 0, 1)
    )
  })

  dimnames(encoded) <- dimnames(genotype_matrix)

  encoded
}

genotypes <- encode_genotypes(genotype_characters)


# ------------------------------------------------------------------------------
# Variant identifiers, taken from the VCF itself.
# ------------------------------------------------------------------------------

vcf_chrom <- as.character(vcf@fix[, "CHROM"])
vcf_pos   <- as.integer(vcf@fix[, "POS"])
vcf_ref   <- stringr::str_to_upper(as.character(vcf@fix[, "REF"]))
vcf_alt   <- stringr::str_to_upper(as.character(vcf@fix[, "ALT"]))

snp_ids     <- paste0(vcf_chrom, "_", vcf_pos)
variant_ids <- paste(vcf_chrom, vcf_pos, vcf_ref, vcf_alt, sep = "_")

rownames(genotypes) <- snp_ids

if (!is.null(depth_matrix)) {
  rownames(depth_matrix) <- snp_ids
}

if (anyDuplicated(snp_ids) > 0) {
  stop(
    "The VCF contains more than one record at the same CHROM and POS. ",
    "The downstream joins use CHROM_POS and would be ambiguous.",
    call. = FALSE
  )
}

if (anyDuplicated(variant_ids) > 0) {
  stop("The VCF contains duplicated CHROM_POS_REF_ALT records.", call. = FALSE)
}


# Samples as rows, variants as columns — the orientation the frequency and
# F_ST calculations below expect.
genotypes_by_sample <- t(genotypes)

stopifnot(identical(rownames(genotypes_by_sample), metadata$SAMPLE))

message(sprintf(
  "\nGenotype matrix: %d samples x %d variants, %.2f%% missing",
  nrow(genotypes_by_sample),
  ncol(genotypes_by_sample),
  100 * mean(is.na(genotypes_by_sample))
))


# ==============================================================================
# 7. PER-VARIANT QUALITY SUMMARIES
# ==============================================================================

variant_table <- tibble::tibble(
  snp_id     = snp_ids,
  variant_id = variant_ids,
  chrom      = vcf_chrom,
  pos        = vcf_pos,
  ref        = vcf_ref,
  alt        = vcf_alt,
  qual       = suppressWarnings(as.numeric(vcf@fix[, "QUAL"]))
)

alt_frequency_overall <- rowMeans(genotypes, na.rm = TRUE)
alt_frequency_overall[is.nan(alt_frequency_overall)] <- NA_real_

mean_depth <- if (is.null(depth_matrix)) {
  rep(NA_real_, nrow(genotypes))
} else {
  depths <- rowMeans(depth_matrix, na.rm = TRUE)
  depths[is.nan(depths)] <- NA_real_
  depths
}

variant_quality <- variant_table %>%
  dplyr::mutate(
    missing_rate          = rowMeans(is.na(genotypes)),
    overall_alt_frequency = alt_frequency_overall,
    overall_maf           = pmin(
      alt_frequency_overall,
      1 - alt_frequency_overall
    ),
    mean_depth            = mean_depth,
    n_called_overall      = rowSums(!is.na(genotypes))
  )


# ==============================================================================
# 8. CHECK THE VCF AND ANNOTATION DESCRIBE THE SAME VARIANTS
# ==============================================================================

matching_report <- variant_table %>%
  dplyr::select(variant_id, vcf_chrom = chrom, vcf_pos = pos) %>%
  dplyr::full_join(
    annotation %>%
      dplyr::select(variant_id, GENE_NAME, EFFECT, LABEL, IS_RESISTANCE),
    by = "variant_id"
  ) %>%
  dplyr::mutate(
    in_vcf        = !is.na(vcf_chrom),
    in_annotation = !is.na(GENE_NAME)
  )

readr::write_tsv(
  matching_report,
  file.path(OUTPUT_DIR, "vcf_annotation_matching.tsv")
)

n_vcf_only        <- sum(matching_report$in_vcf & !matching_report$in_annotation)
n_annotation_only <- sum(!matching_report$in_vcf & matching_report$in_annotation)

message(sprintf("\nVCF variants with no annotation:    %d", n_vcf_only))
message(sprintf("Annotated variants not in the VCF:  %d", n_annotation_only))

if (n_vcf_only > 0 || n_annotation_only > 0) {
  warning(
    "The VCF and the annotation table do not cover the same variants. ",
    "Unannotated variants will be dropped by the figure script. See ",
    "vcf_annotation_matching.tsv.",
    call. = FALSE
  )
}


# ==============================================================================
# 9. ALLELE FREQUENCIES PER PERIOD
# ==============================================================================

group_counts <- table(metadata$group, useNA = "ifany")

for (required_group in c(BASELINE_GROUP, RESURGENCE_GROUP)) {

  if (!required_group %in% names(group_counts)) {
    stop(
      "No samples were found for the '", required_group, "' period.",
      call. = FALSE
    )
  }

  if (group_counts[[required_group]] < 3) {
    stop(
      "Fewer than three samples in the '", required_group, "' period; ",
      "allele frequencies would not be meaningful.",
      call. = FALSE
    )
  }
}

baseline_rows   <- which(metadata$group == BASELINE_GROUP)
resurgence_rows <- which(metadata$group == RESURGENCE_GROUP)

column_frequency <- function(matrix_subset) {
  frequencies <- colMeans(matrix_subset, na.rm = TRUE)
  frequencies[is.nan(frequencies)] <- NA_real_
  frequencies
}

af_baseline <- column_frequency(
  genotypes_by_sample[baseline_rows, , drop = FALSE]
)

af_resurgence <- column_frequency(
  genotypes_by_sample[resurgence_rows, , drop = FALSE]
)

n_called_baseline <- colSums(
  !is.na(genotypes_by_sample[baseline_rows, , drop = FALSE])
)

n_called_resurgence <- colSums(
  !is.na(genotypes_by_sample[resurgence_rows, , drop = FALSE])
)

delta_af     <- af_resurgence - af_baseline
abs_delta_af <- abs(delta_af)


# ==============================================================================
# 10. PER-SNP F_ST
# ==============================================================================
#
# A Weir and Cockerham style estimator for two groups of haploid samples.
# Reproduced from the original analysis without modification.
#
#   p1, p2   alternate-allele frequency in each group
#   n1, n2   number of successfully called samples in each group
#   nc       variance-corrected sample size
#   MST      among-group component
#   MSG      within-group component
#
# Notes on the edge cases, which are deliberate:
#   - fewer than three called samples in either group -> NA
#   - the variant is fixed the same way in both groups -> 0
#   - a non-positive denominator -> NA
#   - results are clamped to [0, 1]
#
# ==============================================================================

compute_per_snp_fst <- function(
  genotype_matrix,
  group_one_rows,
  group_two_rows
) {

  fst_values <- vapply(
    seq_len(ncol(genotype_matrix)),
    function(column_index) {

      group_one <- genotype_matrix[group_one_rows, column_index]
      group_two <- genotype_matrix[group_two_rows, column_index]

      group_one <- group_one[!is.na(group_one)]
      group_two <- group_two[!is.na(group_two)]

      if (length(group_one) < 3 || length(group_two) < 3) {
        return(NA_real_)
      }

      p1 <- mean(group_one)
      p2 <- mean(group_two)

      n1 <- length(group_one)
      n2 <- length(group_two)

      # Fixed identically in both groups: no differentiation.
      if ((p1 == 0 && p2 == 0) || (p1 == 1 && p2 == 1)) {
        return(0)
      }

      p_bar <- (n1 * p1 + n2 * p2) / (n1 + n2)

      nc <- (n1 + n2) - (n1^2 + n2^2) / (n1 + n2)

      mst <- n1 * (p1 - p_bar)^2 + n2 * (p2 - p_bar)^2

      msg <- (n1 * p1 * (1 - p1) + n2 * p2 * (1 - p2)) / (n1 + n2 - 2)

      denominator <- mst + (nc - 1) * msg

      if (is.na(denominator) || denominator <= 0) {
        return(NA_real_)
      }

      max(0, min(1, (mst - msg) / denominator))
    },
    numeric(1)
  )

  tibble::tibble(
    snp_id     = colnames(genotype_matrix),
    fst        = fst_values,
    comparison = paste0(BASELINE_GROUP, " vs ", RESURGENCE_GROUP)
  )
}

fst_results <- compute_per_snp_fst(
  genotypes_by_sample,
  baseline_rows,
  resurgence_rows
)


# ==============================================================================
# 11. ASSEMBLE THE STATISTICS TABLE
# ==============================================================================

frequency_table <- tibble::tibble(
  snp_id                   = colnames(genotypes_by_sample),
  AF_Baseline_Kigali       = as.numeric(af_baseline),
  AF_Resurgence            = as.numeric(af_resurgence),
  N_called_Baseline_Kigali = as.integer(n_called_baseline),
  N_called_Resurgence      = as.integer(n_called_resurgence),
  delta_AF                 = as.numeric(delta_af),
  abs_delta_AF             = as.numeric(abs_delta_af)
)

annotation_for_join <- annotation %>%
  dplyr::transmute(
    variant_id,
    gene             = GENE_NAME,
    effect           = EFFECT,
    protein_change   = PROTEIN_CHANGE,
    is_resistance    = IS_RESISTANCE,
    snp_label        = LABEL,
    biological_label = dplyr::coalesce(
      if ("BIOLOGICAL_LABEL" %in% names(annotation)) BIOLOGICAL_LABEL
      else NA_character_,
      LABEL
    ),
    genomic_label    = if ("GENOMIC_LABEL" %in% names(annotation)) {
      GENOMIC_LABEL
    } else {
      paste0(GENE_NAME, ":", POS, "_", REF, ">", ALT)
    },
    hgvs_c           = if ("SNPEFF_HGVS_C" %in% names(annotation)) {
      SNPEFF_HGVS_C
    } else {
      NA_character_
    },
    hgvs_p           = if ("SNPEFF_HGVS_P" %in% names(annotation)) {
      SNPEFF_HGVS_P
    } else {
      NA_character_
    },
    snpeff_impact    = if ("SNPEFF_IMPACT" %in% names(annotation)) {
      SNPEFF_IMPACT
    } else {
      NA_character_
    }
  )

selection_statistics <- variant_table %>%
  dplyr::left_join(frequency_table, by = "snp_id") %>%
  dplyr::left_join(fst_results, by = "snp_id") %>%
  dplyr::left_join(
    dplyr::select(
      variant_quality,
      snp_id, missing_rate, overall_alt_frequency,
      overall_maf, mean_depth, n_called_overall
    ),
    by = "snp_id"
  ) %>%
  dplyr::left_join(annotation_for_join, by = "variant_id")


# ==============================================================================
# 12. EMPIRICAL 95TH-PERCENTILE THRESHOLDS
# ==============================================================================
#
# There is no null model here. A variant is called an outlier if it sits in
# the top 5% of the observed distribution of that statistic across the
# variants analysed. This is an empirical ranking, not a significance test,
# and the figure legend should describe it that way.
# ==============================================================================

fst_threshold <- as.numeric(stats::quantile(
  selection_statistics$fst,
  probs = 0.95, na.rm = TRUE, names = FALSE, type = 7
))

delta_threshold <- as.numeric(stats::quantile(
  selection_statistics$abs_delta_AF,
  probs = 0.95, na.rm = TRUE, names = FALSE, type = 7
))

for (threshold_check in list(
  list(name = "F_ST", value = fst_threshold),
  list(name = "|ΔAF|", value = delta_threshold)
)) {
  if (length(threshold_check$value) != 1 ||
      !is.finite(threshold_check$value)) {
    stop(
      "The ", threshold_check$name, " threshold did not resolve to a single ",
      "finite number. Are there enough non-missing values?",
      call. = FALSE
    )
  }
}

message(sprintf("\nEmpirical 95th-percentile thresholds:"))
message(sprintf("  F_ST   %.6f", fst_threshold))
message(sprintf("  |ΔAF|  %.6f", delta_threshold))

selection_statistics <- selection_statistics %>%
  dplyr::mutate(
    fst_threshold_95       = fst_threshold,
    abs_delta_threshold_95 = delta_threshold,

    fst_outlier   = !is.na(fst) & fst >= fst_threshold,
    delta_outlier = !is.na(abs_delta_AF) & abs_delta_AF >= delta_threshold,

    double_outlier = fst_outlier & delta_outlier,

    direction = dplyr::case_when(
      is.na(delta_AF) ~ NA_character_,
      delta_AF > 0    ~ "Increased in Resurgence",
      delta_AF < 0    ~ "Decreased in Resurgence",
      TRUE            ~ "No change"
    ),

    resistance_or_outlier =
      tidyr::replace_na(is_resistance, FALSE) | fst_outlier | delta_outlier
  ) %>%
  dplyr::select(
    variant_id, snp_id,
    chrom, pos, ref, alt,
    gene, effect, snpeff_impact,
    hgvs_c, hgvs_p, protein_change,
    snp_label, biological_label, genomic_label,
    is_resistance,
    qual, mean_depth, missing_rate,
    overall_alt_frequency, overall_maf,
    N_called_Baseline_Kigali, N_called_Resurgence, n_called_overall,
    AF_Baseline_Kigali, AF_Resurgence,
    delta_AF, abs_delta_AF, direction,
    fst,
    fst_threshold_95, abs_delta_threshold_95,
    fst_outlier, delta_outlier, double_outlier, resistance_or_outlier,
    comparison
  ) %>%
  dplyr::arrange(chrom, pos, ref, alt)


# ==============================================================================
# 13. WRITE THE OUTPUT TABLES
# ==============================================================================

write_table <- function(data, filename) {
  readr::write_tsv(data, file.path(OUTPUT_DIR, filename))
  invisible(data)
}

write_table(
  selection_statistics,
  "selection_statistics_all_variants.tsv"
)

write_table(
  tibble::tibble(
    Statistic = c("Per-SNP F_ST", "Absolute allele-frequency change"),
    Quantile  = c(0.95, 0.95),
    Threshold = c(fst_threshold, delta_threshold),
    Definition = c(
      sprintf(
        "Variants with F_ST >= %.6f are in the top 5%% of the observed F_ST distribution.",
        fst_threshold
      ),
      sprintf(
        "Variants with |delta AF| >= %.6f are in the top 5%% of the observed |delta AF| distribution.",
        delta_threshold
      )
    )
  ),
  "empirical_thresholds.tsv"
)

write_table(
  selection_statistics %>%
    dplyr::filter(fst_outlier) %>%
    dplyr::arrange(dplyr::desc(fst)),
  "fst_outliers.tsv"
)

write_table(
  selection_statistics %>%
    dplyr::filter(delta_outlier) %>%
    dplyr::arrange(dplyr::desc(abs_delta_AF)),
  "delta_af_outliers.tsv"
)

double_outliers <- write_table(
  selection_statistics %>%
    dplyr::filter(double_outlier) %>%
    dplyr::arrange(dplyr::desc(fst), dplyr::desc(abs_delta_AF)),
  "double_outliers.tsv"
)

resistance_frequencies <- write_table(
  selection_statistics %>%
    dplyr::filter(tidyr::replace_na(is_resistance, FALSE)) %>%
    dplyr::select(
      gene, snp_label, biological_label, variant_id, snp_id,
      N_called_Baseline_Kigali, N_called_Resurgence,
      AF_Baseline_Kigali, AF_Resurgence,
      delta_AF, abs_delta_AF, direction,
      fst, fst_outlier, delta_outlier, double_outlier
    ) %>%
    dplyr::arrange(gene, snp_label),
  "resistance_marker_frequencies.tsv"
)

gene_summary <- write_table(
  selection_statistics %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      N_variants            = dplyr::n(),
      N_with_fst            = sum(!is.na(fst)),
      N_with_delta_AF       = sum(!is.na(abs_delta_AF)),
      N_resistance_markers  = sum(tidyr::replace_na(is_resistance, FALSE)),
      N_fst_outliers        = sum(fst_outlier, na.rm = TRUE),
      N_delta_outliers      = sum(delta_outlier, na.rm = TRUE),
      N_double_outliers     = sum(double_outlier, na.rm = TRUE),
      max_fst               = suppressWarnings(max(fst, na.rm = TRUE)),
      max_abs_delta_AF      = suppressWarnings(max(abs_delta_AF, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      dplyr::across(
        c(max_fst, max_abs_delta_AF),
        ~ ifelse(is.finite(.x), .x, NA_real_)
      )
    ) %>%
    dplyr::arrange(gene),
  "summary_by_gene.tsv"
)

write_table(
  dplyr::count(metadata, group, name = "N_samples"),
  "sample_counts.tsv"
)

write_table(
  dplyr::filter(selection_statistics, is.na(gene)),
  "variants_without_annotation.tsv"
)


# ==============================================================================
# 14. SUMMARY TO THE CONSOLE
# ==============================================================================

message("\n", strrep("=", 72))
message("SELECTION-SCAN STATISTICS")
message(strrep("=", 72))

message(sprintf("Samples retained:           %d", nrow(metadata)))
message(sprintf("  %-24s  %d", BASELINE_GROUP, group_counts[[BASELINE_GROUP]]))
message(sprintf("  %-24s  %d", RESURGENCE_GROUP, group_counts[[RESURGENCE_GROUP]]))
message(sprintf("Variants analysed:          %d", nrow(selection_statistics)))
message(sprintf(
  "  with annotation:          %d",
  sum(!is.na(selection_statistics$gene))
))
message(sprintf(
  "F_ST outliers:              %d",
  sum(selection_statistics$fst_outlier, na.rm = TRUE)
))
message(sprintf(
  "|ΔAF| outliers:             %d",
  sum(selection_statistics$delta_outlier, na.rm = TRUE)
))
message(sprintf(
  "Outliers on both:           %d",
  sum(selection_statistics$double_outlier, na.rm = TRUE)
))
message(sprintf(
  "Resistance markers present: %d",
  sum(tidyr::replace_na(selection_statistics$is_resistance, FALSE))
))

message("\nPer gene:")
print(as.data.frame(gene_summary), row.names = FALSE)

message("\nResistance-marker frequency change:")
print(as.data.frame(resistance_frequencies), row.names = FALSE)

if (nrow(double_outliers) > 0) {
  message("\nOutliers on both statistics:")
  print(
    as.data.frame(dplyr::select(
      double_outliers, gene, snp_label, fst, delta_AF, direction
    )),
    row.names = FALSE
  )
}


# ==============================================================================
# 15. FINAL CHECKS
# ==============================================================================

if (nrow(selection_statistics) != nrow(vcf@fix)) {
  warning(
    "The statistics table does not have exactly one row per VCF record.",
    call. = FALSE
  )
}

if (dplyr::n_distinct(selection_statistics$variant_id) !=
    nrow(selection_statistics)) {
  warning(
    "The statistics table contains duplicated variant IDs.",
    call. = FALSE
  )
}


# ==============================================================================
# 16. SAVE INTERMEDIATES AND SESSION INFORMATION
# ==============================================================================
#
# The genotype matrix and metadata are saved so the PCoA and other analyses
# can reuse them without re-reading the VCF.
# ==============================================================================

saveRDS(
  genotypes_by_sample,
  file.path(OUTPUT_DIR, "genotypes_samples_by_variants.rds")
)

saveRDS(
  metadata,
  file.path(OUTPUT_DIR, "selection_scan_metadata.rds")
)

saveRDS(
  selection_statistics,
  file.path(OUTPUT_DIR, "selection_statistics_all_variants.rds")
)

capture.output(
  sessionInfo(),
  file = file.path(OUTPUT_DIR, "sessionInfo.txt")
)

message("\n", strrep("=", 72))
message("DONE")
message(strrep("=", 72))
message("Main table:")
message(sprintf(
  "  %s",
  file.path(OUTPUT_DIR, "selection_statistics_all_variants.tsv")
))
message("\nNext step — draw the figure:")
message("  Rscript scripts/06_figures/figure05_selection_scan_panels.R")
message(strrep("=", 72))
