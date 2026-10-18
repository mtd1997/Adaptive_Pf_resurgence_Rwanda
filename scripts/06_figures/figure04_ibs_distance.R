################################################################################
# FIGURE 4 — THE IBS DISTANCE MATRIX
################################################################################
#
# WHAT IT DOES
#   Computes the pairwise genetic distance matrix that Figure 4 and
#   Supplementary Figure S4 are built from, and saves it as an .rds file.
#
# WHAT THE DISTANCE IS
#   For every pair of samples, the proportion of genotypes that differ,
#   counted only over sites where BOTH samples were successfully called:
#
#       d(i,j) = 1 - (matching sites / sites called in both)
#
#   Missing data therefore reduces the denominator for that pair rather than
#   counting as a mismatch. Each pair is compared over whatever sites it has
#   in common, which is why the measure is described as pairwise-complete.
#
#   The original analysis called this "IBS distance", and the file name keeps
#   that label. Worth knowing when writing the Methods: the original script
#   also defined functions named `p_distance` and `hamming_dist`, and all
#   three compute exactly this same quantity. For any pair, every comparable
#   site either matches or does not, so
#
#       1 - matches/n   is identically   mismatches/n
#
#   The three names describe one measure. The Methods should describe the
#   computation above rather than naming a package's IBS implementation,
#   which is a different thing.
#
# WHY THIS SCRIPT EXISTS
#   The matrix was originally produced by a block inside a long exploratory
#   script, mixed in with PCA, trees and several abandoned distance measures.
#   This isolates just the computation that fed the published figure.
#
# THE SPEED CHANGE, AND WHY IT IS SAFE
#   The original used two nested R loops over every pair: for 1,178 samples
#   that is about 693,000 iterations, which is why it was slow.
#
#   This version computes the same thing with matrix products, which is
#   orders of magnitude faster. That is only worth doing if it gives the same
#   answer, so the script PROVES it:
#
#     1. It recomputes a random subset with the original nested-loop code and
#        compares. This always runs.
#     2. If you point --verify-against at the original IBS_distance_matrix.rds,
#        it compares every value against it.
#
#   Both checks report the largest absolute difference. Anything above
#   1e-10 is treated as a failure and stops the script. Floating-point
#   arithmetic means the difference will not be exactly zero; it should be
#   around 1e-16.
#
# INPUT
#   The population VCF used for the population-genetic analyses,
#   all_1179_withDhpsDhfr.vcf.gz. Set POPGEN_VCF in config.R.
#
# OUTPUTS
#   data/processed/ibs_distance_matrix.rds       the dist object
#   data/processed/ibs_distance_matrix.tsv.gz    the same, as a table
#   results/tables/figure04/ibs_distance_report.txt
#
# HOW TO RUN
#   From the top of the repository:
#
#     Rscript scripts/06_figures/figure04_ibs_distance.R
#
#   Optionally, to check against the original:
#
#     Rscript scripts/06_figures/figure04_ibs_distance.R \
#         --verify-against /path/to/IBS_distance_matrix.rds
#
################################################################################


# ==============================================================================
# 1. SETUP
# ==============================================================================

options(error = NULL)

suppressPackageStartupMessages({
  library(vcfR)
  library(adegenet)
})


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
    "Could not find config.R. Run this from the top of the repository.",
    call. = FALSE
  )
}

PROJECT_ROOT <- find_project_root()

source(file.path(PROJECT_ROOT, "config.R"))


arguments <- commandArgs(trailingOnly = TRUE)

get_argument <- function(flag, default = NULL) {

  position <- match(flag, arguments)

  if (is.na(position)) {
    return(default)
  }

  if (position == length(arguments)) {
    stop("The ", flag, " flag needs a value after it.", call. = FALSE)
  }

  arguments[[position + 1]]
}

VCF_FILE <- get_argument(
  "--vcf",
  if (exists("POPGEN_VCF", inherits = TRUE)) {
    get("POPGEN_VCF", inherits = TRUE)
  } else {
    NULL
  }
)

VERIFY_AGAINST <- get_argument("--verify-against", NULL)

# How many samples to recompute with the original loop, as a check.
VERIFY_SUBSET_SIZE <- as.integer(get_argument("--verify-n", "60"))

TOLERANCE <- 1e-10

if (is.null(VCF_FILE)) {
  stop(
    "No VCF given.\n\n",
    "Set POPGEN_VCF in config.R, or pass:\n",
    "  Rscript scripts/06_figures/figure04_ibs_distance.R \\\n",
    "      --vcf /path/to/all_1179_withDhpsDhfr.vcf.gz\n",
    call. = FALSE
  )
}

if (!file.exists(VCF_FILE)) {
  stop("VCF not found:\n  ", VCF_FILE, call. = FALSE)
}

PROCESSED_DIR <- file.path(PROJECT_ROOT, "data", "processed")
TABLE_DIR     <- file.path(PROJECT_ROOT, "results", "tables", "figure04")

dir.create(PROCESSED_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR,     recursive = TRUE, showWarnings = FALSE)

OUTPUT_RDS <- file.path(PROCESSED_DIR, "ibs_distance_matrix.rds")
OUTPUT_TSV <- file.path(PROCESSED_DIR, "ibs_distance_matrix.tsv.gz")
REPORT     <- file.path(TABLE_DIR, "ibs_distance_report.txt")


report_lines <- character(0)

add_to_report <- function(...) {
  line <- paste0(...)
  report_lines <<- c(report_lines, line)
  message(line)
}


# ==============================================================================
# 2. THE GENOTYPE MATRIX
# ==============================================================================
#
# Reproduced from the original analysis exactly:
#
#   vcfR2genlight()                 VCF to a genlight object
#   drop the sample named "3D7"     the laboratory reference control
#   tab(..., NA.method = "asis")    allele counts, missing data kept as NA
#
# NA.method = "asis" matters. The alternative, "mean", would fill missing
# genotypes with the site mean, which would make every pair comparable and
# change the distances. Keeping NA is what makes the measure
# pairwise-complete.
# ==============================================================================

add_to_report(strrep("=", 72))
add_to_report("IBS DISTANCE MATRIX")
add_to_report(strrep("=", 72))
add_to_report("VCF:       ", VCF_FILE)
add_to_report("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
add_to_report("")

vcf <- vcfR::read.vcfR(VCF_FILE, verbose = FALSE)

add_to_report(sprintf(
  "VCF: %d variants, %d samples", nrow(vcf@fix), ncol(vcf@gt) - 1L
))

genlight <- vcfR::vcfR2genlight(vcf)

REFERENCE_SAMPLES <- c("3D7", "Pf3D7", "PF3D7")

reference_present <- intersect(REFERENCE_SAMPLES, adegenet::indNames(genlight))

if (length(reference_present) > 0) {

  genlight <- genlight[
    !adegenet::indNames(genlight) %in% reference_present,
  ]

  add_to_report(sprintf(
    "Dropped the reference control: %s",
    paste(reference_present, collapse = ", ")
  ))
}

genotypes <- adegenet::tab(genlight, NA.method = "asis")

n_samples <- nrow(genotypes)
n_sites   <- ncol(genotypes)

add_to_report(sprintf(
  "Genotype matrix: %d samples x %d sites", n_samples, n_sites
))

add_to_report(sprintf(
  "Missing genotypes: %.2f%%", 100 * mean(is.na(genotypes))
))

observed_values <- sort(unique(as.vector(genotypes)))
observed_values <- observed_values[!is.na(observed_values)]

add_to_report(sprintf(
  "Genotype codes present: %s", paste(observed_values, collapse = ", ")
))

if (n_samples < 3) {
  stop("Fewer than three samples; a distance matrix is not meaningful.",
       call. = FALSE)
}


# ==============================================================================
# 3. THE DISTANCE, COMPUTED WITH MATRIX PRODUCTS
# ==============================================================================
#
# For each genotype code v, build an indicator matrix that is TRUE where a
# sample has that code and was called. Then
#
#   tcrossprod(indicator_v)   counts, for every pair, the sites where BOTH
#                             samples carry code v
#
# Summing over the codes gives the number of matching sites for every pair.
# tcrossprod on the "was called" indicator gives the number of jointly called
# sites. The distance follows directly.
#
# This is the same arithmetic as the original nested loops, rearranged so the
# work happens inside a few matrix multiplications.
# ==============================================================================

compute_ibs_distance_fast <- function(genotype_matrix) {

  called <- !is.na(genotype_matrix)

  storage.mode(called) <- "double"

  # Sites called in both samples, for every pair.
  n_comparable <- tcrossprod(called)

  codes <- sort(unique(as.vector(genotype_matrix)))
  codes <- codes[!is.na(codes)]

  n_matching <- matrix(0, nrow(genotype_matrix), nrow(genotype_matrix))

  for (code in codes) {

    indicator <- !is.na(genotype_matrix) & genotype_matrix == code

    indicator[is.na(indicator)] <- FALSE

    storage.mode(indicator) <- "double"

    n_matching <- n_matching + tcrossprod(indicator)
  }

  # A pair with no jointly called site has no defined distance. The original
  # code left those as NA, so they are left as NA here.
  distances <- ifelse(
    n_comparable > 0,
    1 - (n_matching / n_comparable),
    NA_real_
  )

  dimnames(distances) <- list(
    rownames(genotype_matrix),
    rownames(genotype_matrix)
  )

  diag(distances) <- 0

  distances
}


# ------------------------------------------------------------------------------
# The original implementation, kept verbatim so it can be used as the
# reference for the check below. It is never used to produce the output.
# ------------------------------------------------------------------------------

compute_ibs_distance_original <- function(mat) {

  n <- nrow(mat)
  d <- matrix(NA, n, n)

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {

      valid <- !is.na(mat[i, ]) & !is.na(mat[j, ])

      if (sum(valid) > 0) {
        same <- sum(mat[i, valid] == mat[j, valid])
        d[i, j] <- 1 - (same / sum(valid))
        d[j, i] <- d[i, j]
      }
    }
  }

  rownames(d) <- rownames(mat)
  colnames(d) <- rownames(mat)

  d
}


add_to_report("")
add_to_report("Computing distances with matrix products...")

elapsed <- system.time({
  distance_matrix <- compute_ibs_distance_fast(genotypes)
})

add_to_report(sprintf("  took %.1f seconds", elapsed[["elapsed"]]))

undefined_pairs <- sum(is.na(distance_matrix)) / 2

if (undefined_pairs > 0) {
  add_to_report(sprintf(
    "  %d pair(s) share no called site and have no defined distance (NA)",
    undefined_pairs
  ))
}

add_to_report(sprintf(
  "  distance range: %.6f to %.6f",
  min(distance_matrix, na.rm = TRUE),
  max(distance_matrix, na.rm = TRUE)
))


# ==============================================================================
# 4. CHECK 1 — AGAINST THE ORIGINAL LOOP, ON A SUBSET
# ==============================================================================

add_to_report("")
add_to_report(strrep("-", 72))
add_to_report("CHECK 1: the original nested-loop code, on a random subset")
add_to_report(strrep("-", 72))

subset_size <- min(VERIFY_SUBSET_SIZE, n_samples)

set.seed(12345)

subset_rows <- sort(sample.int(n_samples, subset_size))

genotypes_subset <- genotypes[subset_rows, , drop = FALSE]

reference_subset <- compute_ibs_distance_original(genotypes_subset)
fast_subset      <- compute_ibs_distance_fast(genotypes_subset)

# The original leaves the diagonal as NA; the fast version sets it to 0.
# Compare the off-diagonal values only.
diag(reference_subset) <- 0

comparable <- !is.na(reference_subset) & !is.na(fast_subset)

largest_difference <- if (any(comparable)) {
  max(abs(reference_subset[comparable] - fast_subset[comparable]))
} else {
  NA_real_
}

add_to_report(sprintf("  samples compared:        %d", subset_size))
add_to_report(sprintf("  pairs compared:          %d", sum(comparable) / 2))
add_to_report(sprintf("  largest difference:      %.3e", largest_difference))

na_disagreements <- sum(
  is.na(reference_subset) != is.na(fast_subset)
)

add_to_report(sprintf("  NA pattern differences:  %d", na_disagreements))

if (is.na(largest_difference) || largest_difference > TOLERANCE) {
  writeLines(report_lines, REPORT)
  stop(
    "The fast computation does not match the original implementation ",
    "(largest difference ", format(largest_difference), ", tolerance ",
    format(TOLERANCE), "). The output has NOT been written. See:\n  ", REPORT,
    call. = FALSE
  )
}

add_to_report("  PASS — the two implementations agree to floating-point precision.")


# ==============================================================================
# 5. CHECK 2 — AGAINST THE ORIGINAL SAVED MATRIX
# ==============================================================================

add_to_report("")
add_to_report(strrep("-", 72))
add_to_report("CHECK 2: the matrix saved by the original analysis")
add_to_report(strrep("-", 72))

if (is.null(VERIFY_AGAINST)) {

  add_to_report("  Skipped. To run it, pass:")
  add_to_report("    --verify-against /path/to/IBS_distance_matrix.rds")
  add_to_report("")
  add_to_report("  Worth doing once: it confirms this script reproduces the")
  add_to_report("  exact matrix the published figure was built from.")

} else if (!file.exists(VERIFY_AGAINST)) {

  stop("--verify-against file not found:\n  ", VERIFY_AGAINST, call. = FALSE)

} else {

  original <- readRDS(VERIFY_AGAINST)

  original_matrix <- as.matrix(original)

  add_to_report(sprintf("  file:                    %s", VERIFY_AGAINST))
  add_to_report(sprintf(
    "  original dimensions:     %d x %d",
    nrow(original_matrix), ncol(original_matrix)
  ))
  add_to_report(sprintf(
    "  recomputed dimensions:   %d x %d",
    nrow(distance_matrix), ncol(distance_matrix)
  ))

  shared_samples <- intersect(
    rownames(original_matrix), rownames(distance_matrix)
  )

  only_original    <- setdiff(rownames(original_matrix), rownames(distance_matrix))
  only_recomputed  <- setdiff(rownames(distance_matrix), rownames(original_matrix))

  add_to_report(sprintf("  samples in both:         %d", length(shared_samples)))
  add_to_report(sprintf("  only in the original:    %d", length(only_original)))
  add_to_report(sprintf("  only in the recomputed:  %d", length(only_recomputed)))

  if (length(only_original) > 0) {
    add_to_report("    only in the original: ",
                  paste(utils::head(only_original, 10), collapse = ", "))
  }

  if (length(only_recomputed) > 0) {
    add_to_report("    only in the recomputed: ",
                  paste(utils::head(only_recomputed, 10), collapse = ", "))
  }

  if (length(shared_samples) < 3) {

    add_to_report("")
    add_to_report("  Too few shared samples to compare. Is this the same VCF?")

  } else {

    original_shared   <- original_matrix[shared_samples, shared_samples]
    recomputed_shared <- distance_matrix[shared_samples, shared_samples]

    diag(original_shared) <- 0

    comparable_full <- !is.na(original_shared) & !is.na(recomputed_shared)

    difference_full <- max(
      abs(original_shared[comparable_full] - recomputed_shared[comparable_full])
    )

    add_to_report(sprintf(
      "  pairs compared:          %d", sum(comparable_full) / 2
    ))
    add_to_report(sprintf(
      "  largest difference:      %.3e", difference_full
    ))

    if (difference_full > TOLERANCE) {

      add_to_report("")
      add_to_report("  FAIL — the recomputed matrix differs from the original.")
      add_to_report("")
      add_to_report("  The most likely cause is a different input VCF. The")
      add_to_report("  published matrix came from all_1179_withDhpsDhfr.vcf.gz.")
      add_to_report("  Check POPGEN_VCF in config.R before going further: the")
      add_to_report("  PCoA and PERMANOVA are built on this matrix.")

      writeLines(report_lines, REPORT)

      stop(
        "The recomputed matrix does not match the original (largest ",
        "difference ", format(difference_full), "). The output has NOT been ",
        "written. See:\n  ", REPORT,
        call. = FALSE
      )
    }

    add_to_report("")
    add_to_report("  PASS — this script reproduces the published matrix exactly.")
  }
}


# ==============================================================================
# 6. SAVE
# ==============================================================================

distance_object <- stats::as.dist(distance_matrix)

saveRDS(distance_object, OUTPUT_RDS)

connection <- gzfile(OUTPUT_TSV, "w")

utils::write.table(
  distance_matrix,
  file = connection,
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

close(connection)

add_to_report("")
add_to_report(strrep("=", 72))
add_to_report("WRITTEN")
add_to_report(strrep("=", 72))
add_to_report("  ", OUTPUT_RDS)
add_to_report("  ", OUTPUT_TSV)
add_to_report("")
add_to_report("Next step:")
add_to_report("  Rscript scripts/06_figures/figure04_ibs_pcoa.R")
add_to_report(strrep("=", 72))

writeLines(report_lines, REPORT)

capture.output(
  sessionInfo(),
  file = file.path(TABLE_DIR, "ibs_distance_sessionInfo.txt")
)

message("\nReport written to:\n  ", REPORT)
