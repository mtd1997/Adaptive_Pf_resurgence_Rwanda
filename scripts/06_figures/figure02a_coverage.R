################################################################################
# FIGURE 2A — SEQUENCING COVERAGE ACROSS AMPLICON TARGETS
################################################################################
#
# WHAT IT DOES
#   Draws the distribution of mean sequencing depth for each of the seven
#   amplicons: a boxplot per gene with the individual samples shown as points.
#
# INPUT
#   A long-format coverage table with three columns:
#
#     sample    sample name, with a ".targets" suffix
#     gene      amplicon name
#     coverage  mean depth, with a COMMA as the decimal separator
#
#   This is the output of the per-amplicon coverage step in stage 02.
#
# OUTPUTS  (under results/)
#   figures/figure02/figure02a_coverage.svg / .pdf / .png
#   tables/figure02/figure02a_source_data.tsv       every point plotted
#   tables/figure02/figure02a_summary_by_gene.tsv   medians, ranges, counts
#   tables/figure02/figure02a_excluded_samples.tsv  what was left out, and why
#
# -----------------------------------------------------------------------------
# THREE THINGS THIS SCRIPT FIXES, AND WHY THEY MATTERED
# -----------------------------------------------------------------------------
#
# 1. IT REPORTS EXACTLY WHICH SAMPLES THE FIGURE COVERS
#
#    The coverage table holds every sample that was sequenced and aligned,
#    which is more than reach the published analyses: a few yield no variants
#    after VCF filtering, and a few more are excluded at the metadata check.
#
#    Both sets are legitimate. Coverage is a property of the sequencing,
#    measured before any downstream filtering, so a sequencing-quality figure
#    should show every sequenced sample — and that is the default here.
#
#    What the original script did not do was say how many samples it drew, or
#    how that number related to the analysed cohort. This one reconciles the
#    two and reports both, so the figure's sample count is a known number
#    rather than an assumed one, and the legend can state it correctly.
#
#    Pass --analysis-samples-only to restrict the figure to the analysed set.
#
#    Note what this script deliberately does NOT do: guess from sample names.
#    An earlier version dropped anything that looked like a control or a
#    placeholder, and that guess removed dozens of real samples whose names
#    were merely plain. The metadata is the only sound definition of which
#    samples were analysed.
#
# 2. IT SETS A SEED FOR THE POINT JITTER
#
#    geom_quasirandom() spreads overlapping points using random numbers. With
#    no seed, every run produces a visually different figure. The seed below
#    makes the figure reproducible.
#
# 3. IT REPORTS ZERO-COVERAGE POINTS INSTEAD OF LOSING THEM
#
#    The y axis is log scaled, and log10(0) is undefined, so any amplicon
#    with exactly 0x coverage is silently dropped by ggplot with a warning
#    that is easy to overlook. An amplicon that produced no reads is a real
#    result. This script counts those points, reports them, and lists them in
#    the excluded-samples table so they are accounted for rather than absent.
#
# -----------------------------------------------------------------------------
#
# HOW TO RUN
#   From the top of the repository:
#
#     Rscript scripts/06_figures/figure02a_coverage.R \
#         --coverage /path/to/coverage_all_samples.tsv
#
#   Options:
#     --metadata <file>          sample metadata; defaults to the repository
#                                copy. Used for the reconciliation report.
#     --analysis-samples-only    restrict the figure to the analysed samples
#                                instead of every sequenced one.
#
################################################################################


# ==============================================================================
# 1. SETUP
# ==============================================================================

options(error = NULL)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(ggbeeswarm)
  library(scales)
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


# ------------------------------------------------------------------------------
# Arguments
# ------------------------------------------------------------------------------

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

COVERAGE_FILE <- get_argument("--coverage")

METADATA_FILE <- get_argument(
  "--metadata",
  file.path(PROJECT_ROOT, "data", "metadata", "samples_metadata.tsv")
)

# By default every sample in the coverage table is plotted. Coverage is a
# property of the sequencing, measured before any downstream filtering, so a
# sequencing-quality figure should show everything that was sequenced.
#
# --analysis-samples-only restricts the figure to the samples that survived
# into the analysis, for the cases where that is what is wanted.
RESTRICT_TO_ANALYSIS <- "--analysis-samples-only" %in% arguments

if (is.null(COVERAGE_FILE)) {
  stop(
    "No coverage file given.\n\n",
    "Usage:\n",
    "  Rscript scripts/06_figures/figure02a_coverage.R \\\n",
    "      --coverage /path/to/coverage_all_samples.tsv\n",
    call. = FALSE
  )
}

if (!file.exists(COVERAGE_FILE)) {
  stop("Coverage file not found: ", COVERAGE_FILE, call. = FALSE)
}

FIGURE_DIR <- file.path(PROJECT_ROOT, "results", "figures", "figure02")
TABLE_DIR  <- file.path(PROJECT_ROOT, "results", "tables",  "figure02")

dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR,  recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# 2. GENE ORDER AND COLOURS
# ==============================================================================
#
# This order is the one used in the published figure. It is not the order used
# in the selection-scan figure, which follows chromosome position — the two
# figures deliberately differ, and changing either to match the other would
# alter a published figure.
# ==============================================================================

GENE_ORDER <- c("CRT", "DHFR", "DHPS", "MDR1", "K13", "CSP", "MSP2")

# There is deliberately no list of sample names to exclude here.
#
# An earlier version of this script guessed at which names looked like
# controls or placeholders and dropped them. That guess was wrong: names such
# as "sample12" are real study samples that happen to be plainly named, and
# removing them would have quietly deleted dozens of genuine measurements
# from the figure.
#
# The sample metadata is the only defensible definition of which samples were
# analysed, so it does the filtering when filtering is wanted. A name pattern
# is not evidence of anything.


# ==============================================================================
# 3. READ THE COVERAGE TABLE
# ==============================================================================
#
# Read as text. The coverage values use a comma as the decimal separator, so
# letting readr parse them as numbers would either fail or, worse, read
# "1250,18" as something unintended.
# ==============================================================================

raw_coverage <- readr::read_tsv(
  COVERAGE_FILE,
  col_types = readr::cols(.default = readr::col_character()),
  show_col_types = FALSE
)

required_columns <- c("sample", "gene", "coverage")

missing_columns <- setdiff(required_columns, names(raw_coverage))

if (length(missing_columns) > 0) {
  stop(
    "The coverage file is missing: ",
    paste(missing_columns, collapse = ", "),
    "\nExpected three tab-separated columns: sample, gene, coverage.",
    call. = FALSE
  )
}

message(sprintf("Read %d rows from:\n  %s\n", nrow(raw_coverage), COVERAGE_FILE))


# ==============================================================================
# 4. CLEAN
# ==============================================================================

report_lines <- character(0)

add_to_report <- function(...) {
  line <- paste0(...)
  report_lines <<- c(report_lines, line)
  message(line)
}

add_to_report(strrep("=", 72))
add_to_report("FIGURE 2A — COVERAGE")
add_to_report(strrep("=", 72))
add_to_report("Coverage file: ", COVERAGE_FILE)
add_to_report("Generated:     ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
add_to_report("")
add_to_report("Rows read: ", nrow(raw_coverage))


# ------------------------------------------------------------------------------
# Drop the placeholder rows.
#
# The coverage table contains one row per sample with gene = "0" and no
# coverage value. These are an artefact of whatever produced the file — most
# likely a trailing blank line in the BED input — and are not a real amplicon.
# ------------------------------------------------------------------------------

artefact_rows <- raw_coverage %>%
  dplyr::filter(gene == "0" | is.na(coverage) | stringr::str_trim(coverage) == "")

coverage_data <- raw_coverage %>%
  dplyr::filter(
    gene != "0",
    !is.na(coverage),
    stringr::str_trim(coverage) != ""
  )

add_to_report(sprintf(
  "  dropped %d placeholder row(s) with gene = \"0\" or no coverage value",
  nrow(artefact_rows)
))


# ------------------------------------------------------------------------------
# Parse the coverage values and tidy the sample names.
# ------------------------------------------------------------------------------

coverage_data <- coverage_data %>%
  dplyr::mutate(
    sample = stringr::str_remove(stringr::str_trim(sample), "\\.targets$"),
    gene   = stringr::str_trim(gene),

    # Comma decimal separator -> period, then to numeric.
    coverage_text = stringr::str_replace(stringr::str_trim(coverage), ",", "."),
    coverage      = suppressWarnings(as.numeric(coverage_text))
  )

unparseable <- dplyr::filter(coverage_data, is.na(coverage))

if (nrow(unparseable) > 0) {

  add_to_report("")
  add_to_report(sprintf(
    "  WARNING: %d coverage value(s) could not be read as a number:",
    nrow(unparseable)
  ))

  for (value in utils::head(unique(unparseable$coverage_text), 5)) {
    add_to_report(sprintf("    \"%s\"", value))
  }

  coverage_data <- dplyr::filter(coverage_data, !is.na(coverage))
}

coverage_data <- dplyr::select(coverage_data, -coverage_text)


# ------------------------------------------------------------------------------
# Check the genes are the expected seven.
# ------------------------------------------------------------------------------

genes_found <- sort(unique(coverage_data$gene))

genes_unexpected <- setdiff(genes_found, GENE_ORDER)
genes_absent     <- setdiff(GENE_ORDER, genes_found)

if (length(genes_unexpected) > 0) {
  add_to_report("")
  add_to_report(
    "  WARNING: gene name(s) not in GENE_ORDER, which will not be drawn: ",
    paste(genes_unexpected, collapse = ", ")
  )
}

if (length(genes_absent) > 0) {
  add_to_report("")
  add_to_report(
    "  WARNING: expected gene(s) absent from the coverage file: ",
    paste(genes_absent, collapse = ", ")
  )
}

add_to_report("")
add_to_report("Samples in the coverage file: ", dplyr::n_distinct(coverage_data$sample))

# ==============================================================================
# 5. RECONCILE AGAINST THE SAMPLE METADATA
# ==============================================================================
#
# The coverage table holds every sample that was sequenced and aligned. Fewer
# samples than that reach the published analyses, because some are dropped
# downstream: a few yield no variants after VCF filtering, and a few more are
# excluded once the metadata is checked.
#
# Both sets are legitimate, and they answer different questions:
#
#   every sequenced sample     how well did the sequencing work
#   the analysed samples only  what did the analysed cohort look like
#
# Figure 2a is a sequencing-quality figure, so the default is every sequenced
# sample. The reconciliation is reported either way, so the figure's sample
# count is always a known number rather than an assumed one.
# ==============================================================================

all_samples <- sort(unique(coverage_data$sample))

analysis_samples <- NULL

if (file.exists(METADATA_FILE)) {

  metadata <- readr::read_tsv(
    METADATA_FILE,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE
  )

  if (!"SAMPLE" %in% names(metadata)) {
    stop(
      "The metadata has no SAMPLE column: ", METADATA_FILE,
      call. = FALSE
    )
  }

  analysis_samples <- unique(stringr::str_trim(metadata$SAMPLE))

  in_coverage_only <- setdiff(all_samples, analysis_samples)
  in_metadata_only <- setdiff(analysis_samples, all_samples)

  add_to_report("")
  add_to_report("Metadata: ", METADATA_FILE)
  add_to_report("")
  add_to_report("Reconciliation:")
  add_to_report(sprintf(
    "  sequenced (in the coverage table):  %d", length(all_samples)
  ))
  add_to_report(sprintf(
    "  analysed  (in the metadata):        %d", length(analysis_samples)
  ))
  add_to_report(sprintf(
    "  in both:                            %d",
    length(intersect(all_samples, analysis_samples))
  ))
  add_to_report(sprintf(
    "  sequenced but not analysed:         %d", length(in_coverage_only)
  ))
  add_to_report(sprintf(
    "  analysed but no coverage measured:  %d", length(in_metadata_only)
  ))

  if (length(in_coverage_only) > 0) {

    add_to_report("")
    add_to_report(
      "Sequenced but not analysed. These are expected: samples that yielded "
    )
    add_to_report(
      "no variants after filtering, laboratory controls, and samples excluded "
    )
    add_to_report("at the metadata check. Listed in full in the exclusions table.")
    add_to_report("")

    for (name in utils::head(sort(in_coverage_only), 10)) {
      add_to_report("  ", name)
    }

    if (length(in_coverage_only) > 10) {
      add_to_report(sprintf("  ... and %d more", length(in_coverage_only) - 10))
    }
  }

  # This one is not expected and is worth stopping to look at: a sample that
  # the analysis used but for which no coverage was ever measured means the
  # coverage table is out of date with respect to the alignments.
  if (length(in_metadata_only) > 0) {

    add_to_report("")
    add_to_report(strrep("!", 72))
    add_to_report(sprintf(
      "%d analysed sample(s) have no coverage measurement.",
      length(in_metadata_only)
    ))
    add_to_report("")
    add_to_report("The figure therefore does not cover every analysed sample.")
    add_to_report("Re-measure coverage from the current alignments:")
    add_to_report("  bash scripts/02_alignment/compute_amplicon_coverage.sh")
    add_to_report("")

    for (name in utils::head(sort(in_metadata_only), 10)) {
      add_to_report("  ", name)
    }

    if (length(in_metadata_only) > 10) {
      add_to_report(sprintf("  ... and %d more", length(in_metadata_only) - 10))
    }

    add_to_report(strrep("!", 72))
  }

} else {

  add_to_report("")
  add_to_report("No sample metadata found at:")
  add_to_report("  ", METADATA_FILE)
  add_to_report("")
  add_to_report("The figure will show every sample in the coverage table, and")
  add_to_report("no reconciliation against the analysed set can be reported.")

  if (RESTRICT_TO_ANALYSIS) {
    stop(
      "--analysis-samples-only was given but the metadata was not found:\n  ",
      METADATA_FILE,
      "\n\nBuild it first:\n",
      "  Rscript scripts/00_metadata/prepare_metadata.R --input <your file>",
      call. = FALSE
    )
  }
}


# ------------------------------------------------------------------------------
# Decide what to plot.
# ------------------------------------------------------------------------------

exclusions <- tibble::tibble(
  sample = character(0),
  reason = character(0)
)

if (RESTRICT_TO_ANALYSIS && !is.null(analysis_samples)) {

  excluded_names <- setdiff(all_samples, analysis_samples)

  exclusions <- tibble::tibble(
    sample = excluded_names,
    reason = "not in the analysis metadata"
  )

  plot_data <- dplyr::filter(coverage_data, sample %in% analysis_samples)

  add_to_report("")
  add_to_report(
    "--analysis-samples-only given: restricted to the analysed samples."
  )

} else {

  plot_data <- coverage_data

  add_to_report("")
  add_to_report(
    "Plotting every sequenced sample. Pass --analysis-samples-only to "
  )
  add_to_report("restrict the figure to the analysed set instead.")
}

plot_data <- plot_data %>%
  dplyr::filter(gene %in% GENE_ORDER) %>%
  dplyr::mutate(gene = factor(gene, levels = GENE_ORDER))

add_to_report("")
add_to_report("Samples plotted:  ", dplyr::n_distinct(plot_data$sample))
add_to_report("Points plotted:   ", nrow(plot_data))

if (nrow(plot_data) == 0) {
  stop(
    "No data left to plot. Check that the sample names in the coverage file ",
    "and the metadata use the same convention.",
    call. = FALSE
  )
}


# ==============================================================================
# 6. ZERO-COVERAGE POINTS
# ==============================================================================
#
# The y axis is log scaled, so points at exactly 0x cannot be drawn. Rather
# than let ggplot drop them with a warning, they are identified here, counted,
# and recorded in the excluded-samples table.
# ==============================================================================

zero_coverage <- dplyr::filter(plot_data, coverage <= 0)

if (nrow(zero_coverage) > 0) {

  add_to_report("")
  add_to_report(strrep("-", 72))
  add_to_report(sprintf(
    "%d amplicon(s) across %d sample(s) have 0x coverage.",
    nrow(zero_coverage),
    dplyr::n_distinct(zero_coverage$sample)
  ))
  add_to_report("")
  add_to_report("These cannot appear on a logarithmic axis and are therefore")
  add_to_report("absent from the figure. They are real results — an amplicon")
  add_to_report("that produced no reads — and are listed in")
  add_to_report("figure02a_excluded_samples.tsv so they are accounted for.")
  add_to_report("")

  zero_by_gene <- dplyr::count(zero_coverage, gene, name = "n")

  for (index in seq_len(nrow(zero_by_gene))) {
    add_to_report(sprintf(
      "  %-8s %d sample(s) at 0x",
      as.character(zero_by_gene$gene[[index]]),
      zero_by_gene$n[[index]]
    ))
  }

  add_to_report(strrep("-", 72))

  exclusions <- dplyr::bind_rows(
    exclusions,
    zero_coverage %>%
      dplyr::transmute(
        sample,
        reason = paste0("0x coverage at ", gene, "; cannot be shown on a log axis")
      )
  )
}

# Points at 0x are removed explicitly, so the count of what is drawn is known
# rather than being whatever survived the scale transformation.
plot_data_drawn <- dplyr::filter(plot_data, coverage > 0)


# ==============================================================================
# 7. SUMMARY BY GENE
# ==============================================================================

summary_by_gene <- plot_data %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    N_samples        = dplyr::n(),
    median_coverage  = stats::median(coverage, na.rm = TRUE),
    mean_coverage    = mean(coverage, na.rm = TRUE),
    lower_quartile   = stats::quantile(coverage, 0.25, na.rm = TRUE, names = FALSE),
    upper_quartile   = stats::quantile(coverage, 0.75, na.rm = TRUE, names = FALSE),
    minimum_coverage = min(coverage, na.rm = TRUE),
    maximum_coverage = max(coverage, na.rm = TRUE),
    N_below_50x      = sum(coverage < 50, na.rm = TRUE),
    N_below_100x     = sum(coverage < 100, na.rm = TRUE),
    N_at_zero        = sum(coverage <= 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(gene)

add_to_report("")
add_to_report(strrep("-", 72))
add_to_report("COVERAGE BY GENE")
add_to_report(strrep("-", 72))
add_to_report(sprintf(
  "%-8s %7s %10s %10s %10s %7s",
  "gene", "n", "median", "min", "max", "<50x"
))

for (index in seq_len(nrow(summary_by_gene))) {
  add_to_report(sprintf(
    "%-8s %7d %10.1f %10.2f %10.1f %7d",
    as.character(summary_by_gene$gene[[index]]),
    summary_by_gene$N_samples[[index]],
    summary_by_gene$median_coverage[[index]],
    summary_by_gene$minimum_coverage[[index]],
    summary_by_gene$maximum_coverage[[index]],
    summary_by_gene$N_below_50x[[index]]
  ))
}


# ==============================================================================
# 8. DRAW
# ==============================================================================
#
# The styling is carried over from the published figure unchanged. The one
# addition is the seed on geom_quasirandom, which makes the point positions
# reproducible.
# ==============================================================================

JITTER_SEED <- 12345

set.seed(JITTER_SEED)

figure02a <- ggplot2::ggplot(
  plot_data_drawn,
  ggplot2::aes(x = gene, y = coverage)
) +

  ggplot2::geom_boxplot(
    outlier.shape = NA,
    fill      = "grey90",
    colour    = "black",
    linewidth = 0.8
  ) +

  ggbeeswarm::geom_quasirandom(
    ggplot2::aes(colour = gene),
    width       = 0.25,
    alpha       = 0.6,
    size        = 1.6,
    show.legend = FALSE
  ) +

  ggplot2::scale_y_log10(
    labels = scales::comma,
    breaks = c(100, 1e3, 1e4, 1e5)
  ) +

  ggplot2::labs(
    title = "Distribution of sequencing coverage depth across amplicon targets",
    x     = NULL,
    y     = "Coverage depth (log10)"
  ) +

  ggplot2::theme_bw(base_size = 14) +

  ggplot2::theme(
    axis.title  = ggplot2::element_text(
      size = 14, colour = "black", face = "bold"
    ),
    axis.text   = ggplot2::element_text(
      size = 14, colour = "black", face = "bold"
    ),
    axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
    plot.title  = ggplot2::element_text(
      face = "bold", size = 20, colour = "black", hjust = 0.5
    ),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    panel.border     = ggplot2::element_blank(),
    axis.line        = ggplot2::element_line(colour = "black")
  )


# ------------------------------------------------------------------------------
# Save. Dimensions are those used for the published figure.
# ------------------------------------------------------------------------------

FIGURE_WIDTH  <- 11
FIGURE_HEIGHT <- 7

for (extension in c("svg", "pdf", "png")) {

  ggplot2::ggsave(
    filename = file.path(
      FIGURE_DIR, paste0("figure02a_coverage.", extension)
    ),
    plot   = figure02a,
    width  = FIGURE_WIDTH,
    height = FIGURE_HEIGHT,
    units  = "in",
    dpi    = 600
  )
}


# ==============================================================================
# 9. SOURCE DATA AND REPORT
# ==============================================================================

readr::write_tsv(
  plot_data %>%
    dplyr::select(sample, gene, coverage) %>%
    dplyr::arrange(gene, sample),
  file.path(TABLE_DIR, "figure02a_source_data.tsv")
)

readr::write_tsv(
  summary_by_gene,
  file.path(TABLE_DIR, "figure02a_summary_by_gene.tsv")
)

readr::write_tsv(
  dplyr::arrange(exclusions, reason, sample),
  file.path(TABLE_DIR, "figure02a_excluded_samples.tsv")
)

add_to_report("")
add_to_report(strrep("=", 72))
add_to_report("WRITTEN")
add_to_report(strrep("=", 72))
add_to_report("  figures/figure02/figure02a_coverage.svg / .pdf / .png")
add_to_report("  tables/figure02/figure02a_source_data.tsv")
add_to_report("  tables/figure02/figure02a_summary_by_gene.tsv")
add_to_report("  tables/figure02/figure02a_excluded_samples.tsv")
add_to_report("")
add_to_report("Jitter seed: ", JITTER_SEED, " (change it and the points move)")

writeLines(
  report_lines,
  file.path(TABLE_DIR, "figure02a_report.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(TABLE_DIR, "figure02a_sessionInfo.txt")
)

message("\nReport also written to:\n  ",
        file.path(TABLE_DIR, "figure02a_report.txt"))
