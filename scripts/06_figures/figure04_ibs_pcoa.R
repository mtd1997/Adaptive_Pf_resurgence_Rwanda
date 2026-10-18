################################################################################
# FIGURE 4 AND SUPPLEMENTARY FIGURE S4 — PCoA AND PERMANOVA
################################################################################
#
# WHAT IT DOES
#   Takes the IBS distance matrix and produces three analyses:
#
#     Kigali     Baseline-Kigali versus Resurgence   (temporal)
#     Resurgence Kigali versus Hotspot               (spatial)
#     Resurgence coloured by kelch13 genotype
#
#   For each it runs a principal coordinates analysis on four axes and a
#   PERMANOVA, then writes:
#
#     one single panel  -> Figure 4a, 4b, 4c
#     one 2x2 composite -> Supplementary Figure S4a, S4b, S4c
#
#   Figure 4d is a distance tree and is not produced here.
#
# INPUTS
#   data/processed/ibs_distance_matrix.rds   from figure04_ibs_distance.R
#   data/metadata/samples_metadata.tsv
#
# OUTPUTS
#   results/figures/figure04/   the panels, named by published panel letter
#   results/tables/figure04/    PERMANOVA results, PCoA coordinates,
#                               eigenvalues, sample counts
#
# -----------------------------------------------------------------------------
# TWO THINGS TO READ BEFORE TRUSTING THE LABELS
# -----------------------------------------------------------------------------
#
# 1. ON THE kelch13 PANEL, THE PERMANOVA IS NOT ON kelch13
#
#    Panel 4c and S4c colour the samples by kelch13 genotype, but the
#    PERMANOVA printed on them tests REGION, because that is what the
#    original script passed to adonis2. The published figure shows R2 = 0.037
#    on that panel, which is identical to the region panel's value -- visible
#    proof that it is the region test.
#
#    This script reproduces the published figure faithfully, so the label on
#    the kelch13 panels still shows the region test. It does NOT silently
#    change a published figure.
#
#    It does, however, additionally compute the PERMANOVA on kelch13 genotype
#    and write it to the results table, so the number exists if the legend is
#    corrected. See permanova_results.tsv, where every test is listed with
#    what it was actually run on.
#
# 2. "p < 0.001" IS NOT A CORRECT WAY TO REPORT THIS
#
#    adonis2 with 999 permutations estimates p as
#        (number of permutations at least as extreme + 1) / (999 + 1)
#    so the smallest value it can return is 1/1000 = 0.001. The result cannot
#    be less than that, and writing "p < 0.001" claims more than the test can
#    support.
#
#    The published single panels already say "p = 0.001", which is right. One
#    of the original composite functions said "p < ", which is not. This
#    script writes "p = " everywhere and states the permutation count in the
#    report.
# -----------------------------------------------------------------------------
#
# HOW TO RUN
#   From the top of the repository:
#
#     Rscript scripts/06_figures/figure04_ibs_pcoa.R
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
  library(patchwork)
  library(vegan)
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

DISTANCE_FILE <- file.path(
  PROJECT_ROOT, "data", "processed", "ibs_distance_matrix.rds"
)

# ------------------------------------------------------------------------------
# WHICH METADATA, AND WHY IT MATTERS FOR THIS FIGURE
#
# The default is the repository's own metadata. For Figure 4 specifically that
# gives one sample fewer than the published figure, because the published
# figure was built from an earlier roster that still contained RW_SN94 -- a
# Kigali baseline sample later excluded from the 1,176-sample analysis set.
#
# The effect is confined to the Kigali panel:
#
#   with RW_SN94     N = 805, R2 = 0.0535  -> prints as 0.054 (published)
#   without RW_SN94  N = 804, R2 = 0.0527  -> prints as 0.053
#
# Pass --metadata with the older file to reproduce the published panel
# exactly. Use the default to match the sample set the paper states. The
# report records which was used either way, so the figure's N is never
# ambiguous.
# ------------------------------------------------------------------------------

arguments_early <- commandArgs(trailingOnly = TRUE)

metadata_flag_position <- match("--metadata", arguments_early)

FIGURE04_METADATA <- file.path(
  PROJECT_ROOT, "data", "metadata", "samples_metadata_figure04.tsv"
)

DEFAULT_METADATA <- file.path(
  PROJECT_ROOT, "data", "metadata", "samples_metadata.tsv"
)

METADATA_FILE <- if (!is.na(metadata_flag_position) &&
                     metadata_flag_position < length(arguments_early)) {

  arguments_early[[metadata_flag_position + 1]]

} else if (file.exists(FIGURE04_METADATA)) {

  # The Figure 4 roster, when it has been prepared. This is the default so
  # that running the script with no arguments reproduces the published
  # figure rather than coming out one sample short.
  FIGURE04_METADATA

} else {

  DEFAULT_METADATA
}

FIGURE_DIR <- file.path(PROJECT_ROOT, "results", "figures", "figure04")
TABLE_DIR  <- file.path(PROJECT_ROOT, "results", "tables",  "figure04")

dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR,  recursive = TRUE, showWarnings = FALSE)

for (required in c(DISTANCE_FILE, METADATA_FILE)) {
  if (!file.exists(required)) {
    stop(
      "Missing input:\n  ", required,
      "\n\nRun these first:\n",
      "  Rscript scripts/06_figures/figure04_ibs_distance.R\n",
      "  Rscript scripts/00_metadata/prepare_metadata.R --input <your file>",
      call. = FALSE
    )
  }
}

PERMUTATIONS <- 999

# adonis2 and stat_ellipse both draw on the random number generator, so the
# seed is fixed to make the p-values and ellipses reproducible between runs.
SEED <- 12345


report_lines <- character(0)

add_to_report <- function(...) {
  line <- paste0(...)
  report_lines <<- c(report_lines, line)
  message(line)
}


# ==============================================================================
# 2. WHAT EACH PANEL IS
# ==============================================================================
#
# This table is the single place the published panel letters are recorded. If
# a letter is wrong, change it here and every output file is renamed to match.
#
#   analysis       internal name for the subset and comparison
#   subset         which samples
#   permanova_on   the metadata column adonis2 is run on
#   colour_by      the metadata column the points are coloured by
#   shape_by       the metadata column the point shapes come from, or NA
#   main_panel     the published Figure 4 panel letter
#   main_axes      which two PCoA axes that single panel shows
#   supp_panel     the published Supplementary Figure S4 panel letter
#
# The axis choices are not all the same: the Kigali panel shows axes 1 and 3,
# the region panel axes 1 and 2, and the kelch13 panel axes 2 and 3. Those are
# the pairs used in the published figure, read off the panel axis labels.
# ==============================================================================

PANEL_MAP <- tibble::tribble(
  ~analysis,            ~subset,       ~permanova_on, ~colour_by,    ~shape_by, ~main_panel, ~main_axis_x, ~main_axis_y, ~supp_panel,
  "kigali_period",      "kigali",      "period",      "period",      NA,        "4a",        1L,           3L,           "S4a",
  "resurgence_region",  "resurgence",  "region",      "region",      NA,        "4b",        1L,           2L,           "S4b",
  "resurgence_k13",     "resurgence",  "region",      "k13_status",  "region",  "4c",        2L,           3L,           "S4c"
)

# The 2x2 composite always shows these four axis pairs, in this arrangement:
#
#     1 vs 2   |   1 vs 3
#     1 vs 4   |   2 vs 3
#
COMPOSITE_AXES <- list(
  c(1, 2), c(1, 3), c(1, 4), c(2, 3)
)

N_AXES <- 4L

PERIOD_COLOURS <- c(
  "Baseline-Kigali" = "#1f77b4",
  "Resurgence"      = "#d62728"
)

REGION_COLOURS <- c(
  "Kigali"  = "#6B7900",
  "Hotspot" = "#D16100"
)

K13_COLOURS <- c(
  "WT"           = "#2ecc71",
  "K189T_only"   = "#d62f7d",
  "R561H+K189T"  = "#ff7f0e",
  "R561H_only"   = "#f4b400",
  "Other_K13"    = "#6a5acd"
)

PERIOD_SHAPES <- c("Baseline-Kigali" = 16, "Resurgence" = 17)
REGION_SHAPES <- c("Kigali" = 16, "Hotspot" = 17)

# On the kelch13 panels an ellipse is drawn around this genotype only, to
# show that the double mutant occupies a distinct part of the space.
ELLIPSE_GENOTYPE <- "R561H+K189T"

palette_for <- function(column) {
  switch(
    column,
    "period"     = PERIOD_COLOURS,
    "region"     = REGION_COLOURS,
    "k13_status" = K13_COLOURS,
    stop("No palette defined for '", column, "'.", call. = FALSE)
  )
}

shapes_for <- function(column) {
  switch(
    column,
    "period" = PERIOD_SHAPES,
    "region" = REGION_SHAPES,
    stop("No shapes defined for '", column, "'.", call. = FALSE)
  )
}

legend_title_for <- function(column) {
  switch(
    column,
    "period"     = "Period",
    "region"     = "Region",
    "k13_status" = "K13 genotype",
    column
  )
}


# ==============================================================================
# 3. READ AND ALIGN
# ==============================================================================

add_to_report(strrep("=", 72))
add_to_report("FIGURE 4 AND SUPPLEMENTARY FIGURE S4 — PCoA AND PERMANOVA")
add_to_report(strrep("=", 72))
add_to_report("Distance: ", DISTANCE_FILE)
add_to_report("Metadata: ", METADATA_FILE)
add_to_report("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
add_to_report("")

distance_full <- readRDS(DISTANCE_FILE)

distance_samples <- attr(distance_full, "Labels")

if (is.null(distance_samples)) {
  stop(
    "The distance object has no sample labels, so it cannot be matched to ",
    "the metadata.",
    call. = FALSE
  )
}

metadata <- readr::read_tsv(
  METADATA_FILE,
  col_types = readr::cols(.default = readr::col_character()),
  show_col_types = FALSE
)


# ------------------------------------------------------------------------------
# Accept either capitalisation.
# ------------------------------------------------------------------------------

resolve_column <- function(data, candidates, description) {

  normalised_available <- tolower(gsub("[^a-z0-9]", "", tolower(names(data))))
  normalised_wanted    <- tolower(gsub("[^a-z0-9]", "", tolower(candidates)))

  position <- match(normalised_wanted, normalised_available)
  position <- position[!is.na(position)]

  if (length(position) == 0) {
    stop(
      "No column found for ", description, ". Looked for: ",
      paste(candidates, collapse = ", "),
      "\nThe metadata has: ", paste(names(data), collapse = ", "),
      call. = FALSE
    )
  }

  names(data)[[position[[1]]]]
}

column_map <- c(
  SAMPLE     = resolve_column(metadata, c("SAMPLE", "sample"), "sample ID"),
  period     = resolve_column(metadata, c("period", "PERIOD"), "period"),
  region     = resolve_column(metadata, c("region", "REGION"), "region"),
  k13_status = resolve_column(
    metadata, c("k13_status", "K13_status", "K13_STATUS"), "kelch13 status"
  )
)

metadata <- metadata %>%
  dplyr::select(dplyr::all_of(unname(column_map)))

names(metadata) <- names(column_map)

metadata <- metadata %>%
  dplyr::mutate(
    dplyr::across(dplyr::everything(), ~ stringr::str_trim(.x)),

    # Period labels differ between metadata files. The published wording is
    # "Baseline-Kigali" and "Resurgence"; earlier files use "Pre-Resurgence"
    # and "Post-Resurgence" for the same two groups. recode() leaves anything
    # it does not match untouched, so a file already using the published
    # wording passes through unchanged.
    #
    # This matters more than it looks: without it, a file using the older
    # labels silently yields an empty Resurgence subset.
    period = dplyr::recode(
      period,
      "Pre-Resurgence"  = "Baseline-Kigali",
      "Post-Resurgence" = "Resurgence"
    ),

    # Samples with no recorded kelch13 mutation are wild type. The original
    # analysis made this explicit the same way.
    k13_status = tidyr::replace_na(k13_status, "WT")
  ) %>%
  dplyr::filter(!is.na(SAMPLE), SAMPLE != "") %>%
  dplyr::distinct(SAMPLE, .keep_all = TRUE)


# ------------------------------------------------------------------------------
# Keep only samples present in both, in the distance matrix's order.
# ------------------------------------------------------------------------------

shared_samples <- intersect(distance_samples, metadata$SAMPLE)

add_to_report("Sample reconciliation:")
add_to_report(sprintf("  in the distance matrix:  %d", length(distance_samples)))
add_to_report(sprintf("  in the metadata:         %d", nrow(metadata)))
add_to_report(sprintf("  in both:                 %d", length(shared_samples)))
add_to_report(sprintf(
  "  distance matrix only:    %d",
  length(setdiff(distance_samples, metadata$SAMPLE))
))
add_to_report(sprintf(
  "  metadata only:           %d",
  length(setdiff(metadata$SAMPLE, distance_samples))
))

if (length(shared_samples) < 10) {
  stop(
    "Only ", length(shared_samples), " sample(s) are shared between the ",
    "distance matrix and the metadata. Check that both use the same sample ",
    "naming convention.",
    call. = FALSE
  )
}

metadata <- metadata %>%
  dplyr::filter(SAMPLE %in% shared_samples) %>%
  dplyr::arrange(match(SAMPLE, distance_samples))

distance_matrix_full <- as.matrix(distance_full)[metadata$SAMPLE, metadata$SAMPLE]

add_to_report("")
add_to_report("Samples per period and region:")

counts_table <- table(metadata$period, metadata$region, useNA = "ifany")

for (period_name in rownames(counts_table)) {
  add_to_report(sprintf(
    "  %-18s %s",
    period_name,
    paste(
      sprintf("%s=%d", colnames(counts_table), counts_table[period_name, ]),
      collapse = "  "
    )
  ))
}

add_to_report("")
add_to_report("kelch13 genotypes:")

k13_counts <- sort(table(metadata$k13_status), decreasing = TRUE)

for (genotype in names(k13_counts)) {
  add_to_report(sprintf("  %-16s %d", genotype, k13_counts[[genotype]]))
}

unexpected_genotypes <- setdiff(names(k13_counts), names(K13_COLOURS))

if (length(unexpected_genotypes) > 0) {
  add_to_report("")
  add_to_report(
    "  WARNING: genotype(s) with no colour defined, which will be drawn grey: ",
    paste(unexpected_genotypes, collapse = ", ")
  )
}


# ==============================================================================
# 4. SUBSETS
# ==============================================================================
#
# Kigali     every Kigali sample, both periods -- the temporal comparison
# Resurgence every Resurgence sample, both regions -- the spatial comparison
# ==============================================================================

build_subset <- function(name) {

  rows <- switch(
    name,
    "kigali"     = metadata$region == "Kigali",
    "resurgence" = metadata$period == "Resurgence",
    stop("Unknown subset '", name, "'.", call. = FALSE)
  )

  rows[is.na(rows)] <- FALSE

  subset_metadata <- metadata[rows, , drop = FALSE]

  if (nrow(subset_metadata) < 3) {
    stop(
      "The '", name, "' subset has only ", nrow(subset_metadata),
      " sample(s).\n\n",
      "The values actually present are:\n",
      "  period: ", paste(sort(unique(metadata$period)), collapse = ", "), "\n",
      "  region: ", paste(sort(unique(metadata$region)), collapse = ", "), "\n\n",
      "The subsets are built by matching period == \"Resurgence\" and ",
      "region == \"Kigali\". If the labels above use different wording, add ",
      "it to the recode() in section 3.",
      call. = FALSE
    )
  }

  subset_distance <- stats::as.dist(
    distance_matrix_full[subset_metadata$SAMPLE, subset_metadata$SAMPLE]
  )

  list(metadata = subset_metadata, distance = subset_distance)
}

subsets <- list(
  kigali     = build_subset("kigali"),
  resurgence = build_subset("resurgence")
)

add_to_report("")
add_to_report("Subsets:")
for (name in names(subsets)) {
  add_to_report(sprintf(
    "  %-12s %d samples", name, nrow(subsets[[name]]$metadata)
  ))
}


# ==============================================================================
# 5. PCoA AND PERMANOVA
# ==============================================================================
#
# cmdscale is classical multidimensional scaling, which is principal
# coordinates analysis. No correction is applied for negative eigenvalues --
# that matches the original analysis, and the percentages below are computed
# over the positive eigenvalues only.
# ==============================================================================

run_pcoa <- function(distance_object, n_axes = N_AXES) {

  pcoa <- stats::cmdscale(distance_object, k = n_axes, eig = TRUE)

  # ---------------------------------------------------------------------------
  # THE DENOMINATOR HERE IS A DELIBERATE CHOICE, AND IT IS NOT THE OBVIOUS ONE
  #
  # The axis percentages are each eigenvalue over the sum of ALL eigenvalues,
  # which is what the original analysis did and therefore what the published
  # figure shows.
  #
  # This distance is not Euclidean, so cmdscale returns some negative
  # eigenvalues. They make sum(pcoa$eig) smaller than the sum of the positive
  # ones alone, which inflates every percentage -- by roughly a factor of two
  # for this data. Dividing instead by the sum of the positive eigenvalues
  # would be the more defensible statistic, and gives axis 1 of the Kigali
  # analysis as 23.4% rather than the published 46.3%.
  #
  # The published value is kept. If the Methods describe these as "variance
  # explained", it is worth saying they are percentages of the signed
  # eigenvalue total, because a reader who recomputes them the other way will
  # get half the number and think something is wrong.
  # ---------------------------------------------------------------------------

  variance_percent <- 100 * pcoa$eig / sum(pcoa$eig)

  coordinates <- as.data.frame(pcoa$points)

  names(coordinates) <- paste0("Axis", seq_len(ncol(coordinates)))

  coordinates$SAMPLE <- rownames(pcoa$points)

  list(
    coordinates = coordinates,
    variance    = variance_percent,
    eigenvalues = pcoa$eig
  )
}

run_permanova <- function(distance_object, subset_metadata, group_column) {

  values <- subset_metadata[[group_column]]

  keep <- !is.na(values)

  if (sum(keep) < 3) {
    return(NULL)
  }

  kept_metadata <- subset_metadata[keep, , drop = FALSE]

  kept_distance <- stats::as.dist(
    as.matrix(distance_object)[kept_metadata$SAMPLE, kept_metadata$SAMPLE]
  )

  kept_metadata[[group_column]] <- droplevels(
    factor(kept_metadata[[group_column]])
  )

  if (nlevels(kept_metadata[[group_column]]) < 2) {
    return(NULL)
  }

  set.seed(SEED)

  formula <- stats::as.formula(paste("kept_distance ~", group_column))

  result <- vegan::adonis2(
    formula,
    data = kept_metadata,
    permutations = PERMUTATIONS
  )

  list(
    N        = nrow(kept_metadata),
    levels   = nlevels(kept_metadata[[group_column]]),
    Df       = result$Df[1],
    SumOfSqs = result$SumOfSqs[1],
    R2       = result$R2[1],
    F        = result$F[1],
    p_value  = result$`Pr(>F)`[1]
  )
}

# The smallest p-value adonis2 can return is 1/(permutations + 1), so it is
# reported with "=" rather than "<". See the note at the top of this script.
format_permanova_label <- function(statistics) {

  # The superscript 2 is written as an escape rather than pasted in, so this
  # file stays pure ASCII. R turns "²" into the character at parse time,
  # which means the label is correct regardless of the file's encoding or the
  # locale the script is run in.
  sprintf(
    "PERMANOVA: p = %s | R\u00b2 = %.3f",
    formatC(statistics$p_value, format = "f", digits = 3),
    statistics$R2
  )
}


# ==============================================================================
# 6. DRAWING
# ==============================================================================

axis_label <- function(pcoa, axis) {
  sprintf("PCoA%d (%.1f%%)", axis, pcoa$variance[axis])
}

make_panel <- function(
  pcoa,
  subset_metadata,
  colour_by,
  shape_by,
  axis_x,
  axis_y,
  statistics_label
) {

  plot_data <- pcoa$coordinates %>%
    dplyr::left_join(subset_metadata, by = "SAMPLE")

  x_name <- paste0("Axis", axis_x)
  y_name <- paste0("Axis", axis_y)

  colour_values <- palette_for(colour_by)

  plot_object <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data[[x_name]], y = .data[[y_name]])
  )

  if (is.na(shape_by)) {

    plot_object <- plot_object +
      ggplot2::geom_point(
        ggplot2::aes(colour = .data[[colour_by]]),
        size = 3, alpha = 0.9
      )

  } else {

    plot_object <- plot_object +
      ggplot2::geom_point(
        ggplot2::aes(
          colour = .data[[colour_by]],
          shape  = .data[[shape_by]]
        ),
        size = 2.4, alpha = 0.9
      ) +
      ggplot2::scale_shape_manual(
        values = shapes_for(shape_by),
        name = legend_title_for(shape_by),
        drop = FALSE
      )
  }

  # Ellipses. For a group comparison, one per group. For the kelch13 panels,
  # only around the double mutant, as in the published figure.
  if (colour_by == "k13_status") {

    ellipse_data <- dplyr::filter(
      plot_data, .data[[colour_by]] == ELLIPSE_GENOTYPE
    )

    if (nrow(ellipse_data) >= 4) {

      plot_object <- plot_object +
        ggplot2::stat_ellipse(
          data = ellipse_data,
          ggplot2::aes(x = .data[[x_name]], y = .data[[y_name]]),
          inherit.aes = FALSE,
          type = "norm", level = 0.95,
          colour = unname(colour_values[[ELLIPSE_GENOTYPE]]),
          linewidth = 1.2
        )
    }

  } else {

    group_sizes <- table(plot_data[[colour_by]])

    if (all(group_sizes >= 4)) {

      plot_object <- plot_object +
        ggplot2::stat_ellipse(
          ggplot2::aes(colour = .data[[colour_by]]),
          type = "norm", level = 0.95, linewidth = 1.2
        )
    }
  }

  plot_object +

    ggplot2::geom_hline(yintercept = 0, colour = "grey25", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey25", linewidth = 0.5) +

    ggplot2::scale_colour_manual(
      values = colour_values,
      name = legend_title_for(colour_by),
      drop = FALSE,
      na.value = "grey70"
    ) +

    ggplot2::annotate(
      "text", x = Inf, y = Inf, label = statistics_label,
      hjust = 1.05, vjust = 1.5, size = 4.6, fontface = "bold"
    ) +

    ggplot2::coord_cartesian(clip = "off") +

    ggplot2::labs(
      x = axis_label(pcoa, axis_x),
      y = axis_label(pcoa, axis_y)
    ) +

    ggplot2::theme_classic(base_size = 16) +

    ggplot2::theme(
      axis.title   = ggplot2::element_text(face = "bold", colour = "black"),
      axis.text    = ggplot2::element_text(face = "bold", colour = "black"),
      legend.title = ggplot2::element_text(face = "bold", colour = "black"),
      legend.text  = ggplot2::element_text(colour = "black")
    )
}


save_figure <- function(plot_object, basename, width, height) {

  # Only SVG and PNG. R's PDF devices cannot write the superscript 2 in the
  # R-squared label and would drop it silently; the same problem the Figure 5
  # script documents. Convert to PDF from the SVG if a PDF is needed.
  for (extension in c("svg", "png")) {

    ggplot2::ggsave(
      filename  = file.path(FIGURE_DIR, paste0(basename, ".", extension)),
      plot      = plot_object,
      width     = width,
      height    = height,
      units     = "in",
      dpi       = 600,
      bg        = "white",
      limitsize = FALSE
    )
  }
}


# ==============================================================================
# 7. RUN EVERY ANALYSIS
# ==============================================================================

permanova_rows   <- list()
coordinate_rows  <- list()
eigenvalue_rows  <- list()

for (index in seq_len(nrow(PANEL_MAP))) {

  specification <- PANEL_MAP[index, ]

  subset_data <- subsets[[specification$subset]]

  add_to_report("")
  add_to_report(strrep("-", 72))
  add_to_report(sprintf(
    "%s / %s  —  %s",
    specification$main_panel, specification$supp_panel, specification$analysis
  ))
  add_to_report(strrep("-", 72))

  pcoa <- run_pcoa(subset_data$distance)

  add_to_report(sprintf(
    "  variance explained: %s",
    paste(
      sprintf("axis%d %.1f%%", seq_len(N_AXES), pcoa$variance[seq_len(N_AXES)]),
      collapse = "  "
    )
  ))

  # The test whose result is printed on the panel, reproducing the published
  # figure.
  displayed <- run_permanova(
    subset_data$distance, subset_data$metadata, specification$permanova_on
  )

  if (is.null(displayed)) {
    stop(
      "The PERMANOVA on '", specification$permanova_on, "' could not be run ",
      "for the '", specification$subset, "' subset.",
      call. = FALSE
    )
  }

  add_to_report(sprintf(
    "  PERMANOVA on %-10s N = %d   R2 = %.4f   p = %.3f   (shown on the panel)",
    specification$permanova_on, displayed$N, displayed$R2, displayed$p_value
  ))

  permanova_rows[[length(permanova_rows) + 1]] <- tibble::tibble(
    published_panel  = paste(specification$main_panel, specification$supp_panel,
                             sep = " / "),
    analysis         = specification$analysis,
    subset           = specification$subset,
    tested_variable  = specification$permanova_on,
    shown_on_panel   = TRUE,
    N                = displayed$N,
    Df               = displayed$Df,
    SumOfSqs         = displayed$SumOfSqs,
    R2               = displayed$R2,
    F                = displayed$F,
    p_value          = displayed$p_value,
    permutations     = PERMUTATIONS
  )

  # Where the points are coloured by kelch13 genotype but the test shown is
  # not on kelch13, also run the kelch13 test so the number exists. It is
  # written to the results table and is NOT put on the figure.
  if (specification$colour_by == "k13_status" &&
      specification$permanova_on != "k13_status") {

    k13_test <- run_permanova(
      subset_data$distance, subset_data$metadata, "k13_status"
    )

    if (!is.null(k13_test)) {

      add_to_report(sprintf(
        "  PERMANOVA on %-10s N = %d   R2 = %.4f   p = %.3f   (NOT shown; see below)",
        "k13_status", k13_test$N, k13_test$R2, k13_test$p_value
      ))

      permanova_rows[[length(permanova_rows) + 1]] <- tibble::tibble(
        published_panel  = paste(specification$main_panel,
                                 specification$supp_panel, sep = " / "),
        analysis         = specification$analysis,
        subset           = specification$subset,
        tested_variable  = "k13_status",
        shown_on_panel   = FALSE,
        N                = k13_test$N,
        Df               = k13_test$Df,
        SumOfSqs         = k13_test$SumOfSqs,
        R2               = k13_test$R2,
        F                = k13_test$F,
        p_value          = k13_test$p_value,
        permutations     = PERMUTATIONS
      )
    }
  }

  statistics_label <- format_permanova_label(displayed)


  # ----------------------------------------------------------------------------
  # The single main-text panel.
  # ----------------------------------------------------------------------------

  main_panel <- make_panel(
    pcoa, subset_data$metadata,
    colour_by = specification$colour_by,
    shape_by  = specification$shape_by,
    axis_x    = specification$main_axis_x,
    axis_y    = specification$main_axis_y,
    statistics_label = statistics_label
  )

  main_basename <- sprintf(
    "figure%s_%s_axis%d_vs_axis%d",
    specification$main_panel, specification$analysis,
    specification$main_axis_x, specification$main_axis_y
  )

  save_figure(main_panel, main_basename, width = 9, height = 7)

  add_to_report(sprintf(
    "  Figure %s: axes %d and %d  ->  %s.svg",
    specification$main_panel,
    specification$main_axis_x, specification$main_axis_y,
    main_basename
  ))


  # ----------------------------------------------------------------------------
  # The 2x2 supplementary composite.
  # ----------------------------------------------------------------------------

  composite_panels <- lapply(
    COMPOSITE_AXES,
    function(pair) {
      make_panel(
        pcoa, subset_data$metadata,
        colour_by = specification$colour_by,
        shape_by  = specification$shape_by,
        axis_x    = pair[[1]],
        axis_y    = pair[[2]],
        statistics_label = statistics_label
      )
    }
  )

  composite <- (composite_panels[[1]] | composite_panels[[2]]) /
               (composite_panels[[3]] | composite_panels[[4]]) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "right")

  supp_basename <- sprintf(
    "figure%s_%s_four_panels",
    specification$supp_panel, specification$analysis
  )

  save_figure(composite, supp_basename, width = 16, height = 11)

  add_to_report(sprintf(
    "  Figure %s: four axis pairs  ->  %s.svg",
    specification$supp_panel, supp_basename
  ))


  # ----------------------------------------------------------------------------
  # Source data.
  # ----------------------------------------------------------------------------

  coordinate_rows[[length(coordinate_rows) + 1]] <-
    pcoa$coordinates %>%
    dplyr::left_join(subset_data$metadata, by = "SAMPLE") %>%
    dplyr::mutate(
      published_panel = paste(specification$main_panel,
                              specification$supp_panel, sep = " / "),
      analysis = specification$analysis,
      .before = 1
    )

  eigenvalue_rows[[length(eigenvalue_rows) + 1]] <- tibble::tibble(
    analysis            = specification$analysis,
    axis                = seq_along(pcoa$eigenvalues),
    eigenvalue          = pcoa$eigenvalues,
    percent_of_positive = pcoa$variance
  )
}


# ==============================================================================
# 8. WRITE THE TABLES
# ==============================================================================

permanova_results <- dplyr::bind_rows(permanova_rows)

readr::write_tsv(
  permanova_results,
  file.path(TABLE_DIR, "permanova_results.tsv")
)

readr::write_tsv(
  dplyr::bind_rows(coordinate_rows),
  file.path(TABLE_DIR, "pcoa_coordinates.tsv")
)

readr::write_tsv(
  dplyr::bind_rows(eigenvalue_rows),
  file.path(TABLE_DIR, "pcoa_eigenvalues.tsv")
)

readr::write_tsv(
  dplyr::count(metadata, period, region, k13_status, name = "n"),
  file.path(TABLE_DIR, "figure04_sample_counts.tsv")
)


# ==============================================================================
# 9. SUMMARY
# ==============================================================================

add_to_report("")
add_to_report(strrep("=", 72))
add_to_report("PERMANOVA RESULTS")
add_to_report(strrep("=", 72))
add_to_report(sprintf(
  "%-12s %-14s %6s %8s %8s  %s",
  "panel", "tested on", "N", "R2", "p", "on the figure?"
))

for (index in seq_len(nrow(permanova_results))) {

  row <- permanova_results[index, ]

  add_to_report(sprintf(
    "%-12s %-14s %6d %8.4f %8.3f  %s",
    row$published_panel, row$tested_variable, row$N, row$R2, row$p_value,
    if (row$shown_on_panel) "yes" else "no"
  ))
}

add_to_report("")
add_to_report(sprintf(
  "p-values come from %d permutations, so the smallest attainable value is",
  PERMUTATIONS
))
add_to_report(sprintf(
  "%.3f. Report it as \"p = %.3f\", not \"p < %.3f\".",
  1 / (PERMUTATIONS + 1), 1 / (PERMUTATIONS + 1), 1 / (PERMUTATIONS + 1)
))

not_shown <- dplyr::filter(permanova_results, !shown_on_panel)

if (nrow(not_shown) > 0) {

  add_to_report("")
  add_to_report(strrep("!", 72))
  add_to_report("A PANEL COLOURED BY kelch13 GENOTYPE SHOWS A DIFFERENT TEST")
  add_to_report(strrep("!", 72))

  for (index in seq_len(nrow(not_shown))) {

    row <- not_shown[index, ]

    shown <- permanova_results %>%
      dplyr::filter(
        published_panel == row$published_panel, shown_on_panel
      )

    add_to_report("")
    add_to_report(sprintf("Panel %s", row$published_panel))
    add_to_report(sprintf(
      "  printed on the figure:  %s,  R2 = %.4f",
      shown$tested_variable[[1]], shown$R2[[1]]
    ))
    add_to_report(sprintf(
      "  kelch13 genotype test:  k13_status,  R2 = %.4f,  p = %.3f",
      row$R2, row$p_value
    ))
  }

  add_to_report("")
  add_to_report("This reproduces the published figure, which shows the first")
  add_to_report("value. The kelch13 test is given here so the number exists if")
  add_to_report("the legend is corrected. Check what the legend claims was")
  add_to_report("tested before the proofs are returned.")
  add_to_report(strrep("!", 72))
}

add_to_report("")
add_to_report("Figures:  ", FIGURE_DIR)
add_to_report("Tables:   ", TABLE_DIR)
add_to_report("")
add_to_report("Figure 4d is a distance tree and is not produced by this script.")

writeLines(report_lines, file.path(TABLE_DIR, "figure04_report.txt"))

capture.output(
  sessionInfo(),
  file = file.path(TABLE_DIR, "figure04_sessionInfo.txt")
)

message("\nReport written to:\n  ",
        file.path(TABLE_DIR, "figure04_report.txt"))
