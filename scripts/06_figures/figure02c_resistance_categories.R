# Figure 2c — composite antimalarial resistance backgrounds
#
# Inputs:  Manuscript_data_29-07-2026.tsv
#          sample appears once for each mutation it carries. Set MANUSCRIPT_DATA in
#          config.R.
# Outputs: results/figures/figure02/figure02c_resistance_categories.svg / .png
#          results/tables/figure02/figure02c_category_prevalence.tsv
#          results/tables/figure02/figure02c_sample_level_categories.tsv
#          results/tables/figure02/figure02c_denominator_check.tsv
#          results/tables/figure02/figure02c_report.txt
#
# Usage:   Rscript scripts/06_figures/figure02c_resistance_categories.R

options(error = NULL)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(ggplot2)
})


find_project_root <- function() {

  candidate <- normalizePath(getwd(), mustWork = TRUE)

  for (level in seq_len(5)) {
    if (file.exists(file.path(candidate, "config.R"))) {
      return(candidate)
    }
    parent <- dirname(candidate)
    if (identical(parent, candidate)) break
    candidate <- parent
  }

  stop("Could not find config.R. Run this from the top of the repository.",
       call. = FALSE)
}

PROJECT_ROOT <- find_project_root()

source(file.path(PROJECT_ROOT, "config.R"))


arguments <- commandArgs(trailingOnly = TRUE)

get_argument <- function(flag, default = NULL) {
  position <- match(flag, arguments)
  if (is.na(position)) return(default)
  if (position == length(arguments)) {
    stop("The ", flag, " flag needs a value after it.", call. = FALSE)
  }
  arguments[[position + 1]]
}

INPUT_FILE <- get_argument(
  "--input",
  if (exists("MANUSCRIPT_DATA", inherits = TRUE)) {
    get("MANUSCRIPT_DATA", inherits = TRUE)
  } else {
    NULL
  }
)

if (is.null(INPUT_FILE)) {
  stop(
    "No input file given.\n\n",
    "Set MANUSCRIPT_DATA in config.R, or pass:\n",
    "  --input /path/to/Manuscript_data_29-07-2026.tsv",
    call. = FALSE
  )
}

if (!file.exists(INPUT_FILE)) {
  stop("Input file not found:\n  ", INPUT_FILE, call. = FALSE)
}

FIGURE_DIR <- file.path(PROJECT_ROOT, "results", "figures", "figure02")
TABLE_DIR  <- file.path(PROJECT_ROOT, "results", "tables",  "figure02")

dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR,  recursive = TRUE, showWarnings = FALSE)


report_lines <- character(0)

add_to_report <- function(...) {
  line <- paste0(...)
  report_lines <<- c(report_lines, line)
  message(line)
}


# ---- the markers and the categories ------------------------------------

K13_VALIDATED <- c(
  "F446I", "N458Y", "C469Y", "M476I", "Y493H", "R539T", "I543T",
  "P553L", "R561H", "P574L", "C580Y", "R622I", "A675V"
)

K13_CANDIDATE <- c(
  "P441L", "G449A", "C469F", "A481V", "R515K", "P527H",
  "N537I", "N537D", "G538V", "V568G"
)

K13_ALL <- c(K13_VALIDATED, K13_CANDIDATE)

# Single mutations the categories are built from, as GENE:CHANGE.
REQUIRED_MUTATIONS <- c(
  "DHFR:N51I", "DHFR:C59R", "DHFR:S108N", "DHFR:I164L",
  "DHPS:K540E", "DHPS:A581G",
  "MDR1:Y184F",
  "CRT:K76T",
  paste0("K13:", K13_ALL)
)

# Each category, the label on the axis, what it means, the published fill,
# and the count printed on the published panel.
CATEGORY_DEFINITIONS <- tibble::tribble(
  ~category,              ~short_label,     ~definition,                                                   ~fill,      ~published_n,
  "SP_triple_mutant",     "SP triple",      "DHFR N51I + C59R + S108N",                                    "#D9D9D9",  1129L,
  "SP_quadruple_mutant",  "SP quadruple",   "DHFR N51I + C59R + S108N + DHPS K540E",                       "#4A6FD4",  1101L,
  "high_SP_resistance",   "High-level SP",  "DHFR N51I + C59R + S108N + I164L + DHPS K540E + A581G",       "#C6D2F2",   140L,
  "CQ_R",                 "CQ-R",           "CRT K76T",                                                    "#36476B",   128L,
  "partner_drug",         "Partner-drug",   "MDR1 Y184F",                                                  "#4FC3C7",   752L,
  "ART_R",                "ART-R",          "K13 validated or candidate ART-R mutation",                   "#F57C20",   527L,
  "ACT",                  "ACT",            "MDR1 Y184F + K13 ART-R mutation",                             "#A8DADC",   386L,
  "Multi_resistant",      "Multi-drug",     "SP triple + DHPS K540E + MDR1 Y184F + CRT K76T + K13 ART-R",  "#F5B820",    24L
)

# The denominator printed on the published panel.
PUBLISHED_N_TOTAL <- 1176L


# ---- read --------------------------------------------------------------

add_to_report(strrep("=", 72))
add_to_report("FIGURE 2C — COMPOSITE RESISTANCE BACKGROUNDS")
add_to_report(strrep("=", 72))
add_to_report("Input:     ", INPUT_FILE)
add_to_report("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
add_to_report("")

raw <- readr::read_tsv(
  INPUT_FILE,
  col_types = readr::cols(.default = readr::col_character()),
  show_col_types = FALSE
)

required_columns <- c("SAMPLE", "GENE_NAME", "HGVS_PROTEIN")

missing_columns <- setdiff(required_columns, names(raw))

if (length(missing_columns) > 0) {
  stop(
    "The input is missing: ", paste(missing_columns, collapse = ", "),
    "\nIt should have one row per sample and mutation.",
    call. = FALSE
  )
}

add_to_report(sprintf("Rows read: %d", nrow(raw)))

all_samples <- raw %>%
  dplyr::filter(!is.na(SAMPLE), stringr::str_trim(SAMPLE) != "") %>%
  dplyr::mutate(SAMPLE = stringr::str_trim(SAMPLE)) %>%
  dplyr::distinct(SAMPLE) %>%
  dplyr::pull(SAMPLE)

add_to_report(sprintf("Distinct samples in the file: %d", length(all_samples)))


# ---- the sample-by-mutation matrix -------------------------------------

mutation_rows <- raw %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(required_columns), ~ stringr::str_trim(.x)
    )
  ) %>%
  dplyr::filter(
    !is.na(SAMPLE),       SAMPLE != "",
    !is.na(GENE_NAME),    GENE_NAME != "",
    !is.na(HGVS_PROTEIN), HGVS_PROTEIN != "", HGVS_PROTEIN != "."
  ) %>%
  dplyr::mutate(mutation = paste(GENE_NAME, HGVS_PROTEIN, sep = ":")) %>%
  dplyr::distinct(SAMPLE, mutation)

add_to_report(sprintf(
  "Sample-mutation pairs: %d across %d distinct mutations",
  nrow(mutation_rows), dplyr::n_distinct(mutation_rows$mutation)
))

variants_wide <- mutation_rows %>%
  dplyr::mutate(present = 1L) %>%
  tidyr::pivot_wider(
    names_from = mutation,
    values_from = present,
    values_fill = 0L
  )


# A mutation absent from every sample produces no column, so add it as zero.
# Without this, referring to it below would fail rather than evaluate as
# "nobody has it".

absent_mutations <- setdiff(REQUIRED_MUTATIONS, names(variants_wide))

if (length(absent_mutations) > 0) {

  variants_wide[absent_mutations] <- 0L

  add_to_report("")
  add_to_report(sprintf(
    "%d required mutation(s) are absent from the data and were added as zero:",
    length(absent_mutations)
  ))

  for (mutation in absent_mutations) {
    add_to_report("  ", mutation)
  }
}


# ---- the denominator ---------------------------------------------------

n_with_mutations <- nrow(variants_wide)
n_all_samples    <- length(all_samples)

samples_without_mutations <- setdiff(all_samples, variants_wide$SAMPLE)

add_to_report("")
add_to_report(strrep("-", 72))
add_to_report("DENOMINATOR")
add_to_report(strrep("-", 72))
add_to_report(sprintf("  samples in the file:            %d", n_all_samples))
add_to_report(sprintf("  samples with >= 1 mutation:     %d", n_with_mutations))
add_to_report(sprintf("  samples with no mutation:       %d",
                      length(samples_without_mutations)))

readr::write_tsv(
  tibble::tibble(
    measure = c(
      "samples in the file",
      "samples with at least one mutation",
      "samples with no mutation"
    ),
    n = c(n_all_samples, n_with_mutations, length(samples_without_mutations))
  ),
  file.path(TABLE_DIR, "figure02c_denominator_check.tsv")
)

add_to_report(sprintf("  published panel says:           %d",
                      PUBLISHED_N_TOTAL))

if (length(samples_without_mutations) == 0) {

  add_to_report("")
  add_to_report("  Every sample carries at least one mutation, so the two")
  add_to_report("  denominators are the same and the distinction does not")
  add_to_report("  affect the figure.")

  if (n_with_mutations != PUBLISHED_N_TOTAL) {
    add_to_report("")
    add_to_report(strrep("!", 72))
    add_to_report(sprintf(
      "  But the denominator is %d, not the published %d.",
      n_with_mutations, PUBLISHED_N_TOTAL
    ))
    add_to_report("  Every percentage below will differ from the paper.")
    add_to_report("  Check that the input file is the one the analysis used.")
    add_to_report(strrep("!", 72))
  }

} else {

  add_to_report("")
  add_to_report(strrep("!", 72))
  add_to_report(sprintf(
    "%d sample(s) carry no mutation and are absent from the matrix.",
    length(samples_without_mutations)
  ))
  add_to_report("")
  add_to_report("Prevalence below is therefore out of the mutated samples, not")
  add_to_report("out of all sequenced samples -- which is what the original")
  add_to_report("code did, so the published figure is reproduced. But the")
  add_to_report("prevalences are higher than they would be over the full set.")
  add_to_report("")
  add_to_report("Check that the figure legend states the denominator. The bars")
  add_to_report("print n/N, so a reader can see it either way.")
  add_to_report(strrep("!", 72))
}

# The original denominator, kept so the published figure is reproduced.
N_TOTAL <- n_with_mutations

if (N_TOTAL == 0) {
  stop("No samples carry any mutation; there is nothing to plot.",
       call. = FALSE)
}


# ---- the composite categories ------------------------------------------

k13_columns <- paste0("K13:", K13_ALL)

categories <- variants_wide %>%
  dplyr::mutate(

    # Any validated or candidate kelch13 marker.
    K13_ART = as.integer(
      rowSums(dplyr::across(dplyr::all_of(k13_columns)), na.rm = TRUE) > 0
    ),

    SP_triple_mutant = as.integer(
      `DHFR:N51I` == 1 & `DHFR:C59R` == 1 & `DHFR:S108N` == 1
    ),

    SP_quadruple_mutant = as.integer(
      `DHFR:N51I` == 1 & `DHFR:C59R` == 1 & `DHFR:S108N` == 1 &
        `DHPS:K540E` == 1
    ),

    high_SP_resistance = as.integer(
      `DHFR:N51I` == 1 & `DHFR:C59R` == 1 & `DHFR:S108N` == 1 &
        `DHFR:I164L` == 1 & `DHPS:K540E` == 1 & `DHPS:A581G` == 1
    ),

    CQ_R = as.integer(`CRT:K76T` == 1),

    partner_drug = as.integer(`MDR1:Y184F` == 1),

    ART_R = as.integer(K13_ART == 1),

    ACT = as.integer(`MDR1:Y184F` == 1 & K13_ART == 1),

    Multi_resistant = as.integer(
      `DHFR:N51I` == 1 & `DHFR:C59R` == 1 & `DHFR:S108N` == 1 &
        `DHPS:K540E` == 1 & `MDR1:Y184F` == 1 & `CRT:K76T` == 1 &
        K13_ART == 1
    )
  )

missing_categories <- setdiff(CATEGORY_DEFINITIONS$category, names(categories))

if (length(missing_categories) > 0) {
  stop(
    "These categories were not created: ",
    paste(missing_categories, collapse = ", "),
    call. = FALSE
  )
}


# A consistency check the original did not do. Some categories are nested: a
# quadruple mutant must also be a triple mutant, and anything multi-resistant
# must be a quadruple mutant. If a count rises where it must fall, a
# definition is wrong.

NESTED_PAIRS <- list(
  c("SP_quadruple_mutant", "SP_triple_mutant"),
  c("high_SP_resistance",  "SP_quadruple_mutant"),
  c("Multi_resistant",     "SP_quadruple_mutant"),
  c("ACT",                 "ART_R"),
  c("ACT",                 "partner_drug"),
  c("Multi_resistant",     "ART_R"),
  c("Multi_resistant",     "CQ_R")
)

nesting_problems <- character(0)

for (pair in NESTED_PAIRS) {

  narrower <- pair[[1]]
  broader  <- pair[[2]]

  violations <- sum(
    categories[[narrower]] == 1 & categories[[broader]] == 0
  )

  if (violations > 0) {
    nesting_problems <- c(
      nesting_problems,
      sprintf(
        "%d sample(s) are %s but not %s", violations, narrower, broader
      )
    )
  }
}

add_to_report("")
add_to_report("Nesting check (a narrower category must imply the broader one):")

if (length(nesting_problems) == 0) {
  add_to_report("  PASS — every nested relationship holds.")
} else {
  for (problem in nesting_problems) {
    add_to_report("  FAIL — ", problem)
  }
  warning(
    "Category nesting is violated; a definition in section 6 is wrong. ",
    "See the report.",
    call. = FALSE
  )
}


# ---- prevalence --------------------------------------------------------

prevalence <- categories %>%
  dplyr::select(dplyr::all_of(CATEGORY_DEFINITIONS$category)) %>%
  dplyr::summarise(
    dplyr::across(dplyr::everything(), ~ sum(.x, na.rm = TRUE))
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::everything(), names_to = "category", values_to = "n"
  ) %>%
  dplyr::left_join(CATEGORY_DEFINITIONS, by = "category") %>%
  dplyr::mutate(
    N_total    = N_TOTAL,
    prevalence = 100 * n / N_TOTAL,
    # The published panel labels each bar with the bare count and states the
    # denominator once, above the plot.
    label      = format(n, big.mark = "", trim = TRUE),
    category     = factor(category,    levels = CATEGORY_DEFINITIONS$category),
    short_label  = factor(short_label, levels = CATEGORY_DEFINITIONS$short_label),
    definition   = factor(definition,  levels = CATEGORY_DEFINITIONS$definition)
  ) %>%
  dplyr::arrange(category)

add_to_report("")
add_to_report(strrep("-", 72))
add_to_report("PREVALENCE")
add_to_report(strrep("-", 72))
add_to_report(sprintf("%-16s %8s %8s %12s",
                      "category", "n", "percent", "published"))

for (index in seq_len(nrow(prevalence))) {
  add_to_report(sprintf(
    "%-16s %8d %7.1f%% %12d",
    as.character(prevalence$short_label[[index]]),
    prevalence$n[[index]],
    prevalence$prevalence[[index]],
    prevalence$published_n[[index]]
  ))
}


# Check against the published panel.

count_mismatches <- prevalence %>%
  dplyr::filter(n != published_n)

add_to_report("")

if (nrow(count_mismatches) == 0 && N_TOTAL == PUBLISHED_N_TOTAL) {

  add_to_report("  Every count matches published Figure 2c exactly.")

} else {

  add_to_report(strrep("!", 72))
  add_to_report("  These counts do NOT match published Figure 2c:")
  add_to_report("")

  for (index in seq_len(nrow(count_mismatches))) {
    add_to_report(sprintf(
      "    %-16s computed %6d, published %6d",
      as.character(count_mismatches$short_label[[index]]),
      count_mismatches$n[[index]],
      count_mismatches$published_n[[index]]
    ))
  }

  add_to_report("")
  add_to_report("  Check that the input file is the one the analysis used,")
  add_to_report("  and that the marker lists in section 2 are unchanged.")
  add_to_report(strrep("!", 72))
}


# ---- the figure --------------------------------------------------------

figure <- ggplot2::ggplot(
  prevalence,
  ggplot2::aes(x = short_label, y = prevalence, fill = definition)
) +

  ggplot2::geom_col(width = 0.7) +

  ggplot2::geom_text(
    ggplot2::aes(label = label), vjust = -0.5, size = 3.2, colour = "black"
  ) +

  # The denominator, stated once above the plot as in the published panel.
  ggplot2::annotate(
    "text",
    x = (nrow(prevalence) + 1) / 2, y = 108,
    label = paste0("N = ", N_TOTAL),
    size = 4, colour = "black"
  ) +

  ggplot2::scale_fill_manual(
    values = setNames(
      CATEGORY_DEFINITIONS$fill, CATEGORY_DEFINITIONS$definition
    ),
    drop = FALSE
  ) +

  ggplot2::scale_y_continuous(
    breaks = seq(0, 100, by = 25),
    expand = ggplot2::expansion(mult = c(0, 0.08))
  ) +

  ggplot2::coord_cartesian(ylim = c(0, 110), clip = "off") +

  ggplot2::labs(
    x = NULL,
    y = "Prevalence (%)",
    fill = "Mutation composition"
  ) +

  ggplot2::theme_classic(base_size = 12) +

  ggplot2::theme(
    axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1,
                                            colour = "black"),
    axis.text.y     = ggplot2::element_text(colour = "black"),
    axis.title.y    = ggplot2::element_text(face = "bold"),
    axis.line.x     = ggplot2::element_blank(),
    axis.ticks.x    = ggplot2::element_blank(),
    legend.position = "right",
    legend.title    = ggplot2::element_text(face = "bold"),
    plot.margin     = ggplot2::margin(20, 10, 10, 10)
  )

for (extension in c("svg", "png")) {

  ggplot2::ggsave(
    filename = file.path(
      FIGURE_DIR, paste0("figure02c_resistance_categories.", extension)
    ),
    plot = figure,
    width = 11, height = 7, units = "in",
    dpi = 600, bg = "white"
  )
}


# ---- tables ------------------------------------------------------------

# The Source Data for this panel. `fill` is styling and `label` is the text
# drawn on the bar, so neither belongs in a data table; published_n is kept
# so the comparison against the paper travels with the numbers.
readr::write_tsv(
  prevalence %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.factor), as.character)) %>%
    dplyr::select(
      category, short_label, definition, n, N_total, prevalence, published_n
    ),
  file.path(TABLE_DIR, "figure02c_category_prevalence.tsv")
)

readr::write_tsv(
  categories %>%
    dplyr::select(
      SAMPLE, K13_ART, dplyr::all_of(CATEGORY_DEFINITIONS$category)
    ),
  file.path(TABLE_DIR, "figure02c_sample_level_categories.tsv")
)

add_to_report("")
add_to_report(strrep("=", 72))
add_to_report("WRITTEN")
add_to_report(strrep("=", 72))
add_to_report("  figures/figure02/figure02c_resistance_categories.svg / .png")
add_to_report("  tables/figure02/figure02c_category_prevalence.tsv")
add_to_report("  tables/figure02/figure02c_sample_level_categories.tsv")
add_to_report("  tables/figure02/figure02c_denominator_check.tsv")

writeLines(report_lines, file.path(TABLE_DIR, "figure02c_report.txt"))

capture.output(
  sessionInfo(),
  file = file.path(TABLE_DIR, "figure02c_sessionInfo.txt")
)
