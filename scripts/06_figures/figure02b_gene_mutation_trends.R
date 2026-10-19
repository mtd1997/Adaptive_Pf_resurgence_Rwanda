################################################################################
# FIGURE 2B — PER-GENE MUTATION PREVALENCE BY YEAR
################################################################################
#
# WHAT IT DOES
#   For each resistance gene and each sampling year, counts the samples
#   carrying at least one non-synonymous mutation in that gene, and plots
#   that count as a percentage of the samples sequenced in that year.
#
#   The unit is the sample, not the mutation: a sample carrying three DHFR
#   mutations counts once towards DHFR. This is a prevalence of mutated
#   samples, not a mutation frequency.
#
# INPUT
#   Manuscript_data_29-07-2026.tsv -- one row per sample and mutation. Set
#   MANUSCRIPT_DATA in config.R.
#
# OUTPUTS
#   results/figures/figure02/figure02b_gene_mutation_trends.svg / .png
#   results/tables/figure02/figure02b_gene_prevalence_by_year.tsv
#   results/tables/figure02/figure02b_samples_per_year.tsv
#   results/tables/figure02/figure02b_report.txt
#
#   The prevalence table is the Source Data for this panel.
#
# -----------------------------------------------------------------------------
# WHAT THE PUBLISHED PANEL SHOWS, AND WHAT THIS REPRODUCES
# -----------------------------------------------------------------------------
#   Published Figure 2b carries the per-year denominators under the axis:
#
#       2021 N =  30      2024 N = 234
#       2022 N = 196      2025 N = 685
#       2023 N =  31
#
#   which total 1,176 -- the full analysis set. The script checks its own
#   totals against those numbers and says plainly whether they match, so a
#   drift in the input shows up as a message rather than as a quietly
#   different figure.
#
#   One difference from the published panel is cosmetic and cannot be made
#   in ggplot without hand-placing text: in the published version the gene
#   names sit at the right-hand end of each line rather than in a legend
#   box. That was done when the figure was assembled. The legend here
#   carries the same information.
# -----------------------------------------------------------------------------
#
# -----------------------------------------------------------------------------
# TWO THINGS WORTH KNOWING BEFORE READING THE NUMBERS
# -----------------------------------------------------------------------------
#   1. The period label comes from the year, not from the metadata column.
#      The original code overrode whatever PERIOD the input carried:
#      2021-2023 became Baseline-Kigali and 2024-2025 became Resurgence.
#      That is kept here, because it matches the sampling design -- the
#      2021-2023 samples were all collected at King Faisal Hospital and were
#      treated as the Kigali baseline. The script reports any sample whose
#      metadata period disagrees with its year, so a mismatch is visible
#      rather than silently rewritten.
#
#   2. A gene-year combination with no mutated sample must show zero, not
#      vanish. Counting alone drops those combinations, which would make a
#      gene's line skip a year and join across the gap -- reading as though
#      prevalence had been interpolated rather than observed as zero. The
#      full gene-by-year grid is therefore built first and the counts joined
#      onto it.
# -----------------------------------------------------------------------------
#
# HOW TO RUN
#     Rscript scripts/06_figures/figure02b_gene_mutation_trends.R
#
#   Optional:
#     --input PATH    use a different input file
#     --genes LIST    comma-separated gene list, default DHFR,DHPS,K13,MDR1,CRT
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

TARGET_GENES <- strsplit(
  get_argument("--genes", "DHFR,DHPS,K13,MDR1,CRT"), ",", fixed = TRUE
)[[1]]

TARGET_GENES <- stringr::str_trim(TARGET_GENES)
TARGET_GENES <- TARGET_GENES[TARGET_GENES != ""]

FIGURE_DIR <- file.path(PROJECT_ROOT, "results", "figures", "figure02")
TABLE_DIR  <- file.path(PROJECT_ROOT, "results", "tables",  "figure02")

dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR,  recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------------
# Published styling. These are the values the published panel was drawn with,
# kept as named constants so the figure cannot drift when the data changes.
# ------------------------------------------------------------------------------

GENE_COLOURS <- c(
  DHFR = "#6690FF",
  DHPS = "#F5A000",
  K13  = "#FF5A00",
  MDR1 = "#DE277D",
  CRT  = "#304D87"
)

GENE_SHAPES <- c(
  DHFR = 22,   # square
  DHPS = 24,   # triangle up
  K13  = 25,   # triangle down
  MDR1 = 23,   # diamond
  CRT  = 21    # circle
)

# scale_*_manual stops if a factor level has no value, so any gene requested
# through --genes that is not in the tables above is given a neutral colour
# and a plain circle rather than being allowed to kill the plot at the end.
unstyled_genes <- setdiff(TARGET_GENES, names(GENE_COLOURS))

if (length(unstyled_genes) > 0) {
  GENE_COLOURS[unstyled_genes] <- "grey50"
  GENE_SHAPES[unstyled_genes]  <- 21
}

PERIOD_BANDS <- tibble::tribble(
  ~period,            ~xmin,    ~xmax,    ~fill,      ~label_x, ~label,
  "Baseline-Kigali",  2020.87,  2023.5,   "#EAF6F8",  2022.0,   "Baseline-Kigali",
  "Resurgence",       2023.5,   2025.10,  "#FCEEEF",  2024.3,   "Resurgence"
)

PERIOD_DIVIDER <- 2023.5

X_LIMITS <- c(2020.87, 2025.10)

PERIOD_FROM_YEAR <- tibble::tribble(
  ~year, ~period,
  2021L, "Baseline-Kigali",
  2022L, "Baseline-Kigali",
  2023L, "Baseline-Kigali",
  2024L, "Resurgence",
  2025L, "Resurgence"
)

# The per-year denominators printed on the published panel. Used only as a
# check; nothing is drawn from them.
PUBLISHED_TOTALS <- tibble::tribble(
  ~year, ~published_n,
  2021L,  30L,
  2022L, 196L,
  2023L,  31L,
  2024L, 234L,
  2025L, 685L
)


report_lines <- character(0)

add_to_report <- function(...) {
  line <- paste0(...)
  report_lines <<- c(report_lines, line)
  message(line)
}


# ==============================================================================
# 2. READ
# ==============================================================================

add_to_report(strrep("=", 72))
add_to_report("FIGURE 2B — PER-GENE MUTATION PREVALENCE BY YEAR")
add_to_report(strrep("=", 72))
add_to_report("Input:     ", INPUT_FILE)
add_to_report("Genes:     ", paste(TARGET_GENES, collapse = ", "))
add_to_report("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
add_to_report("")

raw <- readr::read_tsv(
  INPUT_FILE,
  col_types = readr::cols(.default = readr::col_character()),
  show_col_types = FALSE
)

required_columns <- c("SAMPLE", "YEAR", "GENE_NAME", "HGVS_PROTEIN")

missing_columns <- setdiff(required_columns, names(raw))

if (length(missing_columns) > 0) {
  stop(
    "The input is missing: ", paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

has_period_column <- "PERIOD" %in% names(raw)

clean <- raw %>%
  dplyr::mutate(
    SAMPLE    = stringr::str_trim(SAMPLE),
    GENE_NAME = stringr::str_trim(GENE_NAME),
    year      = suppressWarnings(as.integer(stringr::str_trim(YEAR)))
  ) %>%
  dplyr::filter(!is.na(SAMPLE), SAMPLE != "", !is.na(year))

dropped <- nrow(raw) - nrow(clean)

if (dropped > 0) {
  add_to_report(sprintf(
    "%d row(s) dropped for a missing sample name or an unreadable year.",
    dropped
  ))
}

add_to_report(sprintf("Rows used: %d", nrow(clean)))


# ==============================================================================
# 3. PERIOD FROM YEAR
# ==============================================================================

unknown_years <- setdiff(unique(clean$year), PERIOD_FROM_YEAR$year)

if (length(unknown_years) > 0) {
  stop(
    "These years have no period assigned: ",
    paste(sort(unknown_years), collapse = ", "),
    "\n\nAdd them to PERIOD_FROM_YEAR near the top of this script, so the\n",
    "assignment stays written down in one place rather than being guessed.",
    call. = FALSE
  )
}

# The input may already carry a PERIOD column. The lower-case `period` used
# from here on comes from the year, as in the original analysis; the original
# column is left untouched and only used for the check below.
clean <- clean %>%
  dplyr::select(-dplyr::any_of("period")) %>%
  dplyr::left_join(PERIOD_FROM_YEAR, by = "year")


if (has_period_column) {

  comparison <- tibble::tibble(
    SAMPLE = stringr::str_trim(raw$SAMPLE),
    year   = suppressWarnings(as.integer(stringr::str_trim(raw$YEAR))),
    metadata_period = stringr::str_trim(as.character(raw$PERIOD))
  ) %>%
    dplyr::filter(!is.na(year), !is.na(SAMPLE), SAMPLE != "") %>%
    dplyr::left_join(PERIOD_FROM_YEAR, by = "year") %>%
    dplyr::filter(
      !is.na(metadata_period), metadata_period != "",
      metadata_period != period
    ) %>%
    dplyr::distinct(SAMPLE, year, metadata_period, period)

  add_to_report("")

  if (nrow(comparison) == 0) {
    add_to_report("Period check: the PERIOD column agrees with the year for")
    add_to_report("every sample.")
  } else {
    add_to_report(sprintf(
      "Period check: %d sample(s) have a PERIOD that disagrees with their",
      nrow(comparison)
    ))
    add_to_report("year. The year is used, as in the original analysis.")
    add_to_report("")
    add_to_report(sprintf("  %-14s %6s  %-18s %-18s",
                          "sample", "year", "PERIOD says", "year implies"))
    for (index in seq_len(min(nrow(comparison), 20))) {
      add_to_report(sprintf(
        "  %-14s %6d  %-18s %-18s",
        comparison$SAMPLE[[index]], comparison$year[[index]],
        comparison$metadata_period[[index]], comparison$period[[index]]
      ))
    }
    if (nrow(comparison) > 20) {
      add_to_report(sprintf("  ... and %d more", nrow(comparison) - 20))
    }

    readr::write_tsv(
      comparison,
      file.path(TABLE_DIR, "figure02b_period_disagreements.tsv")
    )
  }
}


# ==============================================================================
# 4. DENOMINATOR: SAMPLES PER YEAR
# ==============================================================================

year_totals <- clean %>%
  dplyr::distinct(SAMPLE, year, period) %>%
  dplyr::count(year, period, name = "n_total") %>%
  dplyr::arrange(year)

total_sample_years <- sum(year_totals$n_total)
n_distinct_samples <- dplyr::n_distinct(clean$SAMPLE)

add_to_report("")
add_to_report(strrep("-", 72))
add_to_report("SAMPLES PER YEAR")
add_to_report(strrep("-", 72))
add_to_report(sprintf("  %6s  %-18s %8s %12s",
                      "year", "period", "n", "published"))

check <- year_totals %>%
  dplyr::left_join(PUBLISHED_TOTALS, by = "year")

for (index in seq_len(nrow(check))) {

  published <- check$published_n[[index]]

  add_to_report(sprintf(
    "  %6d  %-18s %8d %12s",
    check$year[[index]],
    check$period[[index]],
    check$n_total[[index]],
    if (is.na(published)) "-" else format(published)
  ))
}

add_to_report(sprintf("  %6s  %-18s %8d %12d",
                      "", "total", total_sample_years,
                      sum(PUBLISHED_TOTALS$published_n)))

mismatched <- check %>%
  dplyr::filter(!is.na(published_n), n_total != published_n)

add_to_report("")

if (nrow(mismatched) == 0 &&
    total_sample_years == sum(PUBLISHED_TOTALS$published_n)) {

  add_to_report("  Every per-year total matches the published panel.")

} else {

  add_to_report(strrep("!", 72))
  add_to_report("  The per-year totals do NOT match the published panel.")
  add_to_report("  The figure below will differ from Figure 2b as published.")
  add_to_report("  Check that the input file is the one the analysis used.")
  add_to_report(strrep("!", 72))
}

if (total_sample_years != n_distinct_samples) {
  add_to_report("")
  add_to_report(sprintf(
    "NOTE: %d distinct sample names but %d sample-years, so %d sample(s)",
    n_distinct_samples, total_sample_years,
    total_sample_years - n_distinct_samples
  ))
  add_to_report("appear under more than one year. Those samples are counted in")
  add_to_report("each of their years. Listed in figure02b_sample_year_duplicates.tsv.")

  readr::write_tsv(
    clean %>%
      dplyr::distinct(SAMPLE, year) %>%
      dplyr::add_count(SAMPLE, name = "n_years") %>%
      dplyr::filter(n_years > 1) %>%
      dplyr::arrange(SAMPLE, year),
    file.path(TABLE_DIR, "figure02b_sample_year_duplicates.tsv")
  )
}

readr::write_tsv(
  year_totals, file.path(TABLE_DIR, "figure02b_samples_per_year.tsv")
)


# ==============================================================================
# 5. NUMERATOR: MUTATED SAMPLES PER GENE PER YEAR
# ==============================================================================

genes_absent <- setdiff(TARGET_GENES, unique(clean$GENE_NAME))

if (length(genes_absent) > 0) {
  add_to_report("")
  add_to_report(
    "These requested genes do not appear in the input at all: ",
    paste(genes_absent, collapse = ", ")
  )
  add_to_report("They will be drawn at zero across every year.")
}

if (length(genes_absent) == length(TARGET_GENES)) {
  stop(
    "None of the requested genes appear in the input.\n",
    "Genes present: ",
    paste(sort(unique(clean$GENE_NAME)), collapse = ", "),
    call. = FALSE
  )
}

mutated <- clean %>%
  dplyr::filter(
    GENE_NAME %in% TARGET_GENES,
    !is.na(HGVS_PROTEIN),
    stringr::str_trim(HGVS_PROTEIN) != "",
    stringr::str_trim(HGVS_PROTEIN) != "."
  ) %>%
  dplyr::distinct(SAMPLE, year, GENE_NAME) %>%
  dplyr::count(year, GENE_NAME, name = "n_mutated")


# The full grid, so a zero stays a plotted zero.
prevalence <- tidyr::crossing(
  year      = sort(unique(clean$year)),
  GENE_NAME = TARGET_GENES
) %>%
  dplyr::left_join(mutated,     by = c("year", "GENE_NAME")) %>%
  dplyr::left_join(year_totals, by = "year") %>%
  dplyr::mutate(
    n_mutated  = tidyr::replace_na(n_mutated, 0L),
    prevalence = 100 * n_mutated / n_total,
    GENE_NAME  = factor(GENE_NAME, levels = TARGET_GENES)
  ) %>%
  dplyr::arrange(GENE_NAME, year)

if (any(is.na(prevalence$n_total))) {
  stop("A year has no denominator; the join in section 5 failed.",
       call. = FALSE)
}

add_to_report("")
add_to_report(strrep("-", 72))
add_to_report("PREVALENCE (mutated samples / samples sequenced that year)")
add_to_report(strrep("-", 72))

for (gene in TARGET_GENES) {

  gene_rows <- prevalence %>% dplyr::filter(GENE_NAME == gene)

  add_to_report("")
  add_to_report("  ", gene)

  for (index in seq_len(nrow(gene_rows))) {
    add_to_report(sprintf(
      "    %d   %5d/%-5d  %6.1f%%",
      gene_rows$year[[index]],
      gene_rows$n_mutated[[index]],
      gene_rows$n_total[[index]],
      gene_rows$prevalence[[index]]
    ))
  }
}


# ==============================================================================
# 6. THE FIGURE
# ==============================================================================

year_axis <- year_totals %>%
  dplyr::mutate(label = sprintf("%d\n(N = %d)", year, n_total))

axis_label_lookup <- setNames(year_axis$label, year_axis$year)

figure <- ggplot2::ggplot(
  prevalence,
  ggplot2::aes(
    x = year, y = prevalence,
    group = GENE_NAME, colour = GENE_NAME,
    shape = GENE_NAME, fill = GENE_NAME
  )
) +

  # Period bands, drawn first so everything else sits on top.
  ggplot2::geom_rect(
    data = PERIOD_BANDS,
    mapping = ggplot2::aes(
      xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf
    ),
    fill = PERIOD_BANDS$fill,
    alpha = 0.80,
    inherit.aes = FALSE
  ) +

  ggplot2::geom_hline(
    yintercept = seq(0, 100, by = 20),
    colour = "grey72", linewidth = 0.55
  ) +

  ggplot2::geom_vline(
    xintercept = PERIOD_DIVIDER,
    linetype = "dashed", colour = "#454545", linewidth = 0.8
  ) +

  ggplot2::geom_line(linewidth = 0.9) +

  ggplot2::geom_point(size = 4.8, stroke = 1.2, alpha = 0.75) +

  ggplot2::geom_text(
    data = PERIOD_BANDS,
    mapping = ggplot2::aes(x = label_x, y = 115, label = label),
    fontface = "bold", size = 7, colour = "#444444",
    inherit.aes = FALSE
  ) +

  ggplot2::scale_colour_manual(
    values = GENE_COLOURS, breaks = TARGET_GENES, name = NULL
  ) +

  ggplot2::scale_fill_manual(
    values = GENE_COLOURS, breaks = TARGET_GENES, name = NULL
  ) +

  ggplot2::scale_shape_manual(
    values = GENE_SHAPES, breaks = TARGET_GENES, name = NULL
  ) +

  ggplot2::scale_x_continuous(
    breaks = year_axis$year,
    labels = function(breaks) axis_label_lookup[as.character(breaks)],
    limits = X_LIMITS,
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +

  ggplot2::scale_y_continuous(
    breaks = seq(0, 100, by = 20),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +

  ggplot2::coord_cartesian(ylim = c(0, 121), clip = "off") +

  ggplot2::labs(x = NULL, y = "Mutated samples (%)") +

  ggplot2::guides(
    colour = ggplot2::guide_legend(
      override.aes = list(size = 5, linewidth = 1)
    ),
    fill  = "none",
    shape = "none"
  ) +

  ggplot2::theme_classic(base_size = 16) +

  ggplot2::theme(
    axis.title.y = ggplot2::element_text(
      size = 22, face = "bold", margin = ggplot2::margin(r = 15)
    ),
    axis.text.y  = ggplot2::element_text(size = 17, colour = "#405785"),
    axis.text.x  = ggplot2::element_text(
      size = 16, colour = "#222222", lineheight = 0.85,
      margin = ggplot2::margin(t = 10)
    ),
    axis.ticks.x = ggplot2::element_line(colour = "#405785", linewidth = 1),
    axis.ticks.y = ggplot2::element_blank(),
    axis.ticks.length.x = ggplot2::unit(0.25, "cm"),
    axis.line.x  = ggplot2::element_line(colour = "#405785", linewidth = 1.2),
    axis.line.y  = ggplot2::element_blank(),
    panel.grid   = ggplot2::element_blank(),
    legend.position  = "right",
    legend.text      = ggplot2::element_text(size = 18, face = "bold"),
    legend.key.width  = ggplot2::unit(1.2, "cm"),
    legend.key.height = ggplot2::unit(1.2, "cm"),
    plot.margin = ggplot2::margin(t = 25, r = 30, b = 15, l = 15)
  )

for (extension in c("svg", "png")) {

  ggplot2::ggsave(
    filename = file.path(
      FIGURE_DIR, paste0("figure02b_gene_mutation_trends.", extension)
    ),
    plot = figure,
    width = 12, height = 8, units = "in",
    dpi = 600, bg = "white"
  )
}


# ==============================================================================
# 7. TABLES
# ==============================================================================

readr::write_tsv(
  prevalence %>%
    dplyr::mutate(GENE_NAME = as.character(GENE_NAME)) %>%
    dplyr::select(
      gene = GENE_NAME, year, period, n_mutated, n_total, prevalence
    ),
  file.path(TABLE_DIR, "figure02b_gene_prevalence_by_year.tsv")
)

add_to_report("")
add_to_report(strrep("=", 72))
add_to_report("WRITTEN")
add_to_report(strrep("=", 72))
add_to_report("  figures/figure02/figure02b_gene_mutation_trends.svg / .png")
add_to_report("  tables/figure02/figure02b_gene_prevalence_by_year.tsv")
add_to_report("  tables/figure02/figure02b_samples_per_year.tsv")

writeLines(report_lines, file.path(TABLE_DIR, "figure02b_report.txt"))

capture.output(
  sessionInfo(),
  file = file.path(TABLE_DIR, "figure02b_sessionInfo.txt")
)
