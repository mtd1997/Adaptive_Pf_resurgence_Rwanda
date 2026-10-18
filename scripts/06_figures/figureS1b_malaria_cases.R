################################################################################
# SUPPLEMENTARY FIGURE S1, PANEL B — MONTHLY MALARIA CASES, 2023 VS 2024
################################################################################
#
# WHAT IT DOES
#   Reads the reported monthly malaria case counts for four Kigali health
#   facilities, tidies them, and draws the 2023 versus 2024 comparison.
#
#   It also writes the comparison tables: monthly change, and matched-month
#   annual totals.
#
# -----------------------------------------------------------------------------
# THE BUG THIS FIXES
# -----------------------------------------------------------------------------
#   The original script was in two halves that disagreed about facility names.
#
#   The first half wrote them short:
#       "KFH", "COR-UNUM HC", "Gikondo HC", "Remera"
#
#   The second half then set the factor levels to the long forms:
#       "King Faisal Hospital", "COR-UNUM Health Center", ...
#
#   Setting factor levels that do not match the data turns every value into
#   NA. facet_wrap() would then draw a single panel labelled NA instead of
#   four facilities, and no error would be raised.
#
#   The published figure has four facets, so the file on disk had been written
#   by an earlier version that used the long names. But running the script end
#   to end would have broken it.
#
#   The fix here is not to correct one half to match the other. It is to
#   define the names ONCE, in FACILITY_NAMES below, and use that single
#   definition for both the tidying and the plotting. They cannot drift apart
#   again.
# -----------------------------------------------------------------------------
#
# INPUT
#   An Excel workbook of monthly counts, laid out as:
#
#     column A  month
#     columns B,C   King Faisal Hospital     2024, 2023
#     columns D,E   COR-UNUM Health Center   2024, 2023
#     columns F,G   Gikondo Health Center    2024, 2023
#     columns H,I   Remera Health Center     2024, 2023
#
#   with the twelve monthly rows starting at row 4. That layout is asserted
#   and checked below rather than assumed, because a shifted row would
#   otherwise be read as data.
#
#   Set MALARIA_CASES_FILE in config.R, or pass --input.
#
# OUTPUTS
#   results/figures/figureS1/figureS1b_malaria_cases.svg
#   results/figures/panels/figureS1b_malaria_cases.rds   for the assembly
#   results/tables/figureS1/figureS1b_monthly_cases_long.tsv
#   results/tables/figureS1/figureS1b_monthly_change.tsv
#   results/tables/figureS1/figureS1b_matched_month_totals.tsv
#   results/tables/figureS1/figureS1b_data_accounting.tsv
#
# A MISSING MONTH
#   COR-UNUM Health Center has no December 2024 figure in the source
#   workbook. It is kept as NA and never converted to zero, because zero
#   would mean "no cases" rather than "not reported".
#
#   The annual comparison therefore uses only months present in BOTH years
#   for each facility, so a facility with eleven months is not compared
#   against another facility's twelve. The number of months compared is
#   reported per facility.
#
# HOW TO RUN
#   From the top of the repository:
#
#     Rscript scripts/06_figures/figureS1b_malaria_cases.R
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
  library(readxl)
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

INPUT_FILE <- get_argument(
  "--input",
  if (exists("MALARIA_CASES_FILE", inherits = TRUE)) {
    get("MALARIA_CASES_FILE", inherits = TRUE)
  } else {
    NULL
  }
)

SHEET_NAME <- get_argument("--sheet", "Feuil1")

if (is.null(INPUT_FILE)) {
  stop(
    "No input workbook given.\n\n",
    "Either set MALARIA_CASES_FILE in config.R, or pass:\n",
    "  Rscript scripts/06_figures/figureS1b_malaria_cases.R \\\n",
    "      --input \"/path/to/Malaria cases Kigali 2023vs2024.xlsx\"\n",
    call. = FALSE
  )
}

if (!file.exists(INPUT_FILE)) {
  stop("Input workbook not found:\n  ", INPUT_FILE, call. = FALSE)
}

FIGURES1_DIR <- file.path(PROJECT_ROOT, "results", "figures", "figureS1")
PANEL_DIR    <- file.path(PROJECT_ROOT, "results", "figures", "panels")
TABLE_DIR    <- file.path(PROJECT_ROOT, "results", "tables",  "figureS1")

for (directory in c(FIGURES1_DIR, PANEL_DIR, TABLE_DIR)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}


# ==============================================================================
# 2. THE SINGLE DEFINITION OF THE FACILITIES
# ==============================================================================
#
# This is the one place facility names are written down. The workbook column
# positions, the tidy table, the facet order and the factor levels all come
# from here, so they cannot disagree.
#
# `column_2024` and `column_2023` are 1-based column positions in the sheet.
# ==============================================================================

FACILITY_NAMES <- tibble::tribble(
  ~facility,                  ~column_2024, ~column_2023,
  "King Faisal Hospital",     2L,           3L,
  "COR-UNUM Health Center",   4L,           5L,
  "Gikondo Health Center",    6L,           7L,
  "Remera Health Center",     8L,           9L
)

MONTH_COLUMN   <- 1L
FIRST_DATA_ROW <- 4L
MONTH_COUNT    <- 12L

YEAR_COLOURS <- c(
  "2023" = "#0072B2",
  "2024" = "#D55E00"
)


# ==============================================================================
# 3. READ THE WORKBOOK
# ==============================================================================

raw <- readxl::read_excel(
  INPUT_FILE,
  sheet = SHEET_NAME,
  col_names = FALSE,
  .name_repair = "minimal"
)

message(sprintf(
  "Read %d rows and %d columns from:\n  %s\n  sheet: %s\n",
  nrow(raw), ncol(raw), INPUT_FILE, SHEET_NAME
))

required_columns <- max(c(MONTH_COLUMN, FACILITY_NAMES$column_2024,
                          FACILITY_NAMES$column_2023))

if (ncol(raw) < required_columns) {
  stop(
    "The sheet has ", ncol(raw), " column(s) but the layout needs at least ",
    required_columns, ". Check the sheet name and the layout described at ",
    "the top of this script.",
    call. = FALSE
  )
}

last_data_row <- FIRST_DATA_ROW + MONTH_COUNT - 1L

if (nrow(raw) < last_data_row) {
  stop(
    "The sheet has ", nrow(raw), " row(s) but the twelve monthly rows are ",
    "expected at rows ", FIRST_DATA_ROW, " to ", last_data_row, ".",
    call. = FALSE
  )
}

monthly_rows <- raw[FIRST_DATA_ROW:last_data_row, ]


# ------------------------------------------------------------------------------
# Verify the month column really holds month names before trusting the row
# offset. A shifted sheet would otherwise be read silently as data.
# ------------------------------------------------------------------------------

month_values <- stringr::str_trim(as.character(monthly_rows[[MONTH_COLUMN]]))

# "Febrary" is a typo in the source workbook.
month_values <- dplyr::recode(month_values, "Febrary" = "February")

unrecognised_months <- setdiff(month_values, month.name)

if (length(unrecognised_months) > 0) {
  stop(
    "Rows ", FIRST_DATA_ROW, "-", last_data_row, " of column ", MONTH_COLUMN,
    " should hold the twelve month names, but these were not recognised: ",
    paste(unrecognised_months, collapse = ", "),
    "\n\nThe sheet layout may have changed. Check FIRST_DATA_ROW near the ",
    "top of this script.",
    call. = FALSE
  )
}


# ==============================================================================
# 4. TIDY
# ==============================================================================

read_facility <- function(index) {

  row <- FACILITY_NAMES[index, ]

  tibble::tibble(
    month    = month_values,
    facility = row$facility,
    `2024`   = suppressWarnings(
      as.numeric(as.character(monthly_rows[[row$column_2024]]))
    ),
    `2023`   = suppressWarnings(
      as.numeric(as.character(monthly_rows[[row$column_2023]]))
    )
  ) %>%
    tidyr::pivot_longer(
      cols = c("2023", "2024"),
      names_to = "year",
      values_to = "cases"
    )
}

malaria_long <- dplyr::bind_rows(
  lapply(seq_len(nrow(FACILITY_NAMES)), read_facility)
) %>%
  dplyr::mutate(
    month = factor(month, levels = month.name, ordered = TRUE),
    year  = factor(year, levels = c("2023", "2024")),
    facility = factor(facility, levels = FACILITY_NAMES$facility)
  ) %>%
  dplyr::arrange(facility, year, month)


# ------------------------------------------------------------------------------
# The factor levels come from FACILITY_NAMES, so by construction nothing can
# fail to match. Assert it anyway, because this is the exact failure the
# original script had and a silent recurrence would be hard to spot.
# ------------------------------------------------------------------------------

if (any(is.na(malaria_long$facility))) {
  stop(
    "Some facility values did not match FACILITY_NAMES. This is the bug the ",
    "header describes; the names have drifted apart again.",
    call. = FALSE
  )
}

if (any(is.na(malaria_long$month))) {
  stop("Some month values did not match the month names.", call. = FALSE)
}

message("Facilities read:")
print(table(malaria_long$facility, malaria_long$year))


# ------------------------------------------------------------------------------
# Report the missing observations rather than filling them.
# ------------------------------------------------------------------------------

missing_observations <- dplyr::filter(malaria_long, is.na(cases))

if (nrow(missing_observations) > 0) {

  message(sprintf(
    "\n%d month(s) have no reported figure and are kept as NA:",
    nrow(missing_observations)
  ))

  print(
    as.data.frame(
      dplyr::select(missing_observations, facility, year, month)
    ),
    row.names = FALSE
  )

  message(
    "These are not converted to zero: zero would assert no cases, which is ",
    "a different claim from not reported."
  )
}


# ==============================================================================
# 5. THE FIGURE
# ==============================================================================

malaria_cases_panel <- ggplot2::ggplot(
  malaria_long,
  ggplot2::aes(x = month, y = cases, group = year, colour = year)
) +

  ggplot2::geom_line(linewidth = 1.3, na.rm = TRUE) +

  ggplot2::geom_point(size = 3.5, na.rm = TRUE) +

  ggplot2::facet_wrap(~ facility, ncol = 2, scales = "fixed") +

  ggplot2::scale_colour_manual(values = YEAR_COLOURS) +

  ggplot2::labs(
    x = "Month",
    y = "Number of reported malaria cases",
    colour = "Year"
  ) +

  ggplot2::theme_classic(base_size = 16) +

  ggplot2::theme(
    axis.title      = ggplot2::element_text(size = 17, face = "bold"),
    axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y     = ggplot2::element_text(size = 13),
    strip.text      = ggplot2::element_text(size = 16, face = "bold"),
    legend.title    = ggplot2::element_text(size = 15, face = "bold"),
    legend.text     = ggplot2::element_text(size = 14),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  filename = file.path(FIGURES1_DIR, "figureS1b_malaria_cases.svg"),
  plot = malaria_cases_panel,
  width = 12, height = 8, units = "in",
  bg = "white"
)

saveRDS(
  list(panel = malaria_cases_panel),
  file.path(PANEL_DIR, "figureS1b_malaria_cases.rds")
)


# ==============================================================================
# 6. COMPARISON TABLES
# ==============================================================================

wide_by_year <- malaria_long %>%
  dplyr::mutate(year = as.character(year)) %>%
  tidyr::pivot_wider(
    names_from = year,
    values_from = cases,
    names_prefix = "cases_"
  )

monthly_change <- wide_by_year %>%
  dplyr::mutate(
    absolute_change = cases_2024 - cases_2023,
    percent_change  = 100 * (cases_2024 - cases_2023) / cases_2023,
    fold_change     = cases_2024 / cases_2023
  ) %>%
  dplyr::arrange(facility, month)


# ------------------------------------------------------------------------------
# Matched months only, for the annual totals.
# ------------------------------------------------------------------------------

matched_months <- dplyr::filter(
  wide_by_year,
  !is.na(cases_2023), !is.na(cases_2024)
)

annual_change <- matched_months %>%
  dplyr::group_by(facility) %>%
  dplyr::summarise(
    months_compared = dplyr::n(),
    cases_2023      = sum(cases_2023),
    cases_2024      = sum(cases_2024),
    absolute_change = sum(cases_2024) - sum(cases_2023),
    percent_change  = 100 * (sum(cases_2024) - sum(cases_2023)) / sum(cases_2023),
    fold_change     = sum(cases_2024) / sum(cases_2023),
    .groups = "drop"
  )

data_accounting <- malaria_long %>%
  dplyr::group_by(facility, year) %>%
  dplyr::summarise(
    months_with_data = sum(!is.na(cases)),
    months_missing   = sum(is.na(cases)),
    total_cases      = sum(cases, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_tsv(
  malaria_long,
  file.path(TABLE_DIR, "figureS1b_monthly_cases_long.tsv")
)

readr::write_tsv(
  monthly_change,
  file.path(TABLE_DIR, "figureS1b_monthly_change.tsv")
)

readr::write_tsv(
  annual_change,
  file.path(TABLE_DIR, "figureS1b_matched_month_totals.tsv")
)

readr::write_tsv(
  data_accounting,
  file.path(TABLE_DIR, "figureS1b_data_accounting.tsv")
)

capture.output(
  sessionInfo(),
  file = file.path(TABLE_DIR, "figureS1b_sessionInfo.txt")
)


# ==============================================================================
# 7. SUMMARY
# ==============================================================================

message("\n", strrep("=", 72))
message("SUPPLEMENTARY FIGURE S1, PANEL B WRITTEN")
message(strrep("=", 72))
message(sprintf("  %s", file.path(FIGURES1_DIR, "figureS1b_malaria_cases.svg")))

message("\nMatched-month comparison, 2023 vs 2024:")
print(
  as.data.frame(
    dplyr::mutate(
      annual_change,
      percent_change = round(percent_change, 1),
      fold_change    = round(fold_change, 2)
    )
  ),
  row.names = FALSE
)

message("\nNext step — assemble Supplementary Figure S1:")
message("  Rscript scripts/06_figures/figureS1_assemble.R")
message(strrep("=", 72))
