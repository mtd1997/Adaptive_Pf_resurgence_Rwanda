################################################################################
# PREPARE THE SAMPLE METADATA
################################################################################
#
# WHAT IT DOES
#   Takes the working metadata spreadsheet and produces the two metadata
#   files the rest of the pipeline uses:
#
#     data/metadata/samples_metadata.tsv
#         De-identified. Committed to git and published with the repository.
#         Every downstream script reads this one.
#
#     data/metadata/samples_metadata_internal.tsv
#         Everything, including exact age and health centre. Stays on your
#         machine. Listed in .gitignore so it cannot be committed by accident.
#
#   It also writes a validation report so you can see exactly what changed.
#
# WHY TWO FILES
#   Exact age, district, sampling year and sex together can single out one
#   patient in a small district, even with the name removed. That is a
#   quasi-identifier, and publishing it is a disclosure risk rather than a
#   theoretical one.
#
#   Age is therefore banded in the public file and the exact value removed.
#   Health centre is removed as well: health centre plus sampling date plus
#   age narrows a patient down even further than district does.
#
#   Nothing in the manuscript needs either at full resolution. The age bands
#   reproduce the demographic table, and the published analyses stratify by
#   district and period, not by health centre.
#
#   If a reviewer asks for the unbanded data, it is available on request
#   under a data-sharing agreement. That is the standard arrangement and is
#   what the Data Availability statement should say.
#
# WHAT IT DOES NOT DO
#   It does not drop samples, change a period or region label, or alter any
#   value other than those listed above. Sample filtering happened upstream
#   and is not revisited here.
#
# HOW TO RUN
#   From the top of the repository:
#
#     Rscript scripts/00_metadata/prepare_metadata.R --input <your file.tsv>
#
#   Add --dry-run to see the report without writing anything. Do that first.
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
# Command-line arguments
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

INPUT_FILE <- get_argument("--input")

DRY_RUN <- "--dry-run" %in% arguments

if (is.null(INPUT_FILE)) {
  stop(
    "No input file given.\n\n",
    "Usage:\n",
    "  Rscript scripts/00_metadata/prepare_metadata.R \\\n",
    "      --input /path/to/metadata_05-08-2026.tsv [--dry-run]\n",
    call. = FALSE
  )
}

if (!file.exists(INPUT_FILE)) {
  stop("Input file not found: ", INPUT_FILE, call. = FALSE)
}

OUTPUT_DIR <- file.path(PROJECT_ROOT, "data", "metadata")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

PUBLIC_FILE   <- file.path(OUTPUT_DIR, "samples_metadata.tsv")
INTERNAL_FILE <- file.path(OUTPUT_DIR, "samples_metadata_internal.tsv")
REPORT_FILE   <- file.path(OUTPUT_DIR, "metadata_preparation_report.txt")


# ==============================================================================
# 2. READ THE INPUT
# ==============================================================================
#
# Read everything as text first. Letting readr guess column types has caused
# trouble before: a sample ID like 0012 loses its leading zero if it is read
# as a number, and that silently breaks every join against the VCF.
# ==============================================================================

raw_metadata <- readr::read_tsv(
  INPUT_FILE,
  col_types = readr::cols(.default = readr::col_character()),
  show_col_types = FALSE
)

message(sprintf(
  "Read %d rows and %d columns from:\n  %s\n",
  nrow(raw_metadata), ncol(raw_metadata), INPUT_FILE
))

message("Columns found:")
message(paste0("  ", paste(names(raw_metadata), collapse = ", ")), "\n")


# ==============================================================================
# 3. MAP THE COLUMNS
# ==============================================================================
#
# Metadata files in this project have used different capitalisation and
# spellings over time. Rather than hard-code one spelling, each piece of
# information is looked up by a list of names it has gone by.
#
# Matching ignores case, spaces, hyphens and underscores, so "Health_Center",
# "health center" and "HealthCentre" all resolve to the same column.
#
# If a column in your file is not found, add its name to the relevant list.
# ==============================================================================

COLUMN_ALIASES <- list(
  sample        = c("SAMPLE", "sample_id", "sampleid", "ID", "Sample_Name"),
  period        = c("PERIOD", "period", "sampling_period", "Period"),
  region        = c("REGION", "region", "Region"),
  district      = c("DISTRICT", "district", "District"),
  health_centre = c("Health_Center", "Health_Centre", "health_center",
                    "HealthCenter", "HC", "site", "Site"),
  year          = c("YEAR", "year", "Year", "collection_year",
                    "sampling_year"),
  age           = c("AGE", "age", "Age", "age_years"),
  sex           = c("GENDER", "gender", "Gender", "SEX", "sex", "Sex"),
  k13_status    = c("K13_status", "k13_status", "K13_STATUS", "kelch13_status")
)

normalise_name <- function(names_vector) {
  tolower(gsub("[^a-z0-9]", "", tolower(names_vector)))
}

find_column <- function(data, aliases) {

  available_normalised <- normalise_name(names(data))
  wanted_normalised    <- normalise_name(aliases)

  position <- match(wanted_normalised, available_normalised)
  position <- position[!is.na(position)]

  if (length(position) == 0) {
    return(NA_character_)
  }

  names(data)[[position[[1]]]]
}

resolved_columns <- vapply(
  COLUMN_ALIASES,
  function(aliases) find_column(raw_metadata, aliases),
  character(1)
)

message("Column mapping:")
for (field in names(resolved_columns)) {
  message(sprintf(
    "  %-14s -> %s",
    field,
    ifelse(is.na(resolved_columns[[field]]), "NOT FOUND", resolved_columns[[field]])
  ))
}
message("")


# ------------------------------------------------------------------------------
# Only the sample identifier is strictly required. Everything else is
# reported as missing and simply left out, so the script still runs on a
# reduced file rather than refusing.
# ------------------------------------------------------------------------------

if (is.na(resolved_columns[["sample"]])) {
  stop(
    "No sample-identifier column was found. Add the name your file uses to ",
    "COLUMN_ALIASES$sample near the top of this script.",
    call. = FALSE
  )
}

missing_fields <- names(resolved_columns)[is.na(resolved_columns)]

if (length(missing_fields) > 0) {
  warning(
    "These fields were not found and will be absent from the output: ",
    paste(missing_fields, collapse = ", "),
    call. = FALSE
  )
}


# ------------------------------------------------------------------------------
# Rename to the canonical names used from here on.
# ------------------------------------------------------------------------------

present_columns <- resolved_columns[!is.na(resolved_columns)]

metadata <- raw_metadata %>%
  dplyr::select(dplyr::all_of(unname(present_columns)))

names(metadata) <- names(present_columns)


# ==============================================================================
# 4. TIDY WHITESPACE AND BLANKS
# ==============================================================================
#
# Trailing spaces in a sample ID are invisible on screen and break joins, so
# every text field is trimmed, and blank strings become proper NA values.
# ==============================================================================

metadata <- metadata %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::everything(),
      ~ dplyr::na_if(stringr::str_trim(.x), "")
    )
  )


# ==============================================================================
# 5. CHECK THE SAMPLE IDENTIFIERS
# ==============================================================================

report_lines <- character(0)

add_to_report <- function(...) {
  line <- paste0(...)
  report_lines <<- c(report_lines, line)
  message(line)
}

add_to_report(strrep("=", 72))
add_to_report("METADATA PREPARATION REPORT")
add_to_report(strrep("=", 72))
add_to_report("Input:     ", INPUT_FILE)
add_to_report("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
add_to_report("")

add_to_report("Rows read: ", nrow(metadata))

blank_samples <- sum(is.na(metadata$sample))

if (blank_samples > 0) {
  stop(
    blank_samples, " rows have no sample identifier. Fix these in the source ",
    "spreadsheet before continuing; dropping them here would hide a problem.",
    call. = FALSE
  )
}

duplicate_samples <- metadata %>%
  dplyr::count(sample, name = "n") %>%
  dplyr::filter(n > 1)

if (nrow(duplicate_samples) > 0) {

  add_to_report("")
  add_to_report("DUPLICATED SAMPLE IDENTIFIERS: ", nrow(duplicate_samples))

  for (index in seq_len(min(10, nrow(duplicate_samples)))) {
    add_to_report(sprintf(
      "  %s appears %d times",
      duplicate_samples$sample[[index]],
      duplicate_samples$n[[index]]
    ))
  }

  stop(
    nrow(duplicate_samples), " sample identifiers appear more than once. ",
    "Every downstream join assumes one row per sample. Resolve these in the ",
    "source spreadsheet.",
    call. = FALSE
  )
}

add_to_report("Unique sample identifiers: ", dplyr::n_distinct(metadata$sample))


# ==============================================================================
# 6. STANDARDISE THE PERIOD AND REGION LABELS
# ==============================================================================
#
# The manuscript uses "Baseline-Kigali" and "Resurgence". Older files use
# "Pre-Resurgence" and "Post-Resurgence". Both are accepted. The original
# value is kept in period_original so nothing is rewritten silently.
# ==============================================================================

if ("period" %in% names(metadata)) {

  metadata <- metadata %>%
    dplyr::mutate(
      period_original = period,
      period = dplyr::recode(
        period,
        "Pre-Resurgence"  = "Baseline-Kigali",
        "Post-Resurgence" = "Resurgence"
      )
    )

  add_to_report("")
  add_to_report("Samples per period:")

  period_counts <- dplyr::count(metadata, period, name = "n")

  for (index in seq_len(nrow(period_counts))) {
    add_to_report(sprintf(
      "  %-20s %d",
      period_counts$period[[index]],
      period_counts$n[[index]]
    ))
  }

  relabelled <- sum(
    metadata$period != metadata$period_original,
    na.rm = TRUE
  )

  if (relabelled > 0) {
    add_to_report(sprintf(
      "  (%d rows relabelled to the published wording; the original value is ",
      relabelled
    ))
    add_to_report("   kept in the period_original column)")
  }
}

if ("region" %in% names(metadata)) {

  add_to_report("")
  add_to_report("Samples per region:")

  region_counts <- dplyr::count(metadata, region, name = "n")

  for (index in seq_len(nrow(region_counts))) {
    add_to_report(sprintf(
      "  %-20s %d",
      region_counts$region[[index]],
      region_counts$n[[index]]
    ))
  }
}


# ==============================================================================
# 7. BAND THE AGE
# ==============================================================================
#
# Bands follow the groups used in malaria epidemiology and reproduce the
# demographic table in the supplement:
#
#   <5    5-14    15-24    25-34    35-49    50+
#
# Ages that are not numbers, or are outside 0-120, become NA and are counted
# in the report rather than quietly bundled into a band.
# ==============================================================================

AGE_BREAKS <- c(0, 5, 15, 25, 35, 50, Inf)
AGE_LABELS <- c("<5", "5-14", "15-24", "25-34", "35-49", "50+")

if ("age" %in% names(metadata)) {

  metadata <- metadata %>%
    dplyr::mutate(
      age_numeric = suppressWarnings(as.numeric(age)),

      age_numeric = dplyr::if_else(
        !is.na(age_numeric) & age_numeric >= 0 & age_numeric <= 120,
        age_numeric,
        NA_real_
      ),

      age_band = cut(
        age_numeric,
        breaks = AGE_BREAKS,
        labels = AGE_LABELS,
        right  = FALSE,
        include.lowest = TRUE
      ),

      age_band = as.character(age_band)
    )

  unparseable_ages <- metadata %>%
    dplyr::filter(!is.na(age), is.na(age_numeric))

  add_to_report("")
  add_to_report("Age banding:")

  band_counts <- metadata %>%
    dplyr::count(age_band, name = "n") %>%
    dplyr::arrange(match(age_band, AGE_LABELS))

  for (index in seq_len(nrow(band_counts))) {
    add_to_report(sprintf(
      "  %-10s %d",
      ifelse(is.na(band_counts$age_band[[index]]), "missing",
             band_counts$age_band[[index]]),
      band_counts$n[[index]]
    ))
  }

  if (nrow(unparseable_ages) > 0) {

    add_to_report("")
    add_to_report(sprintf(
      "  %d age value(s) could not be read as a number and became NA:",
      nrow(unparseable_ages)
    ))

    distinct_bad <- unique(unparseable_ages$age)

    for (value in utils::head(distinct_bad, 10)) {
      add_to_report(sprintf("    \"%s\"", value))
    }

    add_to_report("  Check these in the source spreadsheet.")
  }
}


# ==============================================================================
# 8. BUILD THE TWO OUTPUT TABLES
# ==============================================================================
#
# The internal table is everything. The public table is the internal one with
# the two identifying columns removed.
#
# Columns are listed explicitly rather than excluded by name, so a new column
# added to the source spreadsheet later cannot leak into the public file
# without someone deciding to add it here.
# ==============================================================================

INTERNAL_COLUMNS <- c(
  "sample", "period", "period_original", "region", "district",
  "health_centre", "year", "age", "age_band", "sex", "k13_status"
)

PUBLIC_COLUMNS <- c(
  "sample", "period", "region", "district",
  "year", "age_band", "sex", "k13_status"
)

internal_metadata <- metadata %>%
  dplyr::select(dplyr::any_of(INTERNAL_COLUMNS))

public_metadata <- metadata %>%
  dplyr::select(dplyr::any_of(PUBLIC_COLUMNS))


# ------------------------------------------------------------------------------
# Rename `sample` to SAMPLE on output, because every downstream script and
# the VCF headers use that spelling.
# ------------------------------------------------------------------------------

names(internal_metadata)[names(internal_metadata) == "sample"] <- "SAMPLE"
names(public_metadata)[names(public_metadata) == "sample"]     <- "SAMPLE"


# ------------------------------------------------------------------------------
# A direct check that nothing identifying survived into the public table.
# This is the safety net: if someone edits PUBLIC_COLUMNS carelessly later,
# the script stops rather than writing the file.
# ------------------------------------------------------------------------------

FORBIDDEN_IN_PUBLIC <- c("age", "health_centre", "period_original")

leaked_columns <- intersect(FORBIDDEN_IN_PUBLIC, names(public_metadata))

if (length(leaked_columns) > 0) {
  stop(
    "These columns must not appear in the public metadata file but do: ",
    paste(leaked_columns, collapse = ", "),
    "\nCheck PUBLIC_COLUMNS in section 8.",
    call. = FALSE
  )
}

add_to_report("")
add_to_report(strrep("-", 72))
add_to_report("OUTPUT")
add_to_report(strrep("-", 72))
add_to_report("")
add_to_report("Public file (committed to git, published with the repository):")
add_to_report("  ", basename(PUBLIC_FILE))
add_to_report("  columns: ", paste(names(public_metadata), collapse = ", "))
add_to_report("")
add_to_report("Internal file (gitignored, stays on this machine):")
add_to_report("  ", basename(INTERNAL_FILE))
add_to_report("  columns: ", paste(names(internal_metadata), collapse = ", "))
add_to_report("")
add_to_report("Removed from the public file:")
add_to_report("  age            exact value, replaced by age_band")
add_to_report("  health_centre  identifying in combination with date and age")
add_to_report("  period_original  working label, superseded by period")


# ==============================================================================
# 9. WRITE, OR REPORT ONLY
# ==============================================================================

add_to_report("")
add_to_report(strrep("=", 72))

if (DRY_RUN) {

  add_to_report("DRY RUN — no files were written.")
  add_to_report("Re-run without --dry-run to write them.")
  add_to_report(strrep("=", 72))

  writeLines(report_lines, REPORT_FILE)

  message("\nReport written to:\n  ", REPORT_FILE)

} else {

  readr::write_tsv(public_metadata, PUBLIC_FILE, na = "")
  readr::write_tsv(internal_metadata, INTERNAL_FILE, na = "")

  add_to_report("Files written.")
  add_to_report(strrep("=", 72))

  writeLines(report_lines, REPORT_FILE)

  message("\nWritten:")
  message("  ", PUBLIC_FILE)
  message("  ", INTERNAL_FILE)
  message("  ", REPORT_FILE)

  message(
    "\nBefore committing, confirm the internal file is ignored:\n",
    "  git check-ignore -v data/metadata/samples_metadata_internal.tsv"
  )
}
