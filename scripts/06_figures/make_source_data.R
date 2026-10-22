# Source Data workbook - one sheet per published panel
#
# Inputs:  the tables under results/tables/, written by the figure scripts
# Outputs: results/Source_Data.xlsx
#          results/source_data_report.txt
#
# Usage:   Rscript scripts/06_figures/make_source_data.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("openxlsx is needed: install.packages(\"openxlsx\")", call. = FALSE)
}

PROJECT_ROOT <- normalizePath(".", mustWork = FALSE)
source(file.path(PROJECT_ROOT, "config.R"))

OUTPUT_FILE <- file.path(RESULTS_DIR, "Source_Data.xlsx")
REPORT_FILE <- file.path(RESULTS_DIR, "source_data_report.txt")

# ---- what goes on each sheet ----
# One row per sheet. Several files on one sheet are stacked, each under its
# own sub-heading. Edit the sheet names to match the panel letters in the
# final manuscript.

SHEETS <- tibble::tribble(
  ~sheet,          ~caption,                                                        ~files,

  "Fig 1",         "Figure 1. Samples sequenced per district, and the WHO 2025 global malaria ranks annotated on the map.",
                   "figure01_samples_by_district.tsv,who_malaria_rank_2025.tsv,figure01_sample_level.tsv",

  "Fig 2a",        "Figure 2a. Sequencing coverage across amplicon targets.",
                   "figure02a_source_data.tsv,figure02a_summary_by_gene.tsv",

  "Fig 2b",        "Figure 2b. Per-gene mutation prevalence by year.",
                   "figure02b_gene_prevalence_by_year.tsv,figure02b_samples_per_year.tsv",

  "Fig 2c",        "Figure 2c. Composite antimalarial resistance backgrounds.",
                   "figure02c_category_prevalence.tsv,figure02c_sample_level_categories.tsv",

  "Fig 3a",        "Figure 3a. K13 mutation frequencies along the protein.",
                   "figure03a_mutation_frequencies.tsv,figure03a_domains.tsv",

  "Fig 3b",        "Figure 3b. K13 mutation categories across four sample groups.",
                   "figure03b_category_prevalence.tsv,figure03b_sample_level_categories.tsv",

  "Fig 3d",        "Figure 3d. K13 genotypes by district.",
                   "figure03d_district_genotypes.tsv",

  "Fig 4a-c",      "Figure 4a-c. Principal coordinates of IBS distance. Panels a-c plot the same coordinates coloured by period, region and kelch13 genotype respectively; the axis pair shown is given in the figure legend.",
                   "pcoa_coordinates.tsv,pcoa_eigenvalues.tsv",

  "Fig 4 PERMANOVA", "Figure 4a-c. PERMANOVA results.",
                   "permanova_results.tsv",

  "Fig 4d",        "Figure 4d. Neighbour-joining tree of the Kigali samples. The tree itself is provided as a Newick file in the code repository.",
                   "figure04d_tip_metadata.tsv",

  "Fig 5a-c",      "Figure 5a-c. Selection scan. One row per variant, carrying the FST plotted in panel a, the absolute allele-frequency change in panel b and the signed change in panel c.",
                   "figure05_plotting_table.tsv,empirical_thresholds.tsv,selection_statistics_all_variants.tsv",

  "Supp Fig 1a",   "Supplementary Figure 1a. Baseline samples by district.",
                   "figureS1a_samples_by_district.tsv",

  "Supp Fig 1b",   "Supplementary Figure 1b. Monthly malaria cases, 2023 and 2024.",
                   "figureS1b_monthly_cases_long.tsv,figureS1b_monthly_change.tsv,figureS1b_matched_month_totals.tsv",

  "Supp Fig 3",    "Supplementary Figure 3. K13 category prevalence by group.",
                   "figureS3_category_prevalence.tsv,figureS3_group_sizes.tsv,figureS3_sample_level_categories.tsv",

  "Supp Fig 3 stats", "Supplementary Figure 3. Binomial tests, Fisher exact tests and binomial generalised linear models.",
                   "figureS3_binomial_tests.tsv,figureS3_fisher_tests.tsv,figureS3_logistic_models.tsv",

  "Supp Fig 4",    "Supplementary Figure 4. The same principal coordinates as Figure 4, shown as all pairwise comparisons of the first four axes.",
                   "pcoa_coordinates.tsv,pcoa_eigenvalues.tsv"
)

report <- character()
add_to_report <- function(...) {
  line <- paste0(...)
  report <<- c(report, line)
  message(line)
}

# Searched across results/, not just results/tables/: the Figure 5 panels
# script writes its tables beside the figures.
find_table <- function(name) {
  hits <- list.files(RESULTS_DIR, pattern = paste0("^", name, "$"),
                     recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0) return(NA_character_)
  hits[[1]]
}

SUBCAPTIONS <- c(
  figure01_samples_by_district.tsv       = "Samples per district",
  who_malaria_rank_2025.tsv              = "WHO 2025 global malaria rank",
  figure02a_source_data.tsv              = "Per-sample coverage by gene",
  figure02a_summary_by_gene.tsv          = "Coverage summary by gene",
  figure02b_gene_prevalence_by_year.tsv  = "Mutation prevalence by gene and year",
  figure02b_samples_per_year.tsv         = "Samples per year",
  figure03a_mutation_frequencies.tsv     = "Mutation frequencies",
  figure03a_domains.tsv                  = "Protein domain boundaries",
  pcoa_coordinates.tsv                   = "Principal coordinates",
  pcoa_eigenvalues.tsv                   = "Eigenvalues and variance explained",
  selection_statistics_all_variants.tsv  = "All variants, annotated",
  figure05_plotting_table.tsv            = "Values plotted in panels a to c",
  empirical_thresholds.tsv               = "Empirical 95th-percentile thresholds",
  figureS1b_monthly_cases_long.tsv       = "Monthly cases",
  figureS1b_monthly_change.tsv           = "Change between years",
  figureS1b_matched_month_totals.tsv     = "Matched-month totals",
  figureS3_category_prevalence.tsv       = "Category prevalence by group",
  figureS3_group_sizes.tsv               = "Group sizes",
  figureS3_binomial_tests.tsv            = "Exact binomial tests",
  figureS3_fisher_tests.tsv              = "Fisher exact tests",
  figureS3_logistic_models.tsv           = "Binomial generalised linear models",
  figure01_sample_level.tsv              = "Per-sample detail",
  figure02c_sample_level_categories.tsv  = "Per-sample detail",
  figure03b_sample_level_categories.tsv  = "Per-sample detail",
  figureS3_sample_level_categories.tsv   = "Per-sample detail",
  figure02c_category_prevalence.tsv      = "Category prevalence",
  figure03b_category_prevalence.tsv      = "Category prevalence by group"
)

subcaption <- function(name) {
  if (name %in% names(SUBCAPTIONS)) SUBCAPTIONS[[name]] else name
}

# openxlsx declares a drawing object on every sheet even when none exists,
# which makes strict readers reject the file.
remove_dangling_drawings <- function(path) {
  work <- file.path(tempdir(), "source_data_repack")
  unlink(work, recursive = TRUE)
  dir.create(work, recursive = TRUE)
  utils::unzip(path, exdir = work)

  rels <- list.files(file.path(work, "xl", "worksheets", "_rels"),
                     full.names = TRUE)
  for (file in rels) {
    text <- paste(readLines(file, warn = FALSE), collapse = "")
    text <- gsub("<Relationship[^>]*[Dd]rawing[^>]*/>", "", text)
    writeLines(text, file)
  }

  rebuilt <- file.path(tempdir(), "Source_Data_rebuilt.xlsx")
  if (file.exists(rebuilt)) file.remove(rebuilt)
  previous <- setwd(work)
  status <- try(utils::zip(rebuilt,
                           list.files(".", recursive = TRUE, all.files = TRUE),
                           flags = "-qr9X"), silent = TRUE)
  setwd(previous)

  if (!inherits(status, "try-error") && file.exists(rebuilt)) {
    file.copy(rebuilt, path, overwrite = TRUE)
    TRUE
  } else {
    FALSE
  }
}

# ---- build ----

workbook <- openxlsx::createWorkbook()
heading <- openxlsx::createStyle(textDecoration = "bold")
used <- character()
missing <- character()

for (index in seq_len(nrow(SHEETS))) {

  sheet   <- SHEETS$sheet[[index]]
  caption <- SHEETS$caption[[index]]
  wanted  <- trimws(strsplit(SHEETS$files[[index]], ",")[[1]])

  found <- vapply(wanted, find_table, character(1))
  absent <- wanted[is.na(found)]
  if (length(absent) > 0) {
    missing <- c(missing, paste0(sheet, ": ", absent))
  }
  found <- found[!is.na(found)]
  if (length(found) == 0) {
    add_to_report("skipped ", sheet, " - no table found")
    next
  }

  openxlsx::addWorksheet(workbook, sheet)
  openxlsx::writeData(workbook, sheet, caption, startRow = 1)
  openxlsx::addStyle(workbook, sheet, heading, rows = 1, cols = 1)

  row <- 3
  for (path in found) {
    table <- readr::read_tsv(path, show_col_types = FALSE,
                             progress = FALSE)
    if (length(found) > 1) {
      openxlsx::writeData(workbook, sheet, subcaption(basename(path)),
                          startRow = row)
      openxlsx::addStyle(workbook, sheet, heading, rows = row, cols = 1)
      row <- row + 1
    }
    openxlsx::writeData(workbook, sheet, table, startRow = row,
                        headerStyle = heading)
    row <- row + nrow(table) + 3
    used <- c(used, normalizePath(path))
  }

  openxlsx::setColWidths(workbook, sheet, cols = 1:12, widths = "auto")
  add_to_report("wrote ", sheet, "  (", length(found), " table(s))")
}

openxlsx::saveWorkbook(workbook, OUTPUT_FILE, overwrite = TRUE)

if (!remove_dangling_drawings(OUTPUT_FILE)) {
  add_to_report("")
  add_to_report("Could not repack the workbook; it may prompt Excel to repair it.")
}

# ---- what was left out ----

all_tables <- list.files(RESULTS_DIR, pattern = "\\.tsv$",
                         recursive = TRUE, full.names = TRUE)
unused <- setdiff(normalizePath(all_tables), unique(used))

add_to_report("")
add_to_report(strrep("-", 72))
add_to_report("Sheets written: ", length(openxlsx::sheets(workbook)))
add_to_report("Workbook: ", OUTPUT_FILE)

if (length(missing) > 0) {
  add_to_report("")
  add_to_report("Listed but not found - the figure script may not have run:")
  for (line in missing) add_to_report("  ", line)
}

if (length(unused) > 0) {
  add_to_report("")
  add_to_report("Present under results/ but on no sheet:")
  for (path in unused) {
    add_to_report("  ", sub(paste0(RESULTS_DIR, "/"), "", path, fixed = TRUE))
  }
  add_to_report("")
  add_to_report("Check whether any of these underlies a published graph.")
}

writeLines(report, REPORT_FILE)
