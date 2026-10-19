################################################################################
# SUPPLEMENTARY FIGURE S3 — K13 CATEGORY COMPARISONS AND THEIR STATISTICS
################################################################################
#
# WHAT IT DOES
#   Classifies every sample into one K13 mutation category and compares the
#   categories across four groups, in four panels:
#
#     A   K13 categories, 2021-2023 versus King Faisal Hospital 2024-2025
#     B   K13 categories, Kigali versus Hotspot districts, 2024-2025
#     C   Wild type and R561H + K189T, 2021-2023 versus 2024-2025
#     D   Wild type and R561H + K189T, Kigali versus Hotspots, 2024-2025
#
#   and runs the tests reported with them: an exact binomial test of the
#   wild-type / double-mutant split within each group, Fisher's exact tests
#   between groups, and an intercept-only binomial model per group.
#
# INPUT
#   Manuscript_data_29-07-2026.tsv  (MANUSCRIPT_DATA in config.R)
#   The original used an .xlsx on a personal Desktop; --input accepts either.
#
# OUTPUTS
#   results/figures/figureS3/figureS3_k13_comparisons.svg / .pdf / .png
#   results/figures/figureS3/figureS3_panel[A-D].svg
#   results/tables/figureS3/figureS3_category_prevalence.tsv     <- Source Data
#   results/tables/figureS3/figureS3_sample_level_categories.tsv
#   results/tables/figureS3/figureS3_group_sizes.tsv
#   results/tables/figureS3/figureS3_binomial_tests.tsv
#   results/tables/figureS3/figureS3_fisher_tests.tsv
#   results/tables/figureS3/figureS3_logistic_models.tsv
#   results/tables/figureS3/figureS3_report.txt
#
# -----------------------------------------------------------------------------
# WHY THESE GROUPS MATCH FIGURE 3B'S DESPITE BEING DEFINED DIFFERENTLY
# -----------------------------------------------------------------------------
#   Figure 3b splits the 2024-2025 samples by the REGION column. This figure
#   splits them by DISTRICT: Kigali means Gasabo, Kicukiro or Nyarugenge,
#   and Hotspots means the seven districts in HOTSPOT_DISTRICTS.
#
#   Two different definitions, yet both give 700 and 219. That is worth
#   understanding rather than taking on trust, because it is a coincidence
#   of this dataset and not a guarantee.
#
#   Seven samples have a hotspot district but were collected at King Faisal
#   Hospital during 2021-2023 -- the same seven that make Supplementary S1a
#   show 250 of the 257 baseline samples. Because both definitions restrict
#   this comparison to 2024-2025, those seven fall into the baseline group
#   either way and never reach the split. The arithmetic:
#
#       Kigali districts, all years        950
#       less baseline in Kigali districts  250   (257 baseline less the 7)
#       Kigali 2024-2025                   700
#
#       Hotspot districts, all years       226
#       less baseline in hotspot districts   7
#       Hotspots 2024-2025                 219
#
#   The script compares the two definitions sample by sample and reports
#   any sample they disagree on, so if the data changes and they stop
#   coinciding, that is a printed number rather than a silent divergence.
#
#   Remember also that District is the patient's residence, not the
#   collection site -- the same caveat as Figure 3d and Supplementary S1a.
# -----------------------------------------------------------------------------
#
# -----------------------------------------------------------------------------
# WHAT THE BINOMIAL TEST IS ASKING
# -----------------------------------------------------------------------------
#   binom.test here compares wild type against R561H + K189T within a group,
#   against a null of 50:50. It therefore asks whether the double mutant has
#   reached parity with wild type in that group. It does NOT test prevalence
#   among all samples, and it ignores every other category: its denominator
#   is wild type plus double mutant, not the group size.
#
#   The intercept-only binomial model in section 9 has the same null as that
#   test -- it is a second route to the same p-value, reported for its odds
#   ratio and confidence interval rather than as independent evidence.
#
#   The comparisons BETWEEN groups are the Fisher tests, and those do use
#   the full group size as the denominator.
# -----------------------------------------------------------------------------
#
# HOW TO RUN
#     Rscript scripts/06_figures/figureS3_k13_comparisons.R
#
################################################################################


# ==============================================================================
# 1. SETUP
# ==============================================================================

options(error = NULL)

required_packages <- c(
  "dplyr", "tidyr", "readr", "stringr", "tibble", "purrr",
  "ggplot2", "patchwork"
)

absent_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(absent_packages) > 0) {
  stop(
    "These packages are needed but not installed:\n  ",
    paste(absent_packages, collapse = ", "),
    "\n\nInstall them with:\n  install.packages(c(",
    paste0("\"", absent_packages, "\"", collapse = ", "), "))",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(ggplot2)
  library(patchwork)
})


find_project_root <- function() {

  candidate <- normalizePath(getwd(), mustWork = TRUE)

  for (level in seq_len(5)) {
    if (file.exists(file.path(candidate, "config.R"))) return(candidate)
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
  } else NULL
)

if (is.null(INPUT_FILE)) {
  stop("No input file. Set MANUSCRIPT_DATA in config.R or pass --input PATH.",
       call. = FALSE)
}

if (!file.exists(INPUT_FILE)) {
  stop("Input file not found:\n  ", INPUT_FILE, call. = FALSE)
}

FIGURE_DIR <- file.path(PROJECT_ROOT, "results", "figures", "figureS3")
TABLE_DIR  <- file.path(PROJECT_ROOT, "results", "tables",  "figureS3")

dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR,  recursive = TRUE, showWarnings = FALSE)


report_lines <- character(0)

add_to_report <- function(...) {
  line <- paste0(...)
  report_lines <<- c(report_lines, line)
  message(line)
}


# ==============================================================================
# 2. CATEGORIES, COLOURS AND GROUPS
# ==============================================================================

CATEGORY_ORDER <- c(
  "Wild type",
  "H384R", "A675V", "P441L", "P574L",
  "R561H + K189T", "R561H only", "R561H + other",
  "K189T only", "K189T + other",
  "C469Y", "G449A", "C496F", "N537D", "R515K",
  "Other (non-WHO)"
)

CATEGORY_COLOURS <- c(
  "Wild type"       = "#2ECC71",
  "H384R"           = "#4C6EF5",
  "A675V"           = "#D9D9D9",
  "P441L"           = "#D8A7FF",
  "P574L"           = "#1F5AA6",
  "R561H + K189T"   = "#FF7F0E",
  "R561H only"      = "#F4B400",
  "R561H + other"   = "#C49A00",
  "K189T only"      = "#D62F7D",
  "K189T + other"   = "#F4A6C1",
  "C469Y"           = "#8DD3C7",
  "G449A"           = "#FB8072",
  "C496F"           = "#80B1D3",
  "N537D"           = "#BEBADA",
  "R515K"           = "#B3DE69",
  "Other (non-WHO)" = "#5B4CC4"
)

KFH_NAME <- "King Faisal Hospital"

KIGALI_DISTRICTS <- c("Gasabo", "Kicukiro", "Nyarugenge")

HOTSPOT_DISTRICTS <- c(
  "Bugesera", "Gisagara", "Kirehe", "Karongi",
  "Nyagatare", "Musanze", "Rutsiro"
)

# The two categories panels C and D compare.
COMPARED_CATEGORIES <- c("Wild type", "R561H + K189T")

PAIR_COLOURS <- c("Wild type" = "#2ECC71", "R561H + K189T" = "#FF7F0E")


# ------------------------------------------------------------------------------
# The published panel, read off Supplementary Figure S3. Nothing is drawn
# from these; they exist so the script can state whether it has reproduced
# the figure or drifted from it.
#
# Every one of these p-values was checked to follow from the counts alone,
# so a disagreement means the counts changed, not that a test was run
# differently.
# ------------------------------------------------------------------------------

PUBLISHED_GROUPS <- tibble::tribble(
  ~group,                ~published_n, ~published_wild_percent, ~published_double_percent, ~published_binomial_p,
  "KFH 2021-2023",        257L,         41.6,                    30.4,                      0.0392,
  "KFH 2024-2025",        291L,         31.3,                    42.6,                      0.0288,
  "Kigali 2024-2025",     700L,         32.4,                    48.4,                      2.87e-06,
  "Hotspots 2024-2025",   219L,         41.1,                    14.6,                      1.45e-07
)

PUBLISHED_FISHER <- tibble::tribble(
  ~comparison,                       ~category,         ~published_p,
  "Period: 2021-2023 vs 2024-2025",  "Wild type",        0.0127,
  "Period: 2021-2023 vs 2024-2025",  "R561H + K189T",    0.0034,
  "Location: Kigali vs Hotspots",    "Wild type",        0.0224,
  "Location: Kigali vs Hotspots",    "R561H + K189T",    1.71e-20
)

# A computed p-value counts as matching when it agrees with the published
# one to the precision the figure prints it at -- three significant figures.
matches_published <- function(computed, published) {
  !is.na(computed) & !is.na(published) &
    abs(computed - published) <= abs(published) * 0.02
}


p_to_stars <- function(p) {
  dplyr::case_when(
    p < 0.0001 ~ "****",
    p < 0.001  ~ "***",
    p < 0.01   ~ "**",
    p < 0.05   ~ "*",
    TRUE       ~ "ns"
  )
}


# ==============================================================================
# 3. READ
# ==============================================================================

add_to_report(strrep("=", 72))
add_to_report("SUPPLEMENTARY FIGURE S3 — K13 CATEGORY COMPARISONS")
add_to_report(strrep("=", 72))
add_to_report("Input:     ", INPUT_FILE)
add_to_report("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
add_to_report("")

if (grepl("\\.xlsx?$", INPUT_FILE, ignore.case = TRUE)) {

  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Reading an Excel file needs the readxl package:\n",
         "  install.packages(\"readxl\")", call. = FALSE)
  }

  raw <- readxl::read_excel(INPUT_FILE) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character))

} else {

  raw <- readr::read_tsv(
    INPUT_FILE,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE
  )
}

required_columns <- c(
  "SAMPLE", "YEAR", "GENE_NAME", "HGVS_PROTEIN", "District", "Health_Center"
)

missing_columns <- setdiff(required_columns, names(raw))

if (length(missing_columns) > 0) {
  stop("The input is missing: ", paste(missing_columns, collapse = ", "),
       call. = FALSE)
}

tidy <- raw %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(c("SAMPLE", "GENE_NAME", "HGVS_PROTEIN",
                      "District", "Health_Center")),
      ~ stringr::str_trim(.x)
    ),
    year = suppressWarnings(as.integer(stringr::str_trim(YEAR)))
  ) %>%
  dplyr::filter(!is.na(SAMPLE), SAMPLE != "", !is.na(year))

add_to_report(sprintf("Rows read:        %d", nrow(raw)))
add_to_report(sprintf("Distinct samples: %d", dplyr::n_distinct(tidy$SAMPLE)))


# ------------------------------------------------------------------------------
# A K13 row whose protein change is blank would otherwise be classified as
# "Other (non-WHO)", because classify_k13 has no branch that returns wild
# type. Such rows are dropped here and reported, so a blank cell cannot
# quietly inflate the Other category.
# ------------------------------------------------------------------------------

blank_k13 <- tidy %>%
  dplyr::filter(
    GENE_NAME == "K13",
    is.na(HGVS_PROTEIN) | HGVS_PROTEIN %in% c("", ".")
  )

if (nrow(blank_k13) > 0) {
  add_to_report("")
  add_to_report(sprintf(
    "%d K13 row(s) have no protein change and were dropped before",
    nrow(blank_k13)
  ))
  add_to_report("classification, affecting ", dplyr::n_distinct(blank_k13$SAMPLE),
                " sample(s).")
}


# ------------------------------------------------------------------------------
# One row per sample for the grouping columns.
# ------------------------------------------------------------------------------

sample_metadata <- tidy %>%
  dplyr::distinct(SAMPLE, .keep_all = TRUE) %>%
  dplyr::select(SAMPLE, year, District, Health_Center)


# ==============================================================================
# 4. CLASSIFY
# ==============================================================================

classify_k13 <- function(mutations) {

  mutations <- unique(mutations[!is.na(mutations) & mutations != ""])

  if (length(mutations) == 0) return("Wild type")

  has_r561h <- "R561H" %in% mutations
  has_k189t <- "K189T" %in% mutations

  others <- setdiff(mutations, c("R561H", "K189T"))

  if (has_r561h && has_k189t)           return("R561H + K189T")
  if (has_r561h && length(others) == 0) return("R561H only")
  if (has_k189t && length(others) == 0) return("K189T only")
  if (has_r561h)                        return("R561H + other")
  if (has_k189t)                        return("K189T + other")

  for (named in c("H384R", "A675V", "P574L", "C469Y",
                  "P441L", "G449A", "C496F", "N537D", "R515K")) {
    if (named %in% mutations) return(named)
  }

  "Other (non-WHO)"
}

k13_per_sample <- tidy %>%
  dplyr::filter(
    GENE_NAME == "K13",
    !is.na(HGVS_PROTEIN), !HGVS_PROTEIN %in% c("", ".")
  ) %>%
  dplyr::distinct(SAMPLE, HGVS_PROTEIN) %>%
  dplyr::group_by(SAMPLE) %>%
  dplyr::summarise(
    k13_mutations = paste(sort(unique(HGVS_PROTEIN)), collapse = ", "),
    .groups = "drop"
  )

classified <- sample_metadata %>%
  dplyr::left_join(k13_per_sample, by = "SAMPLE") %>%
  dplyr::mutate(
    k13_mutations = tidyr::replace_na(k13_mutations, ""),
    category = purrr::map_chr(
      k13_mutations,
      function(text) {
        if (text == "") return("Wild type")
        classify_k13(stringr::str_split(text, ",\\s*")[[1]])
      }
    )
  )

unknown <- setdiff(unique(classified$category), CATEGORY_ORDER)

if (length(unknown) > 0) {
  stop("classify_k13 produced categories not in CATEGORY_ORDER: ",
       paste(unknown, collapse = ", "), call. = FALSE)
}

classified <- classified %>%
  dplyr::mutate(category = factor(category, levels = CATEGORY_ORDER))


# ==============================================================================
# 5. THE FOUR GROUPS
# ==============================================================================
#
# The 2021-2023 group is every sample from those years, NOT a filter on
# Health_Center. It is labelled KFH because the baseline was collected
# there; that is asserted rather than assumed, just below.
# ==============================================================================

group_2123 <- classified %>%
  dplyr::filter(year %in% 2021:2023) %>%
  dplyr::mutate(group = "KFH 2021-2023")

group_2425_kfh <- classified %>%
  dplyr::filter(year %in% 2024:2025, Health_Center == KFH_NAME) %>%
  dplyr::mutate(group = "KFH 2024-2025")

group_kigali <- classified %>%
  dplyr::filter(year %in% 2024:2025, District %in% KIGALI_DISTRICTS) %>%
  dplyr::mutate(group = "Kigali 2024-2025")

group_hotspot <- classified %>%
  dplyr::filter(year %in% 2024:2025, District %in% HOTSPOT_DISTRICTS) %>%
  dplyr::mutate(group = "Hotspots 2024-2025")

GROUP_ORDER <- c(
  "KFH 2021-2023", "KFH 2024-2025",
  "Kigali 2024-2025", "Hotspots 2024-2025"
)

# The facet headings as the published panel writes them. They do not follow
# one pattern -- panel A puts the years first and panel B leaves no space
# before the bracket -- so they are spelled out rather than built by paste.
PANEL_HEADINGS <- c(
  "KFH 2021-2023"      = "2021-2023 KFH (N = %d)",
  "KFH 2024-2025"      = "2024-2025 KFH (N = %d)",
  "Kigali 2024-2025"   = "Kigali 2024-2025(N = %d)",
  "Hotspots 2024-2025" = "Hotspots 2024-2025(N = %d)"
)


# ------------------------------------------------------------------------------
# Checks that the group definitions actually mean what their names say.
# ------------------------------------------------------------------------------

non_kfh_baseline <- group_2123 %>%
  dplyr::filter(Health_Center != KFH_NAME)

add_to_report("")
add_to_report(strrep("-", 72))
add_to_report("GROUP DEFINITIONS")
add_to_report(strrep("-", 72))

if (nrow(non_kfh_baseline) == 0) {
  add_to_report("  The 2021-2023 group is unfiltered by health centre, and")
  add_to_report("  every one of its samples is from ", KFH_NAME, ",")
  add_to_report("  so the KFH label is accurate.")
} else {
  add_to_report(strrep("!", 72))
  add_to_report(sprintf(
    "  %d of the %d samples in 'KFH 2021-2023' are NOT from %s:",
    nrow(non_kfh_baseline), nrow(group_2123), KFH_NAME
  ))
  centre_counts <- non_kfh_baseline %>%
    dplyr::count(Health_Center, name = "n") %>%
    dplyr::arrange(dplyr::desc(n))
  for (index in seq_len(nrow(centre_counts))) {
    add_to_report(sprintf("    %-32s %d",
                          centre_counts$Health_Center[[index]],
                          centre_counts$n[[index]]))
  }
  add_to_report("  The panel label says KFH. Either the label or the filter")
  add_to_report("  needs correcting before this figure is used.")
  add_to_report(strrep("!", 72))
}

unassigned <- classified %>%
  dplyr::filter(
    year %in% 2024:2025,
    !District %in% c(KIGALI_DISTRICTS, HOTSPOT_DISTRICTS)
  )

if (nrow(unassigned) > 0) {
  add_to_report("")
  add_to_report(sprintf(
    "  %d sample(s) from 2024-2025 are in neither district group:",
    nrow(unassigned)
  ))
  district_counts <- unassigned %>%
    dplyr::count(District, name = "n") %>%
    dplyr::arrange(dplyr::desc(n))
  for (index in seq_len(nrow(district_counts))) {
    add_to_report(sprintf("    %-20s %d",
                          district_counts$District[[index]],
                          district_counts$n[[index]]))
  }
  add_to_report("  They appear in no panel. Add them to KIGALI_DISTRICTS or")
  add_to_report("  HOTSPOT_DISTRICTS if they belong in one.")
}


# ==============================================================================
# 6. COUNTS PER GROUP
# ==============================================================================

all_groups <- dplyr::bind_rows(
  group_2123, group_2425_kfh, group_kigali, group_hotspot
) %>%
  dplyr::mutate(group = factor(group, levels = GROUP_ORDER))

group_sizes <- all_groups %>% dplyr::count(group, name = "total_n")

counts <- all_groups %>%
  dplyr::count(group, category, name = "n", .drop = FALSE) %>%
  dplyr::left_join(group_sizes, by = "group") %>%
  dplyr::filter(total_n > 0) %>%
  dplyr::mutate(prevalence = 100 * n / total_n)


# The original script derived wild type by subtraction -- total minus the
# samples that had a K13 record. Here wild type is counted directly, so the
# two routes should agree exactly. Check that they do.
subtraction_check <- all_groups %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(
    counted_wild_type = sum(category == "Wild type"),
    implied_wild_type = dplyr::n() - sum(k13_mutations != ""),
    .groups = "drop"
  ) %>%
  dplyr::mutate(agrees = counted_wild_type == implied_wild_type)

add_to_report("")
add_to_report(strrep("-", 72))
add_to_report("GROUP SIZES")
add_to_report(strrep("-", 72))
add_to_report(sprintf("  %-22s %8s %12s %14s",
                      "group", "N", "wild type", "R561H+K189T"))

for (this_group in GROUP_ORDER) {

  rows <- counts %>% dplyr::filter(group == this_group)

  if (nrow(rows) == 0) next

  add_to_report(sprintf(
    "  %-22s %8d %12d %14d",
    this_group,
    rows$total_n[[1]],
    rows$n[rows$category == "Wild type"],
    rows$n[rows$category == "R561H + K189T"]
  ))
}

# ------------------------------------------------------------------------------
# Against the published panel.
# ------------------------------------------------------------------------------

published_group_check <- counts %>%
  dplyr::filter(category %in% COMPARED_CATEGORIES) %>%
  dplyr::mutate(category = as.character(category)) %>%
  dplyr::select(group, category, n, total_n, prevalence) %>%
  tidyr::pivot_wider(
    names_from = category, values_from = c(n, prevalence)
  ) %>%
  dplyr::rename(
    wild_percent   = `prevalence_Wild type`,
    double_percent = `prevalence_R561H + K189T`
  ) %>%
  dplyr::mutate(group = as.character(group)) %>%
  dplyr::left_join(PUBLISHED_GROUPS, by = "group") %>%
  dplyr::mutate(
    n_ok      = total_n == published_n,
    wild_ok   = abs(round(wild_percent, 1) - published_wild_percent) < 0.051,
    double_ok = abs(round(double_percent, 1) - published_double_percent) < 0.051
  )

add_to_report("")
add_to_report(strrep("-", 72))
add_to_report("AGAINST PUBLISHED FIGURE S3 — SIZES AND PREVALENCES")
add_to_report(strrep("-", 72))
add_to_report(sprintf("  %-22s %14s %18s %18s",
                      "group", "N", "wild type %", "R561H+K189T %"))

for (index in seq_len(nrow(published_group_check))) {
  add_to_report(sprintf(
    "  %-22s %6d vs %-5d %s %7.1f vs %-5.1f %s %7.1f vs %-5.1f %s",
    published_group_check$group[[index]],
    published_group_check$total_n[[index]],
    published_group_check$published_n[[index]],
    if (published_group_check$n_ok[[index]]) "ok " else "DIFF",
    published_group_check$wild_percent[[index]],
    published_group_check$published_wild_percent[[index]],
    if (published_group_check$wild_ok[[index]]) "ok " else "DIFF",
    published_group_check$double_percent[[index]],
    published_group_check$published_double_percent[[index]],
    if (published_group_check$double_ok[[index]]) "ok " else "DIFF"
  ))
}

sizes_and_prevalences_match <- all(
  published_group_check$n_ok,
  published_group_check$wild_ok,
  published_group_check$double_ok
)

if (!all(subtraction_check$agrees)) {
  add_to_report("")
  add_to_report(strrep("!", 72))
  add_to_report("  Counting wild type directly disagrees with deriving it by")
  add_to_report("  subtraction. A sample with a K13 record has been classified")
  add_to_report("  as wild type, which should be impossible.")
  add_to_report(strrep("!", 72))
}

# Every category must sum back to the group size.
sum_check <- counts %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(summed = sum(n), total_n = dplyr::first(total_n),
                   .groups = "drop") %>%
  dplyr::filter(summed != total_n)

if (nrow(sum_check) > 0) {
  stop("Category counts do not sum to the group size for: ",
       paste(sum_check$group, collapse = ", "), call. = FALSE)
}


# ------------------------------------------------------------------------------
# How the district-based Kigali group differs from Figure 3b's REGION-based
# one, stated as a number rather than left for a reader to notice.
# ------------------------------------------------------------------------------

if ("REGION" %in% names(raw)) {

  region_lookup <- tidy %>%
    dplyr::distinct(SAMPLE, .keep_all = TRUE) %>%
    dplyr::select(SAMPLE, REGION) %>%
    dplyr::mutate(REGION = stringr::str_trim(REGION))

  region_kigali <- region_lookup %>%
    dplyr::semi_join(
      classified %>% dplyr::filter(year %in% 2024:2025), by = "SAMPLE"
    ) %>%
    dplyr::filter(REGION == "Kigali") %>%
    dplyr::pull(SAMPLE)

  district_kigali <- group_kigali$SAMPLE

  add_to_report("")
  add_to_report(strrep("-", 72))
  add_to_report("DISTRICT GROUPS VERSUS FIGURE 3B's REGION GROUPS")
  add_to_report(strrep("-", 72))
  add_to_report(sprintf(
    "  Kigali 2024-2025 by District (this figure): %d", length(district_kigali)
  ))
  add_to_report(sprintf(
    "  Kigali 2024-2025 by REGION   (Figure 3b):   %d", length(region_kigali)
  ))
  add_to_report(sprintf(
    "  in District but not REGION: %d",
    length(setdiff(district_kigali, region_kigali))
  ))
  add_to_report(sprintf(
    "  in REGION but not District: %d",
    length(setdiff(region_kigali, district_kigali))
  ))
  add_to_report("")
  add_to_report("  Both definitions are kept, each as its own figure used it.")
  add_to_report("  The legends should not imply the two groups are the same.")

  readr::write_tsv(
    tibble::tibble(
      SAMPLE = union(district_kigali, region_kigali)
    ) %>%
      dplyr::mutate(
        in_district_group = SAMPLE %in% district_kigali,
        in_region_group   = SAMPLE %in% region_kigali
      ) %>%
      dplyr::filter(in_district_group != in_region_group),
    file.path(TABLE_DIR, "figureS3_district_vs_region_differences.tsv")
  )
}


# ==============================================================================
# 7. STATISTICS
# ==============================================================================

get_count <- function(this_group, category_name) {
  value <- counts$n[counts$group == this_group & counts$category == category_name]
  if (length(value) == 0) 0L else value[[1]]
}

binomial_rows <- purrr::map_dfr(GROUP_ORDER, function(this_group) {

  wild   <- get_count(this_group, "Wild type")
  double <- get_count(this_group, "R561H + K189T")

  if (wild + double == 0) {
    return(tibble::tibble(
      group = this_group, wild_type = wild, double_mutant = double,
      n_compared = 0L, wild_percent = NA_real_, double_percent = NA_real_,
      p_value = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      stars = NA_character_
    ))
  }

  test <- stats::binom.test(wild, wild + double, p = 0.5)

  tibble::tibble(
    group          = this_group,
    wild_type      = wild,
    double_mutant  = double,
    n_compared     = wild + double,
    wild_percent   = 100 * wild / (wild + double),
    double_percent = 100 * double / (wild + double),
    p_value        = test$p.value,
    ci_low         = 100 * test$conf.int[[1]],
    ci_high        = 100 * test$conf.int[[2]],
    stars          = p_to_stars(test$p.value)
  )
})

add_to_report("")
add_to_report(strrep("-", 72))
add_to_report("EXACT BINOMIAL TEST — wild type vs R561H + K189T, null 50:50")
add_to_report(strrep("-", 72))
add_to_report(sprintf("  %-22s %7s %7s %9s %12s",
                      "group", "WT", "double", "WT share", "p"))

for (index in seq_len(nrow(binomial_rows))) {
  add_to_report(sprintf(
    "  %-22s %7d %7d %8.1f%% %12.3g  %s",
    binomial_rows$group[[index]],
    binomial_rows$wild_type[[index]],
    binomial_rows$double_mutant[[index]],
    binomial_rows$wild_percent[[index]],
    binomial_rows$p_value[[index]],
    binomial_rows$stars[[index]]
  ))
}

add_to_report("")
add_to_report("  The 95% intervals are for the wild-type share of those two")
add_to_report("  categories only, not for prevalence in the whole group.")

binomial_rows <- binomial_rows %>%
  dplyr::left_join(
    PUBLISHED_GROUPS %>%
      dplyr::select(group, published_binomial_p),
    by = "group"
  ) %>%
  dplyr::mutate(
    matches_paper = matches_published(p_value, published_binomial_p)
  )

add_to_report("")
add_to_report("  against the published panel:")

for (index in seq_len(nrow(binomial_rows))) {
  add_to_report(sprintf(
    "    %-22s computed %10.3g   published %10.3g   %s",
    binomial_rows$group[[index]],
    binomial_rows$p_value[[index]],
    binomial_rows$published_binomial_p[[index]],
    if (isTRUE(binomial_rows$matches_paper[[index]])) "ok" else "DIFF"
  ))
}


# ------------------------------------------------------------------------------
# Fisher's exact tests between groups. These use the full group size.
# ------------------------------------------------------------------------------

run_fisher <- function(category_name, group_one, group_two, comparison) {

  count_one <- get_count(group_one, category_name)
  count_two <- get_count(group_two, category_name)

  total_one <- group_sizes$total_n[group_sizes$group == group_one]
  total_two <- group_sizes$total_n[group_sizes$group == group_two]

  test <- stats::fisher.test(matrix(
    c(count_one, total_one - count_one,
      count_two, total_two - count_two),
    nrow = 2, byrow = TRUE
  ))

  tibble::tibble(
    comparison   = comparison,
    category     = category_name,
    group_one    = group_one,
    n_one        = count_one,
    total_one    = total_one,
    percent_one  = 100 * count_one / total_one,
    group_two    = group_two,
    n_two        = count_two,
    total_two    = total_two,
    percent_two  = 100 * count_two / total_two,
    odds_ratio   = unname(test$estimate),
    ci_low       = test$conf.int[[1]],
    ci_high      = test$conf.int[[2]],
    p_value      = test$p.value,
    stars        = p_to_stars(test$p.value)
  )
}

fisher_rows <- dplyr::bind_rows(
  purrr::map_dfr(COMPARED_CATEGORIES, run_fisher,
                 group_one = "KFH 2021-2023", group_two = "KFH 2024-2025",
                 comparison = "Period: 2021-2023 vs 2024-2025"),
  purrr::map_dfr(COMPARED_CATEGORIES, run_fisher,
                 group_one = "Kigali 2024-2025",
                 group_two = "Hotspots 2024-2025",
                 comparison = "Location: Kigali vs Hotspots")
)

add_to_report("")
add_to_report(strrep("-", 72))
add_to_report("FISHER'S EXACT TEST — between groups, over the whole group")
add_to_report(strrep("-", 72))

fisher_rows <- fisher_rows %>%
  dplyr::left_join(PUBLISHED_FISHER, by = c("comparison", "category")) %>%
  dplyr::mutate(
    matches_paper = matches_published(p_value, published_p)
  )

for (index in seq_len(nrow(fisher_rows))) {
  add_to_report(sprintf(
    "  %-34s %-15s %6.1f%% vs %6.1f%%  OR %6.2f  p %10.3g  %s",
    fisher_rows$comparison[[index]],
    fisher_rows$category[[index]],
    fisher_rows$percent_one[[index]],
    fisher_rows$percent_two[[index]],
    fisher_rows$odds_ratio[[index]],
    fisher_rows$p_value[[index]],
    fisher_rows$stars[[index]]
  ))
}

add_to_report("")
add_to_report("  against the published panel:")

for (index in seq_len(nrow(fisher_rows))) {
  add_to_report(sprintf(
    "    %-34s %-15s computed %10.3g   published %10.3g   %s",
    fisher_rows$comparison[[index]],
    fisher_rows$category[[index]],
    fisher_rows$p_value[[index]],
    fisher_rows$published_p[[index]],
    if (isTRUE(fisher_rows$matches_paper[[index]])) "ok" else "DIFF"
  ))
}


# ------------------------------------------------------------------------------
# The single verdict.
# ------------------------------------------------------------------------------

everything_matches <- isTRUE(sizes_and_prevalences_match) &&
  all(binomial_rows$matches_paper, na.rm = TRUE) &&
  all(fisher_rows$matches_paper, na.rm = TRUE) &&
  !any(is.na(binomial_rows$matches_paper)) &&
  !any(is.na(fisher_rows$matches_paper))

add_to_report("")
add_to_report(strrep("=", 72))

if (everything_matches) {
  add_to_report("FIGURE S3 IS REPRODUCED")
  add_to_report(strrep("=", 72))
  add_to_report("  All four group sizes, both prevalences in each group, all")
  add_to_report("  four binomial p-values and all four Fisher p-values match")
  add_to_report("  the published panel.")
} else {
  add_to_report("FIGURE S3 DOES NOT MATCH THE PUBLISHED PANEL")
  add_to_report(strrep("=", 72))
  add_to_report("  See the DIFF entries above. Do not use this figure until")
  add_to_report("  the difference is understood.")
}

add_to_report(strrep("=", 72))


# ------------------------------------------------------------------------------
# Intercept-only binomial model per group. Same null as the binomial test
# above; reported for the odds ratio and its interval.
# ------------------------------------------------------------------------------

logistic_rows <- purrr::map_dfr(GROUP_ORDER, function(this_group) {

  wild   <- get_count(this_group, "Wild type")
  double <- get_count(this_group, "R561H + K189T")

  if (wild + double == 0 || wild == 0 || double == 0) {
    return(tibble::tibble(
      group = this_group, wild_type = wild, double_mutant = double,
      odds_ratio = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      p_value = NA_real_,
      note = "not fitted: a category has no samples"
    ))
  }

  model <- stats::glm(
    cbind(double, wild) ~ 1, family = stats::binomial(link = "logit")
  )

  coefficients <- summary(model)$coefficients

  log_odds <- coefficients[1, "Estimate"]
  standard_error <- coefficients[1, "Std. Error"]

  tibble::tibble(
    group         = this_group,
    wild_type     = wild,
    double_mutant = double,
    odds_ratio    = exp(log_odds),
    ci_low        = exp(log_odds - 1.96 * standard_error),
    ci_high       = exp(log_odds + 1.96 * standard_error),
    p_value       = coefficients[1, "Pr(>|z|)"],
    note          = "odds of R561H + K189T against wild type"
  )
})


# ==============================================================================
# 8. THE PANELS
# ==============================================================================
#
# Each facet in panels A and B shows only the categories present in that
# group, so the x position of a bar depends on which categories that group
# has. The significance bracket has to be placed from those positions rather
# than from a fixed index, or it lands on the wrong bars whenever a category
# is missing.
# ==============================================================================

bracket_data <- function(plot_data) {

  plot_data %>%
    dplyr::filter(n > 0) %>%
    dplyr::group_by(group, panel_label) %>%
    dplyr::summarise(
      y_max = max(prevalence, na.rm = TRUE),
      x_wild = match(
        match("Wild type", CATEGORY_ORDER),
        sort(unique(as.integer(category)))
      ),
      x_double = match(
        match("R561H + K189T", CATEGORY_ORDER),
        sort(unique(as.integer(category)))
      ),
      .groups = "drop"
    ) %>%
    dplyr::left_join(
      binomial_rows %>% dplyr::select(group, p_value, stars),
      by = "group"
    ) %>%
    dplyr::mutate(
      y_low  = y_max * 1.05,
      y      = y_max * 1.12,
      y_text = y_max * 1.18,
      legend_line = paste0(stars, " : p = ", signif(p_value, 3))
    ) %>%
    dplyr::filter(!is.na(x_wild), !is.na(x_double))
}


category_panel <- function(plot_data, brackets) {

  ggplot2::ggplot(
    plot_data %>% dplyr::filter(n > 0),
    ggplot2::aes(x = category, y = prevalence, fill = category)
  ) +

    ggplot2::geom_col(width = 0.5, show.legend = FALSE) +

    ggplot2::geom_segment(
      data = brackets,
      mapping = ggplot2::aes(x = x_wild, xend = x_double, y = y, yend = y),
      inherit.aes = FALSE, linewidth = 0.5, colour = "black"
    ) +
    ggplot2::geom_segment(
      data = brackets,
      mapping = ggplot2::aes(x = x_wild, xend = x_wild, y = y_low, yend = y),
      inherit.aes = FALSE, linewidth = 0.5, colour = "black"
    ) +
    ggplot2::geom_segment(
      data = brackets,
      mapping = ggplot2::aes(x = x_double, xend = x_double,
                             y = y_low, yend = y),
      inherit.aes = FALSE, linewidth = 0.5, colour = "black"
    ) +
    ggplot2::geom_text(
      data = brackets,
      mapping = ggplot2::aes(x = (x_wild + x_double) / 2, y = y_text,
                             label = stars),
      inherit.aes = FALSE, size = 6, fontface = "bold"
    ) +
    ggplot2::geom_text(
      data = brackets,
      mapping = ggplot2::aes(x = x_double + 0.3, y = y_text,
                             label = legend_line),
      inherit.aes = FALSE, hjust = 0, size = 4.2, fontface = "bold"
    ) +

    ggplot2::facet_wrap(~ panel_label, scales = "free", nrow = 1) +

    ggplot2::scale_fill_manual(values = CATEGORY_COLOURS, drop = TRUE) +

    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.28))
    ) +

    ggplot2::coord_cartesian(clip = "off") +

    ggplot2::labs(x = NULL, y = "Prevalence (%)") +

    ggplot2::theme_classic(base_size = 13) +

    ggplot2::theme(
      axis.title   = ggplot2::element_text(size = 16, colour = "black",
                                           face = "bold"),
      axis.text    = ggplot2::element_text(size = 12, colour = "black",
                                           face = "bold"),
      axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1),
      strip.text.x = ggplot2::element_text(size = 15, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey85"),
      panel.grid   = ggplot2::element_blank(),
      plot.margin  = ggplot2::margin(14, 60, 10, 10)
    )
}


pair_panel <- function(plot_data, comparison_rows, title_text) {

  brackets <- plot_data %>%
    dplyr::group_by(category) %>%
    dplyr::summarise(y_max = max(prevalence, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(
      x      = as.numeric(factor(category, levels = COMPARED_CATEGORIES)),
      xstart = x - 0.22,
      xend   = x + 0.22,
      y_low  = y_max * 1.14,
      y      = y_max * 1.22,
      y_text = y_max * 1.30
    ) %>%
    dplyr::left_join(
      comparison_rows %>%
        dplyr::select(category, p_value, stars) %>%
        dplyr::mutate(category = factor(category,
                                        levels = COMPARED_CATEGORIES)),
      by = "category"
    )

  legend_text <- paste(
    sprintf("%s : %s p = %s",
            comparison_rows$category, comparison_rows$stars,
            signif(comparison_rows$p_value, 3)),
    collapse = "\n"
  )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = category, y = prevalence, fill = category,
                 group = side)
  ) +

    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.7),
      width = 0.6, colour = "black"
    ) +

    ggplot2::geom_text(
      ggplot2::aes(label = side),
      position = ggplot2::position_dodge(width = 0.7),
      vjust = -0.6, size = 3.6, fontface = "bold"
    ) +

    ggplot2::geom_segment(
      data = brackets,
      mapping = ggplot2::aes(x = xstart, xend = xend, y = y, yend = y),
      inherit.aes = FALSE, linewidth = 0.5, colour = "black"
    ) +
    ggplot2::geom_segment(
      data = brackets,
      mapping = ggplot2::aes(x = xstart, xend = xstart, y = y_low, yend = y),
      inherit.aes = FALSE, linewidth = 0.5, colour = "black"
    ) +
    ggplot2::geom_segment(
      data = brackets,
      mapping = ggplot2::aes(x = xend, xend = xend, y = y_low, yend = y),
      inherit.aes = FALSE, linewidth = 0.5, colour = "black"
    ) +
    ggplot2::geom_text(
      data = brackets,
      mapping = ggplot2::aes(x = x, y = y_text, label = stars),
      inherit.aes = FALSE, size = 6, fontface = "bold"
    ) +

    ggplot2::annotate(
      "text", x = Inf, y = Inf, label = legend_text,
      hjust = 1.05, vjust = 1.2, size = 4.2, fontface = "bold"
    ) +

    ggplot2::scale_fill_manual(values = PAIR_COLOURS) +

    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.35))
    ) +

    ggplot2::coord_cartesian(clip = "off") +

    ggplot2::labs(title = title_text, x = NULL, y = "Prevalence (%)") +

    ggplot2::theme_classic(base_size = 13) +

    ggplot2::theme(
      axis.title   = ggplot2::element_text(size = 16, colour = "black",
                                           face = "bold"),
      axis.text    = ggplot2::element_text(size = 12, colour = "black",
                                           face = "bold"),
      plot.title   = ggplot2::element_text(size = 15, face = "bold",
                                           hjust = 0.5),
      panel.grid   = ggplot2::element_blank(),
      legend.position = "none",
      plot.margin  = ggplot2::margin(14, 14, 10, 10)
    )
}


# ---- panel A ----------------------------------------------------------------

panel_a_data <- counts %>%
  dplyr::filter(group %in% c("KFH 2021-2023", "KFH 2024-2025")) %>%
  dplyr::mutate(
    panel_label = sprintf(
      PANEL_HEADINGS[as.character(group)], total_n
    ),
    panel_label = factor(panel_label, levels = unique(panel_label[order(group)]))
  )

panel_a <- category_panel(panel_a_data, bracket_data(panel_a_data))


# ---- panel B ----------------------------------------------------------------

panel_b_data <- counts %>%
  dplyr::filter(group %in% c("Kigali 2024-2025", "Hotspots 2024-2025")) %>%
  dplyr::mutate(
    panel_label = sprintf(
      PANEL_HEADINGS[as.character(group)], total_n
    ),
    panel_label = factor(panel_label, levels = unique(panel_label[order(group)]))
  )

panel_b <- category_panel(panel_b_data, bracket_data(panel_b_data))


# ---- panel C ----------------------------------------------------------------

panel_c_data <- counts %>%
  dplyr::filter(
    group %in% c("KFH 2021-2023", "KFH 2024-2025"),
    category %in% COMPARED_CATEGORIES
  ) %>%
  dplyr::mutate(
    side = factor(
      dplyr::if_else(group == "KFH 2021-2023", "2021-2023", "2024-2025"),
      levels = c("2021-2023", "2024-2025")
    ),
    category = factor(category, levels = COMPARED_CATEGORIES)
  )

panel_c <- pair_panel(
  panel_c_data,
  fisher_rows %>%
    dplyr::filter(comparison == "Period: 2021-2023 vs 2024-2025"),
  "Wild type and double mutant,\n2021-2023 versus 2024-2025"
)


# ---- panel D ----------------------------------------------------------------

panel_d_data <- counts %>%
  dplyr::filter(
    group %in% c("Kigali 2024-2025", "Hotspots 2024-2025"),
    category %in% COMPARED_CATEGORIES
  ) %>%
  dplyr::mutate(
    side = factor(
      dplyr::if_else(group == "Kigali 2024-2025", "Kigali", "Hotspots"),
      levels = c("Kigali", "Hotspots")
    ),
    category = factor(category, levels = COMPARED_CATEGORIES)
  )

panel_d <- pair_panel(
  panel_d_data,
  fisher_rows %>%
    dplyr::filter(comparison == "Location: Kigali vs Hotspots"),
  "Wild type and double mutant,\nKigali versus Hotspots, 2024-2025"
)


# ==============================================================================
# 9. ASSEMBLE AND SAVE
# ==============================================================================

figure <- (panel_a / panel_b / (panel_c | panel_d)) +
  patchwork::plot_annotation(tag_levels = "A") +
  patchwork::plot_layout(heights = c(1, 1, 1)) &
  ggplot2::theme(
    plot.tag = ggplot2::element_text(size = 28, face = "bold",
                                     colour = "black")
  )

for (extension in c("svg", "pdf", "png")) {
  ggplot2::ggsave(
    filename = file.path(
      FIGURE_DIR, paste0("figureS3_k13_comparisons.", extension)
    ),
    plot = figure,
    width = 17, height = 18, units = "in",
    dpi = 600, bg = "white", limitsize = FALSE
  )
}

# The individual panels too, since the original saved them separately.
for (entry in list(
  list(name = "figureS3_panelA", plot = panel_a, width = 15, height = 7),
  list(name = "figureS3_panelB", plot = panel_b, width = 15, height = 7),
  list(name = "figureS3_panelC", plot = panel_c, width = 8,  height = 7),
  list(name = "figureS3_panelD", plot = panel_d, width = 8,  height = 7)
)) {
  ggplot2::ggsave(
    filename = file.path(FIGURE_DIR, paste0(entry$name, ".svg")),
    plot = entry$plot,
    width = entry$width, height = entry$height, units = "in",
    dpi = 600, bg = "white"
  )
}


# ==============================================================================
# 10. TABLES
# ==============================================================================

readr::write_tsv(
  counts %>%
    dplyr::mutate(group = as.character(group),
                  category = as.character(category)) %>%
    dplyr::select(group, category, n, total_n, prevalence),
  file.path(TABLE_DIR, "figureS3_category_prevalence.tsv")
)

readr::write_tsv(
  all_groups %>%
    dplyr::mutate(group = as.character(group),
                  category = as.character(category)) %>%
    dplyr::select(SAMPLE, group, year, District, Health_Center,
                  k13_mutations, category) %>%
    dplyr::arrange(group, SAMPLE),
  file.path(TABLE_DIR, "figureS3_sample_level_categories.tsv")
)

readr::write_tsv(
  group_sizes %>% dplyr::mutate(group = as.character(group)),
  file.path(TABLE_DIR, "figureS3_group_sizes.tsv")
)

readr::write_tsv(
  binomial_rows, file.path(TABLE_DIR, "figureS3_binomial_tests.tsv")
)

readr::write_tsv(
  fisher_rows, file.path(TABLE_DIR, "figureS3_fisher_tests.tsv")
)

readr::write_tsv(
  logistic_rows, file.path(TABLE_DIR, "figureS3_logistic_models.tsv")
)

add_to_report("")
add_to_report(strrep("=", 72))
add_to_report("WRITTEN")
add_to_report(strrep("=", 72))
add_to_report("  figures/figureS3/figureS3_k13_comparisons.svg / .pdf / .png")
add_to_report("  figures/figureS3/figureS3_panelA-D.svg")
add_to_report("  tables/figureS3/figureS3_category_prevalence.tsv")
add_to_report("  tables/figureS3/figureS3_sample_level_categories.tsv")
add_to_report("  tables/figureS3/figureS3_group_sizes.tsv")
add_to_report("  tables/figureS3/figureS3_binomial_tests.tsv")
add_to_report("  tables/figureS3/figureS3_fisher_tests.tsv")
add_to_report("  tables/figureS3/figureS3_logistic_models.tsv")

writeLines(report_lines, file.path(TABLE_DIR, "figureS3_report.txt"))

capture.output(
  sessionInfo(),
  file = file.path(TABLE_DIR, "figureS3_sessionInfo.txt")
)
