# Figure 3b — K13 mutation categories across four sample groups
#
# Inputs:  Manuscript_data_29-07-2026.tsv  (MANUSCRIPT_DATA in config.R)
# Outputs: results/figures/figure03/figure03b_k13_category_pies.svg / .png
#          results/tables/figure03/figure03b_category_prevalence.tsv   <- Source Data
#          results/tables/figure03/figure03b_sample_level_categories.tsv
#          results/tables/figure03/figure03b_report.txt
#
# Usage:   Rscript scripts/06_figures/figure03b_k13_category_pies.R

options(error = NULL)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(ggplot2)
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

FIGURE_DIR <- file.path(PROJECT_ROOT, "results", "figures", "figure03")
TABLE_DIR  <- file.path(PROJECT_ROOT, "results", "tables",  "figure03")

dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR,  recursive = TRUE, showWarnings = FALSE)

KFH_NAME <- "KFH"


report_lines <- character(0)

add_to_report <- function(...) {
  line <- paste0(...)
  report_lines <<- c(report_lines, line)
  message(line)
}


# ---- categories, colours and the published numbers ---------------------

CATEGORY_ORDER <- c(
  "Wild type",
  "H384R", "A675V", "P441L", "P574L",
  "R561H + K189T", "R561H only", "R561H + other",
  "K189T only", "K189T + other",
  "C469Y", "G449A", "C496F", "N537D", "R515K",
  "Other (non-WHO)"
)

# Read off the published legend
CATEGORY_COLOURS <- c(
  "Wild type"      = "#2ECC71",
  "H384R"          = "#6690FF",
  "A675V"          = "#DCE4F7",
  "P441L"          = "#5B6BBF",
  "P574L"          = "#36476B",
  "R561H + K189T"  = "#FF7F0E",
  "R561H only"     = "#F4B400",
  "R561H + other"  = "#FFE0B2",
  "K189T only"     = "#D62F7D",
  "K189T + other"  = "#F8AFCB",
  "C469Y"          = "#8C6D31",
  "G449A"          = "#B5705A",
  "C496F"          = "#7FA8A0",
  "N537D"          = "#9E9E9E",
  "R515K"          = "#C9B458",
  "Other (non-WHO)" = "#6A5ACD"
)

# Only these two carry a printed percentage on the published panel
LABELLED_CATEGORIES <- c("Wild type", "R561H + K189T")

PUBLISHED_GROUPS <- tibble::tribble(
  ~group,                 ~published_n, ~published_wild_type, ~published_r561h_k189t,
  "Baseline-Kigali",       257L,         41.6,                 30.4,
  "KFH Resurgence",        291L,         31.3,                 42.6,
  "Kigali Resurgence",     700L,         32.4,                 48.4,
  "Hotspot Resurgence",    219L,         41.1,                 14.6
)

GROUP_ORDER <- PUBLISHED_GROUPS$group

GROUP_SUBTITLES <- c(
  "Baseline-Kigali"    = "(2021-2023)",
  "KFH Resurgence"     = "(2024-2025)",
  "Kigali Resurgence"  = "(2024-2025)",
  "Hotspot Resurgence" = "(2024-2025)"
)


# ---- read --------------------------------------------------------------

add_to_report(strrep("=", 72))
add_to_report("FIGURE 3B — K13 CATEGORIES ACROSS FOUR GROUPS")
add_to_report(strrep("=", 72))
add_to_report("Input:     ", INPUT_FILE)
add_to_report("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
add_to_report("")

raw <- readr::read_tsv(
  INPUT_FILE,
  col_types = readr::cols(.default = readr::col_character()),
  show_col_types = FALSE
)

required_columns <- c(
  "SAMPLE", "YEAR", "Site_Group", "District",
  "PERIOD", "REGION", "GENE_NAME", "HGVS_PROTEIN"
)

missing_columns <- setdiff(required_columns, names(raw))

if (length(missing_columns) > 0) {
  stop("The input is missing: ", paste(missing_columns, collapse = ", "),
       call. = FALSE)
}

tidy <- raw %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(c("SAMPLE", "Site_Group", "District",
                      "PERIOD", "REGION", "GENE_NAME", "HGVS_PROTEIN")),
      ~ stringr::str_trim(.x)
    ),
    year = suppressWarnings(as.integer(stringr::str_trim(YEAR)))
  ) %>%
  dplyr::filter(!is.na(SAMPLE), SAMPLE != "")

add_to_report(sprintf("Rows read: %d", nrow(raw)))


# One row per sample

sample_metadata <- tidy %>%
  dplyr::distinct(SAMPLE, .keep_all = TRUE) %>%
  dplyr::select(
    SAMPLE, year, Site_Group, District, PERIOD, REGION
  )

add_to_report(sprintf("Distinct samples: %d", nrow(sample_metadata)))


# Guard against a health-centre name that has drifted

if (!KFH_NAME %in% sample_metadata$Site_Group) {

  stop(
    "No sample has Site_Group == \"", KFH_NAME, "\".\n\n",
    "The values present are:\n  ",
    paste(sort(unique(sample_metadata$Site_Group)), collapse = "\n  "),
    "\n\nUpdate KFH_NAME near the top of this script to the exact spelling.",
    call. = FALSE
  )
}


# ---- classify each sample ----------------------------------------------

classify_k13 <- function(mutations) {

  mutations <- unique(mutations[!is.na(mutations) & mutations != ""])

  if (length(mutations) == 0) return("Wild type")

  has_r561h <- "R561H" %in% mutations
  has_k189t <- "K189T" %in% mutations

  others <- setdiff(mutations, c("R561H", "K189T"))

  if (has_r561h && has_k189t)            return("R561H + K189T")
  if (has_r561h && length(others) == 0)  return("R561H only")
  if (has_k189t && length(others) == 0)  return("K189T only")
  if (has_r561h)                         return("R561H + other")
  if (has_k189t)                         return("K189T + other")

  for (named in c("H384R", "A675V", "P574L", "C469Y",
                  "P441L", "G449A", "C496F", "N537D", "R515K")) {
    if (named %in% mutations) return(named)
  }

  "Other (non-WHO)"
}

k13_per_sample <- tidy %>%
  dplyr::filter(
    GENE_NAME == "K13",
    !is.na(HGVS_PROTEIN), HGVS_PROTEIN != "", HGVS_PROTEIN != "."
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

unknown_categories <- setdiff(unique(classified$category), CATEGORY_ORDER)

if (length(unknown_categories) > 0) {
  stop(
    "classify_k13 produced categories not in CATEGORY_ORDER: ",
    paste(unknown_categories, collapse = ", "),
    call. = FALSE
  )
}

classified <- classified %>%
  dplyr::mutate(category = factor(category, levels = CATEGORY_ORDER))

add_to_report(sprintf(
  "Samples carrying at least one K13 mutation: %d",
  sum(classified$k13_mutations != "")
))


# ---- the four groups ---------------------------------------------------

group_membership <- dplyr::bind_rows(

  classified %>%
    dplyr::filter(
      Site_Group == KFH_NAME,
      PERIOD == "Baseline-Kigali",
      year %in% 2021:2023
    ) %>%
    dplyr::mutate(group = "Baseline-Kigali"),

  classified %>%
    dplyr::filter(
      Site_Group == KFH_NAME,
      PERIOD == "Resurgence",
      year %in% 2024:2025
    ) %>%
    dplyr::mutate(group = "KFH Resurgence"),

  classified %>%
    dplyr::filter(
      PERIOD == "Resurgence",
      REGION == "Kigali",
      year %in% 2024:2025
    ) %>%
    dplyr::mutate(group = "Kigali Resurgence"),

  classified %>%
    dplyr::filter(
      PERIOD == "Resurgence",
      REGION == "Hotspot",
      year %in% 2024:2025
    ) %>%
    dplyr::mutate(group = "Hotspot Resurgence")

) %>%
  dplyr::mutate(group = factor(group, levels = GROUP_ORDER))

group_sizes <- group_membership %>%
  dplyr::count(group, name = "n_total")

add_to_report("")
add_to_report(strrep("-", 72))
add_to_report("GROUP SIZES")
add_to_report(strrep("-", 72))
add_to_report(sprintf("  %-22s %8s %12s", "group", "n", "published"))

size_check <- group_sizes %>%
  dplyr::left_join(
    PUBLISHED_GROUPS %>% dplyr::select(group, published_n),
    by = "group"
  )

for (index in seq_len(nrow(size_check))) {
  add_to_report(sprintf(
    "  %-22s %8d %12d",
    as.character(size_check$group[[index]]),
    size_check$n_total[[index]],
    size_check$published_n[[index]]
  ))
}

size_mismatches <- size_check %>% dplyr::filter(n_total != published_n)

add_to_report("")

if (nrow(size_mismatches) == 0) {
  add_to_report("  All four group sizes match published Figure 3b.")
} else {
  add_to_report(strrep("!", 72))
  add_to_report("  Group sizes do NOT match published Figure 3b.")
  add_to_report("  Check the Site_Group spelling, PERIOD and REGION values.")
  add_to_report(strrep("!", 72))
}


n_baseline  <- size_check$n_total[size_check$group == "Baseline-Kigali"]
n_kigali    <- size_check$n_total[size_check$group == "Kigali Resurgence"]
n_hotspot   <- size_check$n_total[size_check$group == "Hotspot Resurgence"]
n_kfh_resur <- size_check$n_total[size_check$group == "KFH Resurgence"]

kfh_inside_kigali <- group_membership %>%
  dplyr::filter(group == "KFH Resurgence") %>%
  dplyr::pull(SAMPLE) %>%
  setdiff(
    group_membership %>%
      dplyr::filter(group == "Kigali Resurgence") %>%
      dplyr::pull(SAMPLE)
  )

add_to_report("")
add_to_report("Consistency:")
add_to_report(sprintf(
  "  baseline + Kigali + Hotspot = %d + %d + %d = %d  (analysis set is %d)",
  n_baseline, n_kigali, n_hotspot,
  n_baseline + n_kigali + n_hotspot, nrow(sample_metadata)
))

if (length(kfh_inside_kigali) == 0) {
  add_to_report(sprintf(
    "  all %d KFH Resurgence samples are inside Kigali Resurgence",
    n_kfh_resur
  ))
} else {
  add_to_report(sprintf(
    "  WARNING: %d KFH Resurgence sample(s) are NOT in Kigali Resurgence",
    length(kfh_inside_kigali)
  ))
}


# ---- composition of each group -----------------------------------------

prevalence <- group_membership %>%
  dplyr::count(group, category, name = "n", .drop = FALSE) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(
    n_total    = sum(n),
    prevalence = 100 * n / n_total
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n_total > 0)

empty_categories <- prevalence %>%
  dplyr::group_by(category) %>%
  dplyr::summarise(total = sum(n), .groups = "drop") %>%
  dplyr::filter(total == 0) %>%
  dplyr::pull(category) %>%
  as.character()

add_to_report("")

if (length(empty_categories) > 0) {
  add_to_report(sprintf(
    "%d category/categories have no samples in any group and are not drawn:",
    length(empty_categories)
  ))
  add_to_report("  ", paste(empty_categories, collapse = ", "))
}

present_categories <- setdiff(CATEGORY_ORDER, empty_categories)

drawn <- prevalence %>%
  dplyr::filter(category %in% present_categories) %>%
  dplyr::mutate(category = factor(category, levels = present_categories))

add_to_report("")
add_to_report(strrep("-", 72))
add_to_report("COMPOSITION")
add_to_report(strrep("-", 72))

for (this_group in GROUP_ORDER) {

  rows <- drawn %>%
    dplyr::filter(group == this_group, n > 0) %>%
    dplyr::arrange(dplyr::desc(n))

  add_to_report("")
  add_to_report(sprintf(
    "  %s  (N = %d)", this_group, rows$n_total[[1]]
  ))

  for (index in seq_len(nrow(rows))) {
    add_to_report(sprintf(
      "    %-18s %6d  %6.1f%%",
      as.character(rows$category[[index]]),
      rows$n[[index]], rows$prevalence[[index]]
    ))
  }
}


# Check the two percentages printed on the published panel

published_check <- PUBLISHED_GROUPS %>%
  dplyr::left_join(
    drawn %>%
      dplyr::filter(category == "Wild type") %>%
      dplyr::select(group, computed_wild_type = prevalence) %>%
      dplyr::mutate(group = as.character(group)),
    by = "group"
  ) %>%
  dplyr::left_join(
    drawn %>%
      dplyr::filter(category == "R561H + K189T") %>%
      dplyr::select(group, computed_r561h_k189t = prevalence) %>%
      dplyr::mutate(group = as.character(group)),
    by = "group"
  ) %>%
  dplyr::mutate(
    wild_type_matches =
      !is.na(computed_wild_type) &
      abs(round(computed_wild_type, 1) - published_wild_type) < 0.051,
    r561h_k189t_matches =
      !is.na(computed_r561h_k189t) &
      abs(round(computed_r561h_k189t, 1) - published_r561h_k189t) < 0.051
  )

add_to_report("")
add_to_report(strrep("-", 72))
add_to_report("AGAINST THE PUBLISHED PANEL")
add_to_report(strrep("-", 72))
add_to_report(sprintf(
  "  %-20s %-18s %-18s", "group", "Wild type", "R561H + K189T"
))

for (index in seq_len(nrow(published_check))) {
  add_to_report(sprintf(
    "  %-20s %6.1f vs %-6.1f %s   %6.1f vs %-6.1f %s",
    published_check$group[[index]],
    published_check$computed_wild_type[[index]],
    published_check$published_wild_type[[index]],
    if (published_check$wild_type_matches[[index]]) "ok " else "DIFF",
    published_check$computed_r561h_k189t[[index]],
    published_check$published_r561h_k189t[[index]],
    if (published_check$r561h_k189t_matches[[index]]) "ok " else "DIFF"
  ))
}

all_match <- all(published_check$wild_type_matches) &&
  all(published_check$r561h_k189t_matches) &&
  nrow(size_mismatches) == 0

add_to_report("")

if (all_match) {
  add_to_report("  Figure 3b is reproduced: sizes and both labelled")
  add_to_report("  percentages match the published panel in every group.")
} else {
  add_to_report(strrep("!", 72))
  add_to_report("  This does NOT match published Figure 3b. Do not use the")
  add_to_report("  figure until the difference is understood.")
  add_to_report(strrep("!", 72))
}


# ---- the figure --------------------------------------------------------

panel_labels <- drawn %>%
  dplyr::distinct(group, n_total) %>%
  dplyr::mutate(
    facet_label = sprintf(
      "%s\n%s N = %d",
      as.character(group),
      GROUP_SUBTITLES[as.character(group)],
      n_total
    )
  )

facet_lookup <- setNames(panel_labels$facet_label,
                         as.character(panel_labels$group))

plot_data <- drawn %>%
  dplyr::mutate(
    facet_label = factor(
      facet_lookup[as.character(group)],
      levels = facet_lookup[GROUP_ORDER]
    )
  )

slice_positions <- plot_data %>%
  dplyr::group_by(facet_label) %>%
  dplyr::arrange(category, .by_group = TRUE) %>%
  dplyr::mutate(
    cumulative = cumsum(prevalence),
    midpoint   = cumulative - prevalence / 2
  ) %>%
  dplyr::ungroup()

slice_labels <- slice_positions %>%
  dplyr::filter(
    category %in% LABELLED_CATEGORIES,
    prevalence >= 5
  )

figure <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = "", y = prevalence, fill = category)
) +

  # reverse = TRUE makes the slices follow the category order, which the
  # label midpoints below assume.
  ggplot2::geom_col(
    width = 1, colour = "white", linewidth = 0.3,
    position = ggplot2::position_stack(reverse = TRUE)
  ) +

  # direction = 1 sweeps clockwise from twelve o'clock, matching the
  # published panel.
  ggplot2::coord_polar(theta = "y", direction = 1) +

  ggplot2::geom_text(
    data = slice_labels,
    mapping = ggplot2::aes(x = 1.08, y = midpoint, label = sprintf("%.1f", prevalence)),
    inherit.aes = FALSE,
    size = 3.6, fontface = "bold", colour = "black"
  ) +

  ggplot2::facet_wrap(~ facet_label, nrow = 1) +

  ggplot2::scale_fill_manual(
    values = CATEGORY_COLOURS[present_categories],
    drop = FALSE,
    name = "K13 mutation categories"
  ) +

  ggplot2::labs(
    x = NULL, y = NULL,
    title = "Distribution of K13 mutation in Kigali & Hotspots",
    subtitle = "(Prevalence %)"
  ) +

  ggplot2::guides(
    fill = ggplot2::guide_legend(nrow = 2, title.position = "top",
                                 title.hjust = 0.5)
  ) +

  ggplot2::theme_void(base_size = 12) +

  ggplot2::theme(
    plot.title    = ggplot2::element_text(face = "bold", hjust = 0.5,
                                          size = 14),
    plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12),
    strip.text    = ggplot2::element_text(size = 11, lineheight = 1.2,
                                          margin = ggplot2::margin(t = 6, b = 6)),
    legend.position = "bottom",
    legend.title  = ggplot2::element_text(face = "bold"),
    legend.box.background = ggplot2::element_rect(colour = "grey25",
                                                  fill = NA),
    legend.box.margin = ggplot2::margin(6, 6, 6, 6),
    plot.margin   = ggplot2::margin(12, 12, 12, 12)
  )

for (extension in c("svg", "png")) {
  ggplot2::ggsave(
    filename = file.path(
      FIGURE_DIR, paste0("figure03b_k13_category_pies.", extension)
    ),
    plot = figure,
    width = 14, height = 6.5, units = "in",
    dpi = 600, bg = "white"
  )
}


# ---- tables ------------------------------------------------------------

readr::write_tsv(
  prevalence %>%
    dplyr::mutate(
      group = as.character(group), category = as.character(category)
    ) %>%
    dplyr::select(group, category, n, n_total, prevalence),
  file.path(TABLE_DIR, "figure03b_category_prevalence.tsv")
)

readr::write_tsv(
  group_membership %>%
    dplyr::mutate(
      group = as.character(group), category = as.character(category)
    ) %>%
    dplyr::select(
      SAMPLE, group, year, Site_Group, District, PERIOD, REGION,
      k13_mutations, category
    ) %>%
    dplyr::arrange(group, SAMPLE),
  file.path(TABLE_DIR, "figure03b_sample_level_categories.tsv")
)

add_to_report("")
add_to_report(strrep("=", 72))
add_to_report("WRITTEN")
add_to_report(strrep("=", 72))
add_to_report("  figures/figure03/figure03b_k13_category_pies.svg / .png")
add_to_report("  tables/figure03/figure03b_category_prevalence.tsv")
add_to_report("  tables/figure03/figure03b_sample_level_categories.tsv")

writeLines(report_lines, file.path(TABLE_DIR, "figure03b_report.txt"))

capture.output(
  sessionInfo(),
  file = file.path(TABLE_DIR, "figure03b_sessionInfo.txt")
)
