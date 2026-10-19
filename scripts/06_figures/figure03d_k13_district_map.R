# Figure 3d — K13 genotypes by district
#
# Inputs:  Manuscript_data_29-07-2026.tsv        (MANUSCRIPT_DATA in config.R)
#          gadm41_RWA_2.shp and companions       (GADM_DIR in config.R)
# Outputs: results/figures/figure03/figure03d_k13_district_map.svg / .png
#          results/tables/figure03/figure03d_district_genotypes.tsv   <- Source Data
#          results/tables/figure03/figure03d_pie_coordinates.tsv
#          results/tables/figure03/figure03d_report.txt
#
# Usage:   Rscript scripts/06_figures/figure03d_k13_district_map.R

options(error = NULL)

# sf and scatterpie are not part of the tidyverse and are easy to be without.
# Check for them by name so a missing one produces an instruction rather than
# a bare "there is no package called ..." part way down the script.
required_packages <- c(
  "dplyr", "tidyr", "readr", "stringr", "tibble", "ggplot2",
  "sf", "scatterpie"
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
  library(ggplot2)
  library(sf)
  library(scatterpie)
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

SHAPEFILE <- get_argument(
  "--shapefile",
  if (exists("GADM_DIR", inherits = TRUE)) {
    file.path(
      path.expand(get("GADM_DIR", inherits = TRUE)),
      "gadm41_RWA_shp", "gadm41_RWA_2.shp"
    )
  } else NULL
)

if (is.null(SHAPEFILE) || !file.exists(SHAPEFILE)) {
  stop(
    "Rwanda district shapefile not found:\n  ",
    if (is.null(SHAPEFILE)) "(GADM_DIR not set in config.R)" else SHAPEFILE,
    "\n\nSet GADM_DIR in config.R to the folder holding gadm41_RWA_shp/,",
    "\nor pass --shapefile /path/to/gadm41_RWA_2.shp",
    call. = FALSE
  )
}

FIGURE_DIR <- file.path(PROJECT_ROOT, "results", "figures", "figure03")
TABLE_DIR  <- file.path(PROJECT_ROOT, "results", "tables",  "figure03")

dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR,  recursive = TRUE, showWarnings = FALSE)


GENOTYPE_ORDER <- c(
  "Wild type", "K189T only", "R561H + K189T", "R561H only", "Other K13"
)

GENOTYPE_COLOURS <- c(
  "Wild type"     = "#2ECC71",
  "K189T only"    = "#D62F7D",
  "R561H + K189T" = "#FF7F0E",
  "R561H only"    = "#F4B400",
  "Other K13"     = "#6A5ACD"
)

# Pie radius in degrees for the district with the most samples. Everything
# else is scaled down from this by sqrt(n / n_max).
MAXIMUM_PIE_RADIUS <- 0.16
MINIMUM_PIE_RADIUS <- 0.045


# Where a pie is drawn, when its centroid will not do.

DISTRICT_PIE_OFFSETS <- tibble::tribble(
  ~district,      ~offset_longitude, ~offset_latitude,
  "Gasabo",        0.10,              0.07,
  "Nyarugenge",   -0.20,             -0.11,
  "Kicukiro",      0.24,             -0.12
)


report_lines <- character(0)

add_to_report <- function(...) {
  line <- paste0(...)
  report_lines <<- c(report_lines, line)
  message(line)
}


# ---- read the samples --------------------------------------------------

add_to_report(strrep("=", 72))
add_to_report("FIGURE 3D — K13 GENOTYPES BY DISTRICT")
add_to_report(strrep("=", 72))
add_to_report("Samples:   ", INPUT_FILE)
add_to_report("Districts: ", SHAPEFILE)
add_to_report("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
add_to_report("")

raw <- readr::read_tsv(
  INPUT_FILE,
  col_types = readr::cols(.default = readr::col_character()),
  show_col_types = FALSE
)

required_columns <- c("SAMPLE", "District", "GENE_NAME", "HGVS_PROTEIN")

missing_columns <- setdiff(required_columns, names(raw))

if (length(missing_columns) > 0) {
  stop("The input is missing: ", paste(missing_columns, collapse = ", "),
       call. = FALSE)
}

tidy <- raw %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(required_columns), ~ stringr::str_trim(.x)
    )
  ) %>%
  dplyr::filter(!is.na(SAMPLE), SAMPLE != "")


# ---- classify each sample ----------------------------------------------

k13_flags <- tidy %>%
  dplyr::filter(
    GENE_NAME == "K13",
    !is.na(HGVS_PROTEIN), HGVS_PROTEIN != "", HGVS_PROTEIN != "."
  ) %>%
  dplyr::group_by(SAMPLE) %>%
  dplyr::summarise(
    has_r561h  = any(HGVS_PROTEIN == "R561H"),
    has_k189t  = any(HGVS_PROTEIN == "K189T"),
    has_any_k13 = TRUE,
    .groups = "drop"
  )

sample_districts <- tidy %>%
  dplyr::distinct(SAMPLE, District)

# A sample must sit in exactly one district, or it would be counted twice.
duplicated_districts <- sample_districts %>%
  dplyr::count(SAMPLE, name = "n_districts") %>%
  dplyr::filter(n_districts > 1)

if (nrow(duplicated_districts) > 0) {
  stop(
    nrow(duplicated_districts), " sample(s) are listed under more than one ",
    "district, so they would be counted more than once:\n  ",
    paste(utils::head(duplicated_districts$SAMPLE, 10), collapse = ", "),
    call. = FALSE
  )
}

classified <- sample_districts %>%
  dplyr::left_join(k13_flags, by = "SAMPLE") %>%
  dplyr::mutate(
    has_r561h   = tidyr::replace_na(has_r561h, FALSE),
    has_k189t   = tidyr::replace_na(has_k189t, FALSE),
    has_any_k13 = tidyr::replace_na(has_any_k13, FALSE),
    genotype = dplyr::case_when(
      has_r561h &  has_k189t ~ "R561H + K189T",
      has_r561h & !has_k189t ~ "R561H only",
      !has_r561h & has_k189t ~ "K189T only",
      has_any_k13            ~ "Other K13",
      TRUE                   ~ "Wild type"
    ),
    genotype = factor(genotype, levels = GENOTYPE_ORDER)
  )

add_to_report(sprintf("Samples classified: %d", nrow(classified)))
add_to_report("")
add_to_report("Genotype totals across all districts:")

overall <- classified %>%
  dplyr::count(genotype, name = "n", .drop = FALSE) %>%
  dplyr::mutate(percent = 100 * n / sum(n))

for (index in seq_len(nrow(overall))) {
  add_to_report(sprintf(
    "  %-16s %6d  %5.1f%%",
    as.character(overall$genotype[[index]]),
    overall$n[[index]], overall$percent[[index]]
  ))
}


# ---- per-district counts -----------------------------------------------

district_counts <- classified %>%
  dplyr::filter(!is.na(District), District != "") %>%
  dplyr::count(District, genotype, name = "n", .drop = FALSE) %>%
  tidyr::pivot_wider(
    names_from = genotype, values_from = n, values_fill = 0
  ) %>%
  dplyr::mutate(
    n_total = rowSums(dplyr::across(dplyr::all_of(GENOTYPE_ORDER)))
  ) %>%
  dplyr::filter(n_total > 0) %>%
  dplyr::arrange(dplyr::desc(n_total))

missing_district <- classified %>%
  dplyr::filter(is.na(District) | District == "") %>%
  nrow()

if (missing_district > 0) {
  add_to_report("")
  add_to_report(sprintf(
    "%d sample(s) have no district and are not on the map.", missing_district
  ))
}

add_to_report("")
add_to_report(strrep("-", 72))
add_to_report("SAMPLES BY DISTRICT")
add_to_report(strrep("-", 72))

for (index in seq_len(nrow(district_counts))) {
  add_to_report(sprintf(
    "  %-18s %5d",
    district_counts$District[[index]], district_counts$n_total[[index]]
  ))
}

add_to_report(sprintf(
  "  %-18s %5d", "total on map", sum(district_counts$n_total)
))


# ---- read the district boundaries --------------------------------------

districts <- sf::st_read(SHAPEFILE, quiet = TRUE)

if (!"NAME_2" %in% names(districts)) {
  stop(
    "The shapefile has no NAME_2 column, so districts cannot be matched.\n",
    "Columns present: ", paste(names(districts), collapse = ", "),
    call. = FALSE
  )
}


# Names have to agree between the sample table and the shapefile. A silent
# mismatch would simply drop a district's pie, which is why this is checked
# rather than left to the join.

shapefile_names <- sort(unique(as.character(districts$NAME_2)))
sample_names    <- sort(unique(district_counts$District))

unmatched <- setdiff(sample_names, shapefile_names)

add_to_report("")

if (length(unmatched) == 0) {
  add_to_report(sprintf(
    "District names: all %d matched the shapefile.", length(sample_names)
  ))
} else {
  add_to_report(strrep("!", 72))
  add_to_report(sprintf(
    "%d district name(s) in the sample table are not in the shapefile and",
    length(unmatched)
  ))
  add_to_report("their samples will NOT appear on the map:")
  for (name in unmatched) {
    affected <- district_counts$n_total[district_counts$District == name]
    add_to_report(sprintf("  %-20s %d sample(s)", name, affected))
  }
  add_to_report("")
  add_to_report("Shapefile district names are:")
  add_to_report("  ", paste(shapefile_names, collapse = ", "))
  add_to_report(strrep("!", 72))
}


# ---- centroids and pie sizes -------------------------------------------

map_data <- districts %>%
  dplyr::left_join(
    district_counts, by = c("NAME_2" = "District")
  ) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(c(GENOTYPE_ORDER, "n_total")),
      ~ tidyr::replace_na(.x, 0)
    )
  )

centroid_coordinates <- suppressWarnings(
  sf::st_coordinates(sf::st_centroid(sf::st_geometry(map_data)))
)

map_data$longitude <- centroid_coordinates[, 1]
map_data$latitude  <- centroid_coordinates[, 2]

pie_data <- map_data %>%
  sf::st_drop_geometry() %>%
  dplyr::filter(n_total > 0) %>%
  dplyr::mutate(
    # Area, not radius, carries the count.
    radius = pmax(
      MINIMUM_PIE_RADIUS,
      MAXIMUM_PIE_RADIUS * sqrt(n_total / max(n_total))
    )
  ) %>%
  dplyr::select(
    district = NAME_2, longitude, latitude, radius, n_total,
    dplyr::all_of(GENOTYPE_ORDER)
  ) %>%
  dplyr::left_join(DISTRICT_PIE_OFFSETS, by = "district") %>%
  dplyr::mutate(
    offset_longitude = tidyr::replace_na(offset_longitude, 0),
    offset_latitude  = tidyr::replace_na(offset_latitude, 0),
    pie_longitude    = longitude + offset_longitude,
    pie_latitude     = latitude + offset_latitude,
    displaced        = offset_longitude != 0 | offset_latitude != 0
  )

if (nrow(pie_data) == 0) {
  stop("No district has any samples after matching; nothing to draw.",
       call. = FALSE)
}


# An offset that is too small leaves the pies overlapping anyway, and one
# that is too large pushes a pie across the map. Neither is visible from the
# table above, so check the drawn geometry directly.

overlapping_pairs <- character(0)

if (nrow(pie_data) > 1) {

  for (first in seq_len(nrow(pie_data) - 1)) {
    for (second in seq(first + 1, nrow(pie_data))) {

      separation <- sqrt(
        (pie_data$pie_longitude[[first]] - pie_data$pie_longitude[[second]])^2 +
        (pie_data$pie_latitude[[first]]  - pie_data$pie_latitude[[second]])^2
      )

      needed <- pie_data$radius[[first]] + pie_data$radius[[second]]

      if (separation < needed) {
        overlapping_pairs <- c(
          overlapping_pairs,
          sprintf(
            "%s and %s overlap by %.3f degrees",
            pie_data$district[[first]], pie_data$district[[second]],
            needed - separation
          )
        )
      }
    }
  }
}

add_to_report("")

if (length(overlapping_pairs) == 0) {
  add_to_report("Pie placement: no two pies overlap.")
} else {
  add_to_report(sprintf("Pie placement: %d overlapping pair(s):",
                        length(overlapping_pairs)))
  for (pair in overlapping_pairs) add_to_report("  ", pair)
  add_to_report("  Adjust DISTRICT_PIE_OFFSETS or MAXIMUM_PIE_RADIUS.")
}

if (any(pie_data$displaced)) {
  add_to_report("")
  add_to_report("Displaced from their centroid, with a connecting line drawn:")
  for (index in which(pie_data$displaced)) {
    add_to_report(sprintf(
      "  %-14s %+.2f lon  %+.2f lat",
      pie_data$district[[index]],
      pie_data$offset_longitude[[index]],
      pie_data$offset_latitude[[index]]
    ))
  }
}

add_to_report("")
add_to_report(sprintf(
  "Districts drawn: %d, covering %d sample(s).",
  nrow(pie_data), sum(pie_data$n_total)
))


# ---- the figure --------------------------------------------------------

figure <- ggplot2::ggplot() +

  ggplot2::geom_sf(
    data = map_data,
    fill = "#EEF2F7", colour = "grey55", linewidth = 0.3
  ) +

  # Leader lines first, so the pies cover their inner ends.
  ggplot2::geom_segment(
    data = pie_data[pie_data$displaced, ],
    mapping = ggplot2::aes(
      x = longitude, y = latitude, xend = pie_longitude, yend = pie_latitude
    ),
    colour = "grey45", linewidth = 0.4
  ) +

  ggplot2::geom_point(
    data = pie_data[pie_data$displaced, ],
    mapping = ggplot2::aes(x = longitude, y = latitude),
    size = 1.1, colour = "grey25"
  ) +

  scatterpie::geom_scatterpie(
    data = pie_data,
    mapping = ggplot2::aes(x = pie_longitude, y = pie_latitude, r = radius),
    cols = GENOTYPE_ORDER,
    colour = "white", linewidth = 0.25
  ) +

  ggplot2::geom_text(
    data = pie_data,
    mapping = ggplot2::aes(
      x = pie_longitude, y = pie_latitude + radius + 0.035, label = district
    ),
    size = 3.1, fontface = "bold", colour = "black"
  ) +

  ggplot2::scale_fill_manual(
    values = GENOTYPE_COLOURS,
    breaks = GENOTYPE_ORDER,
    name = "K13 mutation category"
  ) +

  ggplot2::coord_sf(expand = TRUE) +

  ggplot2::labs(x = NULL, y = NULL) +

  ggplot2::theme_void(base_size = 12)


# Placing the legend inside the panel changed syntax in ggplot2 3.5.0: the
# old legend.position = c(x, y) is deprecated there and removed later. Both
# spellings are kept so the script works on either side of that change.

legend_inside <- c(0.86, 0.16)

if (utils::packageVersion("ggplot2") >= "3.5.0") {
  figure <- figure + ggplot2::theme(
    legend.position = "inside",
    legend.position.inside = legend_inside
  )
} else {
  figure <- figure + ggplot2::theme(legend.position = legend_inside)
}

figure <- figure + ggplot2::theme(
  legend.title      = ggplot2::element_text(face = "bold", size = 11),
  legend.text       = ggplot2::element_text(size = 10),
  legend.background = ggplot2::element_rect(fill = "white", colour = NA),
  legend.key.size   = ggplot2::unit(0.5, "cm"),
  plot.margin       = ggplot2::margin(10, 10, 10, 10)
)

for (extension in c("svg", "png")) {
  ggplot2::ggsave(
    filename = file.path(
      FIGURE_DIR, paste0("figure03d_k13_district_map.", extension)
    ),
    plot = figure,
    width = 10, height = 10, units = "in",
    dpi = 600, bg = "white"
  )
}


# ---- tables ------------------------------------------------------------

readr::write_tsv(
  district_counts,
  file.path(TABLE_DIR, "figure03d_district_genotypes.tsv")
)

readr::write_tsv(
  pie_data,
  file.path(TABLE_DIR, "figure03d_pie_coordinates.tsv")
)

add_to_report("")
add_to_report(strrep("=", 72))
add_to_report("WRITTEN")
add_to_report(strrep("=", 72))
add_to_report("  figures/figure03/figure03d_k13_district_map.svg / .png")
add_to_report("  tables/figure03/figure03d_district_genotypes.tsv")
add_to_report("  tables/figure03/figure03d_pie_coordinates.tsv")
add_to_report("")
add_to_report("Reminder: District is the patient's residence, not the")
add_to_report("collection site. The figure legend should say so.")

writeLines(report_lines, file.path(TABLE_DIR, "figure03d_report.txt"))

capture.output(
  sessionInfo(),
  file = file.path(TABLE_DIR, "figure03d_sessionInfo.txt")
)
