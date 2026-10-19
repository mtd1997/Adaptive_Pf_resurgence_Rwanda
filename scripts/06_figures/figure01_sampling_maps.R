# Figure 1 — sampling maps
#
# Inputs:  data/metadata/samples_metadata.tsv   sample metadata (district, year,
#          period, region)
#          a Rwanda district shapefile          GADM level 2
#          four neighbouring-country files      GADM level 0: COD, UGA, TZA, BDI
# Outputs: results/figures/figure01/figure01_resurgence_map.svg / .pdf / .png
#          results/figures/figureS1/figureS1a_baseline_map.svg
#          results/figures/panels/figureS1a_baseline_map.rds
#          results/tables/figure01/figure01_samples_by_district.tsv
#          results/tables/figure01/figureS1a_samples_by_district.tsv
#          results/tables/figure01/who_malaria_rank_2025.tsv
#
# Usage:   From the top of the repository:

options(error = NULL)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(sf)
  library(ggspatial)
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


path_or_default <- function(variable_name, default_value) {

  if (exists(variable_name, inherits = TRUE)) {

    value <- get(variable_name, inherits = TRUE)

    if (is.character(value) && length(value) == 1 && nzchar(value)) {
      return(value)
    }
  }

  default_value
}

METADATA_FILE <- path_or_default(
  "METADATA_FILE",
  file.path(PROJECT_ROOT, "data", "metadata", "samples_metadata.tsv")
)

GADM_DIRECTORY <- path_or_default(
  "GADM_DIR",
  path.expand("~/Downloads")
)

FIGURE01_DIR <- file.path(PROJECT_ROOT, "results", "figures", "figure01")
FIGURES1_DIR <- file.path(PROJECT_ROOT, "results", "figures", "figureS1")
PANEL_DIR    <- file.path(PROJECT_ROOT, "results", "figures", "panels")
TABLE_DIR    <- file.path(PROJECT_ROOT, "results", "tables",  "figure01")

for (directory in c(FIGURE01_DIR, FIGURES1_DIR, PANEL_DIR, TABLE_DIR)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}


# ---- constants ---------------------------------------------------------

KIGALI_DISTRICTS <- c("Gasabo", "Kicukiro", "Nyarugenge")

BASELINE_PERIOD   <- "Baseline-Kigali"
RESURGENCE_PERIOD <- "Resurgence"

BASELINE_YEARS   <- 2021:2023
RESURGENCE_YEARS <- 2024:2025

REGION_COLOURS <- c(
  "Kigali"  = "#0072B2",
  "Hotspot" = "#D55E00"
)

ADMINISTRATIVE_COLOURS <- c(
  "Kigali City"     = "#D9EDF2",
  "Other districts" = "white"
)

RANK_COLOUR <- "#B2182B"

# WHO World malaria report 2025, Annex 4H. Global rank by estimated malaria
# cases in 2024.
WHO_RANKS <- tibble::tribble(
  ~country_code, ~map_label,   ~who_rank, ~map_rank_label, ~full_name,
  "COD",         "DR CONGO",   2L,        "2nd",           "Democratic Republic of the Congo",
  "UGA",         "UGANDA",     3L,        "3rd",           "Uganda",
  "TZA",         "TANZANIA",   7L,        "7th",           "United Republic of Tanzania",
  "BDI",         "BURUNDI",    19L,       "19th",          "Burundi"
)


# ---- read the sample metadata ------------------------------------------

if (!file.exists(METADATA_FILE)) {
  stop(
    "The sample metadata was not found:\n  ", METADATA_FILE,
    "\n\nBuild it first:\n",
    "  Rscript scripts/00_metadata/prepare_metadata.R --input <your file>",
    call. = FALSE
  )
}

metadata <- readr::read_tsv(
  METADATA_FILE,
  col_types = readr::cols(.default = readr::col_character()),
  show_col_types = FALSE
)


# Accept either capitalisation for the four columns used here. Metadata files
# in this project have used both over time.

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
  sample   = resolve_column(metadata, c("SAMPLE", "sample"), "sample ID"),
  period   = resolve_column(metadata, c("period", "PERIOD"), "period"),
  region   = resolve_column(metadata, c("region", "REGION"), "region"),
  district = resolve_column(metadata, c("district", "District", "DISTRICT"), "district"),
  year     = resolve_column(metadata, c("year", "YEAR"), "year")
)

sample_metadata <- metadata %>%
  dplyr::select(dplyr::all_of(unname(column_map)))

names(sample_metadata) <- names(column_map)

sample_metadata <- sample_metadata %>%
  dplyr::mutate(
    dplyr::across(dplyr::everything(), ~ stringr::str_trim(.x)),
    year = suppressWarnings(as.integer(year))
  ) %>%
  dplyr::filter(!is.na(sample), sample != "") %>%
  dplyr::distinct(sample, .keep_all = TRUE)

message(sprintf(
  "Metadata: %d samples\n", nrow(sample_metadata)
))

message("Samples per period and region:")
print(table(sample_metadata$period, sample_metadata$region, useNA = "ifany"))


# ---- read the shapefiles -----------------------------------------------

find_gadm_shapefile <- function(country_code, level) {

  pattern <- sprintf("^gadm41_%s_%d\\.shp$", country_code, level)

  matches <- list.files(
    path = GADM_DIRECTORY,
    pattern = pattern,
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(matches) == 0) {
    stop(
      "Could not find ", sprintf("gadm41_%s_%d.shp", country_code, level),
      " under:\n  ", GADM_DIRECTORY,
      "\n\nDownload it from https://gadm.org/download_country.html, unzip it ",
      "there, or set GADM_DIR in config.R to wherever the files are.",
      call. = FALSE
    )
  }

  if (length(matches) > 1) {
    warning(
      "Several copies of gadm41_", country_code, "_", level,
      ".shp were found; using:\n  ", matches[[1]],
      call. = FALSE
    )
  }

  matches[[1]]
}

rwanda_districts <- sf::st_read(
  find_gadm_shapefile("RWA", 2),
  quiet = TRUE
)

# GADM calls the district column NAME_2.
names(rwanda_districts)[names(rwanda_districts) == "NAME_2"] <- "district"

rwanda_districts <- rwanda_districts %>%
  dplyr::mutate(
    administrative_area = dplyr::if_else(
      district %in% KIGALI_DISTRICTS,
      "Kigali City",
      "Other districts"
    )
  )


# Check the district names in the metadata match the shapefile, because a
# mismatch silently drops a district's samples from the map rather than
# failing.

shapefile_districts <- sort(unique(rwanda_districts$district))
metadata_districts  <- sort(unique(stats::na.omit(sample_metadata$district)))

unmatched_districts <- setdiff(metadata_districts, shapefile_districts)

if (length(unmatched_districts) > 0) {
  warning(
    length(unmatched_districts),
    " district name(s) in the metadata are not in the shapefile, so their ",
    "samples will not appear on the maps: ",
    paste(unmatched_districts, collapse = ", "),
    call. = FALSE
  )
}


# A point inside each district, for placing the sample circles and labels.

district_points <- sf::st_point_on_surface(rwanda_districts)

district_coordinates <- sf::st_coordinates(district_points)

district_points <- district_points %>%
  dplyr::mutate(
    longitude = district_coordinates[, 1],
    latitude  = district_coordinates[, 2]
  ) %>%
  sf::st_drop_geometry() %>%
  dplyr::select(district, longitude, latitude)


# ---- count samples per district ----------------------------------------

count_by_district <- function(samples) {

  samples %>%
    dplyr::count(district, name = "sample_count") %>%
    dplyr::left_join(district_points, by = "district")
}

size_breaks <- function(counts, n = 4) {

  breaks <- pretty(counts, n = n)
  breaks[breaks > 0 & breaks <= max(counts, na.rm = TRUE)]
}


baseline_samples <- sample_metadata %>%
  dplyr::filter(
    period == BASELINE_PERIOD,
    region == "Kigali",
    year %in% BASELINE_YEARS
  )

resurgence_samples <- sample_metadata %>%
  dplyr::filter(
    period == RESURGENCE_PERIOD,
    year %in% RESURGENCE_YEARS
  )

baseline_counts   <- count_by_district(baseline_samples)
resurgence_counts <- count_by_district(resurgence_samples)

resurgence_counts <- resurgence_counts %>%
  dplyr::left_join(
    dplyr::distinct(resurgence_samples, district, region),
    by = "district"
  )

baseline_breaks   <- size_breaks(baseline_counts$sample_count, n = 4)
resurgence_breaks <- size_breaks(resurgence_counts$sample_count, n = 5)

message(sprintf(
  "\nBaseline-Kigali 2021-2023: %d samples across %d district(s)",
  sum(baseline_counts$sample_count), nrow(baseline_counts)
))

message(sprintf(
  "Resurgence 2024-2025:      %d samples across %d district(s)",
  sum(resurgence_counts$sample_count), nrow(resurgence_counts)
))

missing_coordinates <- dplyr::filter(
  dplyr::bind_rows(baseline_counts, resurgence_counts),
  is.na(longitude)
)

if (nrow(missing_coordinates) > 0) {
  warning(
    "These districts have samples but no coordinates, so they will not be ",
    "plotted: ",
    paste(unique(missing_coordinates$district), collapse = ", "),
    call. = FALSE
  )
}


# ---- baseline map  (supplementary figure s1, panel a) ------------------

baseline_map <- ggplot2::ggplot() +

  ggplot2::geom_sf(
    data = rwanda_districts,
    fill = "grey95", colour = "grey65", linewidth = 0.4
  ) +

  ggplot2::geom_point(
    data = baseline_counts,
    ggplot2::aes(x = longitude, y = latitude, size = sample_count),
    colour = "#0072B2", alpha = 0.85
  ) +

  ggrepel::geom_text_repel(
    data = baseline_counts,
    ggplot2::aes(x = longitude, y = latitude, label = district),
    size = 4.2, fontface = "bold",
    box.padding = 0.5, point.padding = 0.3,
    max.overlaps = Inf,
    seed = 12345
  ) +

  ggplot2::scale_size_continuous(
    name = "Number of samples",
    range = c(3, 9),
    breaks = baseline_breaks
  ) +

  ggspatial::annotation_scale(location = "bl", width_hint = 0.25) +

  ggspatial::annotation_north_arrow(
    location = "tl",
    which_north = "true",
    style = ggspatial::north_arrow_fancy_orienteering
  ) +

  ggplot2::theme_minimal(base_size = 16) +

  ggplot2::theme(
    legend.title     = ggplot2::element_text(size = 14, face = "bold"),
    legend.text      = ggplot2::element_text(size = 13),
    legend.position  = "right",
    axis.title       = ggplot2::element_blank(),
    panel.grid       = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank()
  )


# ---- neighbouring countries --------------------------------------------

neighbour_shapes <- lapply(
  seq_len(nrow(WHO_RANKS)),
  function(index) {

    row <- WHO_RANKS[index, ]

    shape <- sf::st_read(
      find_gadm_shapefile(row$country_code, 0),
      quiet = TRUE
    )

    shape <- sf::st_transform(shape, sf::st_crs(rwanda_districts))

    shape %>%
      dplyr::mutate(
        country_label = row$map_label,
        rank_label    = row$map_rank_label,
        who_rank      = row$who_rank
      )
  }
)

neighbours <- dplyr::bind_rows(neighbour_shapes)


# Map window: Rwanda's bounding box with a margin, so the neighbours show as
# context rather than in full.

rwanda_bbox <- sf::st_bbox(rwanda_districts)

LONGITUDE_PADDING <- 0.75
LATITUDE_PADDING  <- 0.65

map_xmin <- as.numeric(rwanda_bbox["xmin"]) - LONGITUDE_PADDING
map_xmax <- as.numeric(rwanda_bbox["xmax"]) + LONGITUDE_PADDING
map_ymin <- as.numeric(rwanda_bbox["ymin"]) - LATITUDE_PADDING
map_ymax <- as.numeric(rwanda_bbox["ymax"]) + LATITUDE_PADDING

neighbours_cropped <- suppressWarnings(
  sf::st_crop(
    neighbours,
    xmin = map_xmin, xmax = map_xmax,
    ymin = map_ymin, ymax = map_ymax
  )
)


# Label positions for the neighbours.

country_label_points <- sf::st_point_on_surface(neighbours_cropped)

country_label_coordinates <- sf::st_coordinates(country_label_points)

country_labels <- country_label_points %>%
  dplyr::mutate(
    label_longitude = country_label_coordinates[, 1],
    label_latitude  = country_label_coordinates[, 2],
    label_latitude_adjusted = dplyr::case_when(
      country_label == "UGANDA"  ~ label_latitude + 0.07,
      country_label == "BURUNDI" ~ label_latitude - 0.08,
      TRUE                       ~ label_latitude
    )
  ) %>%
  sf::st_drop_geometry()

rwanda_outline <- sf::st_union(rwanda_districts)

# The WHO rank key sits in the lower right of the surrounding area. Drawn as
# an annotation rather than a legend, because a legend for a text annotation
# would need an invented aesthetic.
key_x <- map_xmax - 0.08
key_y <- map_ymin + 0.12


# ---- resurgence map  (figure 1) ----------------------------------------

resurgence_map <- ggplot2::ggplot() +

  ggplot2::geom_sf(
    data = neighbours_cropped,
    fill = "grey92", colour = "black", linewidth = 0.6
  ) +

  ggplot2::geom_sf(
    data = rwanda_districts,
    ggplot2::aes(fill = administrative_area),
    colour = "grey65", linewidth = 0.4
  ) +

  ggplot2::geom_sf(
    data = rwanda_outline,
    fill = NA, colour = "black", linewidth = 0.6
  ) +

  ggplot2::geom_text(
    data = country_labels,
    ggplot2::aes(
      x = label_longitude, y = label_latitude_adjusted, label = country_label
    ),
    colour = "grey20", size = 4.2, fontface = "bold", vjust = 0.8
  ) +

  ggplot2::geom_text(
    data = country_labels,
    ggplot2::aes(
      x = label_longitude, y = label_latitude_adjusted, label = rank_label
    ),
    colour = RANK_COLOUR, size = 5.2, fontface = "bold", vjust = -0.8
  ) +

  ggplot2::geom_point(
    data = resurgence_counts,
    ggplot2::aes(
      x = longitude, y = latitude, size = sample_count, colour = region
    ),
    alpha = 0.85
  ) +

  ggrepel::geom_text_repel(
    data = resurgence_counts,
    ggplot2::aes(x = longitude, y = latitude, label = district),
    size = 4.2, fontface = "bold",
    box.padding = 0.5, point.padding = 0.3,
    max.overlaps = Inf,
    seed = 12345
  ) +

  ggplot2::annotate(
    "text", x = key_x, y = key_y + 0.12,
    label = "WHO 2025 malaria rank",
    hjust = 1, vjust = 0,
    colour = RANK_COLOUR, size = 4.1, fontface = "bold"
  ) +

  ggplot2::annotate(
    "text", x = key_x, y = key_y,
    label = "Global rank by estimated\nmalaria cases (2024)",
    hjust = 1, vjust = 0,
    colour = "grey20", size = 3.5, lineheight = 1.05
  ) +

  ggplot2::scale_fill_manual(
    name = "Administrative area", values = ADMINISTRATIVE_COLOURS
  ) +

  ggplot2::scale_colour_manual(
    name = "Region", values = REGION_COLOURS
  ) +

  ggplot2::scale_size_continuous(
    name = "Number of samples", range = c(4, 10), breaks = resurgence_breaks
  ) +

  ggplot2::guides(
    fill   = ggplot2::guide_legend(order = 1),
    colour = ggplot2::guide_legend(order = 2),
    size   = ggplot2::guide_legend(order = 3)
  ) +

  ggspatial::annotation_scale(location = "bl", width_hint = 0.25) +

  ggspatial::annotation_north_arrow(
    location = "tl",
    which_north = "true",
    style = ggspatial::north_arrow_fancy_orienteering
  ) +

  ggplot2::coord_sf(
    xlim = c(map_xmin, map_xmax),
    ylim = c(map_ymin, map_ymax),
    expand = FALSE
  ) +

  ggplot2::theme_minimal(base_size = 16) +

  ggplot2::theme(
    legend.title     = ggplot2::element_text(size = 14, face = "bold"),
    legend.text      = ggplot2::element_text(size = 13),
    legend.position  = "right",
    legend.box       = "vertical",
    axis.title       = ggplot2::element_blank(),
    panel.grid       = ggplot2::element_blank(),
    panel.background = ggplot2::element_rect(fill = "white", colour = NA)
  )


# ---- save --------------------------------------------------------------

for (extension in c("svg", "pdf", "png")) {

  ggplot2::ggsave(
    filename = file.path(
      FIGURE01_DIR, paste0("figure01_resurgence_map.", extension)
    ),
    plot = resurgence_map,
    width = 14, height = 10, units = "in",
    dpi = 600, bg = "white", limitsize = FALSE
  )
}

ggplot2::ggsave(
  filename = file.path(FIGURES1_DIR, "figureS1a_baseline_map.svg"),
  plot = baseline_map,
  width = 7.5, height = 8, units = "in",
  bg = "white"
)


# Hand the baseline panel on to the Supplementary Figure S1 assembly.

saveRDS(
  list(
    panel  = baseline_map,
    breaks = baseline_breaks,
    n      = sum(baseline_counts$sample_count)
  ),
  file.path(PANEL_DIR, "figureS1a_baseline_map.rds")
)


# ---- source data -------------------------------------------------------

readr::write_tsv(
  dplyr::arrange(resurgence_counts, dplyr::desc(sample_count)),
  file.path(TABLE_DIR, "figure01_samples_by_district.tsv")
)

readr::write_tsv(
  dplyr::arrange(baseline_counts, dplyr::desc(sample_count)),
  file.path(TABLE_DIR, "figureS1a_samples_by_district.tsv")
)

readr::write_tsv(
  dplyr::transmute(
    WHO_RANKS,
    country = full_name,
    rank_estimated_cases_2024 = who_rank,
    map_label = map_rank_label,
    source = "WHO World malaria report 2025, Annex 4H"
  ),
  file.path(TABLE_DIR, "who_malaria_rank_2025.tsv")
)

readr::write_tsv(
  dplyr::bind_rows(
    dplyr::mutate(baseline_samples, map = "Supplementary Figure S1a"),
    dplyr::mutate(resurgence_samples, map = "Figure 1")
  ),
  file.path(TABLE_DIR, "figure01_sample_level.tsv")
)

capture.output(
  sessionInfo(),
  file = file.path(TABLE_DIR, "sessionInfo.txt")
)


# ---- summary -----------------------------------------------------------

message("\n", strrep("=", 72))
message("SAMPLING MAPS WRITTEN")
message(strrep("=", 72))
message("Figure 1:")
message(sprintf("  %s", file.path(FIGURE01_DIR, "figure01_resurgence_map.svg")))
message("Supplementary Figure S1, panel A:")
message(sprintf("  %s", file.path(FIGURES1_DIR, "figureS1a_baseline_map.svg")))
message("\nPanel saved for the S1 assembly:")
message(sprintf("  %s", file.path(PANEL_DIR, "figureS1a_baseline_map.rds")))

message("\nResurgence samples per district:")
print(
  as.data.frame(
    dplyr::arrange(
      dplyr::select(resurgence_counts, district, region, sample_count),
      dplyr::desc(sample_count)
    )
  ),
  row.names = FALSE
)

message("\nBaseline samples per district:")
print(
  as.data.frame(
    dplyr::arrange(
      dplyr::select(baseline_counts, district, sample_count),
      dplyr::desc(sample_count)
    )
  ),
  row.names = FALSE
)

message("\nNext step for Supplementary Figure S1:")
message("  Rscript scripts/06_figures/figureS1b_malaria_cases.R")
message(strrep("=", 72))
