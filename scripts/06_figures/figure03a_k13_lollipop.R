################################################################################
# FIGURE 3A — K13 MUTATION SPECTRUM ALONG THE PROTEIN
################################################################################
#
# WHAT IT DOES
#   Counts how many samples carry each K13 amino-acid change, places each
#   change at its codon position, and draws the spectrum as a lollipop plot
#   above a track showing the four K13 protein domains.
#
# INPUT
#   Manuscript_data_29-07-2026.tsv  (MANUSCRIPT_DATA in config.R)
#
# OUTPUTS
#   results/figures/figure03/figure03a_k13_lollipop.svg / .png
#   results/tables/figure03/figure03a_mutation_frequencies.tsv   <- Source Data
#   results/tables/figure03/figure03a_domains.tsv
#   results/tables/figure03/figure03a_excluded_singletons.tsv
#   results/tables/figure03/figure03a_report.txt
#
# -----------------------------------------------------------------------------
# THE DENOMINATOR IS NOT THE ANALYSIS SET
# -----------------------------------------------------------------------------
#   Frequency here is the percentage of K13-MUTATED samples carrying each
#   change, not the percentage of all sequenced samples. A mutation found in
#   half of the K13-mutant samples is therefore shown at 50%, even though it
#   is present in a much smaller share of the cohort.
#
#   This is the original definition and it is kept. It is flagged because
#   the published axis is titled only "Mutation Frequency (%)", so the
#   denominator is not visible to a reader of the figure. Both numbers are
#   printed below and both are written to the Source Data table, so the
#   figure legend can state which one it means.
# -----------------------------------------------------------------------------
#
# -----------------------------------------------------------------------------
# SINGLETONS ARE EXCLUDED
# -----------------------------------------------------------------------------
#   Changes seen in exactly one sample are dropped, as in the original. They
#   are written to figure03a_excluded_singletons.tsv rather than discarded,
#   so the figure's omissions are recoverable.
#
#   The original did this with `Mutated_Sample_Count != "1"` -- a count
#   compared against a string. R coerces the number to text, so it happens
#   to give the right answer, but it would also drop a count of 1 written
#   any other way and is fragile. The comparison here is numeric.
# -----------------------------------------------------------------------------
#
# -----------------------------------------------------------------------------
# HOW THE LABELS ARE PLACED
# -----------------------------------------------------------------------------
#   Nearly every change sits below 10%, so the labels, not the lollipops,
#   are what has to be laid out. Three things follow from the labels being
#   rotated upright:
#
#     - A rotated label is tall and narrow. Two labels collide when they are
#       close in CODON position, not when they are close in height. So
#       crowding is resolved by sliding labels apart along the protein and
#       drawing a short leader back to the true position -- not by stacking
#       them higher, which does nothing for two labels at the same codon.
#       K189T and K189N are both at codon 189 and would otherwise be drawn
#       on top of each other.
#
#     - Labels sit low, just above the backbone, because that is where the
#       data is. Pushing them high leaves a band of empty plot between the
#       labels and the two tall lollipops.
#
#     - Heights still step up slightly within a cluster of nearby changes,
#       which is what gives the published panel its staircase look.
#
#   The placement is computed, not repelled: ggrepel is seeded, so the
#   figure would shift between runs of the same script, and it cannot be
#   told to keep a label off the K189T or R561H stem.
# -----------------------------------------------------------------------------
#
# -----------------------------------------------------------------------------
# WHAT IS NOT REPRODUCED
# -----------------------------------------------------------------------------
#   Published Figure 3a has a thin box drawn around the P574L label. Nothing
#   in the analysis distinguishes P574L from its neighbours, so the box
#   appears to be a hand annotation added when the figure was assembled and
#   it is not drawn here.
#
#   The per-mutation colours ARE reproduced, from MUTATION_COLOURS below.
#   The rule behind them is not recorded in any script, so they are listed
#   explicitly rather than derived.
# -----------------------------------------------------------------------------
#
# HOW TO RUN
#     Rscript scripts/06_figures/figure03a_k13_lollipop.R
#
#   Optional:
#     --min-samples N   keep changes seen in at least N samples (default 2)
#
################################################################################


# ==============================================================================
# 1. SETUP
# ==============================================================================

options(error = NULL)

# patchwork is easy to be without; name it rather than let a bare "there is
# no package called ..." surface part way down the script.
required_packages <- c(
  "dplyr", "tidyr", "readr", "stringr", "tibble", "ggplot2", "patchwork"
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

MINIMUM_SAMPLES <- suppressWarnings(
  as.integer(get_argument("--min-samples", "2"))
)

if (is.na(MINIMUM_SAMPLES) || MINIMUM_SAMPLES < 1) {
  stop("--min-samples must be a whole number of 1 or more.", call. = FALSE)
}

FIGURE_DIR <- file.path(PROJECT_ROOT, "results", "figures", "figure03")
TABLE_DIR  <- file.path(PROJECT_ROOT, "results", "tables",  "figure03")

dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR,  recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------------
# K13 protein domains. Boundaries are the ones used in the original; the
# colours and the text colours are those of the published track.
#
# The gaps between domains (342-349 and 438-442) are deliberate and are in
# the original definition: they are linker regions belonging to no domain.
# ------------------------------------------------------------------------------

DOMAINS <- tibble::tribble(
  ~domain,       ~start, ~end,  ~fill,      ~text_colour,
  "N-terminal",     1L,   210L, "#E8468A",  "#FFFFFF",
  "CCC",          211L,   341L, "#8B7FE0",  "#FFFFFF",
  "BTB/POZ",      350L,   437L, "#5B9BE8",  "#20304F",
  "Propeller",    443L,   726L, "#F5B942",  "#20304F"
)

PROTEIN_LENGTH <- 726L

# Colours read off published Figure 3a. Anything not listed is drawn in
# BACKBONE_COLOUR. See the note at the top: the rule behind the blue and
# grey emphasis is not recorded anywhere, so the list is literal.
MUTATION_COLOURS <- c(
  "K189T" = "#E8468A",
  "R561H" = "#F5A623",
  "H384R" = "#4A7FE8",
  "P441L" = "#4A7FE8",
  "P574L" = "#4A7FE8",
  "A675V" = "#AEB6CC"
)

BACKBONE_COLOUR <- "#5B6B8C"

# A change at or above this frequency carries its label at the top of its
# own stem instead of down in the label band.
TALL_THRESHOLD <- 30

# Label band. Heights are where the BOTTOM of the rotated text sits.
LABEL_RUNGS       <- c(8, 16, 24)
LABEL_MIN_SPACING <- 11     # codons between neighbouring labels
LABEL_CLUSTER_GAP <- 15     # closer than this counts as the same cluster
TALL_LABEL_GAP    <- 3      # percent above a tall point for its own label

# Rough height of one character of rotated label text, in percent units, at
# this panel's aspect ratio. Used only to leave headroom above the tallest
# label, so a loose estimate is fine.
TEXT_PERCENT_PER_CHARACTER <- 4.3


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
add_to_report("FIGURE 3A — K13 MUTATION SPECTRUM")
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
  stop("The input is missing: ", paste(missing_columns, collapse = ", "),
       call. = FALSE)
}

tidy <- raw %>%
  dplyr::mutate(
    SAMPLE       = stringr::str_trim(SAMPLE),
    GENE_NAME    = stringr::str_trim(GENE_NAME),
    HGVS_PROTEIN = stringr::str_trim(HGVS_PROTEIN)
  ) %>%
  dplyr::filter(!is.na(SAMPLE), SAMPLE != "")

n_all_samples <- dplyr::n_distinct(tidy$SAMPLE)

k13 <- tidy %>%
  dplyr::filter(
    GENE_NAME == "K13",
    !is.na(HGVS_PROTEIN), HGVS_PROTEIN != "", HGVS_PROTEIN != "."
  ) %>%
  dplyr::distinct(SAMPLE, HGVS_PROTEIN)

n_k13_samples <- dplyr::n_distinct(k13$SAMPLE)

if (n_k13_samples == 0) {
  stop("No K13 mutations found in the input.", call. = FALSE)
}

add_to_report(sprintf("Samples in the file:            %d", n_all_samples))
add_to_report(sprintf("Samples with a K13 mutation:    %d", n_k13_samples))
add_to_report(sprintf("Distinct K13 changes:           %d",
                      dplyr::n_distinct(k13$HGVS_PROTEIN)))
add_to_report("")
add_to_report("Frequencies below are over the K13-mutated samples.")


# ==============================================================================
# 3. COUNT AND POSITION
# ==============================================================================

frequencies <- k13 %>%
  dplyr::count(HGVS_PROTEIN, name = "n_samples") %>%
  dplyr::mutate(
    position = suppressWarnings(
      as.integer(stringr::str_extract(HGVS_PROTEIN, "[0-9]+"))
    ),
    percent_of_k13_mutated = 100 * n_samples / n_k13_samples,
    percent_of_all_samples = 100 * n_samples / n_all_samples
  ) %>%
  dplyr::arrange(dplyr::desc(n_samples))


unplaceable <- frequencies %>% dplyr::filter(is.na(position))

if (nrow(unplaceable) > 0) {
  add_to_report("")
  add_to_report(sprintf(
    "%d change(s) have no readable codon position and cannot be drawn:",
    nrow(unplaceable)
  ))
  for (index in seq_len(nrow(unplaceable))) {
    add_to_report(sprintf("  %-14s  %d sample(s)",
                          unplaceable$HGVS_PROTEIN[[index]],
                          unplaceable$n_samples[[index]]))
  }
}

out_of_range <- frequencies %>%
  dplyr::filter(!is.na(position), position < 1 | position > PROTEIN_LENGTH)

if (nrow(out_of_range) > 0) {
  add_to_report("")
  add_to_report(sprintf(
    "%d change(s) sit outside codons 1-%d and are probably mis-parsed:",
    nrow(out_of_range), PROTEIN_LENGTH
  ))
  for (index in seq_len(nrow(out_of_range))) {
    add_to_report(sprintf("  %-14s  position %d",
                          out_of_range$HGVS_PROTEIN[[index]],
                          out_of_range$position[[index]]))
  }
}


# ==============================================================================
# 4. DROP SINGLETONS
# ==============================================================================

excluded <- frequencies %>%
  dplyr::filter(n_samples < MINIMUM_SAMPLES)

drawn <- frequencies %>%
  dplyr::filter(n_samples >= MINIMUM_SAMPLES, !is.na(position)) %>%
  dplyr::arrange(position)

add_to_report("")
add_to_report(sprintf(
  "Changes seen in fewer than %d sample(s): %d excluded, %d drawn.",
  MINIMUM_SAMPLES, nrow(excluded), nrow(drawn)
))

readr::write_tsv(
  excluded %>% dplyr::arrange(dplyr::desc(n_samples), HGVS_PROTEIN),
  file.path(TABLE_DIR, "figure03a_excluded_singletons.tsv")
)

if (nrow(drawn) == 0) {
  stop("Nothing left to draw after excluding rare changes.", call. = FALSE)
}


# ==============================================================================
# 5. DOMAINS AND COLOURS
# ==============================================================================

assign_domain <- function(position) {

  hit <- DOMAINS$domain[position >= DOMAINS$start & position <= DOMAINS$end]

  if (length(hit) == 0) return(NA_character_)
  if (length(hit) > 1) {
    stop("Position ", position, " falls in more than one domain; the ",
         "boundaries in DOMAINS overlap.", call. = FALSE)
  }

  hit
}

drawn <- drawn %>%
  dplyr::rowwise() %>%
  dplyr::mutate(domain = assign_domain(position)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    colour = dplyr::coalesce(
      unname(MUTATION_COLOURS[HGVS_PROTEIN]), BACKBONE_COLOUR
    ),
    is_tall = percent_of_k13_mutated >= TALL_THRESHOLD
  )

in_linker <- drawn %>% dplyr::filter(is.na(domain))

if (nrow(in_linker) > 0) {
  add_to_report("")
  add_to_report(sprintf(
    "%d change(s) fall in a linker between domains (this is expected):",
    nrow(in_linker)
  ))
  add_to_report("  ", paste(in_linker$HGVS_PROTEIN, collapse = ", "))
}


# ==============================================================================
# 6. REPORT THE SPECTRUM
# ==============================================================================

add_to_report("")
add_to_report(strrep("-", 72))
add_to_report("K13 CHANGES DRAWN, MOST COMMON FIRST")
add_to_report(strrep("-", 72))
add_to_report(sprintf(
  "  %-12s %6s %8s %10s %10s  %s",
  "change", "pos", "n", "% K13mut", "% all", "domain"
))

ranked <- drawn %>% dplyr::arrange(dplyr::desc(n_samples))

for (index in seq_len(nrow(ranked))) {
  add_to_report(sprintf(
    "  %-12s %6d %8d %9.1f%% %9.1f%%  %s",
    ranked$HGVS_PROTEIN[[index]],
    ranked$position[[index]],
    ranked$n_samples[[index]],
    ranked$percent_of_k13_mutated[[index]],
    ranked$percent_of_all_samples[[index]],
    if (is.na(ranked$domain[[index]])) "(linker)" else ranked$domain[[index]]
  ))
}

top_two <- ranked$HGVS_PROTEIN[seq_len(min(2, nrow(ranked)))]

add_to_report("")

if (identical(sort(top_two), sort(c("K189T", "R561H")))) {
  add_to_report("  Top two changes are K189T and R561H, as in published Figure 3a.")
  add_to_report(sprintf(
    "    K189T  %5.1f%%    R561H  %5.1f%%",
    ranked$percent_of_k13_mutated[ranked$HGVS_PROTEIN == "K189T"],
    ranked$percent_of_k13_mutated[ranked$HGVS_PROTEIN == "R561H"]
  ))
} else {
  add_to_report(strrep("!", 72))
  add_to_report("  The two most common changes are ",
                paste(top_two, collapse = " and "),
                ", not K189T and R561H")
  add_to_report("  as shown in published Figure 3a. Check the input file.")
  add_to_report(strrep("!", 72))
}


# ==============================================================================
# 7. LABEL LAYOUT
# ==============================================================================
#
# Tall changes carry their label at the top of their own stem. Everything
# else goes into the label band: slid apart along the protein until no two
# labels are closer than LABEL_MIN_SPACING, then stepped up a rung for each
# further member of a cluster.
#
# Sliding only ever pushes a label to the right of a neighbour that is
# already placed, so the order of labels along the protein is preserved and
# a leader line never crosses another.
# ==============================================================================

drawn <- drawn %>% dplyr::arrange(position)

label_x <- as.numeric(drawn$position)
label_y <- numeric(nrow(drawn))

previous_position <- NA_real_
rung <- 1L
last_placed_x <- -Inf

for (index in seq_len(nrow(drawn))) {

  if (drawn$is_tall[[index]]) {
    label_y[[index]] <-
      drawn$percent_of_k13_mutated[[index]] + TALL_LABEL_GAP
    next
  }

  this_x <- max(as.numeric(drawn$position[[index]]),
                last_placed_x + LABEL_MIN_SPACING)

  rung <- if (is.na(previous_position) ||
              drawn$position[[index]] - previous_position >= LABEL_CLUSTER_GAP) {
    1L
  } else {
    (rung %% length(LABEL_RUNGS)) + 1L
  }

  label_x[[index]] <- this_x
  label_y[[index]] <- LABEL_RUNGS[[rung]]

  last_placed_x     <- this_x
  previous_position <- drawn$position[[index]]
}

drawn$label_x <- label_x
drawn$label_y <- label_y
drawn$label_shift <- drawn$label_x - drawn$position

band <- drawn[!drawn$is_tall, ]

if (nrow(band) > 1) {

  gaps <- diff(band$label_x)
  tight <- which(gaps < LABEL_MIN_SPACING - 1e-9)

  add_to_report("")

  if (length(tight) == 0) {
    add_to_report(sprintf(
      "Label layout: %d labels, all at least %.0f codons apart.",
      nrow(band), LABEL_MIN_SPACING
    ))
    add_to_report(sprintf(
      "  Largest sideways shift: %.0f codons (%s).",
      max(abs(band$label_shift)),
      band$HGVS_PROTEIN[[which.max(abs(band$label_shift))]]
    ))
  } else {
    add_to_report(sprintf("Label layout: %d pair(s) still too close:",
                          length(tight)))
    for (index in tight) {
      add_to_report("  ", band$HGVS_PROTEIN[[index]], " and ",
                    band$HGVS_PROTEIN[[index + 1]])
    }
  }
}


# ==============================================================================
# 8. THE FIGURE
# ==============================================================================

longest_label <- max(nchar(drawn$HGVS_PROTEIN))

y_top <- max(
  90,
  max(drawn$label_y) + longest_label * TEXT_PERCENT_PER_CHARACTER + 3
)

X_BREAKS <- seq(0, 700, by = 100)
X_RANGE  <- c(-12, PROTEIN_LENGTH + 18)


lollipop_panel <- ggplot2::ggplot(drawn) +

  ggplot2::geom_hline(
    yintercept = seq(20, 80, by = 20),
    colour = "grey88", linewidth = 0.4
  ) +

  # The protein backbone. Low-frequency changes sit on it rather than
  # floating above an invisible baseline.
  ggplot2::annotate(
    "segment",
    x = 1, xend = PROTEIN_LENGTH, y = 0, yend = 0,
    colour = BACKBONE_COLOUR, linewidth = 2.4, lineend = "round"
  ) +

  ggplot2::geom_segment(
    ggplot2::aes(
      x = position, xend = position,
      y = 0, yend = percent_of_k13_mutated, colour = colour
    ),
    linewidth = 1.0
  ) +

  ggplot2::geom_point(
    ggplot2::aes(x = position, y = percent_of_k13_mutated, colour = colour),
    size = 3.0
  ) +

  # Leaders, drawn only where the label had to move sideways.
  ggplot2::geom_segment(
    data = drawn[!drawn$is_tall & abs(drawn$label_shift) > 0.5, ],
    mapping = ggplot2::aes(
      x = position, y = percent_of_k13_mutated,
      xend = label_x, yend = label_y
    ),
    colour = "grey72", linewidth = 0.45
  ) +

  ggplot2::geom_text(
    ggplot2::aes(
      x = label_x, y = label_y, label = HGVS_PROTEIN, colour = colour
    ),
    angle = 90, hjust = 0, vjust = 0.5,
    size = 4.2, fontface = "bold"
  ) +

  ggplot2::scale_colour_identity() +

  ggplot2::scale_x_continuous(
    breaks = X_BREAKS, limits = X_RANGE,
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +

  ggplot2::scale_y_continuous(
    breaks = seq(20, 80, by = 20),
    limits = c(-2, y_top),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +

  ggplot2::labs(x = "Amino Acid Position", y = "Mutation Frequency (%)") +

  ggplot2::theme_classic(base_size = 13) +

  ggplot2::theme(
    axis.title.x = ggplot2::element_text(
      face = "bold", size = 17, colour = "#4A4A4A",
      margin = ggplot2::margin(t = 8)
    ),
    axis.title.y = ggplot2::element_text(
      face = "bold", size = 17, colour = "#2B2B2B",
      margin = ggplot2::margin(r = 8)
    ),
    axis.text    = ggplot2::element_text(size = 13, colour = "#6E6E6E"),
    axis.ticks   = ggplot2::element_blank(),
    axis.line    = ggplot2::element_blank(),
    panel.grid   = ggplot2::element_blank(),
    plot.margin  = ggplot2::margin(12, 16, 2, 10)
  )


# ------------------------------------------------------------------------------
# The domain track. Each domain is a block with a pointed right-hand end,
# built as a five-sided polygon so the shape travels with the data rather
# than being drawn by hand.
# ------------------------------------------------------------------------------

ARROW_TIP <- 16   # codons taken up by the pointed end

make_arrow <- function(start, end, domain, fill) {

  shoulder <- max(start, end - ARROW_TIP)

  tibble::tibble(
    domain = domain,
    fill   = fill,
    x = c(start, shoulder, end, shoulder, start),
    y = c(0.18,  0.18,     0.5, 0.82,     0.82)
  )
}

domain_polygons <- do.call(
  rbind,
  lapply(
    seq_len(nrow(DOMAINS)),
    function(index) {
      make_arrow(
        DOMAINS$start[[index]], DOMAINS$end[[index]],
        DOMAINS$domain[[index]], DOMAINS$fill[[index]]
      )
    }
  )
)

domain_labels <- DOMAINS %>%
  dplyr::mutate(centre = (start + end) / 2)

domain_panel <- ggplot2::ggplot() +

  ggplot2::geom_polygon(
    data = domain_polygons,
    mapping = ggplot2::aes(x = x, y = y, group = domain, fill = fill)
  ) +

  ggplot2::geom_text(
    data = domain_labels,
    mapping = ggplot2::aes(x = centre, y = 0.5, label = domain,
                           colour = text_colour),
    fontface = "bold", size = 5.2
  ) +

  ggplot2::scale_fill_identity() +
  ggplot2::scale_colour_identity() +

  ggplot2::scale_x_continuous(
    breaks = X_BREAKS, limits = X_RANGE,
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +

  ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +

  ggplot2::labs(x = NULL, y = "K13 domains") +

  ggplot2::theme_classic(base_size = 13) +

  ggplot2::theme(
    axis.title.y = ggplot2::element_text(
      angle = 0, vjust = 0.5, hjust = 1, face = "bold",
      size = 16, colour = "#2B2B2B",
      margin = ggplot2::margin(r = 8)
    ),
    axis.text    = ggplot2::element_blank(),
    axis.ticks   = ggplot2::element_blank(),
    axis.line    = ggplot2::element_blank(),
    panel.grid   = ggplot2::element_blank(),
    plot.margin  = ggplot2::margin(2, 16, 10, 10)
  )


figure <- lollipop_panel / domain_panel +
  patchwork::plot_layout(heights = c(6, 1.15))

for (extension in c("svg", "png")) {
  ggplot2::ggsave(
    filename = file.path(
      FIGURE_DIR, paste0("figure03a_k13_lollipop.", extension)
    ),
    plot = figure,
    width = 15, height = 6.5, units = "in",
    dpi = 600, bg = "white"
  )
}


# ==============================================================================
# 9. TABLES
# ==============================================================================

readr::write_tsv(
  drawn %>%
    dplyr::select(
      change = HGVS_PROTEIN, position, domain, n_samples,
      percent_of_k13_mutated, percent_of_all_samples
    ) %>%
    dplyr::arrange(position),
  file.path(TABLE_DIR, "figure03a_mutation_frequencies.tsv")
)

readr::write_tsv(
  DOMAINS %>% dplyr::select(domain, start, end),
  file.path(TABLE_DIR, "figure03a_domains.tsv")
)

add_to_report("")
add_to_report(strrep("=", 72))
add_to_report("WRITTEN")
add_to_report(strrep("=", 72))
add_to_report("  figures/figure03/figure03a_k13_lollipop.svg / .png")
add_to_report("  tables/figure03/figure03a_mutation_frequencies.tsv")
add_to_report("  tables/figure03/figure03a_domains.tsv")
add_to_report("  tables/figure03/figure03a_excluded_singletons.tsv")
add_to_report("")
add_to_report(sprintf(
  "Denominator for the drawn percentages: %d K13-mutated samples.",
  n_k13_samples
))
add_to_report(paste0(
  "The table also carries the percentage of all ", n_all_samples, " samples."
))

writeLines(report_lines, file.path(TABLE_DIR, "figure03a_report.txt"))

capture.output(
  sessionInfo(),
  file = file.path(TABLE_DIR, "figure03a_sessionInfo.txt")
)
