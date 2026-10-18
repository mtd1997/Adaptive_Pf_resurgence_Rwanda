################################################################################
# FIGURE 5 — SELECTION-SCAN PANELS
################################################################################
#
# WHAT IT DOES
#   Draws the three selection-scan panels from the statistics table produced
#   by figure05_selection_scan_statistics.R:
#
#     A   per-SNP F_ST between the two periods
#     B   absolute allele-frequency change, |ΔAF|
#     C   signed allele-frequency change, ΔAF
#
#   It computes nothing. Every number it draws — frequencies, F_ST,
#   thresholds, outlier flags, resistance flags, labels — is read from the
#   table. If a value looks wrong, the statistics script is where to look.
#
# IMPORTANT — THIS IS NOT THE PUBLISHED FIGURE
#   The Figure 5 that appears in the paper was assembled by hand in Adobe
#   Illustrator from script-generated SVG panels, to fix legend placement and
#   panel proportions that could not be resolved in ggplot2. This script
#   produces the panels that went into that assembly. It does not reproduce
#   the final composed figure, and no script does.
#
#   This is normal practice for publication figures and is stated plainly
#   rather than left for a reader to discover. The underlying data and
#   statistics are fully reproducible; only the final visual arrangement was
#   manual. See scripts/06_figures/README.md for the full account.
#
#   One of the two SVGs used in the assembly came from an earlier variant of
#   this script that was edited in place and not kept. The surviving version
#   is the one with the hand-built legend and matched panel widths, which is
#   what is published here.
#
# WHY THE LEGEND IS DRAWN BY HAND
#   ggplot2's automatic guides could not produce, at the same time, coloured
#   gene symbols, two distinct point shapes, and a visibly dashed threshold
#   key. Sections 7 and 8 therefore draw the legend as its own small plot and
#   place it in a fixed right-hand column. That column is present in all
#   three rows — empty for panel A — which is what keeps the three data
#   panels exactly the same width.
#
# OUTPUTS  (all under results/figures/figure05/)
#   figure05_selection_scan.svg / .pdf / .png   the three panels together
#   figure05_panelA_fst.svg                     individual panels, for
#   figure05_panelB_abs_delta_af.svg            assembly in a vector editor
#   figure05_panelC_signed_delta_af.svg
#   figure05_plotting_table.tsv                 every value drawn, including
#                                               the x positions — the Source
#                                               Data for Figure 5
#   figure05_labelled_variants.tsv              only the labelled points
#   figure05_gene_accounting.tsv                variants and outliers per gene
#   sessionInfo.txt
#
# HOW TO RUN
#   From the top of the repository:
#
#     Rscript scripts/06_figures/figure05_selection_scan_panels.R
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
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(grid)
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
    "Could not find config.R. Run this script from the top of the ",
    "repository.",
    call. = FALSE
  )
}

PROJECT_ROOT <- find_project_root()

source(file.path(PROJECT_ROOT, "config.R"))

STATISTICS_FILE <- file.path(
  PROJECT_ROOT, "results", "tables", "figure05",
  "selection_statistics_all_variants.tsv"
)

OUTPUT_DIR <- file.path(PROJECT_ROOT, "results", "figures", "figure05")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(STATISTICS_FILE)) {
  stop(
    "The statistics table was not found:\n  ", STATISTICS_FILE,
    "\n\nRun this first:\n",
    "  Rscript scripts/06_figures/figure05_selection_scan_statistics.R",
    call. = FALSE
  )
}


# ==============================================================================
# 2. READ THE STATISTICS TABLE
# ==============================================================================

selection_data <- readr::read_tsv(STATISTICS_FILE, show_col_types = FALSE)

required_columns <- c(
  "variant_id", "snp_id", "chrom", "pos", "ref", "alt",
  "gene", "snp_label", "biological_label", "is_resistance",
  "AF_Baseline_Kigali", "AF_Resurgence",
  "fst", "fst_threshold_95", "fst_outlier",
  "delta_AF", "abs_delta_AF", "abs_delta_threshold_95", "delta_outlier"
)

missing_columns <- setdiff(required_columns, names(selection_data))

if (length(missing_columns) > 0) {
  stop(
    "The statistics table is missing: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

selection_data <- selection_data %>%
  dplyr::mutate(
    pos              = as.numeric(pos),
    gene             = as.character(gene),
    snp_label        = as.character(snp_label),
    biological_label = as.character(biological_label),
    dplyr::across(
      c(is_resistance, fst_outlier, delta_outlier),
      ~ tidyr::replace_na(as.logical(.x), FALSE)
    )
  )


# ==============================================================================
# 3. GENE ORDER, COLOURS AND SHAPES
# ==============================================================================
#
# Genes appear left to right in the order the amplicons were designed.
#
# CRT is deliberately absent: no CRT variant survived filtering into the
# selection-scan VCF. The accounting table records that as an explicit zero
# rather than leaving it to be inferred from a missing panel section.
#
# The palette is the Okabe-Ito set, which stays distinguishable under the
# common forms of colour-vision deficiency.
# ==============================================================================

GENE_ORDER <- c("MSP2", "CSP", "DHFR", "MDR1", "DHPS", "K13")

GENE_COLOURS <- c(
  "MSP2" = "#56B4E9",   # sky blue
  "CSP"  = "#E69F00",   # orange
  "DHFR" = "#009E73",   # bluish green
  "MDR1" = "#D55E00",   # vermillion
  "DHPS" = "#0072B2",   # blue
  "K13"  = "#CC79A7"    # reddish purple
)

VARIANT_SHAPES <- c(
  "Other variant"                 = 16,   # filled circle
  "Resistance-associated variant" = 18    # filled diamond
)

THRESHOLD_COLOUR <- "#D62728"

genes_present <- intersect(GENE_ORDER, unique(selection_data$gene))

genes_unexpected <- setdiff(
  unique(stats::na.omit(selection_data$gene)),
  GENE_ORDER
)

if (length(genes_unexpected) > 0) {
  warning(
    "These genes are in the statistics table but not in GENE_ORDER, so they ",
    "will not be drawn: ",
    paste(genes_unexpected, collapse = ", "),
    call. = FALSE
  )
}

selection_data <- selection_data %>%
  dplyr::mutate(gene = factor(gene, levels = GENE_ORDER)) %>%
  dplyr::filter(!is.na(gene))

message(sprintf(
  "Variants drawn: %d across %d genes (%s)",
  nrow(selection_data),
  length(genes_present),
  paste(genes_present, collapse = ", ")
))


# ==============================================================================
# 4. THRESHOLDS
# ==============================================================================
#
# Read from the table rather than recalculated, so the dashed lines cannot
# drift away from the outlier flags they correspond to.
# ==============================================================================

read_single_threshold <- function(values, description) {

  unique_values <- unique(stats::na.omit(values))

  if (length(unique_values) != 1) {
    stop(
      "Expected exactly one ", description, " threshold in the table but ",
      "found ", length(unique_values), ". The statistics script writes one ",
      "value to every row; more than one means the table has been edited ",
      "or two runs have been concatenated.",
      call. = FALSE
    )
  }

  as.numeric(unique_values)
}

FST_THRESHOLD <- read_single_threshold(
  selection_data$fst_threshold_95, "F_ST"
)

DELTA_THRESHOLD <- read_single_threshold(
  selection_data$abs_delta_threshold_95, "absolute delta-AF"
)

message(sprintf("F_ST threshold:  %.6f", FST_THRESHOLD))
message(sprintf("|ΔAF| threshold: %.6f", DELTA_THRESHOLD))


# ==============================================================================
# 5. X POSITIONS
# ==============================================================================
#
# The seven amplicons sit on different chromosomes, megabases apart, so real
# genomic coordinates cannot share one axis. Each gene is given its own block
# instead: positions within a gene keep their true spacing, and a fixed gap
# separates one gene from the next.
#
# cumulative_pos is therefore a drawing coordinate, not a genomic one. The
# true chromosome and position stay in the `chrom` and `pos` columns and are
# written to the plotting table.
# ==============================================================================

GENE_GAP <- 25000

gene_ranges <- selection_data %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    gene_min = min(pos, na.rm = TRUE),
    gene_max = max(pos, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  dplyr::arrange(gene) %>%
  dplyr::mutate(
    # A gene with a single variant has zero width; give it one base so the
    # cumulative offsets below stay strictly increasing.
    gene_width = pmax(gene_max - gene_min, 1),
    offset     = dplyr::lag(cumsum(gene_width + GENE_GAP), default = 0)
  )

plot_data <- selection_data %>%
  dplyr::left_join(gene_ranges, by = "gene") %>%
  dplyr::mutate(cumulative_pos = (pos - gene_min) + offset)

gene_centres <- plot_data %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    centre = (min(cumulative_pos, na.rm = TRUE) +
              max(cumulative_pos, na.rm = TRUE)) / 2,
    .groups = "drop"
  )

gene_bands <- plot_data %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    xmin = min(cumulative_pos, na.rm = TRUE) - GENE_GAP / 2,
    xmax = max(cumulative_pos, na.rm = TRUE) + GENE_GAP / 2,
    .groups = "drop"
  ) %>%
  dplyr::arrange(gene) %>%
  dplyr::mutate(shade_band = dplyr::row_number() %% 2 == 0)

gene_separators <- gene_bands %>%
  dplyr::arrange(xmin) %>%
  dplyr::slice(-1) %>%
  dplyr::transmute(separator = xmin)

# One x-range shared by all three panels. This is what makes them line up.
SHARED_X_MIN <- min(gene_bands$xmin, na.rm = TRUE)
SHARED_X_MAX <- max(gene_bands$xmax, na.rm = TRUE)


# ==============================================================================
# 6. WHICH POINTS ARE EMPHASISED AND LABELLED
# ==============================================================================
#
# A point is enlarged and labelled if it is either a known resistance marker
# or an empirical outlier on the statistic that panel shows. Resistance
# markers are labelled in bold so they can be told from outliers at a glance.
# ==============================================================================

plot_data <- plot_data %>%
  dplyr::mutate(
    variant_category = factor(
      dplyr::if_else(
        is_resistance,
        "Resistance-associated variant",
        "Other variant"
      ),
      levels = names(VARIANT_SHAPES)
    ),
    emphasise_fst   = fst_outlier   | is_resistance,
    emphasise_delta = delta_outlier | is_resistance,
    label_face      = dplyr::if_else(is_resistance, "bold", "plain")
  )


# ==============================================================================
# 7. SHARED THEME
# ==============================================================================

theme_selection <- ggplot2::theme_classic(base_size = 20) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold", size = 22, colour = "black",
      margin = ggplot2::margin(b = 5)
    ),
    plot.subtitle = ggplot2::element_text(
      size = 16, colour = "grey10", lineheight = 1.08,
      margin = ggplot2::margin(b = 10)
    ),
    plot.tag = ggplot2::element_text(
      face = "bold", size = 34, colour = "black"
    ),
    plot.tag.position = c(-0.055, 1.045),
    axis.title.y = ggplot2::element_text(
      face = "bold", size = 22, colour = "black",
      margin = ggplot2::margin(r = 13)
    ),
    axis.title.x = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      face = "bold", size = 16, colour = "black"
    ),
    axis.text.y = ggplot2::element_text(
      face = "bold", size = 15, colour = "black"
    ),
    axis.ticks.x = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_line(colour = "black", linewidth = 0.6),
    panel.grid = ggplot2::element_blank(),
    legend.position = "none",
    plot.margin = ggplot2::margin(t = 20, r = 16, b = 16, l = 40)
  )


# ==============================================================================
# 8. ONE FUNCTION THAT BUILDS ANY OF THE THREE PANELS
# ==============================================================================
#
# All three panels share their structure: shaded gene bands, separators,
# threshold lines, small background points, enlarged emphasised points, and
# repelled labels. The only differences are which column goes on the y axis
# and whether the threshold is drawn once or mirrored above and below zero.
# ==============================================================================

build_panel <- function(
  y_column,
  emphasise_column,
  threshold,
  y_label,
  title,
  subtitle,
  panel_tag,
  signed = FALSE
) {

  panel_data <- plot_data %>%
    dplyr::filter(
      !is.na(.data[[y_column]]),
      !is.na(cumulative_pos)
    )

  plot_object <- ggplot2::ggplot(
    panel_data,
    ggplot2::aes(x = cumulative_pos, y = .data[[y_column]])
  ) +

    # Alternating grey bands, so the eye can tell one gene from the next.
    ggplot2::geom_rect(
      data = dplyr::filter(gene_bands, shade_band),
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "grey94", alpha = 0.7
    ) +

    ggplot2::geom_vline(
      data = gene_separators,
      ggplot2::aes(xintercept = separator),
      inherit.aes = FALSE,
      colour = "grey82", linewidth = 0.55
    )

  # The empirical threshold. Panel C is signed, so it needs the line
  # mirrored, plus a zero reference.
  if (signed) {

    plot_object <- plot_object +
      ggplot2::geom_hline(
        yintercept = 0, colour = "grey25", linewidth = 0.85
      ) +
      ggplot2::geom_hline(
        yintercept = c(threshold, -threshold),
        colour = THRESHOLD_COLOUR, linetype = "dashed",
        linewidth = 0.9, lineend = "butt"
      )

  } else {

    plot_object <- plot_object +
      ggplot2::geom_hline(
        yintercept = threshold,
        colour = THRESHOLD_COLOUR, linetype = "dashed",
        linewidth = 0.9, lineend = "butt"
      )
  }

  plot_object +

    # Background points: small and semi-transparent, so the emphasised ones
    # read clearly without hiding the overall distribution.
    ggplot2::geom_point(
      data = dplyr::filter(panel_data, !.data[[emphasise_column]]),
      ggplot2::aes(colour = gene, shape = variant_category),
      size = 2.8, alpha = 0.50, stroke = 0.9
    ) +

    ggplot2::geom_point(
      data = dplyr::filter(panel_data, .data[[emphasise_column]]),
      ggplot2::aes(colour = gene, shape = variant_category),
      size = 5.5, alpha = 1, stroke = 1.1
    ) +

    # A fixed seed keeps the label positions identical between runs, so
    # regenerating the figure does not shuffle them.
    ggrepel::geom_text_repel(
      data = dplyr::filter(panel_data, .data[[emphasise_column]]),
      ggplot2::aes(label = snp_label, colour = gene, fontface = label_face),
      size = 4.7,
      max.overlaps = Inf,
      segment.size = 0.50,
      segment.alpha = 0.70,
      box.padding = 0.80,
      point.padding = 0.50,
      min.segment.length = 0,
      show.legend = FALSE,
      seed = 12345
    ) +

    ggplot2::scale_colour_manual(values = GENE_COLOURS, drop = FALSE) +
    ggplot2::scale_shape_manual(values = VARIANT_SHAPES, drop = FALSE) +

    ggplot2::scale_x_continuous(
      breaks = gene_centres$centre,
      labels = as.character(gene_centres$gene),
      limits = c(SHARED_X_MIN, SHARED_X_MAX),
      expand = ggplot2::expansion(mult = c(0.01, 0.01))
    ) +

    # clip = "off" lets the bold panel letter sit outside the plotting area.
    ggplot2::coord_cartesian(clip = "off") +

    ggplot2::labs(
      tag = panel_tag, title = title, subtitle = subtitle,
      x = NULL, y = y_label
    ) +

    theme_selection
}


# ==============================================================================
# 9. THE THREE PANELS
# ==============================================================================

panel_a <- build_panel(
  y_column         = "fst",
  emphasise_column = "emphasise_fst",
  threshold        = FST_THRESHOLD,
  y_label          = expression(bolditalic(F)[ST]),
  title            = paste(
    "Per-SNP differentiation between Baseline-Kigali and Resurgence"
  ),
  subtitle         = sprintf(
    paste(
      "The red dashed line is the empirical 95th percentile of per-SNP",
      "F_ST (%.3f); points above it are the top 5%%."
    ),
    FST_THRESHOLD
  ),
  panel_tag        = "A"
)

panel_b <- build_panel(
  y_column         = "abs_delta_AF",
  emphasise_column = "emphasise_delta",
  threshold        = DELTA_THRESHOLD,
  y_label          = expression(bold("|" * Delta * "AF|")),
  title            = "Absolute allele-frequency change between periods",
  subtitle         = sprintf(
    paste(
      "Absolute allele-frequency change: the magnitude of change",
      "irrespective of direction. The red dashed line is the empirical",
      "95th percentile (%.3f)."
    ),
    DELTA_THRESHOLD
  ),
  panel_tag        = "B"
)

panel_c <- build_panel(
  y_column         = "delta_AF",
  emphasise_column = "emphasise_delta",
  threshold        = DELTA_THRESHOLD,
  y_label          = expression(
    bold(Delta * "AF (Resurgence - Baseline-Kigali)")
  ),
  title            = paste(
    "Signed allele-frequency change from Baseline-Kigali to Resurgence"
  ),
  subtitle         = sprintf(
    paste(
      "Signed change, AF(Resurgence) minus AF(Baseline-Kigali).",
      "Positive values increased, negative values decreased.",
      "Red dashed lines are plus and minus %.3f."
    ),
    DELTA_THRESHOLD
  ),
  panel_tag        = "C",
  signed           = TRUE
)


# ==============================================================================
# 10. THE HAND-BUILT LEGEND
# ==============================================================================
#
# Drawn as its own plot on an invisible 0-1 canvas. Nothing here is mapped
# from the data, so the entries must be kept consistent with sections 3 and 4
# by hand — that is the cost of escaping ggplot2's automatic guides.
# ==============================================================================

legend_gene_entries <- tibble::tibble(
  label  = GENE_ORDER,
  y      = seq(0.90, 0.52, length.out = length(GENE_ORDER)),
  colour = unname(GENE_COLOURS[GENE_ORDER])
)

legend_panel <- ggplot2::ggplot() +

  ggplot2::annotate(
    "text", x = 0.05, y = 0.98, label = "Gene",
    hjust = 0, fontface = "bold", size = 6.0
  ) +

  ggplot2::geom_point(
    data = legend_gene_entries,
    ggplot2::aes(x = 0.10, y = y),
    colour = legend_gene_entries$colour, size = 5.0
  ) +

  ggplot2::geom_text(
    data = legend_gene_entries,
    ggplot2::aes(x = 0.19, y = y, label = label),
    hjust = 0, size = 5.0
  ) +

  ggplot2::annotate(
    "text", x = 0.05, y = 0.43, label = "Variant category",
    hjust = 0, fontface = "bold", size = 6.0
  ) +

  ggplot2::annotate(
    "point", x = 0.10, y = 0.34,
    shape = VARIANT_SHAPES[["Other variant"]], size = 5.0, colour = "black"
  ) +

  ggplot2::annotate(
    "text", x = 0.19, y = 0.34, label = "Other variant",
    hjust = 0, size = 5.0
  ) +

  ggplot2::annotate(
    "point", x = 0.10, y = 0.25,
    shape = VARIANT_SHAPES[["Resistance-associated variant"]],
    size = 6.0, colour = "black"
  ) +

  ggplot2::annotate(
    "text", x = 0.19, y = 0.25, label = "Resistance-associated variant",
    hjust = 0, size = 5.0
  ) +

  ggplot2::annotate(
    "text", x = 0.05, y = 0.15, label = "Threshold",
    hjust = 0, fontface = "bold", size = 6.0
  ) +

  ggplot2::annotate(
    "segment", x = 0.07, xend = 0.31, y = 0.075, yend = 0.075,
    colour = THRESHOLD_COLOUR, linetype = "dashed", linewidth = 1.1
  ) +

  ggplot2::annotate(
    "text", x = 0.36, y = 0.075,
    label = "Empirical 95th-percentile threshold",
    hjust = 0, size = 4.8
  ) +

  ggplot2::xlim(0, 1.35) +
  ggplot2::ylim(0, 1) +
  ggplot2::theme_void() +
  ggplot2::theme(plot.margin = ggplot2::margin(t = 5, r = 10, b = 5, l = 5))


# ------------------------------------------------------------------------------
# Panel C's direction key, kept outside the data area so the boxes cannot
# land on top of a point.
# ------------------------------------------------------------------------------

direction_panel <- ggplot2::ggplot() +

  ggplot2::annotate(
    "label", x = 0.05, y = 0.67,
    label = "Positive values\nIncreased in Resurgence",
    hjust = 0, size = 5.2, fontface = "bold",
    fill = "white", colour = "black", linewidth = 0.45,
    label.padding = grid::unit(0.32, "lines")
  ) +

  ggplot2::annotate(
    "label", x = 0.05, y = 0.35,
    label = "Negative values\nDecreased in Resurgence",
    hjust = 0, size = 5.2, fontface = "bold",
    fill = "white", colour = "black", linewidth = 0.45,
    label.padding = grid::unit(0.32, "lines")
  ) +

  ggplot2::xlim(0, 1.35) +
  ggplot2::ylim(0, 1) +
  ggplot2::theme_void() +
  ggplot2::theme(plot.margin = ggplot2::margin(t = 20, r = 10, b = 5, l = 5))


# ==============================================================================
# 11. ASSEMBLE
# ==============================================================================
#
# A fixed two-column grid. Every row has the same column widths, so the three
# data panels are exactly as wide as one another. Panel A's right-hand cell is
# empty; B's holds the legend; C's holds the direction key.
# ==============================================================================

figure05 <- (
  (panel_a | patchwork::plot_spacer()) /
  (panel_b | legend_panel) /
  (panel_c | direction_panel)
) +
  patchwork::plot_layout(
    widths  = c(5.0, 1.75),
    heights = c(1, 1, 1)
  ) +
  patchwork::plot_annotation(
    title = paste(
      "Selection scan of P. falciparum populations in Rwanda, 2021-2025"
    ),
    subtitle = paste(
      "Colours identify genes; shapes identify variant categories.",
      "The scan contains one MDR1 variant and no CRT variants after",
      "filtering."
    ),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", size = 25, colour = "black",
        margin = ggplot2::margin(b = 6)
      ),
      plot.subtitle = ggplot2::element_text(
        size = 16.5, colour = "grey10",
        margin = ggplot2::margin(b = 14)
      ),
      plot.margin = ggplot2::margin(t = 16, r = 20, b = 12, l = 28)
    )
  )


# ==============================================================================
# 12. SAVE
# ==============================================================================
#
# SVG is the format to take into a vector editor. PDF and PNG are for
# submission and for quick viewing. limitsize = FALSE is needed because the
# figure is larger than ggsave's default sanity limit.
# ==============================================================================

FIGURE_WIDTH  <- 19
FIGURE_HEIGHT <- 16.5

# ------------------------------------------------------------------------------
# A NOTE ON THE GREEK DELTA, BECAUSE IT WILL BITE AGAIN
#
# R's PDF devices cannot write a Greek delta into ordinary text. Both
# grDevices::pdf() and cairo_pdf() fail with
#
#   conversion failure on 'Positive ΔAF' in 'mbcsToSbcs': for Δ (U+0394)
#
# and drop the character silently, leaving a gap in the figure. This happens
# even when capabilities("cairo") reports TRUE, so checking for Cairo does
# not help. SVG and PNG are unaffected.
#
# The fix is not a different device. It is to keep the delta out of plain
# text:
#
#   - Axis titles use plotmath, expression(Delta * "AF ..."), which draws the
#     symbol from the maths font and works on every device. Those are fine.
#   - Subtitles and annotation labels are plain strings, so they spell the
#     quantity out in words instead.
#
# If you add a subtitle later, write "allele-frequency change" rather than
# pasting a delta, or the PDF will lose it without telling you.
# ------------------------------------------------------------------------------

save_figure <- function(plot_object, basename, width, height) {

  for (extension in c("svg", "pdf", "png")) {

    ggplot2::ggsave(
      filename  = file.path(OUTPUT_DIR, paste0(basename, ".", extension)),
      plot      = plot_object,
      width     = width,
      height    = height,
      units     = "in",
      dpi       = 600,
      limitsize = FALSE
    )
  }
}

save_figure(figure05, "figure05_selection_scan", FIGURE_WIDTH, FIGURE_HEIGHT)

# Individual panels as SVG, which is what the manual assembly used.
panel_files <- list(
  figure05_panelA_fst            = panel_a,
  figure05_panelB_abs_delta_af   = panel_b,
  figure05_panelC_signed_delta_af = panel_c
)

for (panel_name in names(panel_files)) {

  ggplot2::ggsave(
    filename = file.path(OUTPUT_DIR, paste0(panel_name, ".svg")),
    plot     = panel_files[[panel_name]],
    width    = 15.5,
    height   = 6.2,
    units    = "in",
    dpi      = 600
  )
}


# ==============================================================================
# 13. SOURCE DATA
# ==============================================================================
#
# Everything the figure draws, in one table. cumulative_pos is included so a
# reader can reconstruct the x axis exactly; chrom and pos give the real
# coordinates.
# ==============================================================================

plotting_table <- plot_data %>%
  dplyr::select(
    variant_id, snp_id, gene, chrom, pos, ref, alt, cumulative_pos,
    snp_label, biological_label,
    is_resistance, variant_category,
    AF_Baseline_Kigali, AF_Resurgence,
    fst, fst_threshold_95, fst_outlier, emphasise_fst,
    delta_AF, abs_delta_AF, abs_delta_threshold_95,
    delta_outlier, emphasise_delta
  )

readr::write_tsv(
  plotting_table,
  file.path(OUTPUT_DIR, "figure05_plotting_table.tsv")
)

readr::write_tsv(
  plotting_table %>%
    dplyr::filter(emphasise_fst | emphasise_delta) %>%
    dplyr::arrange(gene, dplyr::desc(fst), dplyr::desc(abs_delta_AF)),
  file.path(OUTPUT_DIR, "figure05_labelled_variants.tsv")
)


# ------------------------------------------------------------------------------
# Gene accounting, with CRT recorded as an explicit zero.
# ------------------------------------------------------------------------------

gene_accounting <- selection_data %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    N_variants            = dplyr::n(),
    N_fst_outliers        = sum(fst_outlier, na.rm = TRUE),
    N_delta_af_outliers   = sum(delta_outlier, na.rm = TRUE),
    N_resistance_variants = sum(is_resistance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(gene = as.character(gene))

genes_with_no_variants <- setdiff(
  c(GENE_ORDER, "CRT"),
  gene_accounting$gene
)

if (length(genes_with_no_variants) > 0) {

  gene_accounting <- dplyr::bind_rows(
    gene_accounting,
    tibble::tibble(
      gene                  = genes_with_no_variants,
      N_variants            = 0L,
      N_fst_outliers        = 0L,
      N_delta_af_outliers   = 0L,
      N_resistance_variants = 0L
    )
  )
}

readr::write_tsv(
  gene_accounting,
  file.path(OUTPUT_DIR, "figure05_gene_accounting.tsv")
)

capture.output(
  sessionInfo(),
  file = file.path(OUTPUT_DIR, "sessionInfo.txt")
)


# ==============================================================================
# 14. SUMMARY
# ==============================================================================

message("\n", strrep("=", 72))
message("FIGURE 5 PANELS WRITTEN")
message(strrep("=", 72))
message(sprintf("Output directory: %s", OUTPUT_DIR))

message("\nGene accounting:")
print(as.data.frame(gene_accounting), row.names = FALSE)

message("\nLabelled points:")
print(
  as.data.frame(
    plotting_table %>%
      dplyr::filter(emphasise_fst | emphasise_delta) %>%
      dplyr::select(gene, snp_label, is_resistance, fst, delta_AF) %>%
      dplyr::arrange(gene, dplyr::desc(fst))
  ),
  row.names = FALSE
)

message(
  "\nReminder: the published Figure 5 was composed by hand in Adobe ",
  "Illustrator from these panels. This script does not reproduce the final ",
  "composed figure. See scripts/06_figures/README.md."
)
message(strrep("=", 72))
